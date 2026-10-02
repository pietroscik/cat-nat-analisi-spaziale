#!/usr/bin/env python3
"""Curva di eccedenza sismica (EP) e AAL numerico dalla curva di hazard MPS04.

Il benchmark EAL (pricing_model.py) usa CURVE = 4: punto di design RP475 con
MDR ~ ag^2, moltiplicato per l'integrazione oltre il design point (con pendenza
lambda(ag) ~ ag^(-k), k = 3: 1 + k/(k-2) = 4). La matrice estesa porta tre punti
della curva di pericolosita' per comune (ag a RP30, RP72 e RP475), cosi' l'AAL si
calcola per integrazione numerica invece che per ipotesi.

CORREZIONE (v2; diagnosi completa in validazione_assunzioni.py e §8 di
docs/pricing_coerenza.md). La v1 somava al trapezoid sui punti noti una "coda
frequente" con forma chiusa (k/(k-2)) * MDR(ag_30)/30. Quella forma chiusa e'
l'integrale della coda RARA da RP30 a infinito SENZA cap di MDR (non la coda
frequente), quindi ricontava la banda RP30-RP475 gia' integrata dal trapezoid:
il totale "2,30 mld" della v1 e' ritirato. La coda frequente RP<30, inoltre, non
e' stimabile dai tre punti disponibili: con qualunque legge di potenza k >= 2
l'integrale diverge (pendenza osservata del tratto RP30-RP72: 3,70 mediana).
Decomposizione corretta:

1. TRAPEZOID (RP30-RP475): interpolazione lineare di MDR in lambda sui punti
   noti; nessuna assunzione oltre ai dati.
2. CODA RARA (RP475-inf): lambda(ag) = lambda_475 * (ag/ag_475)^(-k), k = 3
   (dichiarata; sensibilita' k = 2,5 / 3,5), MDR con cap 1:
   coda = k/(k-2) * m_475/475 - (2/(k-2)) * lambda_cap,
   con lambda_cap = (1/475) * (ag_475 * sqrt(6))^k, il punto oltre il quale
   l'estrapolazione satura MDR (verificata contro quadratura in §8.1).
3. CODA FREQUENTE (RP<30): NON stimabile dai tre punti (divergenza k >= 2):
   l'AAL numerico e' un LIMITE INFERIORE della sola parte RP>=30, dichiarato.

CURVE_eff: per comune, (trapezoid + coda rara) * 475 / MDR(ag_475) — il
moltiplicatore empirico totale identificabile, contro il 4 assunto del benchmark
(che numericamente vale la sola banda: trapezoid * 475/m_475, mediana 4,53).

Curva EP nazionale (sismico puro, rate NON calibrate): loss attesa a scenario per
ogni RP noto, con banda epistemica 16/84 per centile sul RP475, e la loss
dell'evento 1-in-475 espressa in anni di AAL (e' il numero che dimensiona il
rischio di coda: quanto vale un evento raro in costi annuali).

Output: results/ep_curve.json (deterministico, solo stdlib).
Esecuzione dalla root del repo dopo pricing_model.py:
    python3 scripts/pricing/ep_curve.py
"""
import csv, json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))

MATRICE = os.environ.get('MATRICE', os.path.join(ROOT, 'data', 'Matrice_Modello_Savelli_Final_sismico.csv'))
OUT_JSON = os.environ.get('OUT_JSON', os.path.join(ROOT, 'results', 'ep_curve.json'))

# parametri identici a pricing_model.py (illustrativi, vedi docs/pricing_coerenza.md §2)
T_DESIGN = 475.0
MDR_COEFF = 6.0
CURVE_FACTOR = 4.0
K_RARE = 3.0             # pendenza della coda rara lambda(ag) ~ ag^(-k) (dichiarata)
K_SENS = (2.5, 3.5)


def mdr(ag):
    return min(1.0, MDR_COEFF * ag * ag)


def aal_comune(ag30, ag72, ag475, k=K_RARE):
    """AAL rate sismica: trapezoid sui 3 punti noti + coda rara oltre RP475 (cap 1)."""
    trap, coda = 0.0, 0.0
    if ag30 > 0:
        l30, l72, l475 = 1.0 / 30.0, 1.0 / 72.0, 1.0 / 475.0
        m30, m72, m475 = mdr(ag30), mdr(ag72), mdr(ag475)
        trap = (l30 - l72) * (m30 + m72) / 2.0 + (l72 - l475) * (m72 + m475) / 2.0
    if ag475 > 0 and k > 2.0:
        l475 = 1.0 / T_DESIGN
        m475 = mdr(ag475)
        # lambda oltre il quale l'estrapolazione ag(lambda) supera ag_cap = 1/sqrt(6)
        l_cap = l475 * (ag475 * math.sqrt(MDR_COEFF)) ** k
        coda = k / (k - 2.0) * m475 * l475 - 2.0 / (k - 2.0) * l_cap
    return trap, coda, trap + coda


def main():
    com = []
    with open(MATRICE, newline='', encoding='utf-8') as f:
        for r in csv.DictReader(f):
            com.append({'asset': float(r['asset_tot_EUR']),
                        'ag30': float(r['ag_RP30']), 'ag72': float(r['ag_RP72']),
                        'ag475': float(r['ag_RP475']),
                        'a16': float(r['ag_RP475_16perc']), 'a84': float(r['ag_RP475_84perc'])})
    assert len(com) == 3823, len(com)

    # ---------- AAL numerico (RP>=30): trapezoid + coda rara ----------
    tot_trap = tot_k3 = tot_bmk = 0.0
    tot_sens = {k: 0.0 for k in K_SENS}
    curve_eff = []          # (trapezoid + coda rara k3) * 475 / m475
    curve_eff_trap = []     # solo banda: confronto col 4 assunto
    for c in com:
        trap, coda, tot = aal_comune(c['ag30'], c['ag72'], c['ag475'])
        m475 = mdr(c['ag475'])
        bmk = CURVE_FACTOR * m475 / T_DESIGN
        tot_trap += trap * c['asset']
        tot_k3 += tot * c['asset']
        tot_bmk += bmk * c['asset']
        for k in K_SENS:
            _, _, t_k = aal_comune(c['ag30'], c['ag72'], c['ag475'], k=k)
            tot_sens[k] += t_k * c['asset']
        if m475 > 0:                      # CURVE_eff definito solo con hazard al design point
            curve_eff.append(tot * T_DESIGN / m475)
            curve_eff_trap.append(trap * T_DESIGN / m475)
    curve_eff.sort()
    curve_eff_trap.sort()
    q = lambda p: curve_eff[min(len(curve_eff) - 1, int(p * len(curve_eff)))]
    qt = lambda p: curve_eff_trap[min(len(curve_eff_trap) - 1, int(p * len(curve_eff_trap)))]
    tot_coda = tot_k3 - tot_trap

    # ---------- curva EP nazionale (sismico puro, non calibrato) ----------
    def loss_rp(get_ag):
        return sum(mdr(get_ag(c)) * c['asset'] for c in com)

    ep = {}
    for rp, col in ((30, 'ag30'), (72, 'ag72'), (475, 'ag475')):
        ep['RP%d' % rp] = round(loss_rp(lambda c, col=col: c[col]), 0)
    ep['RP475_band'] = {'p16': round(loss_rp(lambda c: c['a16']), 0),
                        'p84': round(loss_rp(lambda c: c['a84']), 0)}
    anni_aal_475 = ep['RP475'] / tot_k3

    out = {
        'modello': 'AAL sismico numerico dalla curva MPS04 (trapezoid RP30-475 + coda rara RP475-inf) e curva EP nazionale',
        'fonti': 'ag RP30/72/475 (+ percentili 16/84 al RP475) e asset dalla matrice estesa; '
                 'benchmark EAL da results/eal_comuni.csv (pricing_model.py)',
        'correzione_v2': 'la v1 somava al trapezoid una "coda frequente" (k/(k-2))*MDR(ag30)/30: forma chiusa '
                         'che e\' in realta\' la coda RARA da RP30 a infinito senza cap di MDR, e quindi ricontava '
                         'la banda RP30-475 gia\' integrata dal trapezoid (totale v1 "2,30 mld" ritirato). '
                         'Diagnosi numerica in results/validazione_assunzioni.json e §8 di docs/pricing_coerenza.md',
        'assunzioni': {
            'MDR(ag) = min(1, 6*ag^2)': 'identico a pricing_model.py §2',
            'trapezoid sui punti noti': 'interpolazione lineare di MDR in lambda tra RP30/72/475; '
                                        'nessuna assunzione aggiuntiva',
            'coda rara RP475-inf': 'lambda(ag) = lambda_475*(ag/ag_475)^(-k), k=3 (dichiarata; sensibilita\' '
                                   'k=2,5/3,5); MDR con cap 1: chiusa = k/(k-2)*m_475/475 - (2/(k-2))*lambda_cap',
            'coda frequente RP30': 'NON modellata: non stimabile dai tre punti (l\'integrale diverge per k>=2); '
                                  'l\'AAL numerico e\' un limite inferiore della parte RP>=30',
        },
        'aal_numerico_sismico': {
            'trapezoid_RP30_475_EUR': round(tot_trap, 0),
            'coda_rara_RP475_inf_k3_EUR': round(tot_coda, 0),
            'totale_k3_EUR': round(tot_k3, 0),
            'nota_limite': 'limite inferiore dichiarato: manca la coda frequente RP<30, non stimabile '
                           'dai tre punti della matrice (divergenza per k>=2)',
            'sensibilita_coda_rara': {'k2.5_EUR': round(tot_sens[2.5], 0), 'k3.5_EUR': round(tot_sens[3.5], 0)},
            'benchmark_CURVE4_EUR': round(tot_bmk, 0),
            'ratio_numerico_vs_benchmark': round(tot_k3 / tot_bmk, 3),
        },
        'curve_effettiva': {
            'definizione': '(trapezoid + coda rara k3) * 475 / MDR(ag_RP475): moltiplicatore empirico '
                           'totale identificabile (RP>=30) rispetto al solo design point',
            'mediana': round(q(0.5), 3), 'p10': round(q(0.1), 3), 'p90': round(q(0.9), 3),
            'banda_RP30_475_solo_trapezoid': {'mediana': round(qt(0.5), 3),
                                              'p10': round(qt(0.1), 3), 'p90': round(qt(0.9), 3),
                                              'nota': 'la banda RP30-475 da sola: il 4 assunto del benchmark '
                                                      'e\' numericamente vicino al moltiplicatore di banda, '
                                                      'non al totale identificabile'},
            'n_comuni_con_hazard': len(curve_eff),
        },
        'ep_nazionale_sismica': {
            **{k: v for k, v in ep.items() if k != 'RP475_band'},
            'RP475_banda_epistemica_16_84': ep['RP475_band'],
            'evento_1in475_in_anni_di_AAL': round(anni_aal_475, 1),
            'nota': 'loss fisica a scenario su MDR(ag)*asset, rate NON calibrate, solo hazard '
                    'sismico (frana e idraulico sono gia AAL nel benchmark, non hanno RP di scenario)',
        },
    }
    with open(OUT_JSON, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
        f.write('\n')
    print('AAL numerico sismico (RP>=30, limite inferiore): %.0f M EUR (trapezoid %.0f + coda rara %.0f) '
          'vs benchmark CURVE=4: %.0f M EUR (ratio %.3f)'
          % (tot_k3 / 1e6, tot_trap / 1e6, tot_coda / 1e6, tot_bmk / 1e6, tot_k3 / tot_bmk))
    print('sensibilita coda rara: k2.5 %.0f M | k3.5 %.0f M' % (tot_sens[2.5] / 1e6, tot_sens[3.5] / 1e6))
    print('CURVE_eff totale (banda+coda rara): mediana %.2f (p10 %.2f, p90 %.2f) | banda sola: %.2f'
          % (q(0.5), q(0.1), q(0.9), qt(0.5)))
    print('EP nazionale: RP30 %.0f M | RP72 %.0f M | RP475 %.0f M (banda 16/84: %.0f-%.0f)'
          % (ep['RP30'] / 1e6, ep['RP72'] / 1e6, ep['RP475'] / 1e6,
             ep['RP475_band']['p16'] / 1e6, ep['RP475_band']['p84'] / 1e6))
    print('evento 1-in-475 = %.1f anni di AAL numerico' % anni_aal_475)
    print('OK: %s' % os.path.relpath(OUT_JSON, ROOT))


if __name__ == '__main__':
    main()
