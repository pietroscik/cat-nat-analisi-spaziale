#!/usr/bin/env python3
"""Curva di eccedenza sismica (EP) e AAL numerico dalla curva di hazard MPS04.

Il benchmark EAL (pricing_model.py) integra la curva di hazard con un moltiplicatore
costante (CURVE = 4): il punto RP475 da solo, con MDR ~ ag^2, sottostima l'AAL di ~4x.
Qui il fattore viene stimato empiricamente: la matrice estesa porta tre punti della
curva di pericolosita' sismica per comune (ag a RP30, RP72 e RP475), cosi' l'AAL si
calcola per integrazione numerica (trapezoid sui punti noti) invece che per ipotesi.

Tre contributi:

1. AAL NUMERICO SISMICO (3 punti): trapezoid su (lambda, MDR) ai tre punti noti
   RP30/72/475; nessuna assunzione oltre ai punti dati.
2. CODA FREQUENTE (lambda > 1/30): gli eventi piu' frequenti di RP30 non sono in
   matrice; contributo chiuso sotto assunzione dichiarata lambda(ag) = lambda_30 *
   (ag/ag_30)^(-k), k = 3 (pendenza tipica MPS04; sensibilita' k = 2,5 / 3,5).
   Con MDR = 6*ag^2 e k > 2: coda = (k/(k-2)) * MDR(ag_30)/30.
3. CURVA EFFETTIVA: per comune, CURVE_eff = AAL_numerico * 475 / MDR(ag_475) —
   quanto vale il "vero" moltiplicatore rispetto all'ipotesi CURVE = 4 del
   benchmark, e quanto ci si allontana (la geografia del loss ratio e' gia'
   robusta a CURVE x0,5/x2, §7.1 di docs/pricing_coerenza.md).

Curva EP nazionale (sismico puro, rate NON calibrate): loss attesa a scenario per
ogni RP noto, con banda epistemica 16/84 per centile sul RP475, e la loss del
evento 1-in-475 espressa in anni di AAL (e' il numero che dimensiona il rischio di
coda: quanto vale un evento raro in costi annuali).

Output: results/ep_curve.json (deterministico, solo stdlib).
Esecuzione dalla root del repo dopo pricing_model.py:
    python3 scripts/pricing/ep_curve.py
"""
import csv, json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))

MATRICE = os.environ.get('MATRICE', os.path.join(ROOT, 'data', 'Matrice_Modello_Savelli_Final_sismico.csv'))
EAL_CSV = os.environ.get('EAL_CSV', os.path.join(ROOT, 'results', 'eal_comuni.csv'))
OUT_JSON = os.environ.get('OUT_JSON', os.path.join(ROOT, 'results', 'ep_curve.json'))

# parametri identici a pricing_model.py (illustrativi, vedi docs/pricing_coerenza.md §2)
T_DESIGN = 475.0
MDR_COEFF = 6.0
CURVE_FACTOR = 4.0
K_TAIL = 3.0            # pendenza della coda frequente lambda(ag) ~ ag^(-k) (dichiarata)
K_SENS = (2.5, 3.5)


def mdr(ag):
    return min(1.0, MDR_COEFF * ag * ag)


def aal_comune(ag30, ag72, ag475, k=K_TAIL):
    """AAL rate (frazione di asset/anno) sismica: trapezoid sui 3 punti + coda frequente."""
    l30, l72, l475 = 1.0 / 30.0, 1.0 / 72.0, 1.0 / 475.0
    m30, m72, m475 = mdr(ag30), mdr(ag72), mdr(ag475)
    trap = (l30 - l72) * (m30 + m72) / 2.0 + (l72 - l475) * (m72 + m475) / 2.0
    coda = (k / (k - 2.0)) * m30 / 30.0 if ag30 > 0 else 0.0
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

    # ---------- AAL numerico vs benchmark ----------
    tot_trap = tot_k3 = tot_bmk = 0.0
    tot_sens = {k: 0.0 for k in K_SENS}
    curve_eff = []
    curve_eff_trap = []
    for c in com:
        trap, coda, tot = aal_comune(c['ag30'], c['ag72'], c['ag475'])
        m475 = mdr(c['ag475'])
        bmk = CURVE_FACTOR * m475 / T_DESIGN
        tot_trap += trap * c['asset']
        tot_k3 += tot * c['asset']
        tot_bmk += bmk * c['asset']
        for k in K_SENS:
            _, c_k, t_k = aal_comune(c['ag30'], c['ag72'], c['ag475'], k=k)
            tot_sens[k] += t_k * c['asset']
        if m475 > 0:                      # CURVE_eff definito solo con hazard al design point
            curve_eff.append(tot * T_DESIGN / m475)
            curve_eff_trap.append(trap * T_DESIGN / m475)
    curve_eff.sort()
    curve_eff_trap.sort()
    q = lambda p: curve_eff[min(len(curve_eff) - 1, int(p * len(curve_eff)))]
    qt = lambda p: curve_eff_trap[min(len(curve_eff_trap) - 1, int(p * len(curve_eff_trap)))]

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
        'modello': 'AAL sismico numerico dalla curva di hazard MPS04 (3 punti) e curva EP nazionale',
        'fonti': 'ag RP30/72/475 (+ percentili 16/84 al RP475) e asset dalla matrice estesa; '
                 'benchmark EAL da results/eal_comuni.csv (pricing_model.py)',
        'assunzioni': {
            'MDR(ag) = min(1, 6*ag^2)': 'identico a pricing_model.py §2',
            'trapezoid sui punti noti': 'interpolazione lineare di MDR in lambda tra RP30/72/475; '
                                        'nessuna assunzione aggiuntiva',
            'coda frequente lambda>1/30': 'lambda(ag) = lambda_30*(ag/ag_30)^(-k), k=3 (dichiarata); '
                                          'contributo chiuso = (k/(k-2))*MDR(ag_30)/30',
        },
        'aal_numerico_sismico': {
            'trapezoid_3punti_EUR': round(tot_trap, 0),
            'coda_frequente_k3_EUR': round(tot_k3 - tot_trap, 0),
            'totale_k3_EUR': round(tot_k3, 0),
            'sensibilita_coda': {'k2.5_EUR': round(tot_sens[2.5], 0), 'k3.5_EUR': round(tot_sens[3.5], 0)},
            'benchmark_CURVE4_EUR': round(tot_bmk, 0),
            'ratio_numerico_vs_benchmark': round(tot_k3 / tot_bmk, 3),
        },
        'curve_effettiva': {
            'definizione': 'AAL_numerico * 475 / MDR(ag_RP475): moltiplicatore empirico della curva '
                           'di hazard rispetto al solo design point',
            'mediana': round(q(0.5), 3), 'p10': round(q(0.1), 3), 'p90': round(q(0.9), 3),
            'trapezoid_solo_banda_RP30_475': {'mediana': round(qt(0.5), 3),
                                             'p10': round(qt(0.1), 3), 'p90': round(qt(0.9), 3),
                                             'nota': 'senza coda frequente: il confronto diretto '
                                                     'con l\'ipotesi CURVE=4 del benchmark'},
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
    print('AAL numerico sismico: %.0f M EUR (trapezoid %.0f + coda %.0f) vs benchmark CURVE=4: %.0f M EUR (ratio %.3f)'
          % (tot_k3 / 1e6, tot_trap / 1e6, (tot_k3 - tot_trap) / 1e6, tot_bmk / 1e6, tot_k3 / tot_bmk))
    print('CURVE_eff mediana: %.2f (p10 %.2f, p90 %.2f) | solo banda RP30-475: %.2f'
          % (q(0.5), q(0.1), q(0.9), qt(0.5)))
    print('EP nazionale: RP30 %.0f M | RP72 %.0f M | RP475 %.0f M (banda 16/84: %.0f-%.0f)'
          % (ep['RP30'] / 1e6, ep['RP72'] / 1e6, ep['RP475'] / 1e6,
             ep['RP475_band']['p16'] / 1e6, ep['RP475_band']['p84'] / 1e6))
    print('evento 1-in-475 = %.1f anni di AAL numerico' % anni_aal_475)
    print('OK: %s' % os.path.relpath(OUT_JSON, ROOT))


if __name__ == '__main__':
    main()
