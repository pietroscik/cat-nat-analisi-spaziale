#!/usr/bin/env python3
"""Chi paga il premio Cat-Nat: incidenza per classe dimensionale (PMI vs Grandi).

La rate IVASS e' costante entro provincia e si applica per 10.000 EUR di asset esposto,
quindi la ripartizione del premio tra PMI e Grandi imprese di un comune segue
esattamente la ripartizione degli asset. La domanda non e' tanto "quanto" (banale:
come gli asset) quanto l'INCIDENZA: quanto pesa il premio e la loss attesa per
impresa, e sul valore prodotto comunale, quando cambia la struttura dimensionale
del tessuto (quota degli asset in mani Grandi).

Domande (vedi docs/pricing_coerenza.md §5.4):

1. QUOTA NAZIONALE: che parte del premio teorico (3,26 mld EUR) e della EAL
   calibrata (3,26 mld) pagano le PMI (micro+piccola+media) e le Grandi?
2. PREMIO PER IMPRESA: quanto paga in media un'impresa PMI vs una Grande, e
   quanto vale la loss attesa che grava su ciascuna?
3. INCIDENZA SUL VALORE PRODOTTO: premio/EBITDA e premio/ricavi a livello comunale
   (EBITDA e ricavi non sono splittabili per classe in matrice: limite dichiarato).
4. GRADIENTE STRUTTURALE: quartili della quota_grandi (asset in mani Grandi):
   incidenza mediana e premio per impresa per quartile — il Cat-Nat pesa di piu'
   dove il tessuto e' dominato dalle Grandi? (collega al concentration risk di §5.3).

Output: results/chi_paga.json (deterministico, solo stdlib).
Esecuzione dalla root del repo dopo pricing_model.py:
    python3 scripts/pricing/chi_paga.py
"""
import csv, json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))

MATRICE = os.environ.get('MATRICE', os.path.join(ROOT, 'data', 'Matrice_Modello_Savelli_Final_sismico.csv'))
EAL_CSV = os.environ.get('EAL_CSV', os.path.join(ROOT, 'results', 'eal_comuni.csv'))
OUT_JSON = os.environ.get('OUT_JSON', os.path.join(ROOT, 'results', 'chi_paga.json'))


def median(v):
    s = sorted(v)
    return s[len(s) // 2] if s else None


def main():
    tess = {}
    with open(MATRICE, newline='', encoding='utf-8') as f:
        for r in csv.DictReader(f):
            tess[r['PRO_COM']] = {
                'asset_pmi': float(r['asset_PMI_EUR']), 'asset_grandi': float(r['asset_grandi_EUR']),
                'asset_tot': float(r['asset_tot_EUR']),
                'n_pmi': int(r['n_micro']) + int(r['n_piccola']) + int(r['n_media']),
                'n_grandi': int(r['n_grande']),
                'EBITDA': float(r['EBITDA_migl_EUR']) * 1e3,
                'ricavi': float(r['ricavi_migl_EUR']) * 1e3,
                'COMUNE': r['COMUNE'], 'Provincia': r['Provincia']}
    eal = {}
    with open(EAL_CSV, newline='', encoding='utf-8') as f:
        for r in csv.DictReader(f):
            eal[r['PRO_COM']] = {'rate': float(r['rate_ivass']),
                                 'rate_bmk': float(r['rate_bmk_calibrato']),
                                 'eal_cal': float(r['eal_calibrato_EUR'])}

    com = []
    for pc, t in tess.items():
        e = eal.get(pc)
        if e is None:
            continue
        c = dict(pc=pc, **t, **e)
        c['quota_grandi'] = c['asset_grandi'] / c['asset_tot'] if c['asset_tot'] > 0 else 0.0
        c['premio_pmi'] = c['rate'] * c['asset_pmi'] / 1e4
        c['premio_grandi'] = c['rate'] * c['asset_grandi'] / 1e4
        c['premio_tot'] = c['rate'] * c['asset_tot'] / 1e4
        # loss attesa calibrata per classe: rate benchmark pesata per gli asset di classe
        c['eal_pmi'] = c['rate_bmk'] * c['asset_pmi'] / 1e4
        c['eal_grandi'] = c['rate_bmk'] * c['asset_grandi'] / 1e4
        com.append(c)
    assert len(com) == 3823, len(com)

    tot = lambda k: sum(c[k] for c in com)
    tot_premio = tot('premio_tot')
    tot_eal = tot('eal_cal')
    n_pmi, n_grandi = tot('n_pmi'), tot('n_grandi')
    prem_pmi, prem_grandi = tot('premio_pmi'), tot('premio_grandi')
    eal_pmi, eal_grandi = tot('eal_pmi'), tot('eal_grandi')

    nazionale = {
        'premio_teorico_EUR': round(tot_premio, 0),
        'quota_premio_PMI': round(prem_pmi / tot_premio, 4),
        'quota_premio_grandi': round(prem_grandi / tot_premio, 4),
        'premio_medio_impresa_PMI_EUR': round(prem_pmi / n_pmi, 0),
        'premio_medio_impresa_grande_EUR': round(prem_grandi / n_grandi, 0),
        'eal_calibrata_per_classe': {'PMI_EUR': round(eal_pmi, 0), 'grandi_EUR': round(eal_grandi, 0),
                                     'quota_eal_PMI': round(eal_pmi / tot_eal, 4)},
        'n_imprese': {'PMI': n_pmi, 'grandi': n_grandi},
        'incidenza_comunale': {
            'nota': 'EBITDA e ricavi sono totali di comune: l\'incidenza non e\' splittabile per classe',
            'premio_su_EBITDA_mediana': round(median([c['premio_tot'] / c['EBITDA'] for c in com if c['EBITDA'] > 0]), 4),
            'premio_su_ricavi_mediana': round(median([c['premio_tot'] / c['ricavi'] for c in com if c['ricavi'] > 0]), 4),
        },
    }

    # ---------- struttura dimensionale del tessuto ----------
    # 2.564 comuni su 3.823 (67%) non hanno Grandi imprese (quota_grandi = 0): i quartili
    # della quota sarebbero degeneri (Q1 e Q2 interamente a zero, fatto informativo in se').
    # Gruppo A: comuni senza Grandi (premio interamente PMI); Gruppo B: comuni con Grandi,
    # splittati in terzili di quota_grandi per rank (deterministico).
    com.sort(key=lambda c: c['quota_grandi'])
    senza = [c for c in com if c['quota_grandi'] == 0.0]
    con_g = [c for c in com if c['quota_grandi'] > 0.0]
    nb = len(con_g)
    terzili = [[c for i, c in enumerate(con_g) if i * 3 // nb == t] for t in range(3)]

    def stats(group):
        prem = sum(c['premio_tot'] for c in group)
        n_p = sum(c['n_pmi'] for c in group)
        n_g = sum(c['n_grandi'] for c in group)
        return {
            'n_comuni': len(group),
            'quota_grandi_asset_range': [round(min(c['quota_grandi'] for c in group), 4),
                                         round(max(c['quota_grandi'] for c in group), 4)],
            'quota_premio_nazionale': round(prem / tot_premio, 4),
            'premio_medio_PMI_EUR': round(sum(c['premio_pmi'] for c in group) / n_p, 0) if n_p else None,
            'premio_medio_grande_EUR': round(sum(c['premio_grandi'] for c in group) / n_g, 0) if n_g else None,
            'premio_su_EBITDA_mediano': round(median([c['premio_tot'] / c['EBITDA'] for c in group if c['EBITDA'] > 0]), 4),
            'premio_su_ricavi_mediano': round(median([c['premio_tot'] / c['ricavi'] for c in group if c['ricavi'] > 0]), 4),
        }

    out = {
        'modello': 'ripartizione del premio teorico Cat-Nat e incidenza per classe dimensionale (PMI vs Grandi)',
        'fonti': 'asset_PMI/grandi, n_imprese per classe, EBITDA e ricavi dalla matrice estesa; '
                 'rate IVASS e EAL calibrata da results/eal_comuni.csv (pricing_model.py)',
        'nazionale': nazionale,
        'gruppi_struttura_dimensionale': {
            'nota': '67% dei comuni non ha Grandi imprese: i quartili di quota_grandi sarebbero '
                    'degeneri; gruppi = senza Grandi vs terzili di quota_grandi tra i comuni che le hanno',
            'senza_grandi': stats(senza),
            'con_grandi_terzili_quota': [stats(g) for g in terzili],
        },
    }
    with open(OUT_JSON, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
        f.write('\n')
    print('premio teorico %.0f mld EUR: PMI %.1f%% (%.0f mld) | Grandi %.1f%% (%.0f mld)'
          % (tot_premio / 1e9, 100 * prem_pmi / tot_premio, prem_pmi / 1e9,
             100 * prem_grandi / tot_premio, prem_grandi / 1e9))
    print('premio medio per impresa: PMI %.0f EUR vs Grande %.0f EUR (x%.0f)'
          % (prem_pmi / n_pmi, prem_grandi / n_grandi, (prem_grandi / n_grandi) / (prem_pmi / n_pmi)))
    print('incidenza mediana: premio/EBITDA %.4f | premio/ricavi %.4f'
          % (nazionale['incidenza_comunale']['premio_su_EBITDA_mediana'],
             nazionale['incidenza_comunale']['premio_su_ricavi_mediana']))
    print('struttura: %d comuni senza Grandi | %d con Grandi (terzili quota %.0f-%.0f%%)'
          % (len(senza), nb,
             100 * min(c['quota_grandi'] for c in terzili[0]),
             100 * max(c['quota_grandi'] for c in terzili[2])))
    s = out['gruppi_struttura_dimensionale']
    print('  senza Grandi:        premio/EBITDA %.4f | premio medio PMI %.0f EUR'
          % (s['senza_grandi']['premio_su_EBITDA_mediano'], s['senza_grandi']['premio_medio_PMI_EUR']))
    for i, g in enumerate(s['con_grandi_terzili_quota']):
        print('  con Grandi T%d:      premio/EBITDA %.4f | premio medio PMI %.0f vs Grande %.0f EUR'
              % (i + 1, g['premio_su_EBITDA_mediano'], g['premio_medio_PMI_EUR'], g['premio_medio_grande_EUR']))
    print('OK: %s' % os.path.relpath(OUT_JSON, ROOT))


if __name__ == '__main__':
    main()
