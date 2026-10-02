#!/usr/bin/env python3
"""Esposizione delle imprese performanti e peso della loss sul tessuto produttivo.

Indagine attuariale che collega il benchmark EAL (results/eal_comuni.csv, generato da
pricing_model.py) agli indicatori del tessuto produttivo comunale (matrice estesa):
performance delle imprese (ISP_std_medio, l'indicatore composito della tesi),
EBITDA, dipendenti, struttura dimensionale degli asset.

Tre domande (vedi docs/pricing_coerenza.md §5):

1. PESO DELLA LOSS SUL TESSUTO: quanta parte del valore prodotto comunale assorbe
   la loss annuale attesa (EAL calibrata / EBITDA, EAL per addetto, EAL per impresa)?
2. ESPOSIZIONE DELLE IMPRESE PERFORMANTI: i comuni con imprese performanti
   (quartili superiori dell'ISP, cluster High-High) espongono piu' o meno asset
   al rischio rispetto al resto? Confronto delle rate benchmark e del peso della
   loss per classe di performance, e quota di EAL vs quota di EBITDA.
3. MODELLO PREDITTIVO: regressione OLS (SE robusti White) del log peso-EBITDA su
   performance (ISP), hazard (ag, frana, idraulico), dimensione e struttura
   (log asset totale, quota Grandi): la performance predice il peso della loss
   dopo aver controllato il rischio fisico e la struttura del tessuto?

Output: results/esposizione_tessuto.json (deterministico, solo stdlib).
Esecuzione dalla root del repo dopo pricing_model.py:
    python3 scripts/pricing/tessuto_produttivo.py
"""
import csv, json, math, os
from pricing_model import ols, spearman

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))

MATRICE = os.environ.get('MATRICE', os.path.join(ROOT, 'data', 'Matrice_Modello_Savelli_Final_sismico.csv'))
EAL_CSV = os.environ.get('EAL_CSV', os.path.join(ROOT, 'results', 'eal_comuni.csv'))
OUT_JSON = os.environ.get('OUT_JSON', os.path.join(ROOT, 'results', 'esposizione_tessuto.json'))


def pearson(v1, v2):
    n = len(v1)
    m1, m2 = sum(v1) / n, sum(v2) / n
    num = sum((a - m1) * (b - m2) for a, b in zip(v1, v2))
    den = math.sqrt(sum((a - m1) ** 2 for a in v1) * sum((b - m2) ** 2 for b in v2))
    return num / den


def mean(v):
    return sum(v) / len(v) if v else None


def median(v):
    s = sorted(v)
    return s[len(s) // 2] if s else None


def main():
    # join matrice (tessuto produttivo) + eal_comuni (benchmark EAL)
    tess = {}
    with open(MATRICE, newline='', encoding='utf-8') as f:
        for r in csv.DictReader(f):
            tess[r['PRO_COM']] = {
                'EBITDA': float(r['EBITDA_migl_EUR']) * 1e3,          # EUR/anno
                'dipendenti': float(r['dipendenti']),
                'n_imprese': int(r['n_imprese']),
                'asset_tot': float(r['asset_tot_EUR']),
                'quota_grandi': (float(r['asset_grandi_EUR']) / float(r['asset_tot_EUR'])
                                 if float(r['asset_tot_EUR']) > 0 else 0.0),
                'ISP': None if r['ISP_std_medio'].strip() in ('', 'nan', 'NaN')
                       else float(r['ISP_std_medio']),
                'LISA': r['LISA_cluster_moda'],
                'Provincia': r['Provincia'], 'COMUNE': r['COMUNE'],
                'COD_PROV': r['COD_PROV']}
    eal = {}
    with open(EAL_CSV, newline='', encoding='utf-8') as f:
        for r in csv.DictReader(f):
            eal[r['PRO_COM']] = {'eal_cal': float(r['eal_calibrato_EUR']),
                                  'eal_raw': float(r['eal_raw_EUR']),
                                  'rate_bmk': float(r['rate_bmk_calibrato']),
                                  'loss_ratio': float(r['loss_ratio']) if r['loss_ratio'] else None,
                                  'ag': float(r['ag_RP475']),
                                  'frana': float(r['hazard_frana_share']),
                                  'idr': float(r['hazard_idraulico_share'])}

    com = []
    for pc, t in tess.items():
        e = eal.get(pc)
        if e is None:
            continue
        c = dict(pc=pc, **t, **e)
        c['peso_ebitda'] = c['eal_cal'] / c['EBITDA'] if c['EBITDA'] > 0 else None  # frazione annua
        c['eal_addetto'] = c['eal_cal'] / c['dipendenti'] if c['dipendenti'] > 0 else None
        c['eal_impresa'] = c['eal_cal'] / c['n_imprese'] if c['n_imprese'] > 0 else None
        com.append(c)
    n = len(com)
    assert n == 3823, n

    tot_eal = sum(c['eal_cal'] for c in com)
    tot_ebitda = sum(c['EBITDA'] for c in com)
    tot_dip = sum(c['dipendenti'] for c in com)
    tot_imp = sum(c['n_imprese'] for c in com)

    # ---------- 1. peso della loss sul tessuto ----------
    pesi = sorted(c['peso_ebitda'] for c in com if c['peso_ebitda'] is not None)
    q = lambda v, p: v[min(len(v) - 1, int(p * len(v)))]
    nazionali = {
        'eal_tot_EUR': tot_eal, 'ebitda_tot_EUR': tot_ebitda,
        'peso_ebitda_nazionale': tot_eal / tot_ebitda,
        'eal_per_addetto_nazionale': tot_eal / tot_dip,
        'eal_per_impresa_nazionale': tot_eal / tot_imp,
        'peso_ebitda_comuni': {'mediana': q(pesi, 0.5), 'p90': q(pesi, 0.9), 'max': pesi[-1],
                               'n_comuni_ebitda_positivo': len(pesi)}}

    # ---------- 2. esposizione per classe di performance (ISP) ----------
    validi = sorted((c for c in com if c['ISP'] is not None), key=lambda c: c['ISP'])
    nv = len(validi)
    q25, q50, q75 = validi[nv // 4]['ISP'], validi[nv // 2]['ISP'], validi[3 * nv // 4]['ISP']
    def classe(isp):
        if isp < q25: return 'Q1 (meno performanti)'
        if isp < q50: return 'Q2'
        if isp < q75: return 'Q3'
        return 'Q4 (piu\' performanti)'
    quad = {}
    for c in validi:
        quad.setdefault(classe(c['ISP']), []).append(c)
    per_quartile = []
    for k, cc in [('Q1 (meno performanti)', 'Q1'), ('Q2', 'Q2'), ('Q3', 'Q3'), ('Q4 (piu\' performanti)', 'Q4')]:
        v = quad[k]
        a_v = sum(x['asset_tot'] for x in v)
        eal_v = sum(x['eal_cal'] for x in v)
        eb_v = sum(x['EBITDA'] for x in v if x['EBITDA'] > 0)
        rate_pesata = eal_v / a_v * 1e4
        per_quartile.append({
            'classe': k, 'n_comuni': len(v),
            'isp_medio': round(mean([x['ISP'] for x in v]), 3),
            'rate_bmk_pesata_asset': round(rate_pesata, 2),
            'peso_ebitda': round(eal_v / eb_v, 4),
            'quota_EAL_nazionale': round(eal_v / tot_eal, 4),
            'quota_EBITDA_nazionale': round(eb_v / tot_ebitda, 4),
            'intensita_esposizione (quotaEAL/quotaEBITDA)': round((eal_v / tot_eal) / (eb_v / tot_ebitda), 3),
            'ag_medio': round(mean([x['ag'] for x in v]), 4),
            'share_frana_medio': round(mean([x['frana'] for x in v]), 4),
            'share_idr_medio': round(mean([x['idr'] for x in v]), 4)})

    # cluster High-High dell'ISP vs resto
    hh = [c for c in validi if c['LISA'] == 'High-High']
    nohh = [c for c in validi if c['LISA'] != 'High-High']
    cluster = {
        'n_HighHigh': len(hh),
        'isp_medio_HH': round(mean([c['ISP'] for c in hh]), 3),
        'rate_bmk_pesata_HH': round(sum(c['eal_cal'] for c in hh) /
                                    sum(c['asset_tot'] for c in hh) * 1e4, 2),
        'rate_bmk_pesata_resto': round(sum(c['eal_cal'] for c in nohh) /
                                       sum(c['asset_tot'] for c in nohh) * 1e4, 2),
        'peso_ebitda_HH': round(sum(c['eal_cal'] for c in hh) /
                                sum(c['EBITDA'] for c in hh if c['EBITDA'] > 0), 4)}

    # correlazioni ISP
    isps = [c['ISP'] for c in validi]
    correlazioni = {
        'ISP_vs_ag': round(pearson(isps, [c['ag'] for c in validi]), 4),
        'ISP_vs_share_frana': round(pearson(isps, [c['frana'] for c in validi]), 4),
        'ISP_vs_share_idr': round(pearson(isps, [c['idr'] for c in validi]), 4),
        'ISP_vs_rate_bmk': round(pearson(isps, [c['rate_bmk'] for c in validi]), 4),
        'ISP_vs_peso_ebitda_spearman': round(
            spearman([c['ISP'] for c in validi if c['peso_ebitda'] is not None],
                     [c['peso_ebitda'] for c in validi if c['peso_ebitda'] is not None]), 4)}

    # ---------- 3. modello predittivo del peso della loss ----------
    # nota: 3 comuni (Sardegna) hanno EAL calibrata = 0 (hazard nullo su tutti e tre i pericoli
    # modellati) e sono esclusi dal log-peso insieme ai 24 senza ISP e ai 5 con EBITDA <= 0
    reg_sample = [c for c in validi if c['peso_ebitda'] is not None and c['peso_ebitda'] > 0 and c['EBITDA'] > 0]
    X = [[1.0, c['ISP'], math.log1p(c['ag']), c['frana'], c['idr'],
          math.log(c['asset_tot']), c['quota_grandi']] for c in reg_sample]
    y = [math.log(c['peso_ebitda']) for c in reg_sample]
    reg, r2, nreg = ols(X, y, ['interc', 'ISP (performance)', 'log1p(ag)',
                               'share frana', 'share idraulico', 'log(asset totale)',
                               'quota Grandi'])

    # ---------- 4. aggregato provinciale ----------
    prov = {}
    for c in com:
        p = prov.setdefault(c['COD_PROV'], {'name': c['Provincia'], 'eal': 0.0, 'eb': 0.0, 'dip': 0.0})
        p['eal'] += c['eal_cal']
        if c['EBITDA'] > 0:
            p['eb'] += c['EBITDA']
        p['dip'] += c['dipendenti']
    prov_rows = sorted(({'provincia': p['name'], 'peso_ebitda': p['eal'] / p['eb'],
                         'eal_addetto': p['eal'] / p['dip']}
                        for p in prov.values()),
                       key=lambda r: -r['peso_ebitda'])
    for r in prov_rows:
        r['peso_ebitda'] = round(r['peso_ebitda'], 4)
        r['eal_addetto'] = round(r['eal_addetto'], 1)

    out = {
        'modello': 'peso della loss EAL sul tessuto produttivo ed esposizione delle imprese performanti',
        'fonti': 'EAL calibrata da results/eal_comuni.csv (pricing_model.py); ISP_std_medio, EBITDA, dipendenti dalla matrice estesa',
        'nazionale': {k: (round(v, 4) if isinstance(v, float) else v) for k, v in nazionali.items()},
        'nazionale_note': {
            'peso_ebitda_nazionale': 'EAL calibrata / EBITDA: quota del valore prodotto assorbita dalla loss attesa annua',
            'n_comuni': n},
        'quartili_ISP': per_quartile,
        'cluster_HighHigh': cluster,
        'correlazioni': correlazioni,
        'regressione_peso_ebitda': {
            'n': nreg, 'R2': round(r2, 4),
            'specifica': 'log(EAL/EBITDA) ~ ISP + log1p(ag) + share_frana + share_idraulico + log(asset_tot) + quota_grandi (OLS, SE White HC1)',
            'coeff': reg},
        'province_peso_ebitda': prov_rows[:15] + [{'...': 'tabella completa ordinata nel campo sopra; prime 15 per peso'}]
        if False else prov_rows,
    }
    with open(OUT_JSON, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=1, ensure_ascii=False)

    # ---------- sintesi console ----------
    print(f"comuni: {n} | EAL calibrata {tot_eal/1e9:.2f} Mld | EBITDA {tot_ebitda/1e9:.1f} Mld")
    print(f"peso nazionale della loss: {tot_eal/tot_ebitda*100:.2f}% dell'EBITDA | "
          f"{tot_eal/tot_dip:.0f} EUR per addetto | {tot_eal/tot_imp:.0f} EUR per impresa")
    for q4 in per_quartile:
        print(f"  {q4['classe']:24s} n={q4['n_comuni']:4d} rate bmk {q4['rate_bmk_pesata_asset']:5.2f} "
              f"peso EBITDA {q4['peso_ebitda']*100:5.2f}%  intensita' {q4['intensita_esposizione (quotaEAL/quotaEBITDA)']:.2f}")
    print(f"cluster High-High (n={cluster['n_HighHigh']}): rate bmk pesata {cluster['rate_bmk_pesata_HH']} "
          f"vs resto {cluster['rate_bmk_pesata_resto']} | peso EBITDA HH {cluster['peso_ebitda_HH']*100:.2f}%")
    print("correlazioni ISP:", correlazioni)
    print(f"regressione log-peso-EBITDA: n={nreg}, R2={r2:.4f}")
    for b in reg:
        print(f"  {b['name']:22s} beta {b['beta']:+.4f}  t {b['t']:+6.2f}  t_robust {b['t_rob']:+6.2f}")
    top3 = prov_rows[:3]
    print("province top per peso EBITDA: " + ", ".join(f"{r['provincia']} {r['peso_ebitda']*100:.1f}%" for r in top3))
    print("OK: esposizione_tessuto.json")


if __name__ == '__main__':
    main()
