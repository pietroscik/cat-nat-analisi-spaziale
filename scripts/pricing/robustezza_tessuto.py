#!/usr/bin/env python3
"""Robustezza del risultato del tessuto produttivo: il beta dell'ISP sopravvive?

Il risultato centrale di §5 (docs/pricing_coerenza.md) — a parita' di hazard,
dimensione e struttura dimensionale, +1 deviazione standard di performance riduce
il peso della loss del ~31% (beta_ISP = -0,375, t_rob = -21,8, n = 3.791) — viene
sottoposto a cinque verifiche, senza toccare i dati e senza selezioni ex post:

1. WINSORIZZAZIONE: log-peso tagliato all'1%/99% (la coda del peso-EBITDA arriva al
   70,7%: il beta dipende dalla coda o e' un gradiente diffuso?);
2. TRIMMING: stesso campione senza l'1% di comuni a peso piu' alto;
3. MISURA ALTERNATIVA: l'ISP (indicatore composito della tesi) sostituito dal ROA
   comunale EBITDA/asset — una misura di performance che non dipende dalla pipeline
   della tesi: se il gradiente e' della performance e non dell'ISP in se', il segno
   e l'ordine di grandezza devono reggere;
4. SLX (spatially lagged X): la specifica base con i ritardi spaziali dei regressori
   (W*ISP, W*log1p(ag), W*share_frana, W*share_idraulico; W KNN k=5): l'OLS resta
   consistente, e il W*ISP dice se il peso della loss di un comune dipende anche
   dalla performance dei comuni vicini (effetto di contesto);
5. MORAN SUI RESIDUI della specifica base: diagnostica spaziale — se i residui sono
   autocorrelati, l'OLS trascura una componente spaziale e il SLX e' la risposta
   dichiarata (non si stima un SAR: si documenta la diagnostica).

Output: results/robustezza_tessuto.json (deterministico, solo stdlib).
Esecuzione dalla root del repo dopo pricing_model.py e tessuto_produttivo.py:
    python3 scripts/pricing/robustezza_tessuto.py
"""
import csv, json, math, os

from pricing_model import ols
from spazializzazione import knn_rows, moran_global

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))

MATRICE = os.environ.get('MATRICE', os.path.join(ROOT, 'data', 'Matrice_Modello_Savelli_Final_sismico.csv'))
EAL_CSV = os.environ.get('EAL_CSV', os.path.join(ROOT, 'results', 'eal_comuni.csv'))
OUT_JSON = os.environ.get('OUT_JSON', os.path.join(ROOT, 'results', 'robustezza_tessuto.json'))


def main():
    tess = {}
    with open(MATRICE, newline='', encoding='utf-8') as f:
        for r in csv.DictReader(f):
            tess[r['PRO_COM']] = {
                'EBITDA': float(r['EBITDA_migl_EUR']) * 1e3,
                'asset': float(r['asset_tot_EUR']),
                'quota_grandi': (float(r['asset_grandi_EUR']) / float(r['asset_tot_EUR'])
                                 if float(r['asset_tot_EUR']) > 0 else 0.0),
                'ISP': None if r['ISP_std_medio'].strip() in ('', 'nan', 'NaN')
                       else float(r['ISP_std_medio'])}
    com = []
    with open(EAL_CSV, newline='', encoding='utf-8') as f:
        for r in csv.DictReader(f):
            t = tess.get(r['PRO_COM'])
            if t is None or t['ISP'] is None or t['EBITDA'] <= 0:
                continue
            eal = float(r['eal_calibrato_EUR'])
            if eal <= 0:
                continue
            c = dict(pc=r['PRO_COM'], **t,
                     eal=eal, ag=float(r['ag_RP475']),
                     frana=float(r['hazard_frana_share']), idr=float(r['hazard_idraulico_share']),
                     lat=float(r['lat']), lon=float(r['long']))
            c['peso'] = eal / t['EBITDA']
            c['roa'] = t['EBITDA'] / t['asset'] if t['asset'] > 0 else None
            com.append(c)
    n = len(com)
    assert n == 3791, n

    # ---------- specifiche ----------
    # X base (identica a tessuto_produttivo.py §3: interc, ISP, log1p(ag), frana, idr, log asset, quota Grandi)
    def x_base(c, perf='ISP'):
        p = c['ISP'] if perf == 'ISP' else c['roa']
        return [1.0, p, math.log1p(c['ag']), c['frana'], c['idr'],
                math.log(c['asset']), c['quota_grandi']]

    y = [math.log(c['peso']) for c in com]

    # winsorizzazione del log-peso all'1%/99% (valori, non osservazioni: gradiente diffuso vs coda)
    ys = sorted(y)
    lo, hi = ys[int(0.01 * n)], ys[int(0.99 * n)]
    y_win = [min(max(v, lo), hi) for v in y]

    # trimming: stesso campione senza l'1% di comuni a peso piu' alto
    peso_ord = sorted(c['peso'] for c in com)
    soglia = peso_ord[int(0.99 * n)]
    trim_idx = [i for i, c in enumerate(com) if c['peso'] <= soglia]

    def run(nome, X, yy, idx, perf='ISP', note=''):
        names = ['interc', 'performance', 'log1p(ag)', 'share frana', 'share idraulico',
                 'log(asset totale)', 'quota Grandi']
        reg, r2, nn = ols([X[i] for i in idx] if idx is not None else X,
                          [yy[i] for i in idx] if idx is not None else yy, names)
        b = reg[1]
        # effetto per deviazione standard della misura di performance (confrontabilita' ISP vs ROA)
        vals = [X[i][1] for i in idx] if idx is not None else [row[1] for row in X]
        m = sum(vals) / len(vals)
        sd = math.sqrt(sum((v - m) ** 2 for v in vals) / len(vals))
        return {'specifica': nome, 'n': nn, 'R2': round(r2, 4),
                'misura': 'ISP (indicatore della tesi)' if perf == 'ISP' else 'ROA comunale = EBITDA/asset',
                'beta_performance': round(b['beta'], 4), 't_rob': round(b['t_rob'], 2),
                'beta_x_deviazione_standard': round(b['beta'] * sd, 4),
                'sd_misura': round(sd, 4),
                'note': note}

    specs = [run('base (tessuto §5.3)', [x_base(c) for c in com], y, None,
                 note='replica di riferimento: beta_ISP = -0,375, t_rob = -21,8'),
             run('winsorizzato log-peso 1%/99%', [x_base(c) for c in com], y_win, None,
                 note='la coda del peso (max 70,7%) non porta il gradiente'),
             run('trim top 1% peso', [x_base(c) for c in com], y, trim_idx,
                 note='senza i comuni col peso piu\' alto; stesso campione del winsorizzato'),
             run('performance = ROA comunale', [x_base(c, perf='roa') for c in com], y, None, perf='roa',
                 note='misura indipendente dalla pipeline della tesi: regge il segno e l\'ordine di grandezza?')]

    # ---------- SLX: ritardi spaziali dei regressori (W KNN k=5) ----------
    cs = math.cos(math.radians(sum(c['lat'] for c in com) / n))
    neigh = knn_rows([(c['lon'] * cs, c['lat']) for c in com])
    wisp, w_ag, w_f, w_i = [], [], [], []
    for i, (idx, w) in enumerate(neigh):
        wisp.append(sum(com[j]['ISP'] for j in idx) * w)
        w_ag.append(sum(math.log1p(com[j]['ag']) for j in idx) * w)
        w_f.append(sum(com[j]['frana'] for j in idx) * w)
        w_i.append(sum(com[j]['idr'] for j in idx) * w)
    X_slx = [x_base(c) + [wisp[i], w_ag[i], w_f[i], w_i[i]] for i, c in enumerate(com)]
    names_slx = ['interc', 'ISP (performance)', 'log1p(ag)', 'share frana', 'share idraulico',
                 'log(asset totale)', 'quota Grandi', 'W*ISP', 'W*log1p(ag)', 'W*share frana', 'W*share idraulico']
    reg_slx, r2_slx, n_slx = ols(X_slx, y, names_slx)
    slx = {'specifica': 'SLX (regressori con ritardo spaziale, W KNN k=5)', 'n': n_slx,
           'R2': round(r2_slx, 4),
           'coeff': [{'name': b['name'], 'beta': round(b['beta'], 4), 't_rob': round(b['t_rob'], 2)}
                     for b in reg_slx]}

    # ---------- Moran sui residui della specifica base ----------
    reg_base, _, _ = ols([x_base(c) for c in com], y)
    resid = [y[i] - sum(x_base(com[i])[a] * reg_base[a]['beta'] for a in range(7)) for i in range(n)]
    moran_res = moran_global(resid, neigh)

    # correlazione ISP-ROA (contesto per la variante 3)
    def pearson(v1, v2):
        m1, m2 = sum(v1) / len(v1), sum(v2) / len(v2)
        num = sum((a - m1) * (b - m2) for a, b in zip(v1, v2))
        den = math.sqrt(sum((a - m1) ** 2 for a in v1) * sum((b - m2) ** 2 for b in v2))
        return num / den
    corr_isp_roa = pearson([c['ISP'] for c in com], [c['roa'] for c in com])

    out = {
        'modello': 'verifiche di robustezza del risultato del tessuto produttivo (beta della performance sul log peso-EBITDA)',
        'fonti': 'results/eal_comuni.csv (pricing_model.py) + matrice estesa; specifica base identica a tessuto_produttivo.py',
        'metodo': 'OLS con SE White HC1; W KNN k=5 sui centroidi (come spazializzazione.py); nessuna selezione ex post',
        'specifiche': specs,
        'slx': slx,
        'moran_residui_base': {**moran_res, 'nota': 'diagnostica spaziale dell\'OLS: residui autocorrelati => '
                                                     'il SLX e\' la risposta dichiarata, il beta_ISP resta la lettura'},
        'correlazione_ISP_ROA': round(corr_isp_roa, 4),
    }
    with open(OUT_JSON, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
        f.write('\n')
    for s in specs:
        print('%-32s n=%d  R2=%.3f  beta_perf=%+.4f (xSD %+.3f)  t_rob=%+.1f  (%s)'
              % (s['specifica'], s['n'], s['R2'], s['beta_performance'], s['beta_x_deviazione_standard'],
                 s['t_rob'], s['misura']))
    print('SLX: W*ISP beta %+.4f t_rob %+.1f | ISP beta %+.4f t_rob %+.1f (R2 %.3f)'
          % (reg_slx[7]['beta'], reg_slx[7]['t_rob'], reg_slx[1]['beta'], reg_slx[1]['t_rob'], r2_slx))
    print('Moran residui base: I=%.4f z=%+.1f p_perm=%.4f | corr ISP-ROA %.3f'
          % (moran_res['I'], moran_res['z_clifford'], moran_res['p_perm_one_sided'], corr_isp_roa))
    print('OK: %s' % os.path.relpath(OUT_JSON, ROOT))


if __name__ == '__main__':
    main()
