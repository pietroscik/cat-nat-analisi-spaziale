#!/usr/bin/env python3
"""Modello di pricing (benchmark di loss attesa) e coerenza asset — Cat-Nat.

Tre analisi (vedi docs/pricing_coerenza.md per la discussione completa):

1. BENCHMARK DI EAL (expected annual loss) a tre hazard per comune:
   - sismico:  EAL_rate(ag) = CURVE * min(1, MDR(ag)) / T,  MDR(ag) = min(1, m_sism*ag^2)
     (T = 475 anni e' il RP del design point; CURVE approssima il contributo
     delle intensita' inferiori che una singola scorciatoia di curva sottostima);
   - frana:    EAL_rate = f_frana * share P3/P4  (frequenza-danno annua su area esposta);
   - idraulico: EAL_rate = f_idr * share P3.
   Parametri ILLUSTRATIVI (ordini di grandezza della letteratura tecnica, dichiarati e
   perturbati nell'analisi di sensibilita'): il benchmark NON e' un modello di pricing
   operativo, misura la coerenza RELATIVA della tariffazione rispetto al rischio modellato.

2. CALIBRAZIONE E LOSS RATIO: un solo scalare c allinea la rate media (pesata per asset)
   del benchmark a quella IVASS (22,8 per 10.000 EUR); loss_ratio = rate_IVASS / rate_bmk
   misura l'adeguatezza relativa (dove la tariffa segue il rischio e dove no).

3. COERENZA ASSET: la tariffa provinciale e' la variabile di rischio "pura" (il premio
   comunale = rate * asset e' meccanico in asset). Regressione OLS a livello provincia
   (n=107) di log(rate) su log1p degli hazard pesati per asset, con SE classici e robusti
   (White HC1); decomposizione della variabilita' del log-premio comunale tra componente
   asset e componente rate.

Output (dalla root del repo, nessuna dipendenza oltre la stdlib):
- results/eal_comuni.csv         ← EAL e loss ratio per comune (3.823 righe)
- results/pricing_benchmark.json ← parametri, calibrazione, tabella province, regressione
- docs/mappa_loss.svg            ← mappa di loss: EAL (EUR) e loss ratio per comune
"""
import csv, json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))

MATRICE = os.environ.get('MATRICE', os.path.join(ROOT, 'data', 'Matrice_Modello_Savelli_Final_sismico.csv'))
OUT_CSV = os.environ.get('OUT_CSV', os.path.join(ROOT, 'results', 'eal_comuni.csv'))
OUT_JSON = os.environ.get('OUT_JSON', os.path.join(ROOT, 'results', 'pricing_benchmark.json'))
OUT_SVG = os.environ.get('OUT_SVG', os.path.join(ROOT, 'docs', 'mappa_loss.svg'))

# ---- parametri del benchmark EAL (illustrativi, vedi docs/pricing_coerenza.md §2) ----
T_DESIGN = 475.0      # RP del design point sismico (10% in 50 anni)
MDR_COEFF = 6.0       # MDR(ag) = min(1, 6*ag^2): ag=0.106 (media naz. pesata) -> 6.7% di danno
                      # al design point; ag=0.25 (Calabria) -> 37%
CURVE_FACTOR = 4.0    # contributo delle intensita' sotto il design point (integrazione curva):
                      # con MDR ~ ag^2 il punto RP475 da solo sottostima l'AAL di ~4 volte
F_FRANA = 0.004       # 0.4%/anno di damage-rate atteso su area P3/P4 esposta
F_IDR = 0.006         # 0.6%/anno di damage-rate atteso su area P3 esposta
                      # -> mix nazionale del benchmark: sismico ~55%, idraulico ~32%, frana ~13%

RATIO_BINS = [0.0, 0.6, 0.8, 1.0, 1.2, 1.5, 99.0]   # classi di loss ratio (tariffa/benchmark)


def eal_rates(ag, share_frana, share_idr):
    """Rate di EAL per 10.000 EUR di asset, per hazard e totale (non calibrata)."""
    r_sism = CURVE_FACTOR * min(1.0, MDR_COEFF * ag * ag) / T_DESIGN * 1e4
    r_fran = F_FRANA * share_frana * 1e4
    r_idr = F_IDR * share_idr * 1e4
    return r_sism, r_fran, r_idr, r_sism + r_fran + r_idr


# ---------- OLS con SE classici e White (HC1) ----------
def ols(X, y):
    n, k = len(y), len(X[0])
    xtx = [[sum(X[i][a] * X[i][b] for i in range(n)) for b in range(k)] for a in range(k)]
    xty = [sum(X[i][a] * y[i] for i in range(n)) for a in range(k)]
    A = [row[:] + [xty[a]] for a, row in enumerate(xtx)]
    for c in range(k):                      # eliminazione gaussiana con pivoting
        p = max(range(c, k), key=lambda r: abs(A[r][c]))
        A[c], A[p] = A[p], A[c]
        for r in range(k):
            if r != c and abs(A[c][c]) > 1e-12:
                f = A[r][c] / A[c][c]
                for cc in range(c, k + 1):
                    A[r][cc] -= f * A[c][cc]
    beta = [A[a][k] / A[a][a] for a in range(k)]
    resid = [y[i] - sum(X[i][a] * beta[a] for a in range(k)) for i in range(n)]
    ssr = sum(e * e for e in resid)
    sst = sum((v - sum(y) / n) ** 2 for v in y)
    # (X'X)^-1 (serve per la covarianza)
    def inv(M):
        m = len(M)
        A = [M[r][:] + [1.0 if r == c else 0.0 for c in range(m)] for r in range(m)]
        for c in range(m):
            p = max(range(c, m), key=lambda r: abs(A[r][c]))
            A[c], A[p] = A[p], A[c]
            for r in range(m):
                if r != c:
                    f = A[r][c] / A[c][c]
                    for cc in range(c, 2 * m):
                        A[r][cc] -= f * A[c][cc]
        return [[A[r][m + c] / A[r][r] for c in range(m)] for r in range(m)]
    xtx_inv = inv(xtx)
    s2 = ssr / (n - k)
    cov_c = [[xtx_inv[a][b] * s2 for b in range(k)] for a in range(k)]
    meat = [[sum(X[i][a] * X[i][b] * resid[i] ** 2 for i in range(n)) for b in range(k)] for a in range(k)]
    hc = n / (n - k)
    cov_w = [[hc * sum(xtx_inv[a][t] * meat[t][u] * xtx_inv[u][b] for t in range(k) for u in range(k))
              for b in range(k)] for a in range(k)]
    names = ['interc', 'log1p(ag aw)', 'log1p(frana aw)', 'log1p(idraulico aw)']
    out = []
    for a in range(k):
        se_c, se_w = math.sqrt(cov_c[a][a]), math.sqrt(cov_w[a][a])
        out.append({'name': names[a], 'beta': beta[a],
                    'se': se_c, 't': beta[a] / se_c,
                    'se_rob': se_w, 't_rob': beta[a] / se_w})
    return out, 1 - ssr / sst, n


def spearman(v1, v2):
    """Spearman = Pearson sui ranghi (concorde sui pareggi)."""
    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for t in range(i, j + 1):
                r[order[t]] = avg
            i = j + 1
        return r
    r1, r2 = ranks(v1), ranks(v2)
    n = len(v1)
    m1, m2 = sum(r1) / n, sum(r2) / n
    num = sum((a - m1) * (b - m2) for a, b in zip(r1, r2))
    den = math.sqrt(sum((a - m1) ** 2 for a in r1) * sum((b - m2) ** 2 for b in r2))
    return num / den


def main():
    rows = []
    with open(MATRICE, newline='', encoding='utf-8') as f:
        for r in csv.DictReader(f):
            rows.append(r)
    n = len(rows)

    # ---------- 1. EAL per comune ----------
    com = []
    for r in rows:
        asset = float(r['asset_tot_EUR'])
        ag = float(r['ag_RP475'])
        fr = float(r['hazard_frana_share'])
        idr = float(r['hazard_idraulico_share'])
        rate_ivass = float(r['premio_10k_prov'])
        rs, rf, ri, rtot = eal_rates(ag, fr, idr)
        com.append({'PRO_COM': r['PRO_COM'], 'COMUNE': r['COMUNE'], 'Provincia': r['Provincia'],
                    'COD_PROV': r['COD_PROV'], 'COD_REG': r['COD_REG'],
                    'lat': float(r['lat']), 'lon': float(r['long']),
                    'asset': asset, 'ag': ag, 'frana': fr, 'idr': idr,
                    'rate_ivass': rate_ivass,
                    'r_sism': rs, 'r_fran': rf, 'r_idr': ri, 'r_bmk_raw': rtot,
                    'eal_raw_eur': asset * rtot / 1e4,
                    'sism_nc': int(r['sismico_non_classificato'])})
    tot_asset = sum(c['asset'] for c in com)

    # ---------- 2. calibrazione sull'aggregato IVASS ----------
    aw_ivass = sum(c['asset'] * c['rate_ivass'] for c in com) / tot_asset
    aw_raw = sum(c['asset'] * c['r_bmk_raw'] for c in com) / tot_asset
    calib = aw_ivass / aw_raw
    for c in com:
        c['r_bmk'] = c['r_bmk_raw'] * calib
        c['eal_eur'] = c['asset'] * c['r_bmk'] / 1e4
        c['loss_ratio'] = c['rate_ivass'] / c['r_bmk'] if c['r_bmk'] > 0 else None
    eal_raw_tot = sum(c['eal_raw_eur'] for c in com)
    prem_tot = sum(c['asset'] * c['rate_ivass'] for c in com) / 1e4

    # ---------- 3. coerenza asset: decomposizione del log-premio ----------
    # log(premio) = log(asset) + log(rate) - log(1e4): identita' esatta
    lg_prem = [math.log(c['asset'] * c['rate_ivass'] / 1e4) for c in com]
    lg_ast = [math.log(c['asset']) for c in com]
    lg_rate = [math.log(c['rate_ivass']) for c in com]
    m = lambda v: sum(v) / len(v)
    var = lambda v, mu=None: sum((x - (mu or m(v))) ** 2 for x in v) / len(v)
    mu = m(lg_prem)
    v_prem = var(lg_prem, mu)
    v_ast, v_rate = var(lg_ast), var(lg_rate)
    ma, mr = m(lg_ast), m(lg_rate)
    cov_ar = sum((a - ma) * (b - mr) for a, b in zip(lg_ast, lg_rate)) / n
    # termini mancanti: log1p(asset*rate/1e4) ~ log(asset)+log(rate)-log(1e4) per asset*rate >> 1e4
    decomp = {'var_log_premio': v_prem, 'var_log_asset': v_ast, 'var_log_rate': v_rate,
              '2cov_asset_rate': 2 * cov_ar,
              'quota_asset': (v_ast + cov_ar) / v_prem, 'quota_rate': (v_rate + cov_ar) / v_prem,
              'cor_log_asset_log_rate': cov_ar / math.sqrt(v_ast * v_rate)}

    # ---------- 4. regressione provinciale (n=107) ----------
    prov = {}
    for c in com:
        p = prov.setdefault(c['COD_PROV'], {'name': c['Provincia'], 'com': []})
        p['com'].append(c)
    prov_rows = []
    for cp, p in sorted(prov.items(), key=lambda kv: kv[1]['name'].lower()):
        cc = p['com']
        ta = sum(c['asset'] for c in cc)
        ag_w = sum(c['ag'] * c['asset'] for c in cc) / ta
        fr_w = sum(c['frana'] * c['asset'] for c in cc) / ta
        idr_w = sum(c['idr'] * c['asset'] for c in cc) / ta
        rate = cc[0]['rate_ivass']
        bmk = sum(c['r_bmk'] * c['asset'] for c in cc) / ta
        prov_rows.append({'provincia': p['name'], 'n': len(cc),
                          'rate_ivass': round(rate, 2), 'bmk_calibrato': round(bmk, 2),
                          'loss_ratio': round(rate / bmk, 3) if bmk > 0 else None,
                          'eal_raw_M EUR': round(sum(c['eal_raw_eur'] for c in cc) / 1e6, 1)})
    Xp = [[1.0, math.log1p(sum(c['ag'] * c['asset'] for c in p['com']) / sum(c['asset'] for c in p['com'])),
           math.log1p(sum(c['frana'] * c['asset'] for c in p['com']) / sum(c['asset'] for c in p['com'])),
           math.log1p(sum(c['idr'] * c['asset'] for c in p['com']) / sum(c['asset'] for c in p['com']))]
          for p in prov.values()]
    yp = [math.log(p['com'][0]['rate_ivass']) for p in prov.values()]
    reg, r2, nreg = ols(Xp, yp)

    # ---------- 5. sensibilita' dei parametri ----------
    base_ratio = [c['loss_ratio'] for c in com]
    base_eal = [c['eal_raw_eur'] for c in com]
    # (varianti definite esplicitamente qui sotto)
    def ratios_with(mdr_c=MDR_COEFF, curve=CURVE_FACTOR, f_f=F_FRANA, f_i=F_IDR):
        out = []
        for c in com:
            rs = curve * min(1.0, mdr_c * c['ag'] ** 2) / T_DESIGN * 1e4
            rt = rs + f_f * c['frana'] * 1e4 + f_i * c['idr'] * 1e4
            out.append(rt)
        aw = sum(c['asset'] * v for c, v in zip(com, out)) / tot_asset
        return [c['rate_ivass'] / (v * calib) if v > 0 else None
                for c, v in zip(com, out)], sum(c['asset'] * v for c, v in zip(com, out)) / 1e4 * 0 + aw
    sens = []
    for lbl, kw in [('MDR sismico x0.5', {'mdr_c': MDR_COEFF / 2}), ('MDR sismico x2', {'mdr_c': MDR_COEFF * 2}),
                    ('fattore curva x0.5', {'curve': CURVE_FACTOR / 2}), ('fattore curva x2', {'curve': CURVE_FACTOR * 2}),
                    ('frequenza frana x0.5', {'f_f': F_FRANA / 2}), ('frequenza frana x2', {'f_f': F_FRANA * 2}),
                    ('frequenza idraulico x0.5', {'f_i': F_IDR / 2}), ('frequenza idraulico x2', {'f_i': F_IDR * 2})]:
        rr, _ = ratios_with(**kw)
        pairs = [(a, b) for a, b in zip(rr, base_ratio) if a is not None and b is not None]
        sens.append({'variante': lbl,
                     'spearman_loss_ratio': round(spearman([a for a, _ in pairs], [b for _, b in pairs]), 4)})
    # stabilita' del ranking: perturba tutti i parametri insieme
    rr_lo, _ = ratios_with(mdr_c=MDR_COEFF / 2, curve=CURVE_FACTOR / 2, f_f=F_FRANA / 2, f_i=F_IDR / 2)
    pairs_all = [(a, b) for a, b in zip(rr_lo, base_ratio) if a is not None and b is not None]
    sens.append({'variante': 'tutti i parametri x0.5', 'spearman_loss_ratio': round(spearman([a for a, _ in pairs_all], [b for _, b in pairs_all]), 4)})

    # ---------- 6. statistiche di sintesi ----------
    ratios = sorted(c['loss_ratio'] for c in com if c['loss_ratio'] is not None)
    eal_sorted = sorted(c['eal_raw_eur'] for c in com)
    q = lambda v, p: v[min(len(v) - 1, int(p * len(v)))]
    sard = [c for c in com if c['COD_REG'] == '20']

    out = {
        'modello': 'benchmark EAL 3-hazard (parametri illustrativi) + calibrazione aggregata IVASS',
        'parametri': {'T_design_anni': T_DESIGN, 'MDR_coeff': MDR_COEFF, 'curve_factor': CURVE_FACTOR,
                      'f_frana': F_FRANA, 'f_idraulico': F_IDR},
        'calibrazione': {'rate_ivass_media_pesata_asset': round(aw_ivass, 2),
                         'rate_bmk_raw_media_pesata_asset': round(aw_raw, 2),
                         'fattore_calibrazione': round(calib, 3)},
        'nazionale': {'n_comuni': n, 'n_province': len(prov),
                      'asset_tot_EUR': tot_asset,
                      'eal_raw_tot_EUR': eal_raw_tot,
                      'premio_teorico_tot_EUR': prem_tot,
                      'eal_raw_su_premio': round(eal_raw_tot / prem_tot, 3),
                      'contributi_eal_raw': {
                          'sismico_pct': round(sum(c['asset'] * c['r_sism'] for c in com) / 1e4 / eal_raw_tot * 100, 1),
                          'frana_pct': round(sum(c['asset'] * c['r_fran'] for c in com) / 1e4 / eal_raw_tot * 100, 1),
                          'idraulico_pct': round(sum(c['asset'] * c['r_idr'] for c in com) / 1e4 / eal_raw_tot * 100, 1)}},
        'loss_ratio': {'mediana': round(q(ratios, 0.5), 3), 'p10': round(q(ratios, 0.1), 3),
                       'p90': round(q(ratios, 0.9), 3), 'min': round(ratios[0], 3), 'max': round(ratios[-1], 2),
                       'distribuzione_classi': {
                           f"{RATIO_BINS[i]:g}-{(f'{RATIO_BINS[i+1]:g}' if RATIO_BINS[i+1] < 90 else 'inf')}":
                           sum(1 for v in ratios if RATIO_BINS[i] <= v < RATIO_BINS[i + 1])
                           for i in range(len(RATIO_BINS) - 1)},
                       'n_sovra_sottoprezzo': {'ratio_gt_1': sum(1 for v in ratios if v > 1),
                                               'ratio_lt_1': sum(1 for v in ratios if v <= 1)}},
        'decomposizione_log_premio': {k: round(v, 4) if isinstance(v, float) else v for k, v in decomp.items()},
        'regressione_provinciale': {'n': nreg, 'R2': round(r2, 4), 'coeff': reg},
        'province': prov_rows,
        'sensibilita': sens,
        'nota_sardegna': {'n_comuni': len(sard),
                          'loss_ratio_mediano': round(sorted(c['loss_ratio'] for c in sard
                                                             if c['loss_ratio'] is not None)[len(sard) // 2], 2) if sard else None,
                          'n_comuni_senza_rischio_modellato': sum(1 for c in com if c['loss_ratio'] is None),
                          'nota': 'benchmark ~ 0 per i 3 hazard modellati: il loss ratio non e\' interpretabile (la tariffa copre un perimetro piu\' ampio)'}
    }
    with open(OUT_JSON, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=1, ensure_ascii=False)

    with open(OUT_CSV, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['PRO_COM', 'COMUNE', 'Provincia', 'COD_PROV', 'lat', 'long',
                    'asset_tot_EUR', 'ag_RP475', 'hazard_frana_share', 'hazard_idraulico_share',
                    'rate_eal_sismico', 'rate_eal_frana', 'rate_eal_idraulico',
                    'rate_bmk_raw', 'rate_bmk_calibrato', 'eal_raw_EUR', 'eal_calibrato_EUR',
                    'rate_ivass', 'loss_ratio'])
        for c in com:
            w.writerow([c['PRO_COM'], c['COMUNE'], c['Provincia'], c['COD_PROV'],
                        c['lat'], c['lon'], c['asset'], c['ag'], c['frana'], c['idr'],
                        round(c['r_sism'], 3), round(c['r_fran'], 3), round(c['r_idr'], 3),
                        round(c['r_bmk_raw'], 3), round(c['r_bmk'], 3),
                        round(c['eal_raw_eur'], 2), round(c['eal_eur'], 2),
                        c['rate_ivass'],
                        round(c['loss_ratio'], 4) if c['loss_ratio'] is not None else ''])

    # ---------- 7. mappa SVG ----------
    draw_svg(com)

    print(f"comuni: {n} | province: {len(prov)}")
    print(f"rate IVASS media pesata: {aw_ivass:.2f} | benchmark raw: {aw_raw:.2f} | calibrazione c={calib:.3f}")
    print(f"EAL raw nazionale: {eal_raw_tot/1e9:.2f} Mld EUR/anno ({eal_raw_tot/prem_tot*100:.1f}% del premio teorico)")
    print(f"contributi EAL raw: {out['nazionale']['contributi_eal_raw']}")
    print(f"loss ratio: mediana {out['loss_ratio']['mediana']}, p10 {out['loss_ratio']['p10']}, p90 {out['loss_ratio']['p90']}")
    print(f"decomposizione log-premio: quota asset {decomp['quota_asset']*100:.1f}%, quota rate {decomp['quota_rate']*100:.1f}%, corr(asset,rate) {decomp['cor_log_asset_log_rate']:.3f}")
    print(f"regressione provinciale R2 = {r2:.3f}")
    for b in reg:
        print(f"  {b['name']:22s} beta {b['beta']:+.4f}  t {b['t']:+6.2f}  t_robust {b['t_rob']:+6.2f}")
    print("OK: pricing_benchmark.json, eal_comuni.csv, mappa_loss.svg")


# ---------- rendering SVG ----------
def draw_svg(com):
    lons = [c['lon'] for c in com]
    lats = [c['lat'] for c in com]
    lon0, lon1, lat0, lat1 = min(lons), max(lons), min(lats), max(lats)

    W_PANEL, H_PANEL = 440, 480
    PAD_X, PAD_TOP = 30, 58
    GAP = 60
    W, H = PAD_X * 2 + W_PANEL * 2 + GAP, H_PANEL + 128
    # proiezione: compensa la contrazione dei meridiani (cos(lat media))
    LAT_M = math.radians((lat0 + lat1) / 2)
    SCALE = min(W_PANEL / ((lon1 - lon0) * math.cos(LAT_M)), H_PANEL / (lat1 - lat0))

    def xy(lon, lat, x0):
        return (x0 + (lon - lon0) * math.cos(LAT_M) * SCALE,
                PAD_TOP + (lat1 - lat) * SCALE)

    # --- pannello A: EAL raw (EUR), 8 classi a quantili, palette calda ---
    eals = sorted(c['eal_raw_eur'] for c in com)
    cuts_e = [eals[int(i / 8 * len(eals))] for i in range(1, 8)]
    PAL_E = ['#f7f0c7', '#f3d99a', '#eec168', '#e8a33c', '#e07b28', '#d4521c', '#bd2b18', '#991b12']
    PAL_R = ['#2b5f9e', '#5b8fc4', '#a9c4df', '#e8e8e3', '#f0c184', '#e08a3c', '#c94b2a', '#9c1f1f']
    bins_r = [0.0, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0, 99.0]

    def fmt_eur(v):
        return f"{v/1e6:.1f} M€" if v >= 1e6 else (f"{v/1e3:.0f} k€" if v >= 1e3 else f"{v:.0f} €")

    s = ['<?xml version="1.0" encoding="UTF-8"?>',
         f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" font-family="Helvetica,Arial,sans-serif">',
         f'<rect width="{W}" height="{H}" fill="#ffffff"/>',
         '<text x="%d" y="26" font-size="17" font-weight="bold" fill="#222">Mappa di loss — EAL attesa per comune e coerenza della tariffazione Cat-Nat</text>' % (W // 2),
         '<text x="%d" y="44" font-size="11" fill="#666">Benchmark EAL a 3 hazard (parametri illustrativi, calibrato sull\u2019aggregato IVASS) · loss ratio = tariffa IVASS / benchmark calibrato · 3.823 comuni</text>' % (W // 2)]

    def panel(x0, title, subtitle, classify, palette, legend_items):
        s.append(f'<text x="{x0 + W_PANEL // 2}" y="{PAD_TOP - 26}" font-size="13" font-weight="bold" fill="#333" text-anchor="middle">{title}</text>')
        s.append(f'<text x="{x0 + W_PANEL // 2}" y="{PAD_TOP - 12}" font-size="10" fill="#777" text-anchor="middle">{subtitle}</text>')
        s.append(f'<rect x="{x0}" y="{PAD_TOP}" width="{W_PANEL}" height="{H_PANEL}" fill="#f4f2ec" stroke="#ccc" stroke-width="1"/>')
        for c in com:
            px, py = xy(c['lon'], c['lat'], x0)
            k = classify(c)
            s.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="2.2" fill="{palette[k]}" fill-opacity="0.92" stroke="none"/>')
        # legenda
        ly = PAD_TOP + H_PANEL + 26
        s.append(f'<text x="{x0}" y="{ly - 8}" font-size="10" fill="#555">classi:</text>')
        lx = x0
        for i, (col, lab) in enumerate(legend_items):
            if lx > x0 + W_PANEL - 88:
                lx = x0
                ly += 20
            s.append(f'<rect x="{lx}" y="{ly}" width="12" height="12" fill="{col}" stroke="#999" stroke-width="0.5"/>')
            s.append(f'<text x="{lx + 16}" y="{ly + 10}" font-size="10" fill="#444">{lab}</text>')
            lx += 16 + 6 * len(lab) + 14

    def legend_eal():
        items = []
        bounds = [0.0] + cuts_e + [eals[-1]]
        for i in range(8):
            items.append((PAL_E[i], f"{fmt_eur(bounds[i])}–{fmt_eur(bounds[i+1])}" if i < 7 else f"≥ {fmt_eur(bounds[7])}"))
        return items

    def legend_ratio():
        items = []
        for i in range(8):
            lo, hi = bins_r[i], bins_r[i + 1]
            items.append((PAL_R[i], (f"{lo:g}–{hi:g}" if i < 7 else f"> {bins_r[7]:g}")))
        return items

    def cls_eal(c):
        v = c['eal_raw_eur']
        for i, cut in enumerate(cuts_e):
            if v < cut:
                return i
        return 7

    def cls_ratio(c):
        v = c['loss_ratio']
        if v is None:
            return 7
        for i in range(8):
            if v < bins_r[i + 1]:
                return i
        return 7

    panel(PAD_X, 'Loss attesa (EAL benchmark, EUR/anno)',
          'più scuro = loss attesa maggiore (pesata per gli asset esposti)',
          cls_eal, PAL_E, legend_eal())
    x1 = PAD_X + W_PANEL + GAP
    panel(x1, 'Loss ratio (coerenza tariffa/rischio)',
          'rosso = tariffa sopra il benchmark, blu = sotto',
          cls_ratio, PAL_R, legend_ratio())

    s.append(f'<text x="{W // 2}" y="{H - 12}" font-size="9.5" fill="#888" text-anchor="middle">Fonti: IVASS (tariffe), INGV MPS04 (ag), ISPRA (PAI P3/P4, P3), AIDA (asset) · generato da scripts/pricing/pricing_model.py · dettagli in docs/pricing_coerenza.md</text>')
    s.append('</svg>')
    with open(OUT_SVG, 'w', encoding='utf-8') as f:
        f.write('\n'.join(s))


if __name__ == '__main__':
    main()
