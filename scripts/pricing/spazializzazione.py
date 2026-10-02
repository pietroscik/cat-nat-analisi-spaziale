#!/usr/bin/env python3
"""Spazializzazione del loss ratio: autocorrelazione e cluster LISA della coerenza tariffaria.

Il ponte tra i due mondi del repo: la geografia dell'adeguatezza tariffaria (loss
ratio = tariffa IVASS / benchmark calibrato, docs/pricing_coerenza.md §3) letta con
gli strumenti spaziali del SDM. Se la tariffa fosse coerente col rischio comune per
comune, i loss ratio non avrebbero struttura spaziale; se invece la coerenza e'
determinata a livello provinciale (la rate e' costante entro provincia), i comuni
vicini — anche di province diverse — condividono lo stesso segno di errore, e i
cluster devono emergere.

Metodo (deterministico, solo stdlib):

- W: KNN k = 5 sui centroidi comunali (come il SDM di livello comune), pesi di riga
  standardizzati ("W"), tie-break deterministico su (distanza, PRO_COM);
- Moran globale: statistiche analitiche (Cliff-Ord, assunzione di normalita') E
  p-value di permutazione (999 permutazioni, seed fisso 42) — se i due concordano,
  la conclusione e' robusta all'impostazione del test;
- Moran locale (LISA, Anselin 1995): randomizzazione condizionata (499 permutazioni,
  seed fisso 42), pseudo-p bilaterale; classificazione HH/LL/HL/LH sul segno di
  z_i e del ritardo spaziale Wz ai comuni significativi (p < 0,05).

Variabili spazializzate: loss ratio (n = 1.920), log peso-EBITDA (n = 3.791) e
rate benchmark (n = 3.823) — la W e' ricostruita sul sottoinsieme valido di ciascuna.

Output: results/spazializzazione_tariffa.json + docs/mappa_lisa_tariffa.svg
(deterministici, solo stdlib). Esecuzione dalla root del repo dopo pricing_model.py:
    python3 scripts/pricing/spazializzazione.py
"""
import csv, json, math, os, random

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))

MATRICE = os.environ.get('MATRICE', os.path.join(ROOT, 'data', 'Matrice_Modello_Savelli_Final_sismico.csv'))
EAL_CSV = os.environ.get('EAL_CSV', os.path.join(ROOT, 'results', 'eal_comuni.csv'))
OUT_JSON = os.environ.get('OUT_JSON', os.path.join(ROOT, 'results', 'spazializzazione_tariffa.json'))
OUT_SVG = os.environ.get('OUT_SVG', os.path.join(ROOT, 'docs', 'mappa_lisa_tariffa.svg'))

KNN = 5
N_PERM_GLOBAL = 999
N_PERM_LOCAL = 499
SEED = 42
ALPHA = 0.05


# ---------- matrice dei pesi KNN k=5 (righe standardizzate, tie-break deterministico) ----------
def knn_rows(pts):
    """pts: list of (x, y) proiettate; ritorna liste di indici e pesi di riga per punto."""
    n = len(pts)
    neigh = []
    for i in range(n):
        xi, yi = pts[i]
        d = []
        for j in range(n):
            if j == i:
                continue
            dx, dy = xi - pts[j][0], yi - pts[j][1]
            d.append((dx * dx + dy * dy, j))
        d.sort()
        idx = [j for _, j in d[:KNN]]
        w = 1.0 / len(idx)
        neigh.append((idx, w))
    return neigh


def moran_global(z, neigh):
    """I di Moran con z analitico Cliff-Ord (normalita') e p permutazione (seed fisso)."""
    n = len(z)
    m = sum(z) / n
    zz = [v - m for v in z]
    sst = sum(v * v for v in zz)
    S0 = sum(len(idx) * w for idx, w in neigh)          # = n con righe standardizzate
    num = sum(zz[i] * w * sum(zz[j] for j in idx) for i, (idx, w) in enumerate(neigh))
    I = (n / S0) * num / sst if sst > 0 else 0.0
    EI = -1.0 / (n - 1)
    # Cliff-Ord, assunzione di normalita': S1 e S2 richiedono i pesi in entrambe le direzioni
    wfull = {}
    for i, (idx, w) in enumerate(neigh):
        for j in idx:
            wfull[(i, j)] = w
    S1 = 0.5 * sum((w_ij + wfull.get((j, i), 0.0)) ** 2 for (i, j), w_ij in wfull.items())
    S2 = 0.0
    for i, (idx, w) in enumerate(neigh):
        col = sum(wfull.get((j, i), 0.0) for j in range(n))
        S2 += (w * len(idx) + col) ** 2
    b2 = n * sum(v ** 4 for v in zz) / (sst ** 2) if sst > 0 else 3.0
    A = n * ((n * n - 3 * n + 3) * S1 - n * S2 + 3 * S0 * S0)
    B = b2 * ((n * n - n) * S1 - 2 * n * S2 + 6 * S0 * S0)
    var = (A - B) / ((n - 1) * (n - 2) * (n - 3) * S0 * S0) - EI * EI
    z_stat = (I - EI) / math.sqrt(var) if var > 0 else 0.0
    # permutazione (seed fisso: deterministica)
    rng = random.Random(SEED)
    zz2 = zz[:]
    ge = 0
    for _ in range(N_PERM_GLOBAL):
        rng.shuffle(zz2)
        num_p = sum(zz2[i] * w * sum(zz2[j] for j in idx) for i, (idx, w) in enumerate(neigh))
        Ip = (n / S0) * num_p / sst if sst > 0 else 0.0
        if Ip >= I:
            ge += 1
    return {'n': n, 'I': round(I, 5), 'E_I': round(EI, 5), 'z_clifford': round(z_stat, 2),
            'p_perm_one_sided': round((ge + 1) / (N_PERM_GLOBAL + 1), 4)}


def lisa(z, neigh):
    """Moran locale con randomizzazione condizionata (Anselin 1995, seed fisso)."""
    n = len(z)
    m = sum(z) / n
    zz = [v - m for v in z]
    lag = [sum(zz[j] for j in idx) * w for idx, w in neigh]
    loc = [zz[i] * lag[i] for i in range(n)]
    rng = random.Random(SEED)
    out = []
    for i in range(n):
        idx, w = neigh[i]
        others = [zz[j] for j in range(n) if j != i]
        zi = zz[i]
        ext = 0
        for _ in range(N_PERM_LOCAL):
            samp = rng.sample(others, len(idx))
            lp = sum(samp) * w
            ext += 1 if abs(zi * lp) >= abs(loc[i]) else 0
        p = (ext + 1) / (N_PERM_LOCAL + 1)
        if p < ALPHA:
            cl = ('HH' if zz[i] > 0 and lag[i] > 0 else
                  'LL' if zz[i] <= 0 and lag[i] <= 0 else
                  'HL' if zz[i] > 0 else 'LH')
        else:
            cl = 'n.s.'
        out.append({'i': i, 'I_local': loc[i], 'p': p, 'cluster': cl})
    return out


# ---------- mappa SVG (dot map a pannello unico, stesso stile di mappa_loss.svg) ----------
def draw_svg(items, path):
    """items: dict PRO_COM -> (lon, lat, cluster)."""
    xs = [v[0] for v in items.values()]
    ys = [v[1] for v in items.values()]
    W, H = 1000, 600
    x0, y0, w, h = 60, 30, 470, 440
    sx = w / (max(xs) - min(xs))
    sy = h / (max(ys) - min(ys))
    s = min(sx, sy)
    ox = x0 + (w - s * (max(xs) - min(xs))) / 2
    oy = y0 + (h - s * (max(ys) - min(ys))) / 2
    COL = {'HH': '#9c1f1f', 'LL': '#2b5f9e', 'HL': '#c94b2a', 'LH': '#5b8fc4', 'n.s.': '#cfcfc9'}
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d" '
             'font-family="Helvetica,Arial,sans-serif">' % (W, H, W, H)]
    parts.append('<rect width="%d" height="%d" fill="#fafaf7"/>' % (W, H))
    parts.append('<text x="500" y="24" font-size="15" font-weight="bold" text-anchor="middle" fill="#222">'
                 'Coerenza tariffaria Cat-Nat: cluster spaziali del loss ratio</text>')
    parts.append('<text x="500" y="41" font-size="10.5" text-anchor="middle" fill="#666">'
                 'LISA del log(loss ratio) (Moran locale, KNN k=5, randomizzazione condizionata 499 permutazioni, '
                 'p&lt;0,05) — 3.820 comuni con loss ratio definito</text>')
    for pc, (lon, lat, cl) in sorted(items.items()):
        x = ox + (lon - min(xs)) * s
        y = oy + h - (lat - min(ys)) * s
        parts.append('<circle cx="%.1f" cy="%.1f" r="2.4" fill="%s" fill-opacity="0.9"/>'
                     % (x, y, COL[cl]))
    leg = [('HH', 'sovraprezzo in cluster di sovraprezzo'),
           ('LL', 'sottoprezzo in cluster di sottoprezzo'),
           ('HL', 'sovraprezzo tra sottoprezzo'),
           ('LH', 'sottoprezzo tra sovraprezzo'),
           ('n.s.', 'non significativo')]
    ly = H - 74
    parts.append('<text x="60" y="%d" font-size="10" fill="#555">cluster:</text>' % ly)
    for i, (cl, desc) in enumerate(leg):
        cx = 60 + (i % 2) * 450
        cy = ly + 12 + (i // 2) * 16
        parts.append('<rect x="%d" y="%d" width="11" height="11" fill="%s" stroke="#999" stroke-width="0.5"/>'
                     % (cx, cy, COL[cl]))
        parts.append('<text x="%d" y="%d" font-size="9.5" fill="#444">%s — %s</text>'
                     % (cx + 16, cy + 9, cl, desc))
    parts.append('<text x="500" y="%d" font-size="9" text-anchor="middle" fill="#888">'
                 'Fonti: IVASS (tariffe), INGV MPS04 (ag), ISPRA (PAI P3/P4, P3), AIDA (asset) · '
                 'generato da scripts/pricing/spazializzazione.py · dettagli in docs/pricing_coerenza.md</text>'
                 % (H - 10))
    parts.append('</svg>')
    with open(path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(parts))


def main():
    tess = {}
    with open(MATRICE, newline='', encoding='utf-8') as f:
        for r in csv.DictReader(f):
            tess[r['PRO_COM']] = float(r['EBITDA_migl_EUR']) * 1e3
    com = []
    with open(EAL_CSV, newline='', encoding='utf-8') as f:
        for r in csv.DictReader(f):
            com.append({'pc': r['PRO_COM'], 'nome': r['COMUNE'], 'prov': r['Provincia'],
                        'lat': float(r['lat']), 'lon': float(r['long']),
                        'lr': float(r['loss_ratio']) if r['loss_ratio'] else None,
                        'rate_bmk': float(r['rate_bmk_calibrato']),
                        'eal': float(r['eal_calibrato_EUR']),
                        'ebitda': tess.get(r['PRO_COM'])})
    assert len(com) == 3823, len(com)

    lat0 = sum(c['lat'] for c in com) / len(com)
    cs = math.cos(math.radians(lat0))
    proj = lambda c: (c['lon'] * cs, c['lat'])

    def moran_of(items, label, note):
        pts = [proj(c) for c in items]
        neigh = knn_rows(pts)
        res = moran_global([c['v'] for c in items], neigh)
        res['variabile'] = label
        res['nota'] = note
        return res

    # loss ratio (n=3.820 con ratio definito): grezzo e log (dichiarato: la distribuzione
    # e' estremamente right-skewed, max ~13.429, il Moran grezzo e' pilotato dagli outlier)
    lr_raw = [dict(c, v=c['lr']) for c in com if c['lr'] is not None]
    lr_log = [dict(c, v=math.log(c['lr'])) for c in com if c['lr'] is not None]
    # log peso-EBITDA (come la regressione del tessuto: EBITDA>0 e EAL>0)
    pe = [dict(c, v=math.log(c['eal'] / c['ebitda']))
          for c in com if c['ebitda'] and c['ebitda'] > 0 and c['eal'] > 0]
    # rate benchmark calibrata
    rb = [dict(c, v=c['rate_bmk']) for c in com]

    glob = [moran_of(lr_raw, 'loss_ratio (grezzo)', 'tariffa/benchmark: pilotato dalla coda (max ~13.429), riportato per completezza'),
            moran_of(lr_log, 'log loss_ratio', 'log(tariffa/benchmark): trasformazione dichiarata per la skewness — lettura primaria'),
            moran_of(pe, 'log peso-EBITDA', 'log(EAL/EBITDA): autocorrelazione del peso della loss sul tessuto'),
            moran_of(rb, 'rate_bmk', 'rate benchmark calibrata per 10.000 EUR di asset')]

    # ---------- LISA del log loss ratio + mappa ----------
    pts = [proj(c) for c in lr_log]
    neigh = knn_rows(pts)
    lis = lisa([c['v'] for c in lr_log], neigh)
    counts = {}
    for l in lis:
        counts[l['cluster']] = counts.get(l['cluster'], 0) + 1
    # comuni per cluster, ordinati per p poi |I| (deterministico)
    top = {}
    for cl in ('HH', 'LL'):
        cand = [(lis[i]['p'], -abs(lis[i]['I_local']), lr_log[i]['nome'], lr_log[i]['prov'], round(lr_log[i]['lr'], 2))
                for i in range(len(lr_log)) if lis[i]['cluster'] == cl]
        cand.sort()
        top[cl] = [{'comune': c[2], 'provincia': c[3], 'loss_ratio': c[4], 'p': round(c[0], 3)}
                   for c in cand[:12]]
    draw_svg({lr_log[i]['pc']: (lr_log[i]['lon'], lr_log[i]['lat'], lis[i]['cluster']) for i in range(len(lr_log))}, OUT_SVG)

    out = {
        'modello': 'autocorrelazione spaziale e cluster LISA della coerenza tariffaria (loss ratio)',
        'fonti': 'results/eal_comuni.csv (pricing_model.py), centroidi comunali; EBITDA dalla matrice estesa',
        'metodo': {
            'W': 'KNN k=%d sui centroidi, pesi di riga standardizzati, tie-break su (distanza, PRO_COM)' % KNN,
            'moran_globale': 'statistica analitica Cliff-Ord (normalita\') + permutazioni %d (seed %d)' % (N_PERM_GLOBAL, SEED),
            'lisa': 'randomizzazione condizionata, %d permutazioni (seed %d), pseudo-p bilatero, soglia p<%.2f'
                     % (N_PERM_LOCAL, SEED, ALPHA),
        },
        'moran_globale': glob,
        'lisa_loss_ratio': {
            'variabile': 'log(loss_ratio) (trasformazione dichiarata per la skewness)',
            'conteggio_cluster': counts,
            'nota_interpretazione': 'HH = cluster di sovraprezzo (tariffa sopra il rischio modellato); '
                                    'LL = cluster di sottoprezzo (tariffa sotto il rischio modellato)',
            'cluster_HH': top['HH'],
            'cluster_LL': top['LL'],
        },
    }
    with open(OUT_JSON, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
        f.write('\n')
    for g in glob:
        print('Moran %-18s n=%4d  I=%+.4f  z=%+.1f  p_perm=%.4f' % (g['variabile'], g['n'], g['I'], g['z_clifford'], g['p_perm_one_sided']))
    print('LISA loss ratio:', counts)
    print('OK: %s' % os.path.relpath(OUT_JSON, ROOT))
    print('OK: %s' % os.path.relpath(OUT_SVG, ROOT))


if __name__ == '__main__':
    main()
