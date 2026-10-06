#!/usr/bin/env python3
"""Quadro grafico dei risultati SDM: sei pannelli + mappe comunali degli hazard.

Come scripts/pricing/grafici.py: ogni pannello legge i valori direttamente dai
JSON di results/ e dalle matrici di data/ (niente numeri hardcoded), la figura
e' una vista, non un'origine di dati, e la CI ne verifica la riproducibilita'
byte per byte insieme al resto.

Output 1: docs/grafici_sdm.svg - quadro a sei pannelli dei risultati SDM:

1. COEFFICIENTI SDM p=4: dot plot con IC 95% (stima +- 1,96 SE) degli otto
   regressori Risk_* e W_Risk_* (da results/FINAL_sismico_k5.json); marker
   pieno = significativo (p<0,05), vuoto = non significativo;
2. EFFETTI LE SAGE-PACE p=4: barre diretti/indiretti/totale per Sismico
   Grandi, Sismico PMI e Frana Grandi (dallo stesso JSON);
3. ROBUSTEZZA DEL BETA SISMICO GRANDI: stima con IC 95% nelle nove
   specifiche (baseline ag RP475, RP30, RP50, Sa01 RP1000/RP2500, k=6/7/8,
   esclusione Sardegna);
4. SELEZIONE k: AIC di SDM/SAR/SEM vs k=3..30 (da results/grid_results.json),
   con k=5 evidenziato;
5. CONFRONTO DI SPECIFICA p=4: LR vs SAR, vs SEM e vs baseline p=2 (scala
   logaritmica) + diagnostica residui (Moran, RESET, BP);
6. RHO PER SPECIFICA: la dipendenza spaziale nelle stesse nove specifiche
   del pannello 3 (nota: noSardegna cambia anche la W, n=3.753).

Output 2: docs/mappe_hazard.svg - dot map comunali (stile di docs/mappa_loss.svg):
ag RP475 (MPS04, in g) e quota di area comunale in frana P3/P4
(hazard_frana_share, 0-1), classi ai quantili, ordinamento per PRO_COM.

Deterministico, solo stdlib. Esecuzione dalla root del repo:
    python3 scripts/grafici_sdm.py
"""
import csv, json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

RESULTS = os.path.join(ROOT, 'results')
FIN_JSON = os.environ.get('FIN_JSON', os.path.join(RESULTS, 'FINAL_sismico_k5.json'))
GRID_JSON = os.environ.get('GRID_JSON', os.path.join(RESULTS, 'grid_results.json'))
MATRICE = os.environ.get('MATRICE', os.path.join(ROOT, 'data', 'Matrice_Modello_Savelli_Final_sismico.csv'))
OUT_QUADRO = os.environ.get('OUT_QUADRO', os.path.join(ROOT, 'docs', 'grafici_sdm.svg'))
OUT_MAPPE = os.environ.get('OUT_MAPPE', os.path.join(ROOT, 'docs', 'mappe_hazard.svg'))

W, H = 1120, 760
PW, PH = 340, 290            # pannello
GX, GY = 20, 18              # gap tra i pannelli
MX, TOP = 20, 66              # margini esterni / inizio griglia
ROW2 = TOP + PH + 40          # seconda riga (spazio per didascalie)

INK, MUT, GRID = '#222', '#777', '#e4e2da'
BLU, ROSSO, VERDE = '#2b5f9e', '#9c1f1f', '#3f7d4e'
ARANCIO = '#b06a1d'


def esc(s):
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


def txt(x, y, s, size=9, fill=INK, anchor='start', bold=False):
    w = ' font-weight="bold"' if bold else ''
    return ('<text x="%.1f" y="%.1f" font-size="%s"%s fill="%s" text-anchor="%s">%s</text>'
            % (x, y, size, w, fill, anchor, esc(s)))


def txt2(x, y, lines, size=8, fill=MUT, anchor='middle'):
    t = ''.join('<tspan x="%.1f" y="%.1f">%s</tspan>' % (x, y, esc(ln)) if i == 0
                else '<tspan x="%.1f" dy="10">%s</tspan>' % (x, esc(ln))
                for i, ln in enumerate(lines))
    return '<text font-size="%s" fill="%s" text-anchor="%s">%s</text>' % (size, fill, anchor, t)


def line(x1, y1, x2, y2, stroke=GRID, dash='', width=1):
    d = ' stroke-dasharray="%s"' % dash if dash else ''
    return ('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="%s"%s/>'
            % (x1, y1, x2, y2, stroke, width, d))


def rect(x, y, w, h, fill, stroke='none'):
    return '<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s" stroke="%s"/>' \
           % (x, y, max(0.0, w), max(0.0, h), fill, stroke)


def panel(x0, y0, num, title, subtitle, labw=40):
    """Cornice del pannello; labw = spazio riservato alle etichette di riga."""
    parts = [rect(x0, y0, PW, PH, '#ffffff', '#d8d6ce')]
    parts.append(txt(x0 + 10, y0 + 20, '%d. %s' % (num, title), 11.5, INK, bold=True))
    parts.append(txt(x0 + 10, y0 + 34, subtitle, 8.6, MUT))
    px = x0 + 34 + labw
    return parts, (px, y0 + 48, PW - 48 - labw, PH - 72)   # (px, py, pw, ph)


def fmt_it(v, dec=2):
    return ('%.*f' % (dec, v)).replace('.', ',')


def header(titolo, sottotitolo):
    return ('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">\n'
            '<rect width="%d" height="%d" fill="#f6f5f0"/>\n'
            % (W, H, W, H, W, H)) + \
        txt(MX, 28, titolo, 15, INK, bold=True) + \
        txt(MX, 44, sottotitolo, 9, MUT)


def footer(fin):
    return txt(W - MX, H - 8,
               'fonte: results/FINAL_sismico_k5.json, results/grid_results.json (SDM ML k=5, n=%d comuni)'
               % fin['n'], 7.5, MUT, anchor='end')


# etichette dei coefficienti (nome JSON -> etichetta breve)
LBL = {
    'Risk_Sismico_Asset_PMI': 'Sismico PMI',
    'Risk_Sismico_Asset_Grandi': 'Sismico Grandi',
    'Risk_Frana_Asset_Grandi': 'Frana Grandi',
    'Risk_Frana_Asset_PMI': 'Frana PMI',
    'W_Risk_Sismico_Asset_PMI': 'W x Sism PMI',
    'W_Risk_Sismico_Asset_Grandi': 'W x Sism Grandi',
    'W_Risk_Frana_Asset_Grandi': 'W x Frana Grandi',
    'W_Risk_Frana_Asset_PMI': 'W x Frana PMI',
}

# specifiche di robustezza: (chiave JSON, etichetta, tick corto)
ROB = [('base', 'baseline ag RP475', 'base'),
       ('RP30', 'ag RP30 (81% in 50)', 'RP30'),
       ('RP50', 'ag RP50 (63% in 50)', 'RP50'),
       ('Sa01_RP1000', 'Sa(0,10s) RP1000', 'Sa1000'),
       ('Sa01_RP2500', 'Sa(0,10s) RP2500', 'Sa2500'),
       ('6', 'k = 6', 'k6'),
       ('7', 'k = 7', 'k7'),
       ('8', 'k = 8', 'k8'),
       ('noSardegna', 'esclusa Sardegna', 'noSard')]


def rob_stats(fin):
    """(etichetta, beta_sismGra, se, rho) per le nove specifiche."""
    out = []
    p4 = fin['sdm_p4']
    for key, lbl, _tick in ROB:
        if key == 'base':
            c = [c for c in p4['coefs'] if c['name'] == 'Risk_Sismico_Asset_Grandi'][0]
            rho = [c for c in p4['coefs'] if c['name'] == 'rho'][0]['est']
            out.append((lbl, c['est'], c['se'], rho))
        elif key.startswith('Sa') or key in ('RP30', 'RP50'):
            r = fin['robust_hazard'][key]
            out.append((lbl, r['b_sismGra'], r['se_sismGra'], r['rho']))
        elif key in ('6', '7', '8'):
            r = fin['robust_k'][key]
            out.append((lbl, r['b_sismGra'], r['se_sismGra'], r['rho']))
        else:
            r = fin['robust_noSardegna']
            out.append((lbl, r['b_sismGra'], r['se_sismGra'], r['rho']))
    return out


# ------------------------------------------------------------------ 1. beta + IC
def p_coef(x0, y0, fin):
    parts, (px, py, pw, ph) = panel(x0, y0, 1, 'Coefficienti SDM (p=4)',
                                    'stima e IC 95% (punto pieno = p < 0,05, vuoto = n.s.)',
                                    labw=88)
    coefs = [c for c in fin['sdm_p4']['coefs'] if c['name'] in LBL]
    coefs.sort(key=lambda c: -c['est'])
    lo = min(c['est'] - 1.96 * c['se'] for c in coefs)
    hi = max(c['est'] + 1.96 * c['se'] for c in coefs)
    m = 0.02
    lo = math.floor((lo - m) / 0.05) * 0.05
    hi = math.ceil((hi + m) / 0.05) * 0.05
    lx = lambda v: px + (v - lo) / (hi - lo) * pw
    n = len(coefs)
    dy = ph / (n + 0.5)
    tick = 0.05
    v = lo
    while v <= hi + 1e-9:
        parts.append(line(lx(v), py, lx(v), py + ph, '#bbb' if abs(v) < 1e-12 else GRID,
                          '' if abs(v) < 1e-12 else '2 3'))
        parts.append(txt(lx(v), py + ph + 11, fmt_it(v, 2), 8.3, MUT, anchor='middle'))
        v += tick
    for i, c in enumerate(coefs):
        yy = py + dy * (i + 0.75)
        ic_lo, ic_hi = c['est'] - 1.96 * c['se'], c['est'] + 1.96 * c['se']
        parts.append(line(lx(ic_lo), yy, lx(ic_hi), yy, BLU, width=1.6))
        for xb in (lx(ic_lo), lx(ic_hi)):
            parts.append(line(xb, yy - 2.6, xb, yy + 2.6, BLU, width=1.6))
        sig = c['p'] is not None and c['p'] < 0.05
        parts.append('<circle cx="%.1f" cy="%.1f" r="3.4" fill="%s" stroke="%s" stroke-width="1.4"/>'
                     % (lx(c['est']), yy, BLU if sig else '#ffffff', BLU))
        parts.append(txt(px - 4, yy + 3, LBL[c['name']], 8.2, INK, anchor='end'))
        parts.append(txt(lx(c['est']) + (5 if c['est'] >= 0 else -5), yy - 4,
                         fmt_it(c['est'], 3), 7.8, MUT, anchor='start' if c['est'] >= 0 else 'end'))
    parts.append(txt(px + pw / 2, py + ph + 24, 'beta (y = log1p premio; X = log1p risk)', 8.4, MUT, anchor='middle'))
    return parts


# ------------------------------------------------------- 2. effetti LeSage-Pace
def p_effetti(x0, y0, fin):
    parts, (px, py, pw, ph) = panel(x0, y0, 2, 'Effetti LeSage-Pace (p=4)',
                                    'diretti, indiretti (spillover) e totali per hazard',
                                    labw=72)
    eff = fin['sdm_p4']['effects']
    groups = [('sismGra', 'Sismico Grandi'), ('sismPMI', 'Sismico PMI'), ('franaGra', 'Frana Grandi')]
    vals = [(eff[k]['direct'], eff[k]['indirect'], eff[k]['total']) for k, _ in groups]
    lo = min(min(v) for v in vals)
    hi = max(max(v) for v in vals)
    lo = math.floor((lo - 0.02) / 0.05) * 0.05
    hi = math.ceil((hi + 0.02) / 0.05) * 0.05
    lx = lambda v: px + (v - lo) / (hi - lo) * pw
    colr = {'d': BLU, 'i': ROSSO, 't': VERDE}
    bh, bv, gs = 10, 3, 12          # altezza barra, gap verticale, gap tra gruppi
    slot = 3 * (bh + bv)
    total_h = 3 * slot + 2 * gs
    ytop = py + ph - total_h
    for gv in (lo, 0.0, hi):
        if lo < gv < hi or gv == 0:
            parts.append(line(lx(gv), ytop, lx(gv), ytop + total_h, '#bbb' if gv == 0 else GRID,
                              '' if gv == 0 else '2 3'))
            parts.append(txt(lx(gv), ytop + total_h + 11, fmt_it(gv, 2), 8.3, MUT, anchor='middle'))
    for gi, (k, lbl) in enumerate(groups):
        gyy = ytop + gi * (slot + gs)
        for bi, (vi, kind) in enumerate(zip(vals[gi], ('d', 'i', 't'))):
            yy = gyy + bi * (bh + bv)
            x1, x2 = lx(min(vi, 0)), lx(max(vi, 0))
            parts.append(rect(min(x1, x2), yy, abs(x2 - x1), bh, colr[kind]))
        parts.append(txt(px - 4, gyy + slot / 2 + 3, lbl, 8.4, INK, anchor='end'))
    lx0, ly0 = px + pw - 148, py + 4
    for i, (kind, name) in enumerate((('d', 'diretti'), ('i', 'indiretti'), ('t', 'totali'))):
        parts.append(rect(lx0 + i * 52, ly0, 8, 8, colr[kind]))
        parts.append(txt(lx0 + i * 52 + 11, ly0 + 7.5, name, 7.8, MUT))
    parts.append(txt(px + pw / 2, py + ph + 24, 'effetti su log1p(premio) di un raddoppio del risk', 8.4, MUT, anchor='middle'))
    return parts


# ------------------------------------------------- 3. robustezza beta sismico
def p_robustezza(x0, y0, fin):
    parts, (px, py, pw, ph) = panel(x0, y0, 3, 'Robustezza: beta Sismico Grandi',
                                    'nove specifiche, stima e IC 95%; linea = baseline',
                                    labw=104)
    stats = rob_stats(fin)
    lo = min(b - 1.96 * se for _, b, se, _ in stats)
    hi = max(b + 1.96 * se for _, b, se, _ in stats)
    lo = math.floor((lo - 0.004) / 0.005) * 0.005
    hi = math.ceil((hi + 0.004) / 0.005) * 0.005
    lx = lambda v: px + (v - lo) / (hi - lo) * pw
    n = len(stats)
    dy = ph / (n + 0.5)
    base = stats[0][1]
    parts.append(line(lx(base), py, lx(base), py + ph, ARANCIO, '3 3'))
    tick = 0.005
    v = lo
    while v <= hi + 1e-9:
        parts.append(line(lx(v), py, lx(v), py + ph, GRID, '2 3'))
        parts.append(txt(lx(v), py + ph + 11, fmt_it(v, 3), 8.3, MUT, anchor='middle'))
        v += tick
    for i, (lbl, b, se, _rho) in enumerate(stats):
        yy = py + dy * (i + 0.75)
        ic_lo, ic_hi = b - 1.96 * se, b + 1.96 * se
        parts.append(line(lx(ic_lo), yy, lx(ic_hi), yy, VERDE, width=1.6))
        for xb in (lx(ic_lo), lx(ic_hi)):
            parts.append(line(xb, yy - 2.6, xb, yy + 2.6, VERDE, width=1.6))
        parts.append('<circle cx="%.1f" cy="%.1f" r="3.2" fill="%s"/>' % (lx(b), yy, VERDE))
        parts.append(txt(px - 4, yy + 3, lbl, 8.2, INK, anchor='end'))
        parts.append(txt(lx(b) + 5, yy - 3, fmt_it(b, 4), 7.8, MUT))
    parts.append(txt(px + pw / 2, py + ph + 24, 'beta Risk Sismico Grandi (log1p risk)', 8.4, MUT, anchor='middle'))
    return parts


# ------------------------------------------------------------- 4. AIC vs k
def p_aic(x0, y0, fin, grid):
    parts, (px, py, pw, ph) = panel(x0, y0, 4, 'Selezione k: AIC',
                                    'SDM vs SAR vs SEM (p=2), k=3..30; k=5 evidenziato', labw=44)
    gs = sorted(grid, key=lambda g: g['k'])
    ks = [g['k'] for g in gs]
    aic = {m: [g[m]['aic'] for g in gs] for m in ('sdm', 'sar', 'sem')}
    allv = [v for m in aic for v in aic[m]]
    lo, hi = min(allv), max(allv)
    pad = (hi - lo) * 0.12
    lo, hi = lo - pad, hi + pad
    lx = lambda k: px + (k - ks[0]) / (ks[-1] - ks[0]) * pw
    ly = lambda v: py + ph - (v - lo) / (hi - lo) * ph
    for v in (lo, (lo + hi) / 2, hi):
        parts.append(line(px, ly(v), px + pw, ly(v)))
        parts.append(txt(px - 5, ly(v) + 3, fmt_it(v, 0), 8.3, MUT, anchor='end'))
    for m, name, col in (('sdm', 'SDM', BLU), ('sar', 'SAR', MUT), ('sem', 'SEM', ROSSO)):
        pts = [(lx(k), ly(v)) for k, v in zip(ks, aic[m])]
        parts.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="%s"/>'
                     % (' '.join('%.1f,%.1f' % p for p in pts), col, 2 if m == 'sdm' else 1.4))
    parts.append(line(lx(5), py, lx(5), py + ph, ARANCIO, '3 3'))
    parts.append(txt(lx(5), py + 10, 'k=5', 8.4, ARANCIO, anchor='middle'))
    for m, name, col in (('sdm', 'SDM', BLU), ('sar', 'SAR', MUT), ('sem', 'SEM', ROSSO)):
        parts.append(txt(px + pw - 4, ly(aic[m][ks.index(5)]) - 5, name, 8.4, col, anchor='end'))
    for k in (3, 10, 20, 30):
        parts.append(txt(lx(k), py + ph + 11, str(k), 8.3, MUT, anchor='middle'))
    parts.append(txt(px + pw / 2, py + ph + 24, 'vicini k (KNN, matrice W)', 8.4, MUT, anchor='middle'))
    parts.append(txt(px - 5, py + 4, 'AIC', 8.3, MUT, anchor='end'))
    return parts


# --------------------------------------------- 5. confronto di specifica p=4
def p_specifica(x0, y0, fin):
    parts, (px, py, pw, ph) = panel(x0, y0, 5, 'Confronto di specifica (p=4)',
                                    'LR vs alternative annidate (scala log) + diagnostica',
                                    labw=88)
    p4 = fin['sdm_p4']
    tests = [('vs baseline p=2', p4['lr_vs_baseline_p2'], p4['df_nested']),
             ('vs SAR', p4['lr_vs_SAR'], 4),
             ('vs SEM', p4['lr_vs_SEM'], 4)]
    lo, hi = 1.0, max(t[1] for t in tests) * 1.6
    lx = lambda v: px + (math.log10(v) - math.log10(lo)) / (math.log10(hi) - math.log10(lo)) * pw
    bh, gap = 16, 13
    for i, (name, v, df) in enumerate(tests):
        yy = py + ph - 34 - (len(tests) - 1 - i) * (bh + gap)
        parts.append(rect(px, yy, lx(v) - px, bh, BLU))
        parts.append(txt(lx(v) + 5, yy + bh - 5, fmt_it(v, 1), 8.6, INK, bold=True))
        parts.append(txt(px - 4, yy + bh / 2 + 2, name, 8.3, INK, anchor='end'))
        parts.append(txt(px - 4, yy + bh / 2 + 11, 'df=%d' % df, 7.2, MUT, anchor='end'))
    parts.append(line(px, py, px, py + ph, '#bbb'))
    mo, re_, bp = p4['moranResid'], p4['reset'], p4['bp']
    parts.append(txt(px, py + ph - 4,
                     'Moran residui I=%s (p=%s) | RESET F=%s | BP LM=%s (df=%s)'
                     % (fmt_it(mo['I'], 4), fmt_it(mo['p'], 3), fmt_it(re_['F'], 2),
                        fmt_it(bp['LM'], 1), bp['df']), 7.8, MUT))
    return parts


# ------------------------------------------------------------- 6. rho
def p_rho(x0, y0, fin):
    parts, (px, py, pw, ph) = panel(x0, y0, 6, 'Dipendenza spaziale (rho)',
                                    'rho stimato nelle nove specifiche del pannello 3', labw=44)
    stats = rob_stats(fin)
    ticks = [t for _, _, t in ROB]
    lo, hi = 0.15, 0.55
    ly = lambda v: py + ph - (v - lo) / (hi - lo) * ph
    for v in (0.2, 0.3, 0.4, 0.5):
        parts.append(line(px, ly(v), px + pw, ly(v)))
        parts.append(txt(px - 5, ly(v) + 3, fmt_it(v, 1), 8.3, MUT, anchor='end'))
    n = len(stats)
    dx = pw / (n - 1)
    for i, (_lbl, _b, _se, rho) in enumerate(stats):
        xx = px + dx * i
        col = ROSSO if ticks[i] == 'noSard' else BLU
        parts.append('<circle cx="%.1f" cy="%.1f" r="3.6" fill="%s"/>' % (xx, ly(rho), col))
        parts.append(txt(xx, py + ph + 10, ticks[i], 7.4, MUT, anchor='middle'))
    parts.append(line(px, py, px, py + ph, '#bbb'))
    parts.append(line(px, py + ph, px + pw, py + ph, '#bbb'))
    parts.append(txt(px + pw - 4, py + 10, 'noSardegna: n=3.753, W ricalcolata', 7.4, ROSSO, anchor='end'))
    return parts


# ============================================================ quadro
def quadro(fin, grid):
    parts = [header('SDM Cat-Nat: quadro dei risultati',
                   'p=4 (Frana + Sismico x PMI/Grandi) da FINAL_sismico_k5.json | '
                   'selezione k e confronto modelli da grid_results.json')]
    x0 = MX
    parts += p_coef(x0, TOP, fin)
    parts += p_effetti(x0 + PW + GX, TOP, fin)
    parts += p_robustezza(x0 + 2 * (PW + GX), TOP, fin)
    parts += p_aic(x0, ROW2, fin, grid)
    parts += p_specifica(x0 + PW + GX, ROW2, fin)
    parts += p_rho(x0 + 2 * (PW + GX), ROW2, fin)
    parts.append(footer(fin))
    return ''.join(parts) + '\n</svg>\n'


# ============================================================ mappe hazard
PALETTE_SIS = ['#f2ecdc', '#e8d5a8', '#dbb06a', '#c9813f', '#a84f26', '#7c2d16']
PALETTE_FRA = ['#eef2ea', '#cfdcc2', '#a8c69a', '#7ba86f', '#4f844f', '#2e6242']
NCLASSI = 6


def quantili(valori, n=NCLASSI):
    """breakpoint quantili deterministici (lista ordinata)."""
    s = sorted(valori)
    return [s[min(len(s) - 1, int(math.ceil(len(s) * i / n)) - 1)] for i in range(1, n)]


def classe(v, breaks):
    for i, b in enumerate(breaks):
        if v <= b:
            return i
    return len(breaks)


def pannello_mappa(x0, y0, w, h, titolo, sub, punti, palette, unit, nota,
                   breaks=None, legenda_classi='classi ai quantili (%s)', legenda_dec=2):
    """punti: lista ordinata (lon, lat, valore); breaks=None -> quantili."""
    parts = [rect(x0, y0, w, h, '#ffffff', '#d8d6ce')]
    parts.append(txt(x0 + 10, y0 + 18, titolo, 11.5, INK, bold=True))
    parts.append(txt(x0 + 10, y0 + 31, sub, 8.4, MUT))
    mx, my, mt = 26, 42, 40
    lonmin = min(p[0] for p in punti); lonmax = max(p[0] for p in punti)
    latmin = min(p[1] for p in punti); latmax = max(p[1] for p in punti)
    lat0 = (latmin + latmax) / 2
    cs = math.cos(math.radians(lat0))
    # proiezione equidistante con correzione cos(lat), come docs/mappa_loss.svg
    span_x = (lonmax - lonmin) * cs
    span_y = (latmax - latmin)
    s = min((w - 2 * mx) / span_x, (h - my - mt) / span_y)
    ox = x0 + (w - span_x * s) / 2
    oy = y0 + my
    vals = [p[2] for p in punti]
    if breaks is None:
        breaks = quantili(vals)
    else:
        breaks = sorted(breaks)
    for i, (lon, lat, v) in enumerate(punti):
        x = ox + (lon - lonmin) * cs * s
        y = oy + (latmax - lat) * s
        parts.append(('<circle cx="%.2f" cy="%.2f" r="1.55" fill="%s"/>' + ('\n' if i % 32 == 31 else ''))
                     % (x, y, palette[classe(v, breaks)]))
    # legenda
    n1 = max(vals)
    edges = [min(vals)] + breaks
    lx0, ly0 = x0 + 12, y0 + h - 54
    parts.append(txt(lx0, ly0 - 6, legenda_classi % unit, 7.6, MUT))
    cw = (w - 24) / NCLASSI
    for i in range(NCLASSI):
        parts.append(rect(lx0 + i * cw, ly0, cw, 10, palette[i], '#ffffff'))
        parts.append(txt(lx0 + i * cw, ly0 + 21, fmt_it(edges[i], legenda_dec), 7.2, MUT))
    parts.append(txt(lx0 + NCLASSI * cw, ly0 + 21, fmt_it(n1, legenda_dec), 7.2, MUT, anchor='end'))
    parts.append(txt(x0 + 12, y0 + h - 8, nota, 7.6, MUT))
    return parts


def mappe():
    rows = []
    with open(MATRICE, newline='', encoding='utf-8') as f:
        for r in csv.DictReader(f):
            rows.append(r)
    sis, fra = [], []
    for r in sorted(rows, key=lambda r: r['PRO_COM']):   # determinismo
        lon, lat = float(r['long']), float(r['lat'])
        ag = float(r['ag_RP475'])
        share = float(r['PAI_area_P3P4_kmq']) / float(r['SUP_kmq'])
        share = max(0.0, min(1.0, share))
        sis.append((lon, lat, ag))
        fra.append((lon, lat, share))
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" width="1040" height="640" viewBox="0 0 1040 640">\n'
             '<rect width="1040" height="640" fill="#f6f5f0"/>\n']
    parts.append(txt(20, 28, 'Hazard comunali: sismico (MPS04) e frana (PAI P3/P4)', 15, INK, bold=True))
    parts.append(txt(20, 44, 'dot map dei 3.823 comuni (data/Matrice_Modello_Savelli_Final_sismico.csv), '
                             'classi ai quantili, proiezione equidistante', 9, MUT))
    parts += pannello_mappa(20, 58, 496, 556,
                            'ag RP 475 anni (accelerazione, g)',
                            'MPS04 INGV, 10% in 50 anni; 70 comuni sardi non classificati (ag = 0)',
                            sis, PALETTE_SIS, 'g',
                            'Sardegna: hazard sismico non classificato nelle mappe MPS04 (ag = 0)')
    parts += pannello_mappa(524, 58, 496, 556,
                            'Quota di area comunale in frana P3/P4',
                            'hazard_frana_share = PAI_area_P3P4_kmq / SUP_kmq (ISPRA, 0-1)',
                            fra, PALETTE_FRA, 'quota area',
                            'Fonte: aree PAI P3/P4 ISPRA; classe 0 = quota nulla o trascurabile',
                            breaks=[0.001, 0.01, 0.05, 0.15, 0.30],
                            legenda_classi='classi fisse (%s)', legenda_dec=3)
    return ''.join(parts) + '\n</svg>\n'


def main():
    with open(FIN_JSON, encoding='utf-8') as f:
        fin = json.load(f)
    with open(GRID_JSON, encoding='utf-8') as f:
        grid = json.load(f)
    with open(OUT_QUADRO, 'w', encoding='utf-8') as f:
        f.write(quadro(fin, grid))
    with open(OUT_MAPPE, 'w', encoding='utf-8') as f:
        f.write(mappe())
    print('OK: %s' % os.path.relpath(OUT_QUADRO, ROOT))
    print('OK: %s' % os.path.relpath(OUT_MAPPE, ROOT))


if __name__ == '__main__':
    main()
