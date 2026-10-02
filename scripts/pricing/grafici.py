#!/usr/bin/env python3
"""Quadro grafico dei risultati del pricing: sei pannelli SVG dai JSON in results/.

Ogni pannello legge i valori direttamente dai JSON generati dagli altri script
(niente numeri hardcoded): la figura e' una vista, non un'origine di dati, e la
CI ne verifica la riproducibilita' byte per byte insieme al resto.

Pannelli (docs/pricing_coerenza.md: §2.1 EP, §5.2 quartili, §5.4 chi paga, §5.3
province, §5.5 robustezza, §6 Moran):

1. CURVA EP SISMICA NAZIONALE (log-log): loss a scenario RP30/72/475 con banda
   epistemica 16/84 al RP475, e le linee di riferimento AAL numerico (1,71 mld,
   limite inferiore RP>=30) e sismico di benchmark CURVE=4 (1,03 mld);
2. PESO DELLA LOSS PER QUARTILE ISP: barre del peso-EBITDA (Q1 3,91% -> Q4 2,06%)
   con l'intensita' di esposizione (quota EAL / quota EBITDA) dentro le barre;
3. CHI PAGA: quote di imprese, premio ed EAL per PMI vs Grandi (le Grandi sono
   il 10,5% delle imprese ma pagano il 63,7% del premio);
4. PROVINCE ESTREME per peso-EBITDA: top e bottom della tabella (107 province);
5. ROBUSTEZZA DEL BETA DELLA PERFORMANCE: le cinque specifiche a confronto
   (base, winsorizzato, trim, SLX, ROA per deviazione standard);
6. MORAN GLOBALE: autocorrelazione delle quattro variabili spazializzate.

Output: docs/grafici_pricing.svg (deterministico, solo stdlib).
Esecuzione dalla root del repo dopo pricing_model.py e gli approfondimenti:
    python3 scripts/pricing/grafici.py
"""
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))

RESULTS = os.path.join(ROOT, 'results')
EP_JSON = os.environ.get('EP_JSON', os.path.join(RESULTS, 'ep_curve.json'))
CHI_JSON = os.environ.get('CHI_JSON', os.path.join(RESULTS, 'chi_paga.json'))
TESS_JSON = os.environ.get('TESS_JSON', os.path.join(RESULTS, 'esposizione_tessuto.json'))
ROB_JSON = os.environ.get('ROB_JSON', os.path.join(RESULTS, 'robustezza_tessuto.json'))
SPA_JSON = os.environ.get('SPA_JSON', os.path.join(RESULTS, 'spazializzazione_tariffa.json'))
OUT_SVG = os.environ.get('OUT_SVG', os.path.join(ROOT, 'docs', 'grafici_pricing.svg'))

W, H = 1120, 760
PW, PH = 340, 290            # pannello
GX, GY = 20, 18              # gap tra i pannelli
MX, TOP = 20, 66              # margini esterni / inizio griglia
ROW2 = TOP + PH + 40          # seconda riga (spazio per didascalie)

INK, MUT, GRID = '#222', '#777', '#e4e2da'
BLU, ROSSO, VERDE = '#2b5f9e', '#9c1f1f', '#3f7d4e'


def txt(x, y, s, size=9, fill=INK, anchor='start', bold=False):
    w = ' font-weight="bold"' if bold else ''
    return ('<text x="%.1f" y="%.1f" font-size="%s"%s fill="%s" text-anchor="%s">%s</text>'
            % (x, y, size, w, fill, anchor, s))


def txt2(x, y, lines, size=8, fill=MUT, anchor='middle'):
    """Testo su piu' righe: gli \\n non renderizzano in SVG, servono tspan."""
    t = ''.join('<tspan x="%.1f" y="%.1f">%s</tspan>' % (x, y, ln) if i == 0
                else '<tspan x="%.1f" dy="10">%s</tspan>' % (x, ln)
                for i, ln in enumerate(lines))
    return '<text font-size="%s" fill="%s" text-anchor="%s">%s</text>' % (size, fill, anchor, t)


def line(x1, y1, x2, y2, stroke=GRID, dash=''):
    d = ' stroke-dasharray="%s"' % dash if dash else ''
    return '<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1"%s/>' \
           % (x1, y1, x2, y2, stroke, d)


def rect(x, y, w, h, fill, stroke='none'):
    return '<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s" stroke="%s"/>' \
           % (x, y, max(0.0, w), max(0.0, h), fill, stroke)


def panel(x0, y0, num, title, subtitle):
    """Cornice del pannello: titolo numerato + sottotitolo; ritorna l'area del tracciato."""
    parts = [rect(x0, y0, PW, PH, '#ffffff', '#d8d6ce')]
    parts.append(txt(x0 + 10, y0 + 20, '%d. %s' % (num, title), 11.5, INK, bold=True))
    parts.append(txt(x0 + 10, y0 + 34, subtitle, 8.6, MUT))
    return parts, (x0 + 34, y0 + 48, PW - 48, PH - 72)   # (px, py, pw, ph)


def fmt_it(v, dec=1):
    return ('%.*f' % (dec, v)).replace('.', ',')


def fmt_migliaia(v):
    return '{:,}'.format(int(v)).replace(',', '.')


# ------------------------------------------------------------------ 1. curva EP
def p_ep(x0, y0, ep):
    parts, (px, py, pw, ph) = panel(x0, y0, 1, 'Curva EP sismica nazionale',
                                    'loss fisica a scenario per periodo di ritorno (log-log, mld EUR)')
    e = ep['ep_nazionale_sismica']
    aal = ep['aal_numerico_sismico']['totale_k3_EUR'] / 1e9
    bmk = ep['aal_numerico_sismico']['benchmark_CURVE4_EUR'] / 1e9
    rps = [(30, e['RP30'] / 1e9), (72, e['RP72'] / 1e9), (475, e['RP475'] / 1e9)]
    band = e['RP475_banda_epistemica_16_84']
    lx = lambda rp: px + (math.log10(rp) - math.log10(20)) / (math.log10(600) - math.log10(20)) * pw
    ly = lambda v: py + ph - (math.log10(v) - math.log10(0.8)) / (math.log10(220) - math.log10(0.8)) * ph
    for t in (1, 10, 100):
        parts.append(line(px, ly(t), px + pw, ly(t)))
        parts.append(txt(px - 5, ly(t) + 3, str(t), 8.5, MUT, anchor='end'))
    for rp in (30, 72, 475):
        parts.append(line(lx(rp), py, lx(rp), py + ph))
        parts.append(txt(lx(rp), py + ph + 12, str(rp), 8.5, MUT, anchor='middle'))
    parts.append(line(px, py, px, py + ph, '#bbb'))
    parts.append(line(px, py + ph, px + pw, py + ph, '#bbb'))
    parts.append(txt(px + pw / 2, py + ph + 23, 'periodo di ritorno (anni)', 8.5, MUT, anchor='middle'))
    pts = [(lx(rp), ly(v)) for rp, v in rps]
    parts.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="2"/>'
                 % (' '.join('%.1f,%.1f' % p for p in pts), ROSSO))
    for (X, Y), (rp, v) in zip(pts, rps):
        parts.append('<circle cx="%.1f" cy="%.1f" r="3.4" fill="%s"/>' % (X, Y, ROSSO))
        parts.append(txt(X, Y - 8, fmt_it(v, 1), 8.6, ROSSO, anchor='middle', bold=True))
    bx = lx(475)
    parts.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="6" stroke-opacity="0.35"/>'
                 % (bx, ly(band['p84'] / 1e9), bx, ly(band['p16'] / 1e9), ROSSO))
    parts.append(txt(bx - 7, ly(band['p16'] / 1e9) + 3, '16-84 perc.', 8, MUT, anchor='end'))
    parts.append(line(px, ly(aal), px + pw, ly(aal), BLU, '4,3'))
    parts.append(txt(px + pw - 3, ly(aal) - 4, 'AAL numerico (RP≥30) %s' % fmt_it(aal, 2), 8.4, BLU, anchor='end'))
    parts.append(line(px, ly(bmk), px + pw, ly(bmk), VERDE, '4,3'))
    parts.append(txt(px + pw - 3, ly(bmk) + 10, 'sismico benchmark (CURVE=4) %s' % fmt_it(bmk, 2), 8.4, VERDE, anchor='end'))
    return parts


# ------------------------------------------------------------ 2. quartili ISP
def p_quartili(x0, y0, tess):
    parts, (px, py, pw, ph) = panel(x0, y0, 2, 'Peso della loss per performance (ISP)',
                                    'EAL/EBITDA per quartile di ISP — il doppio peso dei meno performanti')
    q = tess['quartili_ISP']
    vals = [c['peso_ebitda'] * 100 for c in q]
    inten = [c['intensita_esposizione (quotaEAL/quotaEBITDA)'] for c in q]
    vmax = 4.6
    bw = pw / len(vals) * 0.52
    step = pw / len(vals)
    for t in (1, 2, 3, 4):
        parts.append(line(px, py + ph - t / vmax * ph, px + pw, py + ph - t / vmax * ph))
        parts.append(txt(px - 5, py + ph - t / vmax * ph + 3, str(t), 8.5, MUT, anchor='end'))
    for i, v in enumerate(vals):
        cx = px + step * i + step / 2
        parts.append(rect(cx - bw / 2, py + ph - v / vmax * ph, bw, v / vmax * ph,
                          ROSSO if i == 0 else ('#c94b2a' if i == 1 else ('#e08a3c' if i == 2 else '#f0c184'))))
        parts.append(txt(cx, py + ph - v / vmax * ph - 5, fmt_it(v, 2) + '%', 9, INK, anchor='middle', bold=True))
        parts.append(txt(cx, py + ph - v / vmax * ph + 13, fmt_it(inten[i], 2), 8.2, '#fff', anchor='middle', bold=True))
        parts.append(txt(cx, py + ph + 12, 'Q%d' % (i + 1), 8.6, MUT, anchor='middle'))
    parts.append(line(px, py + ph, px + pw, py + ph, '#bbb'))
    parts.append(txt(px + pw / 2, py + ph + 22,
                     'quartili ISP (Q1 = meno performanti); in barra: intensita\' quotaEAL/quotaEBITDA',
                     7.8, MUT, anchor='middle'))
    return parts


# ----------------------------------------------------------------- 3. chi paga
def p_chi(x0, y0, chi):
    parts, (px, py, pw, ph) = panel(x0, y0, 3, 'Chi paga: PMI vs Grandi',
                                    'quote % di imprese, premio teorico ed EAL calibrata per classe dimensionale')
    naz = chi['nazionale']
    n_pmi, n_g = naz['n_imprese']['PMI'], naz['n_imprese']['grandi']
    n_tot = n_pmi + n_g
    qeal_pmi = naz['eal_calibrata_per_classe']['quota_eal_PMI'] * 100
    pmi = [n_pmi / n_tot * 100, naz['quota_premio_PMI'] * 100, qeal_pmi]
    grandi = [n_g / n_tot * 100, naz['quota_premio_grandi'] * 100, 100 - qeal_pmi]
    labels = ['imprese', 'premio', 'EAL calibrata']
    vmax, step, bw = 100.0, pw / 3, pw / 3 * 0.27
    for t in (25, 50, 75, 100):
        parts.append(line(px, py + ph - t / vmax * ph, px + pw, py + ph - t / vmax * ph))
        parts.append(txt(px - 5, py + ph - t / vmax * ph + 3, '%d' % t, 8.5, MUT, anchor='end'))
    for i, lab in enumerate(labels):
        cx = px + step * i + step / 2
        for j, (vals, col) in enumerate(((pmi, BLU), (grandi, ROSSO))):
            bx = cx - bw - 2 + j * (bw + 4)
            v = vals[i]
            parts.append(rect(bx, py + ph - v / vmax * ph, bw, v / vmax * ph, col))
            parts.append(txt(bx + bw / 2, py + ph - v / vmax * ph - 4, fmt_it(v, 0), 8, col, anchor='middle', bold=True))
        parts.append(txt(cx, py + ph + 12, lab, 8.6, MUT, anchor='middle'))
    parts.append(line(px, py + ph, px + pw, py + ph, '#bbb'))
    parts.append(rect(px, py - 8, 9, 9, BLU))
    parts.append(txt(px + 13, py, 'PMI (%s impr.)' % fmt_migliaia(n_pmi), 8.4, INK))
    parts.append(rect(px + 130, py - 8, 9, 9, ROSSO))
    parts.append(txt(px + 143, py, 'Grandi (%s impr.)' % fmt_migliaia(n_g), 8.4, INK))
    pm, pg = naz['premio_medio_impresa_PMI_EUR'], naz['premio_medio_impresa_grande_EUR']
    parts.append(txt(px + pw, py - 8, 'premio medio: %s vs %s EUR (x%.0f)'
                     % (fmt_migliaia(pm), fmt_migliaia(pg), pg / pm), 8, MUT, anchor='end'))
    return parts


# ------------------------------------------------------------ 4. province estreme
def p_province(x0, y0, tess):
    parts, (px, py, pw, ph) = panel(x0, y0, 4, 'Province estreme: peso della loss',
                                    'EAL/EBITDA per provincia — top e bottom delle 107 province')
    prov = tess['province_peso_ebitda']
    rows = [(r['provincia'], r['peso_ebitda'] * 100, True) for r in prov[:5]] + \
           [(r['provincia'], r['peso_ebitda'] * 100, False) for r in prov[-5:]]
    vmax, rh = 9.0, ph / 10
    for t in (2, 4, 6, 8):
        parts.append(line(px + t / vmax * pw, py, px + t / vmax * pw, py + ph))
        parts.append(txt(px + t / vmax * pw, py + ph + 11, str(t), 7.8, MUT, anchor='middle'))
    for i, (nome, v, is_top) in enumerate(rows):
        yy = py + rh * i + rh / 2
        col = ROSSO if is_top else BLU
        bw_px = v / vmax * pw
        parts.append(rect(px, yy - rh * 0.32, bw_px, rh * 0.64, col))
        parts.append(txt(px + bw_px + 5, yy + 3, fmt_it(v, 1) + '%', 8.2, col))
        if bw_px >= 110:
            parts.append(txt(px + 4, yy + 3, nome, 8.2, '#fff'))
        else:
            parts.append(txt(px + bw_px + 32, yy + 3, nome, 8.2, INK))
    parts.append(txt(px + pw / 2, py + ph + 22, 'peso EAL/EBITDA (%) — rosso: top 5, blu: bottom 5', 7.8, MUT, anchor='middle'))
    return parts


# --------------------------------------------------------------- 5. robustezza
def p_robustezza(x0, y0, rob):
    parts, (px, py, pw, ph) = panel(x0, y0, 5, 'Robustezza del beta della performance',
                                    'coefficiente del log peso-EBITDA per deviazione standard, cinque specifiche')
    base = [(s['specifica'], s['beta_x_deviazione_standard'], s['t_rob']) for s in rob['specifiche']]
    labels = [('base (ISP)', base[0][1], base[0][2]),
              ('winsorizzato 1%/99%', base[1][1], base[1][2]),
              ('trim top 1%', base[2][1], base[2][2]),
              ('SLX (+ W x regressori)', rob['slx']['coeff'][1]['beta'], rob['slx']['coeff'][1]['t_rob']),
              ('ROA comunale (per SD)', base[3][1], base[3][2])]
    xmin, xmax = -0.55, -0.05
    lx = lambda v: px + (v - xmin) / (xmax - xmin) * pw
    rh = ph / len(labels)
    for t in (-0.5, -0.4, -0.3, -0.2, -0.1):
        parts.append(line(lx(t), py, lx(t), py + ph))
        parts.append(txt(lx(t), py + ph + 11, fmt_it(t, 1), 7.8, MUT, anchor='middle'))
    for i, (lab, b, t) in enumerate(labels):
        yy = py + rh * i + rh / 2
        col = VERDE if 'ROA' in lab else BLU
        parts.append(txt(px + pw - 4, yy + 3, lab, 8.2, INK, anchor='end'))
        parts.append(line(lx(xmin), yy, lx(b), yy, '#ddd6c8'))
        parts.append('<circle cx="%.1f" cy="%.1f" r="4.2" fill="%s"/>' % (lx(b), yy, col))
        parts.append(txt(lx(b), yy - 8, '%s (t=%+.1f)' % (fmt_it(b, 3), t), 7.8, col, anchor='middle', bold=True))
    parts.append(txt(px + pw / 2, py + ph + 22,
                     'beta della performance per deviazione standard (negativo = performance riduce il peso)',
                     7.8, MUT, anchor='middle'))
    return parts


# ------------------------------------------------------------------- 6. Moran
def p_moran(x0, y0, spa):
    parts, (px, py, pw, ph) = panel(x0, y0, 6, 'Autocorrelazione spaziale (Moran I)',
                                    'quattro variabili su W KNN k=5; la coerenza tariffaria e\' un fatto spaziale')
    by = {r['variabile']: r for r in spa['moran_globale']}
    rows = [(by['loss_ratio (grezzo)'], ['loss ratio', '(grezzo)']),
            (by['log loss_ratio'], ['log loss', 'ratio']),
            (by['log peso-EBITDA'], ['log', 'peso-EBITDA']),
            (by['rate_bmk'], ['rate', 'benchmark'])]
    vmax = 0.85
    step = pw / len(rows)
    bw = step * 0.5
    for t in (0.2, 0.4, 0.6, 0.8):
        parts.append(line(px, py + ph - t / vmax * ph, px + pw, py + ph - t / vmax * ph))
        parts.append(txt(px - 5, py + ph - t / vmax * ph + 3, fmt_it(t, 1), 8.5, MUT, anchor='end'))
    for i, (r, lab_lines) in enumerate(rows):
        cx = px + step * i + step / 2
        col = ROSSO if r['variabile'] == 'log loss_ratio' else BLU
        parts.append(rect(cx - bw / 2, py + ph - r['I'] / vmax * ph, bw, r['I'] / vmax * ph, col))
        parts.append(txt(cx, py + ph - r['I'] / vmax * ph - 12, fmt_it(r['I'], 3), 8.6, col, anchor='middle', bold=True))
        parts.append(txt(cx, py + ph - r['I'] / vmax * ph - 3, '(z=%+.1f)' % r['z_clifford'], 7.6, col, anchor='middle'))
        parts.append(txt2(cx, py + ph + 12, lab_lines))
    parts.append(line(px, py + ph, px + pw, py + ph, '#bbb'))
    parts.append(txt(px + pw / 2, py + ph + 32,
                     'I = 0,649 sul log loss ratio: i comuni vicini condividono lo stesso errore di tariffazione',
                     7.8, MUT, anchor='middle'))
    return parts


def main():
    ep = json.load(open(EP_JSON, encoding='utf-8'))
    chi = json.load(open(CHI_JSON, encoding='utf-8'))
    tess = json.load(open(TESS_JSON, encoding='utf-8'))
    rob = json.load(open(ROB_JSON, encoding='utf-8'))
    spa = json.load(open(SPA_JSON, encoding='utf-8'))

    parts = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d" '
             'font-family="Helvetica,Arial,sans-serif">' % (W, H, W, H),
             '<rect width="%d" height="%d" fill="#fafaf7"/>' % (W, H),
             txt(W / 2, 26, 'Pricing Cat-Nat: quadro grafico dei risultati', 16, INK, 'middle', bold=True),
             txt(W / 2, 44, 'benchmark EAL vs tariffe IVASS · peso della loss sul tessuto · chi paga · robustezza · spazializzazione',
                 10, MUT, 'middle')]

    x1, x2, x3 = MX, MX + PW + GX, MX + 2 * (PW + GX)
    parts += p_ep(x1, TOP, ep)
    parts += p_quartili(x2, TOP, tess)
    parts += p_chi(x3, TOP, chi)
    parts += p_province(x1, ROW2, tess)
    parts += p_robustezza(x2, ROW2, rob)
    parts += p_moran(x3, ROW2, spa)

    parts.append(txt(W / 2, H - 10,
                     'Fonti: IVASS, INGV MPS04 (ag RP30/72/475), ISPRA (PAI), AIDA · valori letti da results/*.json · '
                     'generato da scripts/pricing/grafici.py · dettagli in docs/pricing_coerenza.md',
                     8.5, MUT, 'middle'))
    parts.append('</svg>')
    with open(OUT_SVG, 'w', encoding='utf-8') as f:
        f.write('\n'.join(parts))
    print('OK: %s (6 pannelli)' % os.path.relpath(OUT_SVG, ROOT))


if __name__ == '__main__':
    main()
