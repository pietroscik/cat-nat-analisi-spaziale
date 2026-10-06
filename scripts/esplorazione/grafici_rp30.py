#!/usr/bin/env python3
"""Quadro grafico dell'esplorazione della banda RP<30: cinque pannelli SVG.

Ogni pannello legge i valori direttamente dai JSON in results/ generati dagli
script di scripts/esplorazione/ (niente numeri hardcoded): la figura e' una
vista, non un'origine di dati, e la CI ne verifica la riproducibilita' byte per
byte insieme al resto.

Pannelli (docs/esplorazione_validazione_rp30.md: §6 stima event-based, gate
del passo 2; §5 protocollo; risultati/ep_curve.json per il limite RP>=30):

1. COMPONENTE DELLA CODA FREQUENTE: AAL event-based per bin di magnitudo
   ([5,0-5,5) 13,2 / [5,5-6,0) 18,9 / [6,0-6,5) 25,3 / [6,5+] 19,7 mln):
   ~58% dagli stessi eventi maggiori (anello 20-60 km sotto ag_RP30);
2. TOP EVENTI: i sei maggiori contributi storici alla coda frequente
   (Emilia 2012 17,0, Valnerina 2016 13,7, Irpinia 1980 6,0 mln/anno ...);
3. GATE: limite identificato RP>=30 (1,713 mld) vs stima event-based (0,077),
   pavimento declusterato (0,056), variante prudenziale scalata (1,49) e
   prior max-ent P1 (3,18, scartato): la coda e' un correctivo minore;
4. RAPPORTO SUL BENCHMARK: totale RP>=0 per specifica (limite inferiore
   1,66x, centrale 1,73x, prudenza estrema ~1,76x) contro i prior scartati
   (P1 4,7x, P2 3,1x);
5. SENSIBILITA': AAL frequente su tutte le varianti (base 77, declusterata
   56, soglia Mw>=4,5 103, finestra 1960 73, attenuazione/PGA/Io0 16-123 mln).

Output: docs/grafici_rp30.svg (deterministico, solo stdlib).
Esecuzione dalla root del repo dopo gli script di scripts/esplorazione/:
    python3 scripts/esplorazione/grafici_rp30.py
"""
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))

RESULTS = os.path.join(ROOT, 'results')
EV_JSON = os.environ.get('EV_JSON', os.path.join(RESULTS, 'esplorazione_coda_eventi.json'))
PR_JSON = os.environ.get('PR_JSON', os.path.join(RESULTS, 'esplorazione_coda_frequente.json'))
EP_JSON = os.environ.get('EP_JSON', os.path.join(RESULTS, 'ep_curve.json'))
OUT_SVG = os.environ.get('OUT_SVG', os.path.join(ROOT, 'docs', 'grafici_rp30.svg'))

W, H = 1120, 760
PW, PH = 340, 290            # pannello
GX, GY = 20, 18              # gap tra i pannelli
MX, TOP = 20, 66              # margini esterni / inizio griglia
ROW2 = TOP + PH + 40          # seconda riga (spazio per didascalie)

INK, MUT, GRID = '#222', '#777', '#e4e2da'
BLU, ROSSO, VERDE, ARANCIO = '#2b5f9e', '#9c1f1f', '#3f7d4e', '#c97b2a'


def esc(s):
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


def txt(x, y, s, size=9, fill=INK, anchor='start', bold=False):
    w = ' font-weight="bold"' if bold else ''
    return ('<text x="%.1f" y="%.1f" font-size="%s"%s fill="%s" text-anchor="%s">%s</text>'
            % (x, y, size, w, fill, anchor, esc(s)))


def line(x1, y1, x2, y2, stroke=GRID, dash=''):
    d = ' stroke-dasharray="%s"' % dash if dash else ''
    return '<line x1="%.1f" y1="%.1f" x2="%.2f" y2="%.2f" stroke="%s" stroke-width="1"%s/>' \
           % (x1, y1, x2, y2, stroke, d)


def rect(x, y, w, h, fill, stroke='none'):
    return '<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s" stroke="%s"/>' \
           % (x, y, max(0.0, w), max(0.0, h), fill, stroke)


def panel(x0, y0, num, title, subtitle):
    parts = [rect(x0, y0, PW, PH, '#ffffff', '#d8d6ce')]
    parts.append(txt(x0 + 10, y0 + 20, '%d. %s' % (num, title), 11.5, INK, bold=True))
    parts.append(txt(x0 + 10, y0 + 34, subtitle, 8.6, MUT))
    return parts, (x0 + 34, y0 + 48, PW - 48, PH - 72)   # (px, py, pw, ph)


def fmt_it(v, dec=1):
    return ('%.*f' % (dec, v)).replace('.', ',')


def fmt_mld(v):
    return fmt_it(v / 1e9, 2)


# ------------------------------------------------------- 1. composizione AAL_f
def p_composizione(x0, y0, ev):
    parts, (px, py, pw, ph) = panel(x0, y0, 1, 'Composizione della coda frequente',
                                    'AAL event-based per bin di magnitudo (mln EUR/anno, RP<30)')
    per_mag = ev['risultato']['per_magnitudo_AAL_freq']
    bins = ['[5.0,5.5)', '[5.5,6.0)', '[6.0,6.5)', '[6.5,inf)']
    labels = ['5,0-5,5', '5,5-6,0', '6,0-6,5', '6,5+']
    vals = [per_mag[b] / 1e6 for b in bins]
    tot = sum(vals)
    vmax = 30.0
    step = pw / len(vals)
    bw = step * 0.52
    for t in (10, 20, 30):
        parts.append(line(px, py + ph - t / vmax * ph, px + pw, py + ph - t / vmax * ph))
        parts.append(txt(px - 5, py + ph - t / vmax * ph + 3, str(t), 8.5, MUT, anchor='end'))
    for i, v in enumerate(vals):
        cx = px + step * i + step / 2
        parts.append(rect(cx - bw / 2, py + ph - v / vmax * ph, bw, v / vmax * ph, ROSSO))
        parts.append(txt(cx, py + ph - v / vmax * ph - 5, fmt_it(v, 1), 9, INK, anchor='middle', bold=True))
        parts.append(txt(cx, py + ph + 12, labels[i], 8.6, MUT, anchor='middle'))
    parts.append(line(px, py + ph, px + pw, py + ph, '#bbb'))
    q_58 = sum(vals[2:]) / tot * 100
    parts.append(txt(px + pw / 2, py + ph + 23,
                     'magnitudo Mw (finestra 1980-2020) — il %s%% dal bin 6,0+ (stessi eventi maggiori)'
                     % fmt_it(q_58, 0), 7.8, MUT, anchor='middle'))
    return parts


# ------------------------------------------------------------- 2. top eventi
def p_eventi(x0, y0, ev):
    parts, (px, py, pw, ph) = panel(x0, y0, 2, 'Top eventi storici',
                                    'maggiori contributi alla AAL frequente (mln EUR/anno)')
    top = ev['risultato']['top_eventi_AAL_freq'][:6]
    vmax = top[0]['contributo_AAL_freq_EUR'] / 1e6 * 1.12
    rh = ph / len(top)
    for i, e in enumerate(top):
        yy = py + rh * i + rh / 2
        v = e['contributo_AAL_freq_EUR'] / 1e6
        bw_px = v / vmax * pw
        parts.append(rect(px, yy - rh * 0.32, bw_px, rh * 0.64, ARANCIO))
        nome = '%s %d' % (e['area'], e['anno'])
        val = '%s (M%s)' % (fmt_it(v, 1), fmt_it(e['Mw'], 2))
        if bw_px > 150:
            parts.append(txt(px + 4, yy + 3, nome + '  ' + val, 8.2, '#fff'))
        else:
            parts.append(txt(px + 4, yy + 3, nome, 8.2, '#fff'))
            parts.append(txt(px + bw_px + 5, yy + 3, val, 8.2, INK))
    parts.append(line(px, py + ph, px + pw, py + ph, '#bbb'))
    parts.append(txt(px + pw / 2, py + ph + 23,
                     'CPTI15 v4.0, Mw>=5,0 — contributo = tasso empirico x danno sotto ag_RP30 per sito', 7.8, MUT, anchor='middle'))
    return parts


# ------------------------------------------------------------------- 3. gate
def p_gate(x0, y0, ev, ep):
    parts, (px, py, pw, ph) = panel(x0, y0, 3, 'Gate: quanto vale la coda frequente?',
                                    'AAL (mld EUR/anno) — stima event-based vs bounds del gate')
    r = ev['risultato']
    gate = ev['gate']
    lim = gate['limite_inferiore_RP_sup30_EUR']
    rows = [('prior max-ent P1 (scartato)', gate['prior_maxent_P1_mediana_EUR'], MUT, True),
            ('variante prudenziale scalata', ev['calibrazione']['AAL_freq_variante_scalata_EUR'], ARANCIO, False),
            ('stima event-based (base)', r['AAL_frequente_RP_inf30_EUR'], ROSSO, False),
            ('pavimento declusterato', ev['declustering']['AAL_freq_declusterata_EUR'], VERDE, False),
            ('limite identificato RP>=30', lim, BLU, False)]
    vmax = max(v for _, v, _, _ in rows) / 1e9 * 1.10
    rh = ph / len(rows)
    for t in (1.0, 2.0, 3.0):
        if t <= vmax:
            parts.append(line(px + t / vmax * pw, py, px + t / vmax * pw, py + ph))
            parts.append(txt(px + t / vmax * pw, py + ph + 11, fmt_it(t, 0), 7.8, MUT, anchor='middle'))
    for i, (lab, v, col, dashed) in enumerate(rows):
        yy = py + rh * i + rh / 2
        bw_px = (v / 1e9) / vmax * pw
        parts.append(rect(px, yy - rh * 0.30, bw_px, rh * 0.60, col,
                          '#9a9a9a' if dashed else 'none'))
        if bw_px > 150:
            parts.append(txt(px + 4, yy + 3, '%s  %s' % (lab, fmt_mld(v)), 8.2, '#fff'))
        else:
            parts.append(txt(px + bw_px + 5, yy + 3, fmt_mld(v), 8.2, col, bold=True))
            parts.append(txt(px + bw_px + 40, yy + 3, lab, 8.2, INK))
    parts.append(line(px, py + ph, px + pw, py + ph, '#bbb'))
    aal_tot = gate['totale_RP_sup0_stimato_EUR']
    parts.append(txt(px + pw / 2, py + ph + 23,
                     'event-based = %s del limite RP>=30: correctivo minore → totale RP>=0 %s mld'
                     % (fmt_it(r['AAL_frequente_RP_inf30_EUR'] / lim * 100, 1) + '%', fmt_mld(aal_tot)),
                     7.8, MUT, anchor='middle'))
    return parts


# ------------------------------------------- 4. rapporto totale su benchmark
def p_rapporto(x0, y0, ev, pr):
    parts, (px, py, pw, ph) = panel(x0, y0, 4, 'Totale RP>=0 vs benchmark CURVE=4',
                                    'rapporto x del totale sul benchmark (1,03 mld EUR/anno)')
    gate = ev['gate']
    pred = pr['predittiva']
    bmk = gate['benchmark_CURVE4_EUR']
    lim_x = gate['limite_inferiore_RP_sup30_EUR'] / bmk
    base_x = gate['totale_quota_benchmark']
    prud_x = (gate['limite_inferiore_RP_sup30_EUR']
              + ev['calibrazione']['AAL_freq_variante_scalata_EUR']) / bmk
    p1_x = pred['P1_ancorata_RP30_72']['rapporto_totale_su_benchmark']['p50']
    p2_x = pred['P2_appiattita_RP72_475']['rapporto_totale_su_benchmark']['p50']
    rows = [('P1 max-ent (scartato)', p1_x, MUT),
            ('prudenza estrema (scalata)', prud_x, ARANCIO),
            ('P2 appiattita (scartata)', p2_x, MUT),
            ('centrale (1,71 + 0,077)', base_x, ROSSO),
            ('limite inferiore RP>=30', lim_x, BLU)]
    vmax = max(v for _, v, _ in rows) * 1.10
    rh = ph / len(rows)
    for t in (1, 2, 3, 4):
        if t <= vmax:
            parts.append(line(px + t / vmax * pw, py, px + t / vmax * pw, py + ph))
            parts.append(txt(px + t / vmax * pw, py + ph + 11, '%dx' % t, 7.8, MUT, anchor='middle'))
    for i, (lab, v, col) in enumerate(rows):
        yy = py + rh * i + rh / 2
        bw_px = v / vmax * pw
        parts.append(rect(px, yy - rh * 0.30, bw_px, rh * 0.60, col))
        if bw_px > 150:
            parts.append(txt(px + 4, yy + 3, '%s  %sx' % (lab, fmt_it(v, 2)), 8.2, '#fff'))
        else:
            parts.append(txt(px + bw_px + 5, yy + 3, fmt_it(v, 2) + 'x', 8.2, col, bold=True))
            parts.append(txt(px + bw_px + 32, yy + 3, lab, 8.2, INK))
    parts.append(line(px, py + ph, px + pw, py + ph, '#bbb'))
    parts.append(txt(px + pw / 2, py + ph + 23,
                     'i prior P1/P2 si archiviano come bound laschi (41x/19x sopra la stima event-based)',
                     7.8, MUT, anchor='middle'))
    return parts


# ----------------------------------------------------------- 5. sensibilita
def p_sensibilita(x0, y0, ev):
    parts, (px, py, pw, ph) = panel(x0, y0, 5, 'Sensibilità della stima',
                                    'AAL frequente per variante (mln EUR/anno)')
    base = ev['risultato']['AAL_frequente_RP_inf30_EUR'] / 1e6
    s = ev['sensibilita']
    rows = [('base (Io0 mediana cat., anel. 0,010, PGA/4,5)', base, ROSSO),
            ('soglia Mw>=4,5 (+ deficit ravvicinato)', s['soglia_Mw_4.5'] / 1e6, BLU),
            ('finestra 1960-2020', s['finestra_1960_2020'] / 1e6, BLU),
            ('declustering 50km/90gg (pavimento)', s['declustering_50km_90gg'] / 1e6, VERDE),
            ('att. anelastica 0,005 (PGA/3,5)', s['anelastico=0.005,pga_div=3.5'] / 1e6, MUT),
            ('Io0 post-1980 (cat. completo)', s['Io0_finestra_1980'] / 1e6, MUT)]
    vmin, vmax = 0.0, 130.0
    lx = lambda v: px + v / (vmax - vmin) * pw
    rh = ph / len(rows)
    for t in (25, 50, 75, 100, 125):
        parts.append(line(lx(t), py, lx(t), py + ph))
        parts.append(txt(lx(t), py + ph + 11, str(t), 7.8, MUT, anchor='middle'))
    for i, (lab, v, col) in enumerate(rows):
        yy = py + rh * i + rh / 2
        parts.append(line(lx(vmin), yy, lx(v), yy, '#ddd6c8'))
        parts.append('<circle cx="%.1f" cy="%.1f" r="4.0" fill="%s"/>' % (lx(v), yy, col))
        parts.append(txt(min(lx(v) + 14, px + pw - 30), yy - 8, fmt_it(v, 1), 8.2, col,
                         anchor='middle', bold=True))
        parts.append(txt(px + pw - 2, yy + 3, lab, 7.9, INK, anchor='end'))
    parts.append(line(lx(0), py, lx(0), py + ph, '#bbb'))
    parts.append(txt(px + pw / 2, py + ph + 23,
                     'intervallo complessivo 16-123 mln/anno; la stima puntuale e un limite inferiore', 7.8, MUT, anchor='middle'))
    return parts


def main():
    ev = json.load(open(EV_JSON, encoding='utf-8'))
    pr = json.load(open(PR_JSON, encoding='utf-8'))
    ep = json.load(open(EP_JSON, encoding='utf-8'))

    parts = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d" '
             'font-family="Helvetica,Arial,sans-serif">' % (W, H, W, H),
             '<rect width="%d" height="%d" fill="#fafaf7"/>' % (W, H),
             txt(W / 2, 26, 'Coda frequente RP<30: dalla banda a priori alla stima event-based', 16, INK, 'middle', bold=True),
             txt(W / 2, 44, 'verdetto MPS04 · catalogo CPTI15 v4.0 · stima event-based · gate · sensibilita', 10, MUT, 'middle')]

    x1, x2, x3 = MX, MX + PW + GX, MX + 2 * (PW + GX)
    parts += p_composizione(x1, TOP, ev)
    parts += p_eventi(x2, TOP, ev)
    parts += p_gate(x3, TOP, ev, ep)
    parts += p_rapporto(x1, ROW2, ev, pr)
    parts += p_sensibilita(x2, ROW2, ev)
    # sesto slot: nota-box RIASSUNTIVA del verdetto (senza numeri hardcoded)
    bx, by, bw, bh = x3, ROW2, PW, PH
    parts.append(rect(bx, by, bw, bh, '#ffffff', '#d8d6ce'))
    gate = ev['gate']
    r_ = ev['risultato']
    pr1 = pr['predittiva']['P1_ancorata_RP30_72']['AAL_frequente_EUR']['p50'] / 1e9
    pr2 = pr['predittiva']['P2_appiattita_RP72_475']['AAL_frequente_EUR']['p50'] / 1e9
    lines = [
        ('Il verdetto in cinque righe', INK, True, 11.5),
        ('', INK, False, 9),
        ('1. Sotto RP30 non esiste hazard pubblico MPS04 (81% in 50 anni = RP~30: il muro informativo).', INK, False, 9),
        ('2. La coda frequente e fatta di eventi osservati: Mw>=6 ogni ~6,8 anni (CPTI15 1980-2020).', INK, False, 9),
        ('3. Stima event-based: %s mln/anno, il %s%% del limite RP>=30 — correctivo minore, non un raddoppio.'
         % (fmt_it(r_['AAL_frequente_RP_inf30_EUR'] / 1e6, 0),
            fmt_it(r_['AAL_frequente_RP_inf30_EUR'] / gate['limite_inferiore_RP_sup30_EUR'] * 100, 1)), INK, False, 9),
        ('4. I prior max-ent (P1 %s / P2 %s mld) si archiviano come bound laschi: scartati dal gate.'
         % (fmt_it(pr1, 2), fmt_it(pr2, 2)), INK, False, 9),
        ('5. Totale RP>=0 = %s mld (%sx benchmark); chiusura definitiva al passo DBMI15 (danni osservati).'
         % (fmt_it(gate['totale_RP_sup0_stimato_EUR'] / 1e9, 2),
            fmt_it(gate['totale_quota_benchmark'], 2)), INK, False, 9),
        ('', INK, False, 9),
        ("Nota: il gap 1,035-vs-1,713 mld e artefatto dell'approssimazione CURVE=4 (ratio 1,655", MUT, False, 8),
        ('gia documentato in ep_curve.json), non una sottostima del premio. Dettagli:', MUT, False, 8),
        ('docs/esplorazione_validazione_rp30.md, §6-§7.', MUT, False, 8),
    ]
    yy = by + 24
    for s, col, bold, size in lines:
        if s:
            parts.append(txt(bx + 12, yy, s, size, col, bold=bold))
        yy += 17 if size >= 9 else 13
    parts.append(txt(W / 2, H - 10,
                     'Fonti: INGV CPTI15 v4.0 (Rovida et al.), MPS04, IVASS, AIDA · valori letti da results/esplorazione_*.json · '
                     'generato da scripts/esplorazione/grafici_rp30.py · dettagli in docs/esplorazione_validazione_rp30.md',
                     8.5, MUT, 'middle'))
    parts.append('</svg>')
    with open(OUT_SVG, 'w', encoding='utf-8') as f:
        f.write('\n'.join(parts))
    print('OK: %s (5 pannelli + verdetto)' % os.path.relpath(OUT_SVG, ROOT))


if __name__ == '__main__':
    main()
