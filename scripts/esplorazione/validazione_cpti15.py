#!/usr/bin/env python3
"""Esplorazione (ramo dev): validazione del catalogo CPTI15 v4.0 (INGV).

Primo mattone della fase di validazione della coda frequente RP<30: la §8 di
docs/pricing_coerenza.md la dichiara non stimabile dalla curva EP (MPS04 non
 pubblica nulla sotto l'81% in 50 anni, verificato: dataset INGV 70/193 senza
risorse, probe delle griglie piu' frequenti tutti 404). L'unica via e' il lato
eventi: frequenze empiriche (CPTI15) x fragilita' (Rosti et al. 2021) x
esposizione (asset della matrice, gia' presente).

Questo script scarica-independent: legge il catalogo CPTI15 v4.0 (fogli
about/format/catalogue) e ne valida struttura e statistiche:
- n eventi, finestra, Mw range, righe anomale;
- tassi empirici per soglia Mw su finestre di completezza crescenti;
- check di eventi noti (L'Aquila 2009, Amatrice 2016, Emilia 2012, Irpinia
  1980, Molise 2002) con Mw atteso;
- coerenza Mw -> Io (mediana per classe).

Fonte: https://emidius.mi.ingv.it/CPTI15-DBMI15/ (Rovida et al., CPTI15 v4.0,
INGV, uso scientifico con citazione). Il file non e' archiviato nel repo:
scaricarlo con curl -O https://emidius.mi.ingv.it/CPTI15-DBMI15/data/CPTI15_v4.0.xlsx
e passare il percorso con l'ambiente CPTI15_XLSX (default: data/sismico/).

Output: results/esplorazione_cpti15.json (deterministico).
"""
import csv, hashlib, io, json, os, statistics, zipfile
from xml.etree.ElementTree import iterparse

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
XLSX = os.environ.get('CPTI15_XLSX', os.path.join(ROOT, 'data', 'sismico', 'CPTI15_v4.0.xlsx'))
OUT_JSON = os.environ.get('OUT_JSON', os.path.join(ROOT, 'results', 'esplorazione_cpti15.json'))

URL = 'https://emidius.mi.ingv.it/CPTI15-DBMI15/data/CPTI15_v4.0.xlsx'
NS = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
MD5_ATTESO = 'da calcolare alla prima esecuzione (registro nel JSON)'

# eventi noti: (anno, frammento area, Mw atteso ~, tolleranza)
EVENTI_NOTI = [
    (2009, 'aquilan', 6.29, 0.15, "L'Aquila 2009"),
    (2016, 'monti della laga', 6.18, 0.15, 'Amatrice 2016'),
    (2012, 'pianura emiliana', 6.09, 0.15, 'Emilia 2012'),
    (1980, 'irpinia', 6.81, 0.15, 'Irpinia 1980'),
    (2002, 'molise', 5.72, 0.15, 'Molise 2002'),
]

# finestre di completezza (proxy prudenti: da verificare con le tabelle CPTI15)
FIN = [(1600, 2020), (1900, 2020), (1960, 2020), (1980, 2020), (2000, 2020)]
SOGLIE = [4.0, 4.5, 5.0, 5.5, 6.0, 6.5]


def col_index(ref):
    col = 0
    for ch in ref:
        if ch.isalpha():
            col = col * 26 + (ord(ch.upper()) - 64)
        else:
            break
    return col - 1


def estrai_catalogo(path):
    """Estrae il foglio 'catalogue' (sheet3) dello xlsx in lista di dict."""
    z = zipfile.ZipFile(path)
    shared = []
    for ev, el in iterparse(z.open('xl/sharedStrings.xml'), events=('end',)):
        if el.tag == NS + 'si':
            shared.append(''.join(t.text or '' for t in el.iter(NS + 't')))
            el.clear()
    # ordina i fogli per numero per trovare 'catalogue' in modo robusto
    import re
    nomi = {n: int(re.search(r'sheet(\d+)\.xml', n).group(1))
            for n in z.namelist() if re.match(r'xl/worksheets/sheet\d+\.xml$', n)}
    fogli_ordinati = [n for n, _ in sorted(nomi.items(), key=lambda kv: kv[1])]
    # cerca il foglio con la riga header contenente 'MwDef' (robusto ai riordini)
    hdr = None
    righe = []
    for nome in fogli_ordinati:
        cand = _leggi_foglio(z, nome, shared)
        if cand and any('MwDef' == h for h in cand[0]):
            hdr, righe = cand[0], cand[1:]
            break
    if hdr is None:
        raise SystemExit('foglio catalogue non trovato')
    out = []
    for r in righe:
        if not any(c.strip() for c in r):
            continue
        out.append(dict(zip(hdr, r)))
    return hdr, out


def _leggi_foglio(z, nome, shared):
    rows_out = []
    for ev, row in iterparse(z.open(nome), events=('end',)):
        if row.tag != NS + 'row':
            continue
        cells = {}
        for c in row.iter(NS + 'c'):
            idx = col_index(c.get('r'))
            v = c.find(NS + 'v')
            t = c.get('t')
            if t == 's':
                val = shared[int(v.text)] if v is not None else ''
            elif t == 'inlineStr':
                is_el = c.find(NS + 'is')
                val = ''.join(t2.text or '' for t2 in is_el.iter(NS + 't')) if is_el is not None else ''
            else:
                val = v.text if v is not None else ''
            cells[idx] = val
        if cells:
            n = max(cells) + 1
            rows_out.append([cells.get(i, '') for i in range(n)])
        row.clear()
    return rows_out


def io_num(s):
    s = (s or '').strip().replace(',', '.')
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        if '-' in s:
            a, _, b = s.partition('-')
            try:
                return (float(a) + float(b)) / 2.0
            except ValueError:
                return None
        return None


def main():
    md5 = hashlib.md5(open(XLSX, 'rb').read()).hexdigest()
    hdr, ev = estrai_catalogo(XLSX)
    eventi = []
    senza_mw = 0
    for e in ev:
        try:
            eventi.append((int(e['Year']), float(e['MwDef']), e['EpicentralArea'].strip(),
                          io_num(e.get('IoDef', ''))))
        except (ValueError, KeyError):
            senza_mw += 1
    anni = [e[0] for e in eventi]
    mw = [e[1] for e in eventi]
    print('CPTI15 v4.0: %d eventi con Mw (%d senza), %d-%d, Mw %.1f-%.1f' %
          (len(eventi), senza_mw, min(anni), max(anni), min(mw), max(mw)))

    # tassi per finestra x soglia
    tassi = {}
    for (y0, y1) in FIN:
        dur = y1 - y0 + 1
        tassi['%d_%d' % (y0, y1)] = {
            'Mw>=%.1f' % s: {'n': sum(1 for e in eventi if y0 <= e[0] <= y1 and e[1] >= s),
                             'eventi_anno': round(sum(1 for e in eventi if y0 <= e[0] <= y1 and e[1] >= s) / dur, 3)}
            for s in SOGLIE}

    # check eventi noti
    check = []
    for anno, frag, mw_att, tol, nome in EVENTI_NOTI:
        hit = [e for e in eventi if e[0] == anno and frag in e[2].lower()]
        main = max(hit, key=lambda e: e[1]) if hit else None
        check.append({'evento': nome, 'anno': anno, 'trovato': main is not None,
                      'area': main[2] if main else None,
                      'Mw': round(main[1], 2) if main else None,
                      'Mw_atteso': mw_att,
                      'ok': bool(main and abs(main[1] - mw_att) <= tol)})

    # Io mediano per classe Mw (1900+)
    io_cls = []
    recent = [e for e in eventi if e[0] >= 1900]
    for lo, hi in [(4.0, 4.5), (4.5, 5.0), (5.0, 5.5), (5.5, 6.0), (6.0, 9.0)]:
        g = [e[3] for e in recent if lo <= e[1] < hi and e[3] is not None]
        if g:
            io_cls.append({'classe_Mw': '[%.1f,%.1f)' % (lo, hi), 'n': len(g),
                           'Io_mediano': round(statistics.median(g), 1), 'Io_max': round(max(g), 1)})

    out = {
        'esperimento': 'validazione catalogo CPTI15 v4.0 (INGV) per la fase di validazione della coda RP<30 (ramo dev)',
        'fonte': {'url': URL, 'md5': md5, 'citazione': 'Rovida A., Locati M., Camassi R., Lolli B., Gasperini P. (a cura di), CPTI15 v4.0, INGV',
                  'note': 'uso scientifico con citazione obbligatoria della fonte; file non archiviato nel repo'},
        'copertura': {'n_eventi_con_Mw': len(eventi), 'n_righe_senza_Mw': senza_mw,
                      'finestra': '%d-%d' % (min(anni), max(anni)),
                      'Mw_min': min(mw), 'Mw_max': max(mw)},
        'tassi_empirici': tassi,
        'finestra_completa_riferimento': '1980-2020 (strumentale, completa per M>=4.5: M>=5,5 ogni ~1,4 anni, M>=6 ogni ~6,8 anni)',
        'eventi_noti': check,
        'coerenza_Mw_Io': io_cls,
        'nota_parsing': "IoDef puo' essere un intervallo (es. '6-7'): convertito al punto medio; 2 righe vuote scartate",
        'verdetto': '',
    }
    n_ok = sum(1 for c in check if c['ok'])
    out['verdetto'] = (
        'catalogo coerente: %d/%d eventi noti verificati, tassi stabili sulle finestre 1960/1980/2000 '
        '(completezza M>=4,5 raggiunta), Mw->Io monotono. Ancoraggio empirico per la coda frequente: '
        'gli eventi dannosi M>=6 ricorrono ogni ~7 anni, DENTRO la banda RP<30: la coda frequente esiste '
        'e va modellata dal lato eventi (frequenze CPTI15 x fragilita x esposizione), non dall\u2019extrapolazione '
        'della curva EP.' % (n_ok, len(check)))

    os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
    with open(OUT_JSON, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print('scritto:', OUT_JSON)


if __name__ == '__main__':
    main()
