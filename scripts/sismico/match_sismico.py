#!/usr/bin/env python3
"""Matching dei centroidi comunali sulle griglie di pericolosita' sismica MPS04 (INGV).

Griglie richieste (vedi docs/sismico_metodologia.md per gli URL di download e la
conversione dei formati; non archiviate nel repo per dimensioni):

- ag RP475 (10%/50 anni): griglia testo italia_ag_002 (passo 0.02 gradi, 104.565 punti,
  colonne: id lon lat ag 16perc 84perc);
- ag RP30 (81% in 50 anni) e ag RP~50 (63% in 50 anni): CSV (ID,Lon,Lat,ag,84perc,16perc;
  16.852 punti) convertiti dagli .xls/.xlsx INGV (il nome file ag_63_RP72.csv e' la
  denominazione d'archivio della griglia 63% in 50 anni);
- Sa(T=0.10s) 50o percentile, RP 475/1000/2500: CSV (ID,Lon,Lat,SA_0.10,...,SA_2.00)
  estratti dagli .xls SA_*.xls INGV (parser BIFF8 incluso: parse_biff.py + recover_sa.py).

Sardegna (COD_REG 20): esclusa dalla classificazione sismica OPCM 3274/2003 ->
ag=0 e flag sismico_non_classificato=1 (evita match spuri a centinaia di km).
Match respinto anche oltre 5 km dal punto griglia piu' vicino (controllo di qualita').

Output: data/sismico/matrice_sismica_ingv.csv (una riga per comune).
Percorsi configurabili via variabili d'ambiente.
"""
import csv, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))

MATRICE = os.environ.get('MATRICE', os.path.join(ROOT, 'data', 'Matrice_Modello_Savelli_Final.csv'))
AG10_TXT = os.environ.get('AG10_TXT', os.path.join(ROOT, 'data', 'sismico', 'griglie', 'italia_ag_002_txt.txt'))
AG81_CSV = os.environ.get('AG81_CSV', os.path.join(ROOT, 'data', 'sismico', 'griglie', 'ag_81_RP30.csv'))
AG63_CSV = os.environ.get('AG63_CSV', os.path.join(ROOT, 'data', 'sismico', 'griglie', 'ag_63_RP72.csv'))
SA475 = os.environ.get('SA475', os.path.join(ROOT, 'data', 'sismico', 'griglie', 'sa0475_50perc.csv'))
SA1000 = os.environ.get('SA1000', os.path.join(ROOT, 'data', 'sismico', 'griglie', 'sa1000_50perc.csv'))
SA2500 = os.environ.get('SA2500', os.path.join(ROOT, 'data', 'sismico', 'griglie', 'sa2500_50perc.csv'))
OUT = os.environ.get('OUT', os.path.join(ROOT, 'data', 'sismico', 'matrice_sismica_ingv.csv'))

KM = 111.32          # km per grado (approssimazione lat)
MAX_DIST_KM = 5.0    # oltre questa distanza il match e' respinto (non classificato)
CELL = 0.2           # celle dello hash spaziale in gradi


def load_grid_txt(path):
    pts = []
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) == 6 and p[0] != 'id':
                pts.append((float(p[1]), float(p[2]), float(p[3]), float(p[4]), float(p[5])))
    return pts


def load_grid_csv(path):
    pts = []
    with open(path, newline='') as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            if any(x == '' for x in row[:6]):
                continue
            pts.append((float(row[1]), float(row[2]), float(row[3]), float(row[4]), float(row[5])))
    return pts


def load_sa(path):
    """SA_*.csv 50 percentile: ID,Lon,Lat,SA_0.10,... -> (lon, lat, sa01, 0, 0)"""
    pts = []
    with open(path, newline='') as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            if any(x == '' for x in row[:4]):
                continue
            pts.append((float(row[1]), float(row[2]), float(row[3]), 0.0, 0.0))
    return pts


def build_hash(pts):
    h = {}
    for i, p in enumerate(pts):
        key = (int(math.floor(p[0] / CELL)), int(math.floor(p[1] / CELL)))
        h.setdefault(key, []).append(i)
    return h


def nearest(lon, lat, pts, h):
    cx, cy = int(math.floor(lon / CELL)), int(math.floor(lat / CELL))
    best, bd = None, 1e18
    for dx in (-2, -1, 0, 1, 2):
        for dy in (-2, -1, 0, 1, 2):
            for i in h.get((cx + dx, cy + dy), []):
                p = pts[i]
                d2 = (p[0] - lon) ** 2 + (p[1] - lat) ** 2
                if d2 < bd:
                    bd, best = d2, i
    if best is None:  # celle vuote: allarga la ricerca ad anelli quadrati successivi
        for r in range(3, 30):
            # righe superiore e inferiore dell'anello (colonne complete)
            for dx in range(-r, r + 1):
                for dy in (-r, r):
                    for i in h.get((cx + dx, cy + dy), []):
                        p = pts[i]
                        d2 = (p[0] - lon) ** 2 + (p[1] - lat) ** 2
                        if d2 < bd:
                            bd, best = d2, i
            # colonne sinistra e destra dell'anello (righe interne, angoli gia' coperti sopra)
            for dx2 in (-r, r):
                for dy2 in range(-r + 1, r):
                    for i in h.get((cx + dx2, cy + dy2), []):
                        p = pts[i]
                        d2 = (p[0] - lon) ** 2 + (p[1] - lat) ** 2
                        if d2 < bd:
                            bd, best = d2, i
            # ogni punto degli anelli successivi dista almeno (r-1)*CELL gradi:
            # se tale bound supera la distanza del migliore trovato, la risposta e' garantita
            if best is not None and (r - 1) * CELL >= math.sqrt(bd):
                break
    if best is None:
        return None, None, None
    return pts[best], math.sqrt(bd), best


def main():
    comuni = []
    with open(MATRICE, newline='', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            if not row.get('PRO_COM'):
                continue
            try:
                lat, lon = float(row['lat']), float(row['long'])
            except (TypeError, ValueError):
                continue
            comuni.append({'PRO_COM': row['PRO_COM'], 'COMUNE': row.get('COMUNE', ''),
                           'Provincia': row.get('Provincia', ''), 'COD_REG': row.get('COD_REG', ''),
                           'lat': lat, 'lon': lon})
    print(f'comuni con coordinate valide: {len(comuni)}', file=sys.stderr)

    g10 = load_grid_txt(AG10_TXT)
    g81 = load_grid_csv(AG81_CSV)
    g63 = load_grid_csv(AG63_CSV)
    gsa475 = load_sa(SA475)
    gsa1000 = load_sa(SA1000)
    gsa2500 = load_sa(SA2500)
    print(f'griglia RP475: {len(g10)} punti | RP30: {len(g81)} | RP50: {len(g63)} | '
          f'SA475: {len(gsa475)} | SA1000: {len(gsa1000)} | SA2500: {len(gsa2500)}', file=sys.stderr)

    h10, h81, h63 = build_hash(g10), build_hash(g81), build_hash(g63)
    hsa475, hsa1000, hsa2500 = build_hash(gsa475), build_hash(gsa1000), build_hash(gsa2500)

    out_rows = []
    n_sardi, n_respinti = 0, 0
    for c in comuni:
        sardegna = (c['COD_REG'] == '20')
        p10, d10, _ = nearest(c['lon'], c['lat'], g10, h10)
        p81, d81, _ = nearest(c['lon'], c['lat'], g81, h81)
        p63, d63, _ = nearest(c['lon'], c['lat'], g63, h63)
        s475, d475, _ = nearest(c['lon'], c['lat'], gsa475, hsa475)
        s1000, _, _ = nearest(c['lon'], c['lat'], gsa1000, hsa1000)
        s2500, _, _ = nearest(c['lon'], c['lat'], gsa2500, hsa2500)
        if p10 is None or p81 is None or p63 is None:
            print('MANCA:', c['PRO_COM'], c['COMUNE'], file=sys.stderr)
            continue
        respinto = (not sardegna) and (d10 * KM > MAX_DIST_KM)
        if sardegna:
            n_sardi += 1
        elif respinto:
            n_respinti += 1
        zero = (sardegna or respinto)
        out_rows.append({
            'PRO_COM': c['PRO_COM'], 'COMUNE': c['COMUNE'], 'Provincia': c['Provincia'],
            'COD_REG': c['COD_REG'], 'lat': c['lat'], 'long': c['lon'],
            'ag_RP475': 0.0 if zero else p10[2],
            'ag_RP475_16perc': 0.0 if zero else p10[3],
            'ag_RP475_84perc': 0.0 if zero else p10[4],
            'ag_RP30': 0.0 if zero else p81[2],
            'ag_RP30_84perc': 0.0 if zero else p81[3],
            'ag_RP30_16perc': 0.0 if zero else p81[4],
            'ag_RP72': 0.0 if zero else p63[2],
            'ag_RP72_84perc': 0.0 if zero else p63[3],
            'ag_RP72_16perc': 0.0 if zero else p63[4],
            'Sa01_RP475': 0.0 if zero else s475[2],
            'Sa01_RP1000': 0.0 if zero else s1000[2],
            'Sa01_RP2500': 0.0 if zero else s2500[2],
            'dist_grid_km_475': d10 * KM, 'dist_grid_km_81': d81 * KM, 'dist_grid_km_63': d63 * KM,
            'sismico_non_classificato': 1 if zero else 0,
        })

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)
    print(f'scritte {len(out_rows)} righe in {OUT} | sardi a zero: {n_sardi} | '
          f'respinti per distanza: {n_respinti}', file=sys.stderr)
    return out_rows


if __name__ == '__main__':
    main()
