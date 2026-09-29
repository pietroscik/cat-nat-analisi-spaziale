#!/usr/bin/env python3
"""Costruisce la matrice v2 con il terzo hazard sismico:
data/Matrice_Modello_Savelli_Final.csv + colonne sismiche INGV (MPS04).

Join su PRO_COM con data/sismico/matrice_sismica_ingv.csv (prodotta da match_sismico.py).
Colonne aggiunte:
- ag_RP475 (+16/84 perc), ag_RP30, ag_RP72, Sa01_RP475, Sa01_RP1000, Sa01_RP2500
- sismico_non_classificato (1 = Sardegna / match respinto, ag=0)
- hazard_sismico = ag_RP475 (accelerazione massima del suolo, in g)
- Risk_Sismico_Asset_PMI   = ag_RP475 * asset_PMI_EUR
- Risk_Sismico_Asset_Grandi = ag_RP475 * asset_grandi_EUR

Nota: per frana/idraulico l'hazard e' una quota di area comunale (share);
per il sismico ag e' un'accelerazione (g): i coefficienti delle due famiglie
di regressori non sono direttamente confrontabili in magnitudine (vedi
docs/sismico_metodologia.md).
"""
import csv, os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))

MATRICE = os.environ.get('MATRICE', os.path.join(ROOT, 'data', 'Matrice_Modello_Savelli_Final.csv'))
SISMICA = os.environ.get('SISMICA', os.path.join(ROOT, 'data', 'sismico', 'matrice_sismica_ingv.csv'))
OUT = os.environ.get('OUT', os.path.join(ROOT, 'data', 'Matrice_Modello_Savelli_Final_sismico.csv'))


def main():
    sism = {}
    with open(SISMICA, newline='', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            sism[row['PRO_COM']] = row
    print('righe sismiche:', len(sism))

    newcols = ['ag_RP475', 'ag_RP475_16perc', 'ag_RP475_84perc',
               'ag_RP30', 'ag_RP30_84perc', 'ag_RP30_16perc',
               'ag_RP72', 'ag_RP72_84perc', 'ag_RP72_16perc',
               'Sa01_RP475', 'Sa01_RP1000', 'Sa01_RP2500',
               'sismico_non_classificato', 'hazard_sismico',
               'Risk_Sismico_Asset_PMI', 'Risk_Sismico_Asset_Grandi']

    n_in = n_join = 0
    with open(MATRICE, newline='', encoding='utf-8') as fin, \
         open(OUT, 'w', newline='', encoding='utf-8') as fout:
        r = csv.DictReader(fin)
        w = csv.DictWriter(fout, fieldnames=r.fieldnames + newcols, extrasaction='raise')
        w.writeheader()
        for row in r:
            n_in += 1
            pc = row['PRO_COM'].strip()
            s = sism.get(pc)
            if s is None:
                raise SystemExit(f'PRO_COM senza match sismico: {pc} {row.get("COMUNE")}')
            n_join += 1
            ag = float(s['ag_RP475'])
            row['ag_RP475'] = ag
            row['ag_RP475_16perc'] = float(s['ag_RP475_16perc'])
            row['ag_RP475_84perc'] = float(s['ag_RP475_84perc'])
            row['ag_RP30'] = float(s['ag_RP30'])
            row['ag_RP30_84perc'] = float(s['ag_RP30_84perc'])
            row['ag_RP30_16perc'] = float(s['ag_RP30_16perc'])
            row['ag_RP72'] = float(s['ag_RP72'])
            row['ag_RP72_84perc'] = float(s['ag_RP72_84perc'])
            row['ag_RP72_16perc'] = float(s['ag_RP72_16perc'])
            row['Sa01_RP475'] = float(s['Sa01_RP475'])
            row['Sa01_RP1000'] = float(s['Sa01_RP1000'])
            row['Sa01_RP2500'] = float(s['Sa01_RP2500'])
            row['sismico_non_classificato'] = int(s['sismico_non_classificato'])
            row['hazard_sismico'] = ag
            row['Risk_Sismico_Asset_PMI'] = ag * float(row['asset_PMI_EUR'])
            row['Risk_Sismico_Asset_Grandi'] = ag * float(row['asset_grandi_EUR'])
            w.writerow(row)
    print(f'lette {n_in} righe, join ok {n_join} -> {OUT}')


if __name__ == '__main__':
    main()
