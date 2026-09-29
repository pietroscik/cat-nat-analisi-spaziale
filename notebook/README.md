# Notebook di riproducibilità

Questa sezione contiene il notebook Jupyter **`riproduce_sdm_comuni.ipynb`**, che riproduce
— da capo e con codice completo — l'intera pipeline dell'analisi SDM definitiva a livello
comune (k = 5) riportata in `results/FINAL_k5.json`.

## Come eseguirlo

```bash
# dipendenza unica: numpy >= 1.24 (tutto il resto è stdlib)
pip install numpy

# esecuzione interattiva
jupyter notebook notebook/riproduce_sdm_comuni.ipynb

# oppure esecuzione completa non interattiva
jupyter nbconvert --to notebook --execute notebook/riproduce_sdm_comuni.ipynb
```

Il notebook legge la matrice dei dati dal percorso relativo `../data/Matrice_Modello_Savelli_Final.csv`:
va quindi eseguito dalla radice del repository (o assicurarsi che il percorso sia raggiungibile).

Tempo di esecuzione: circa 1–2 minuti su un laptop ordinario (la parte costosa è la
verosimiglianza SDM con approssimazione del log-determinante su tracce Monte Carlo).

## Struttura del notebook e corrispondenza con il repository

| Cella | Contenuto | Artefatto del repo corrispondente |
|-------|-----------|----------------------------------|
| C1–C3 | Caricamento CSV, controlli di costruzione (identità asset, formula premio, fix coordinate PRO_COM) e controlli delle transizioni di scala | `data/Matrice_Modello_Savelli_Final.csv`, §1.3 di `docs/capitolo_metodologico_W_k.md` |
| C4 | Matrice W KNN e curva k-distanza (medDistK) | `results/grid_results.json`, `scripts/step1_kdist_moran.py` |
| C5–C6 | Stima ML del SDM (tracce MC T=45, M=120, seed fisso), SE via Hessiana numerica; confronto SAR/SEM con AIC e test LR | `results/FINAL_k5.json`, `results/risultati_definitivi.json` |
| C7 | Effetti diretti/indiretti/totale LeSage–Pace esatti (via soluzione densa) | `results/FINAL_k5.json`, `results/se_definitivi.json` |
| C8 | Diagnostica: Moran sui residui (permutazioni), RESET, Breusch–Pagan | `results/FINAL_k5.json`, `results/diag_results.json` |
| C9 | Robustezza su k = 6, 7, 8 | `results/FINAL_k5.json`, `scripts/sdm_grid.js` |
| C10 | Riepilogo finale con assert automatici | — |

## Corrispondenza numerica e tolleranze

Il notebook contiene **assert automatici** che confrontano ogni risultato con i valori di
riferimento di `results/FINAL_k5.json`. Poiché il log-determinante è approssimato con tracce
Monte Carlo a seme fisso ma implementazione diversa (PCG64 di NumPy vs PRNG degli script
JavaScript), le tolleranze sono quelle dell'incertezza MC (~1e-3 su ρ), ad esempio:

- ρ = 0.4966 ± 0.01 → stima attesa ≈ 0.501
- β_Grandi = 0.14336 ± 0.002 → stima attesa ≈ 0.1433
- logL, AIC e LR entro pochi punti
- effetti, Moran, RESET, BP entro le tolleranze indicate nelle celle

## Percorsi equivalenti per la riproduzione

1. **Script JavaScript** (`scripts/definitivo.js`, `scripts/se_definitivi.js`, `scripts/sdm_grid.js`):
   verificati rigenerare byte-per-byte `results/risultati_definitivi.json` e `results/se_definitivi.json`.
2. **Notebook Python** (questa sezione): stesso algoritmo, stessa matrice, stessi valori di
   riferimento entro la tolleranza MC.
3. **Stima `spreg`** di riferimento in `docs/script_spreg_sdm_catnat.py` (convenzione K diversa:
   K=8 per SDM, K=5 per SAR/SEM).

## Limiti

- L'analisi a **livello impresa** (12 log R, 30.673 imprese) non è riproducibile in questo
  repository: i microdati non sono pubblicabili. Documentazione e riepilogo in
  `docs/log_R_livello_impresa.md`.
- L'**ISP** (Indicatore Sintetico di Performance) è già studiato e documentato nella tesi
  (`pietroscik/tesi-magistrale`); qui è usato solo come variabile già costruita.

## Ambito: modello a due regressori

Il notebook copre il modello **definitivo a due regressori** (rischio frana, k = 5) di
`results/FINAL_k5.json`. L'estensione con il **terzo hazard sismico** (SDM a quattro
regressori, `results/FINAL_sismico_k5.json`) non è nel notebook: è riprodotta dalla pipeline
dedicata `scripts/sismico/` (vedi `docs/sismico_metodologia.md` e
`docs/risultati_sdm_sismico.md`).
