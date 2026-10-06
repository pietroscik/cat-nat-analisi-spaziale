# Notebook di riproducibilità

Questa sezione contiene il notebook Jupyter **`riproduce_sdm_comuni.ipynb`**, che riproduce
— da capo e con codice completo — l'intera pipeline dell'analisi SDM definitiva a livello
comune (k = 5) riportata in `results/FINAL_k5.json`, **e l'estensione al terzo hazard
sismico** (SDM a quattro regressori) riportata in `results/FINAL_sismico_k5.json` (`sdm_p4`).

## Come eseguirlo

```bash
# dipendenza unica: numpy >= 1.24 (tutto il resto è stdlib)
pip install numpy

# esecuzione interattiva
jupyter notebook notebook/riproduce_sdm_comuni.ipynb

# oppure esecuzione completa non interattiva
jupyter nbconvert --to notebook --execute notebook/riproduce_sdm_comuni.ipynb
```

Il notebook legge le matrici dei dati dai percorsi relativi
`../data/Matrice_Modello_Savelli_Final.csv` (modello p=2) e
`../data/Matrice_Modello_Savelli_Final_sismico.csv` (estensione sismica p=4):
va quindi eseguito dalla radice del repository (o assicurarsi che i percorsi siano
raggiungibili).

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
| C11 | Matrice sismica (53 colonne): caricamento, ordine PRO_COM, identità `Risk_Sismico = ag_RP475 x asset`, 70 comuni non classificati, collinearità attese (0,607 / 0,128) | `data/Matrice_Modello_Savelli_Final_sismico.csv`, `results/FINAL_sismico_k5.json` (`corr`) |
| C12 | SDM ML p=4 generalizzato (concentrata su ρ + Newton + Hessiana a 4 angoli): 11 parametri con assert su stime e SE | `results/FINAL_sismico_k5.json` (`sdm_p4`) |
| C13 | Confronto SAR/SEM p=4 (AIC, LR), LR vs baseline p=2 (nested, df=4), effetti LeSage–Pace esatti sui 4 regressori, Moran sui residui | `results/FINAL_sismico_k5.json`, `scripts/sismico/definitivo_sismico.js` |
| C10/C14 | Riepilogo finale con assert automatici (p=2 e p=4) | — |

## Corrispondenza numerica e tolleranze

Il notebook contiene **assert automatici** che confrontano ogni risultato con i valori di
riferimento di `results/FINAL_k5.json`. Poiché il log-determinante è approssimato con tracce
Monte Carlo a seme fisso ma implementazione diversa (PCG64 di NumPy vs PRNG degli script
JavaScript), le tolleranze sono quelle dell'incertezza MC (~1e-3 su ρ), ad esempio:

- ρ = 0.4966 ± 0.01 → stima attesa ≈ 0.501
- β_Grandi = 0.14336 ± 0.002 → stima attesa ≈ 0.1433
- logL, AIC e LR entro pochi punti
- effetti, Moran, RESET, BP entro le tolleranze indicate nelle celle
- estensione sismica p=4 (§7 del notebook): stesse tolleranze di tipo MC sui 11 parametri
  (ρ = 0.37953 ± 0.01, β_Sism_Grandi = 0.12722 ± 0.002, β_Frana_Grandi = 0.02305 ± 0.003),
  logL entro ±5 punti, LR vs SAR/SEM/baseline nei range dichiarati nelle celle C12–C13

## Percorsi equivalenti per la riproduzione

1. **Script JavaScript** (`scripts/definitivo.js`, `scripts/se_definitivi.js`, `scripts/sdm_grid.js`):
   verificati rigenerare byte-per-byte le stime di `results/risultati_definitivi.json` (i cui SE, dopo il
   fix dell'Hessiana, coincidono con `results/se_definitivi.json`).
2. **Notebook Python** (questa sezione): stesso algoritmo, stessa matrice, stessi valori di
   riferimento entro la tolleranza MC.
3. **Stima `spreg`** di riferimento in `scripts/script_spreg_sdm_catnat.py` (convenzione K diversa:
   K=8 per SDM, K=5 per SAR/SEM).

## Limiti

- L'analisi a **livello impresa** (12 log R, 30.673 imprese) non è riproducibile in questo
  repository: i microdati non sono pubblicabili. Documentazione e riepilogo in
  `docs/log_R_livello_impresa.md`.
- L'**ISP** (Indicatore Sintetico di Performance) è già studiato e documentato nella tesi
  (`pietroscik/tesi-magistrale`); qui è usato solo come variabile già costruita.

## Ambito: modello a due regressori + estensione sismica

Il notebook copre il modello **definitivo a due regressori** (rischio frana, k = 5) di
`results/FINAL_k5.json` **e l'estensione con il terzo hazard sismico** (SDM a quattro
regressori di `results/FINAL_sismico_k5.json`, §7: stima p=4, confronto SAR/SEM, effetti,
Moran). Le robustezze della specifica sismica (hazard alternativi RP30/RP50/Sa RP1000/Sa
RP2500, esclusione Sardegna, k = 6–8) restano nella pipeline dedicata `scripts/sismico/`
(vedi `docs/sismico_metodologia.md` e `docs/risultati_sdm_sismico.md`); il quadro grafico
complessivo è `docs/grafici_sdm.svg`.
