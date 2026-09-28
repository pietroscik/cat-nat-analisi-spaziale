# Analisi spaziale Cat-Nat — premi assicurativi e rischio idrogeologico in Italia

Repository dell'elaborazione quantitativa della tesi magistrale sull'assicurazione delle calamità naturali
(Cat-Nat, Decreto-Legge 78/2025 "Misure urgenti per l'inclusione finanziaria e assicurativa"),
a cura di **Pietro Maietta**.

L'oggetto dell'indagine è la relazione spaziale tra il **rischio idrogeologico** (frane, da dati ISPRA
aggregati a livello comunale) e la **dimensione economica delle imprese esposte** (PMI vs Grandi imprese,
da AIDA), in due prospettive complementari: il livello di **impresa** (30.673 unità) e il livello di
**comune** (3.823 unità).

---

## 1. Struttura del repository

```
README.md                              ← questo file (panoramica, risultati, riproducibilità)
docs/
  capitolo_metodologico_W_k.md         ← matrice dati, transizioni di scala, W e k, domanda di ricerca
  risultati_sdm_comuni.md              ← risultati completi del SDM a livello comune (ML, spreg/PySAL)
  script_spreg_sdm_catnat.py           ← script definitivo della stima SDM (Python, spreg/PySAL)
  log_R_livello_impresa.md            ← guida e integrazione dei log R (spatialreg) livello impresa
scripts/
  step1_kdist_moran.py                 ← selezione k via curva k-dist + Moran (livello comune)
  sdm_grid.js, definitivo.js           ← stima ML della SDM (implementazione pura Python/JS)
  se_definitivi.js                     ← errore standard ML (Hessiana numerica, k=5..8)
results/
  grid_results.json                     ← griglia di selezione k (k-dist, Moran, AIC)
  risultati_definitivi.json            ← stime SDM k=5..8: coefficienti, AIC, LR, effetti
  se_definitivi.json                   ← errori standard ML definitivi (k=5..8)
  FINAL_k5.json                         ← quadro riassuntivo del modello definitivo k=5
data/
  Matrice_Modello_Savelli_Final.csv     ← matrice definitiva: 3.823 comuni × 37 colonne
  log_R_livello_impresa/               ← 12 log R (spatialreg) dell'analisi a livello impresa
```

## 2. Dati e transizioni di scala

| Sorgente | Contenuto | Unità |
|---|---|---|
| AIDA (Bureau van Dijk) | bilanci imprese: asset, EBITDA, ricavi, dipendenti, integrazione verticale, ISP (dalla tesi) | impresa |
| ISTAT / ISPRA | confini e centroidi comunali, popolazione, superfici, pericolosità frana (PAI P3/P4) e idraulica (P3) | comune |
| IVASS — elaborazione Cat-Nat | tariffe premi teorici per 110 province, per 10.000 € di asset esposto | provincia |

La **matrice definitiva** (`data/Matrice_Modello_Savelli_Final.csv`, 3.823 comuni × 37 colonne)
è costruita con due transizioni di scala esplicite (dettaglio e formule verificate nel
capitolo metodologico, §1.3):

- **impresa → comune** (bottom-up): ogni impresa AIDA è assegnata al comune via spatial join
  sui poligoni ISTAT; poi conteggi per classe dimensione (Micro/Piccola/Media/Grande; PMI =
  Micro+Piccola+Media), **somme** di asset (PMI/Grandi/totale), EBITDA, ricavi, dipendenti,
  **medie** di integrazione verticale e ISP_std, **moda** del cluster LISA/Gi* d'impresa;
- **provincia → comune** (top-down, IVASS): `Premio_Teorico_Comunale_EUR = premio_10k_prov ×
  asset_tot_EUR / 10.000` (formula verificata su tutti i 3.823 comuni); 70 comuni sardi con
  tariffa ripartita per l'assetto provinciale 2025 (flag `premio_provincia_appross`);
- **rischio incrociato**: `hazard_frana_share = PAI_area_P3P4_kmq / SUP_kmq` e
  `Risk_Frana_Asset_X = hazard_frana_share × asset_X_EUR` (X = PMI, Grandi).

Correzioni di allineamento documentate: fusioni comunali ricodificate (Sovizzo,
Presicce-Acquarica, Figline e Incisa Valdarno, Trapani/Misiliscemi) e centroidi corretti per
**Gazzo (VI), Lucignano (AR), Olgiate Olona (VA), Telese Terme (BN)** — verifica sensibile ai
fini della matrice W.

## 3. Le due analisi

### 3.1 Livello impresa (R, `spatialreg`) — analisi di evoluzione

- Campione: 30.673 imprese georeferenziate (jitter 1e-4 sui duplicati di coordinate);
- Variabile dipendente: **ISP_bn** — Indicatore Sintetico di Performance, indicatore composito
  costruito e validato nella tesi ([repo tesi-magistrale](https://github.com/pietroscik/tesi-magistrale)):
  media pesata di redditività (ROE, EBITDA/vendite, ROI, rotazione del capitale investito) e
  solidità patrimoniale-finanziaria (D/E, D/EBITDA, attivo totale, PFN/EBITDA), con pesi
  settoriali ATECO stimati via LASSO e normalizzazione robusta OrderNorm
  (contesto interpretativo completo in `docs/capitolo_metodologico_W_k.md`, §3.2);
- Matrice W: KNN, **k = 77 nazionale** (criterio combinato curva k-dist + Moran), style "W";
- 20 subset Dimensione (Micro, Piccola, Media, Grande) × Macroarea (Nord-Ovest, Nord-Est, Centro, Sud, Isole),
  con k* ottimali specifici per subset (da 16 a 74);
- Modelli: OLS, SAR(err), SDM(mix), GMM(err), selezione per AIC + diagnostica (Moran residui, BP, AD, RESET);
- Evidenza chiave: autocorrelazione globale positiva ma debole dell'ISP (Moran I = 0,00865, p < 2,2e-16),
  cluster significativi solo nel 12,72% delle imprese; l'effetto spaziale è modesto e instabile
  nei subset con W poco densa (ρ anche negativi estremi nei subset piccoli), mentre ROI è
  costantemente positivo e forte e l'integrazione verticale negativa.

I 12 log originali sono in `data/log_R_livello_impresa/`; la lettura guidata con i risultati di
riferimento (spesso troncati o assenti nei log) è in `docs/log_R_livello_impresa.md`.

### 3.2 Livello comune (Python, `spreg`/PySAL) — analisi definitiva

- Unità: 3.823 comuni; W: KNN row-standardized con **k = 5** (coordinate corrette);
- Specificazione: SDM su `y = log1p(Premio_Teorico_Comunale_EUR)` con regressori
  `log1p(Risk_Frana_Asset_PMI)` e `log1p(Risk_Frana_Asset_Grandi)`;
- Stima: Massima Verosimiglianza (anche confrontata con SAR e SEM via AIC e test LR).

Risultati salienti (dettaglio completo in `docs/risultati_sdm_comuni.md`):

| Parametro | Stima | SE (ML) | z |
|---|---|---|---|
| ρ (SDM) | **0,4966** | 0,0216 | 22,97 |
| β Grandi | **0,1434** | 0,0039 | 37,05 |
| β PMI | 0,0061 | 0,0050 | 1,23 (n.s.) |
| θ PMI (lag spaziale) | **−0,0520** | 0,0063 | −8,33 |
| θ Grandi (lag spaziale) | −0,0374 | 0,0093 | −4,03 |

- LR vs SAR = 126,71***, LR vs SEM = 109,64***, AIC: SDM (12.717,4) < SEM (12.821,1) < SAR (12.838,1);
- Effetti LeSage–Pace: Grandi → diretto +0,151, totale +0,211; PMI → diretto ~0, totale −0,091
  (effetto indiretto negativo dominante);
- Robustezza su k = 6, 7, 8: β_Grandi stabile (0,142–0,144), θ_PMI stabile (−0,042/−0,050), ρ cresce
  con k (0,50→0,61) come atteso da W più densa; segni e significatività mai invertiti;
- Diagnostica: Moran residui I = −0,0486 (p = 0,004), RESET F = 24,4 (forma funzionale da
  approfondire), Breusch–Pagan LM = 124,4 (eteroschedasticità).

**Interpretazione sintetica.** Il premio teorico Cat-Nat a livello comunale cresce con l'esposizione
delle Grandi imprese al rischio frana (elasticità diretta ~0,15, totale ~0,21); l'esposizione delle
PMI non ha effetto proprio rilevante, ma il suo ritardo spaziale è negativo: comuni "vicini" con
elevata esposizione PMI sono associati a premi più bassi, coerente con una capacità attrattiva
delle Grandi imprese nei comuni confinanti. La dipendenza spaziale (ρ ≈ 0,5) conferma che i premi
non sono indipendenti tra comuni contigui, aspetto rilevante per il disegno della tariffazione
Cat-Nat su base provinciale/comunale. Nota di lettura: il premio incorpora per costruzione il
tasso provinciale e gli asset (§2), quindi il modello stima come la tariffazione lascia spazio
a una risposta al rischio locale.

## 4. Riproducibilità

1. Costruire la matrice partendo dai sorgenti (AIDA, ISTAT/ISPRA, IVASS, centroidi ISTAT) con le
   regole di transizione del §2 (formule complete nel capitolo metodologico §1.3); verificare
   3.823 righe × 37 colonne e l'identità `asset_tot = asset_PMI + asset_grandi`.
2. Selezione di k: `scripts/step1_kdist_moran.py` (curva k-dist, Moran vs k) e
   `results/grid_results.json`.
3. Stima SDM ML: `docs/script_spreg_sdm_catnat.py` (via `spreg`) oppure
   `scripts/definitivo.js` + `scripts/se_definitivi.js` (implementazione equivalente in ML pura,
   con Hessiana numerica per gli errori standard).
4. Confronto modelli: SAR, SEM, SDM con AIC e LR; effetti con matrice (I − ρW)^{-1}.
5. Analisi livello impresa: script R (`spatialreg`/`spdep`) i cui output console sono i log in
   `data/log_R_livello_impresa/`; l'ISP dipende dalla pipeline della tesi
   ([tesi-magistrale](https://github.com/pietroscik/tesi-magistrale)); interpretazione guidata
   in `docs/log_R_livello_impresa.md`.

Ambienti: R 4.x con `spdep`, `spatialreg`, `FNN`, `ggplot2`; Python 3 con `spreg`, `libpysal`.
