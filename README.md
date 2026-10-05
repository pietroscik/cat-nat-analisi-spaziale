# Analisi spaziale Cat-Nat — premi assicurativi e rischio idrogeologico in Italia

> ✅ **Stato**: Production-ready — pipeline completa (L. 78/2025), CI attiva, CITATION, risultati riproducibili.

Analisi quantitativa nata dall'applicazione dei modelli spaziali ai premi assicurativi delle
calamità naturali (Cat-Nat, Decreto-Legge 78/2025 "Misure urgenti per l'inclusione finanziaria
e assicurativa"), a cura di **Pietro Maietta**. Non è la tesi magistrale: i dati impresa AIDA
(bilanci, 30.673 unità) e l'indicatore composito di performance (ISP) provengono dalla pipeline
della tesi ([repo tesi-magistrale](https://github.com/pietroscik/tesi-magistrale)); qui vengono
riusi come variabili di esposizione e di contesto produttivo.

L'oggetto dell'indagine è la relazione spaziale tra il **rischio idrogeologico** (frane, da dati ISPRA
aggregati a livello comunale) e la **dimensione economica delle imprese esposte** (PMI vs Grandi imprese,
da AIDA), in due prospettive complementari: il livello di **impresa** (30.673 unità) e il livello di
**comune** (3.823 unità). Il repository include inoltre un'estensione del modello a livello comune con un
**terzo hazard, la pericolosità sismica** (MPS04 INGV), descritta in §3.3.

---

## 1. Struttura del repository

```
README.md                              ← questo file (panoramica, risultati, riproducibilità)
LICENSE                                ← MIT (codice); i dati terzi restano dei rispettivi titolari
CITATION.cff                           ← metadati di citazione (GitHub li espone in sidebar)
requirements.txt                       ← dipendenze Python (notebook: solo numpy; CI: nbconvert)
package.json                           ← script Node della pipeline ML (definitivo, SE, SE robusti, grid, sismico)
run_all.sh                             ← pipeline one-click: catena pricing/validazione (Python) + SE robusti HC1 (Node)
.github/workflows/ci.yml                ← CI: sintassi JS/Python, rigenerazione verificata degli output, notebook
docs/
  capitolo_metodologico_W_k.md         ← matrice dati, transizioni di scala, W e k, domanda di ricerca
  di
zionario_dati.md                    ← dizionario delle 37 (+16 sismiche) colonne delle matrici
  risultati_sdm_comuni.md              ← risultati completi del SDM a livello comune (ML, spreg/PySAL)
  risultati_sdm_sismico.md             ← risultati del SDM esteso con il terzo hazard sismico (p=4)
  sismico_metodologia.md               ← fonti INGV (MPS04), estrazione griglie, matching, pipeline
  pricing_coerenza.md                  ← benchmark EAL, loss ratio vs tariffe IVASS, coerenza asset, esposizione del tessuto
  mappa_loss.svg                       ← mappa di loss: EAL attesa e loss ratio per comune
  mappa_lisa_tariffa.svg                ← cluster LISA della coerenza tariffaria (HH/LL del loss ratio)
  grafici_pricing.svg                   ← quadro grafico dei risultati: EP, quartili ISP, chi paga, robustezza, Moran
  log_R_livello_impresa.md            ← guida e integrazione dei log R (spatialreg) livello impresa
scripts/
  step1_kdist_moran.py                 ← selezione k via curva k-dist + Moran (livello comune)
  sdm_grid.js, definitivo.js           ← stima ML della SDM (implementazione pura Python/JS)
  se_definitivi.js                     ← errore standard ML (Hessiana numerica, k=5..8)
  se_robusti.js                        ← SE robusti HC1 (sandwich Huber–White) del SDM, k=5..8 (riusa definitivo.js)
  script_spreg_sdm_catnat.py           ← stima di riferimento con spreg/PySAL (CSV_PATH = data/…)
  pricing/pricing_model.py             ← benchmark EAL 3-hazard, loss ratio vs IVASS, mappa di loss
  pricing/tessuto_produttivo.py        ← peso della loss sul tessuto: esposizione per performance ISP
  pricing/ep_curve.py                  ← AAL numerico dalla curva MPS04 (RP30/72/475) e curva EP nazionale
  pricing/chi_paga.py                  ← ripartizione del premio per classe dimensionale (PMI vs Grandi)
  pricing/spazializzazione.py          ← Moran e LISA della coerenza tariffaria + mappa dei cluster
  pricing/robustezza_tessuto.py        ← verifiche
 del risultato del tessuto (winsorizzato, ROA, SLX)
  pricing/validazione_assunzioni.py     ← validazione formale delle assunzioni (curve, OLS, W, seed)
  pricing/grafici.py                    ← quadro grafico SVG dei risultati (6 pannelli dai JSON)
  sismico/                             ← pipeline del terzo hazard sismico (estrazione, matching, SDM p=4)
results/
  grid_results.json                     ← griglia di selezione k (k-dist, Moran, AIC)
  risultati_definitivi.json            ← stime SDM k=5..8: coefficienti, AIC, LR, effetti
  se_definitivi.json                   ← errori standard ML definitivi (k=5..8)
  se_robusti_hc1.json                  ← SE robusti HC1 del SDM: confronto SE ML vs HC1 (k=5..8)
  FINAL_k5.json                         ← quadro riassuntivo del modello definitivo k=5
  FINAL_sismico_k5.json                 ← quadro riassuntivo del modello esteso con il sismico (p=4)
  eal_comuni.csv                        ← EAL benchmark e loss ratio per comune (3.823 righe)
  pricing_benchmark.json                ← parametri, calibrazione, tabella province, regressione rate
  esposizione_tessuto.json             ← peso EAL/EBITDA, quartili ISP, regressione peso, province
  ep_curve.json                        ← AAL numerico (trapezoid RP30-475 + coda) e curva EP nazionale
  chi_paga.json                        ← premio PMI vs Grandi: quote, premio medio, incidenze
  spazializzazione_tariffa.json       ← Moran globale e LISA del loss ratio
  robustezza_tessuto.json             ← specifiche di robustezza del beta della performance
  validazione_assunzioni.json          ← validazione delle assunzioni: correzione coda, OLS, Moran, seed
data/
  Matrice_Modello_Savelli_Final.csv     ← matrice definitiva: 3.823 comuni × 37 colonne
  Matrice_Modello_Savelli_Final_sismico.csv ← matrice estesa con le colonne sismiche: 3.823 × 53
  sismico/matrice_sismica_ingv.csv      ← matching comune → griglie INGV (ag, Sa ai vari RP)
  log_R_livello_impresa/               ← 12
 log R (spatialreg) dell'analisi a livello impresa
notebook/
  riproduce_sdm_comuni.ipynb           ← riproduzione end-to-end dell'analisi SDM k=5 (numpy, assert vs FINAL_k5)
  README.md                            ← guida all'esecuzione del notebook e mappa cella → artefatto
```

## 2. Dati e transizioni di scala

| Sorgente | Contenuto | Unità |
|---|---|---|
| AIDA (Bureau van Dijk) | bilanci imprese: asset, EBITDA, ricavi, dipendenti, integrazione verticale, ISP (dalla tesi) | impresa |
| ISTAT / ISPRA | confini e centroidi comunali, popolazione, superfici, pericolosità frana (PAI P3/P4) e idraulica (P3) | comune |
| IVASS — elaborazione Cat-Nat | tariffe premi teorici per 110 province riconciliate (107 presenti in matrice), per 10.000 € di asset esposto | provincia |
| INGV — MPS04 | pericolosità sismica: ag (RP 475/30/50 anni) e Sa(0,10 s) (RP 475/1000/2500 anni), griglie nazionali | punto griglia → comune |

La **matrice definitiva** (`data/Matrice_Modello_Savelli_Final.csv`, 3.823 comuni × 37 colonne)
è costruita con due transizioni di scala esplicite (dettaglio e formule verificate nel
capitolo metodologico, §1.3):

- **impresa → comune** (bottom-up): ogni impresa AIDA è assegnata al comune via spatial join
  sui poligoni ISTAT; poi conteggi per classe dimensione (Micro/Piccola/Media/Grande; PMI =
  Micro+Piccola+Media), **somme** di asset (PMI/Grandi/totale), EBITDA, ricavi, dipendenti,
  **medie** di integrazione verticale e ISP_std, **moda** del cluster LISA/Gi* d'impresa;
- **provincia → comune** (top-down, IVASS): `Premio_Teorico_Comunale_EUR = premio_10k_prov ×
  asset_tot_EUR / 10.000` (formula verificata su tutti i 3.823 comuni, errore relativo < 1e-5); 70 comuni sardi con
  tariffa ripartita per l'assetto provinciale 2025 (flag `premio_provincia_appross`);
- **rischio incrociato**: `hazard_frana_share = PAI_area_P3P4_kmq / SUP_kmq` e
  `Risk_Frana_Asset_X = hazard_frana_share × asset_X_EUR` (X = PMI, Grandi).

L'estensione sismica aggiunge una terza
 transizione **punto-griglia → comune** (matching al
punto INGV più vicino, vedi `docs/sismico_metodologia.md`): `hazard_sismico = ag_RP475` (in g)
e `Risk_Sismico_Asset_X = ag_RP475 × asset_X_EUR`, con la Sardegna non classificata gestita
esplicitamente (ag = 0 + flag). Nota di unità: frana/idraulico sono quote di area (share),
il sismico è un'accelerazione (g) — i coefficienti non sono confrontabili in magnitudine
tra le due famiglie.

Correzioni di allineamento documentate: fusioni comunali ricodificate (Sovizzo,
Presicce-Acquarica, Figline e Incisa Valdarno, Trapani/Misiliscemi) e centroidi corretti per
**Gazzo (VI), Lucignano (AR), Olgiate Olona (VA), Telese Terme (BN)** — verifica sensibile ai
fini della matrice W.

## 3. Le analisi

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
  cluster significativi solo nel 12,72% delle imprese; l'effetto spaziale è modesto 
e instabile
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

- LR vs SAR = 126,71***, LR vs SEM = 109,64***, AIC: SDM (12.717,4) < SEM (12.823,1) < SAR (12.840,1)
  (convenzione uniforme K = 2p+3 per SDM, K = p+3 per SAR/SEM, σ² inclusa);
- Nota di identificazione: `Risk_Frana_Asset_Grandi` = 0 in 3.216 comuni su 3.823 → con
  `log1p(0) = 0` questi comuni formano un gruppo di riferimento e β_Grandi è identificato sui
  607 comuni con esposizione positiva (2.147 per la PMI; dettaglio in
  `docs/risultati_sdm_comuni.md`, §2);
- Effetti LeSage–Pace: Grandi → diretto +0,151, totale +0,211; PMI → diretto ~0, totale −0,091
  (effetto indiretto negativo dominante);
- Robustezza su k = 6, 7, 8: β_Grandi stabile (0,142–0,144), θ_PMI stabile (−0,042/−0,050), ρ cresce
  con k (0,50→0,61) come atteso da W più densa; segni e significatività mai invertiti;
- Diagnostica: Moran residui I = −0,0486 (p = 0,004), RESET F = 24,4 (forma funzionale da
  
approfondare), Breusch–Pagan LM = 124,4 (eteroschedasticità; inferenze confermate con
  SE robusti HC1, §3.5).

**Interpretazione sintetica.** Il premio teorico Cat-Nat a livello comunale cresce con l'esposizione
delle Grandi imprese al rischio frana (elasticità diretta ~0,15, totale ~0,21); l'esposizione delle
PMI non ha effetto proprio rilevante, ma il suo ritardo spaziale è negativo: comuni "vicini" con
elevata esposizione PMI sono associati a premi più bassi, coerente con una capacità attrattiva
delle Grandi imprese nei comuni confinanti. La dipendenza spaziale (ρ ≈ 0,5) conferma che i premi
non sono indipendenti tra comuni contigui, aspetto rilevante per il disegno della tariffazione
Cat-Nat su base provinciale/comunale. Nota di lettura: il premio incorpora per costruzione il
tasso provinciale e gli asset (§2), quindi il modello stima come la tariffazione lascia spazio
a una risposta al rischio locale. **Questa lettura va aggiornata alla luce del terzo hazard
(§3.3): parte dell'effetto "frana" del modello a due regressori riflette in realtà il rischio
sismico.**

### 3.3 Terzo hazard sismico (SDM p = 4)

Estensione del modello a livello comune con la pericolosità sismica INGV (MPS04, ag RP 475,
matching al punto di griglia più vicino; `docs/sismico_metodologia.md`). Regressori aggiunti:
`log1p(Risk_Sismico_Asset_PMI)` e `log1p(Risk_Sismico_Asset_Grandi)` (SDM p = 4, k = 5, n = 3.823;
dettaglio completo in `docs/risultati_sdm_sismico.md`, quadro macchina in
`results/FINAL_sismico_k5.json`):

| Parametro | Stima | SE (ML) | z |
|---|---|---|---|
| ρ (SDM) | **0,3795** | 0,0232 | 16,36 |
| β Sismico Grandi | **0,1272** | 0,0028 | 44,95 |
| β Sismico PMI | **0,1706** | 0,0061 | 28,15 |
| β Frana Grandi | 0,0231 | 0,0041 | 5,58 |
| β Frana PMI | −0,0002 | 0,0042 | n.s. |
| θ Sismico PMI | **−0,1137** | 0,0085 | −13,38 |
| θ Sismico Grandi | −0,0221 | 0,0065 | −3,40 |
| θ Frana Grandi | −0,0187 | 0,0086 | −2,17 |

- AIC: SDM p=4 (**10.396,1**) < SEM p=4 (10.454,3) < 
SAR p=4 (10.595,4) ≪ SDM p=2 (12.717,4);
  LR vs p=2 = 2.329,3*** (df 4);
- Effetti LeSage–Pace: Sismico Grandi → diretto 0,133, totale 0,169; Sismico PMI → diretto 0,171,
  totale 0,092; Frana Grandi → totale 0,007;
- Diagnostica: Moran residui I = −0,0355 (p = 0,004); RESET F = 4,95 (scende da 24,4); BP LM = 470;
- Robustezza: hazard alternativi (RP30/RP50/Sa RP1000/Sa RP2500) → β_Sism_Grandi 0,118–0,136;
  k = 6/7/8 → 0,126–0,127 (ρ 0,41→0,49); esclusione Sardegna → 0,128.

**Sintesi.** Il sismico è il driver dominante del premio teorico Cat-Nat: β grande e
precisissimo, stabile in ogni robustezza. Nel modello a due regressori β_Frana_Grandi
(0,143) era **gonfiato dall'omissione del sismico** (collinearità log1p frana–sismico Grandi
= 0,607): con p = 4 scende a 0,023 e resta significativo. A livello provinciale — dove la
rate dei premi è definita — la correlazione tra log-rate e log1p(ag medio) è **0,838**
(frana 0,175, idraulico 0,109): la tariffazione traccia la pericolosità sismica molto più
di quella idrogeologica. Reset cala da 24,4 a 4,95: gran parte della mancata linearità del
modello a due regressori era dovuta alla variabile omessa.

### 3.4 Modello di pricing e coerenza asset (benchmark EAL)

Estensione attuariale (`docs/pricing_coerenza.md`, quadro macchina in
`results/pricing_benchmark.json`, dati per comune in `results/eal_comuni.csv`): un
benchmark di **loss annuale attesa** a tre hazard (sismico via curva MDR(ag),
idraulico e frana via share di area esposta; parametri illustrativi documentati) è
confrontato con le tariffe IVASS tramite un loss ratio `tariffa / benchmark calibrato`
(calibrazione a un solo scalare sull'aggregato, c = 1,92). Risultati chiave:

- **EAL benchmark nazionale: 1,70 mld €/anno** (52% del premio teorico di 3,26 mld),
  mix sismico 61% / idraulico 28% / frana 11%;
- **loss ratio mediano 1,01** (p10–p90: 0,48–4,13): la tariffa segue il rischio dove il
  sismico domina, mentre le province a forte idraulico/frana sono **sottopre
zzate**
  rispetto al benchmark (Treviso, Udine, Rimini, Forlì-Cesena ~0,5; Valle d'Aosta 0,28)
  e quelle a basso sismico sono sovrapprezzate (Agrigento 4,6; Lecce 3,5) — effetto del
  pavimento tariffario per i perils non modellati;
- **coerenza asset**: il 95,9% della variabilità del log-premio comunale è componente
  asset (esposizione), il 4,1% rate: la variabile di pricing è la rate provinciale;
  la regressione provinciale (n = 107, SE robusti White) dà elasticità **+5,4
  all'accelerazione sismica** (t = 16,7), +1,0 all'idraulico, **frana non
  significativa** — conferma asset-free del SDM: la tariffazione è un fatto sismico;
- la geografia del loss ratio è **robusta ai parametri** (Spearman ≥ 0,96 su ogni
  perturbazione ×0,5/×2);
- mappa di loss a due pannelli (EAL attesa e loss ratio) in `docs/mappa_loss.svg`.

Estensione al tessuto produttivo (`scripts/pricing/tessuto_produttivo.py`, output
`results/esposizione_tessuto.json`, dettagli in §5 di `docs/pricing_coerenza.md`):

- **peso nazionale della loss: 2,54% dell'EBITDA** (3,26 mld € di EAL calibrata su
  128,5 mld € di valore prodotto; 895 € per addetto, 106.408 € per impresa);
- **esposizione per performance (quartili ISP)**: i comuni con imprese meno
  performanti (Q1) assorbono il 3,91% del proprio EBITDA in loss attesa — quasi
  il doppio dei comuni più performanti (Q4: 2,06%); intensità di esposizione
  (quota EAL / quota EBITDA) 1,54 vs 0,81;
- **modello predittivo** (OLS con SE robusti White, n = 3.791 comuni, R² = 0,762):
  a parità di hazard, dimensione e struttura dimensionale, **+1 deviazione standard
  di ISP riduce il peso della loss del ~31%** (β = −0,375, t = −21,8); quota di
  imprese Grandi con un piccolo premio di esposizione (+0,08, t = 2,2 — concentration
  risk comunale); province estreme Vibo Valentia 8,3% / Isernia 8,2% / Avellino
  8,1% contro Monza-Brianza 0,38% / Lecce 0,40%;
- **chi paga** (`results/chi_paga.json`): le Grandi imprese sono il 10,5% delle unità
  ma pagano il *
*63,7% del premio** (premio medio 569.617 € vs 43.884 € della PMI,
  ×13); il 67% dei comuni non ha Grandi imprese;
- **AAL numerico dalla curva MPS04** (`results/ep_curve.json`): la banda RP30–475
  vale 4,53× il design point (moltiplicatore totale identificabile 6,88); l'AAL
  numerico RP≥30 è **1,71 mld €/anno = 1,66× il benchmark**, limite inferiore
  dichiarato (la coda RP<30 non è stimabile dai tre punti); l'evento 1-in-475 vale
  122,9 mld € = **72 anni di AAL** (banda 16/84: 82–154);
- **validazione delle assunzioni** (`results/validazione_assunzioni.json`): diagnostica
  formale di tutto il modello — ha corretto un errore di forma chiusa nella coda della
  v1 (doppio conteggio della banda RP30–475, 2,30 mld ritirati → 1,71) e valida OLS
  (JB/BP/VIF/Cook/RESET), trasformazione log (skewness 61,8→1,07), matrice W (I
  stabile su k=3–10, replica esatta) e seed delle permutazioni;
- **spazializzazione della coerenza tariffaria**
  (`results/spazializzazione_tariffa.json`, `docs/mappa_lisa_tariffa.svg`): il log
  loss ratio ha **Moran I = 0,649** (z = 69) — l'inadeguatezza tariffaria è un fatto
  spaziale, coerente con la rate provinciale; LISA: 630 comuni HH (Sardegna + Nord
  a basso sismico) vs 426 LL (VdA, Udine, Brescia);
- **robustezza** (`results/robustezza_tessuto.json`): il β della performance
  sopravvive a winsorizzazione (−0,361), trimming (−0,359), misura alternativa ROA
  comunale (−0,41 per SD, t = −44,7) e SLX (−0,376, con W×ISP nullo: effetto tutto
  locale).

**Avvertenza**: il benchmark usa parametri fisici illustrativi e modella solo i 3
hazard della matrice: misura coerenza relativa della tariffazione, non è un modello
di pricing operativo.

### 3.5 Robustezza agli errori standard (HC1)

Il Breusch–Pagan del SDM k=5 segnala eteroschedasticità (LM = 124,4, §3.2): gli errori
standard ML della Hessiana numerica potrebbero quindi sbagliare le inferenze. Il problema
è chiuso con la covarianza sandwich Huber–White robusta, con correzione a campio
ne finito
HC1 — lo stesso trattamento già usato per l'OLS del benchmark pricing
(`scripts/pricing/pricing_model.py`), qui esteso al modello spaziale.

`scripts/se_robusti.js` riusa `scripts/definitivo.js` (stessa stima ML, stesse tracce
Monte Carlo, stessa Hessiana a 4 angoli): stime e SE ML replicano esattamente
`results/se_definitivi.json` (scarto relativo 0), quindi SE ML e SE HC1 sono confrontati
a parità di tutto. Risultati in `results/se_robusti_hc1.json` (k = 5..8); a k = 5:

| Parametro | Stima | SE ML | SE HC1 | z HC1 | p HC1 |
|---|---|---|---|---|---|
| β_Grandi | 0,1434 | 0,0039 | 0,0035 | 41,5 | ≈ 0 |
| ρ | 0,4966 | 0,0216 | 0,0200 | 24,8 | ≈ 0 |
| θ_PMI | −0,0520 | 0,0063 | 0,0065 | −8,0 | 9×10⁻¹⁶ |
| θ_Grandi | −0,0374 | 0,0093 | 0,0084 | −4,4 | 9×10⁻⁶ |
| β_PMI | 0,0061 | 0,0050 | 0,0055 | 1,1 | 0,27 |

(p ≈ 0: underflow della normale, p < 10⁻³⁰⁰.) **Nessuna inferenza sostantiva cambia**: i
coefficienti-chiave β_Grandi e ρ restano altamente significativi con SE HC1 anzi più
stretti degli ML (0,89× e 0,93×: la Hessiana era conservativa proprio dove serve
solidità), θ_PMI e θ_Grandi restano significativi (|z| = 8,0 e 4,4) e β_PMI resta non
significativo (p = 0,27), come già dichiarato in §3.2. Il quadro è identico su tutta la
griglia k = 6–8 (`results/se_robusti_hc1.json`).

## 4. Riproducibilità

1. Costruire la matrice partendo dai sorgenti (AIDA, ISTAT/ISPRA, IVASS, centroidi ISTAT) con le
   regole di transizione del §2 (formule complete nel capitolo metodologico §1.3); verificare
   3.823 righe × 37 colonne e l'identità `asset_tot = asset_PMI + asset_grandi`.
2. Selezione di k: `scripts/step1_kdist_moran.py` (curva k-dist, Moran vs k) e
   `results/grid_results.json`.
3. Stima SDM ML: `scripts/script_spreg_sdm_catnat.py` (via `spreg`) oppure
   `scripts/definitivo.js` + `scripts/se_definitivi.js` (implementazione equivalente in ML pura,
   con Hessiana numerica per gli errori standard); errori standard robusti all'eteroschedasticità:
   `scripts/s
e_robusti.js` → `results/se_robusti_hc1.json` (sandwich Huber–White HC1, §3.5).
4. Riproduzione end-to-end: `notebook/riproduce_sdm_comuni.ipynb` — pipeline completa in un unico
   notebook eseguibile (caricamento matrice e controlli, transizioni di scala §1.3, W KNN, stima ML
   del SDM con confronto SAR/SEM, effetti LeSage–Pace, diagnostica, robustezza k=6–8), con assert
   automatici che verificano la corrispondenza con `results/FINAL_k5.json` entro le tolleranze
   Monte Carlo; guida all'esecuzione in `notebook/README.md`.
5. Estensione sismica (terzo hazard): scaricare le griglie INGV (URL in
   `docs/sismico_metodologia.md`, §1), estrarle con `scripts/sismico/recover_sa.py` +
   `parse_biff.py` (i `.xls` SA hanno OLE2 difettoso) e `xlsx_to_csv.py`, quindi
   `scripts/sismico/match_sismico.py` → `data/sismico/matrice_sismica_ingv.csv` (3.823/3.823
   abbinati, 0 respinti), `scripts/sismico/build_matrice_v2.py` → matrice estesa 53 colonne
   (`data/Matrice_Modello_Savelli_Final_sismico.csv`, archiviata e comunque rigenerabile al
   byte: join deterministico dei due file dati) e `node scripts/sismico/definitivo_sismico.js` →
   `results/FINAL_sismico_k5.json`
   (include la baseline p=2 che replica FINAL_k5 e le robustezze RP/k/Sardegna).
6. Confronto modelli: SAR, SEM, SDM con AIC e LR; effetti con matrice (I − ρW)^{-1}.
7. Analisi livello impresa: script R (`spatialreg`/`spdep`) i cui output console sono i log in
   `data/log_R_livello_impresa/`; i dati impresa AIDA e la costruzione dell'ISP provengono
   dalla pipeline della tesi ([repo tesi-magistrale](https://github.com/pietroscik/tesi-magistrale)),
   che qui si riusa sui premi Cat-Nat; interpretazione guidata
   in `docs/log_R_livello_impresa.md`.
8. Benchmark EAL e coerenza asset: `python3 scripts/pricing/pricing_model.py` (solo stdlib,
   deterministico) → `results/pricing_benchmark.json`, `results/eal_comuni.csv`,
   `docs/mappa_loss.svg`; dettagli e limiti in `docs/pricing_coerenza.md`.
9. Peso della l
oss sul tessuto produttivo: `python3 scripts/pricing/tessuto_produttivo.py`
   → `results/esposizione_tessuto.json` (richiede il passo 8; deterministico); quartili di
   esposizione per performance ISP, regressione del peso-EBITDA e tabella provinciale.
10. Approfondimenti attuariali del pricing (tutti deterministici, solo stdlib, richiedono il
   passo 8): `python3 scripts/pricing/ep_curve.py` → `results/ep_curve.json` (AAL numerico
   dalla curva MPS04 e curva EP nazionale); `python3 scripts/pricing/chi_paga.py` →
   `results/chi_paga.json` (ripartizione del premio PMI vs Grandi); `python3
   scripts/pricing/spazializzazione.py` → `results/spazializzazione_tariffa.json` +
   `docs/mappa_lisa_tariffa.svg` (Moran/LISA del loss ratio, permutazioni a seed fisso);
   `python3 scripts/pricing/robustezza_tessuto.py` → `results/robustezza_tessuto.json`
   (winsorizzato, trim, ROA, SLX, Moran sui residui); `python3
   scripts/pricing/validazione_assunzioni.py` → `results/validazione_assunzioni.json`
   (validazione formale delle assunzioni: correzione della coda, diagnostica OLS,
   trasformazione log, matrice W, seed); `python3 scripts/pricing/grafici.py` →
   `docs/grafici_pricing.svg` (quadro grafico a 6 pannelli, legge i JSON dei risultati).
11. SE robusti HC1 del SDM: `node scripts/se_robusti.js` → `results/se_robusti_hc1.json`
   (Node 18+; deterministico: riusa `scripts/definitivo.js` e verifica internamente di
   replicare `results/se_definitivi.json` a scarto relativo 0).

Tutta la catena deterministica (passi 8–11) si esegue con un solo comando:
`bash run_all.sh` — con `--verify` esegue anche la verifica `git diff` sugli output
rigenerati (la stessa della CI), con `--notebook` il notebook end-to-end e con
`--riferimento` la selezione k e la stima spreg/PySAL (§2–3).

Ambienti: R 4.x con `spdep`, `spatialreg`, `FNN`, `ggplot2`; Python 3 con `spreg`, `libpysal`
(pipeline sismica: solo stdlib); Node.js per gli script ML; notebook: Python 3 con sola
dipendenza `numpy
` (≥ 1.24). A ogni push la GitHub Actions (`.github/workflows/ci.yml`)
verifica la sintassi di tutti gli script, rigenera e controlla al byte gli output della
catena deterministica (pricing, validazione, SE robusti HC1) e esegue il notebook per
intero: la riproducibilità è parte del repository, non una dichiarazione.

## 5. Licenza e citazione

Codice e documentazione sono rilasciati sotto licenza **MIT** (`LICENSE`, Copyright © 2026
Pietro Maietta). I **dati** restano dei rispettivi titolari e sono soggetti alle loro
condizioni: AIDA/Bureau van Dijk (asset aziendali), ISTAT (confini e demografia), ISPRA
(pericolosità idrogeologica), INGV (pericolosità sismica MPS04), IVASS (tariffe dei premi
teorici). Il repo ne archivia solo elaborazioni aggregate a livello comunale per
riproducibilità scientifica; per usi commerciali dei dati fare riferimento ai titolari.

Se usi questo lavoro, cita come da `CITATION.cff`:

> Pietro Maietta (2026). *Cat-Nat — analisi spaziale dei premi assicurativi e del rischio
> idrogeologico e sismico in Italia (3.823 comuni)*. GitHub repository,
> https://github.com/pietroscik/cat-nat-analisi-spaziale
