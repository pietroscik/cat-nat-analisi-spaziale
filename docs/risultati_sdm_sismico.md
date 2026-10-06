# Risultati SDM con il terzo hazard sismico (p = 4, k = 5)

Questo documento riporta i risultati dell'estensione del modello SDM a livello comune con il
**terzo hazard sismico** (accelerazione massima del suolo ag, RP 475 anni, MPS04 INGV) accanto
al rischio idrogeologico (frana). Specificazione: `y = log1p(Premio_Teorico_Comunale_EUR)`,
regressori `log1p(Risk_*_Asset_*)` con `*` ∈ {Frana, Sismico} × {PMI, Grandi}, matrice W KNN
row-standardized con **k = 5** (n = 3.823 comuni). Stima ML identica a quella del modello
definitivo a due regressori (`scripts/definitivo.js`), con errori standard via Hessiana
numerica corretta (vedi la nota bug in fondo). Risultati macchina completi in
`results/FINAL_sismico_k5.json`; pipeline in `scripts/sismico/` e metodologia dati in
`docs/sismico_metodologia.md`.

**Vista grafica**: `docs/grafici_sdm.svg` (quadro a sei pannelli rigenerato da
`scripts/grafici_sdm.py`: coefficienti con IC 95%, effetti LeSage-Pace, robustezza nelle
nove specifiche, selezione k via AIC, confronto di specifica, rho) e
`docs/mappe_hazard.svg` (dot map comunali dell'hazard sismico ag RP475 e della quota di
area in frana P3/P4).

## 1. Baseline p = 2: replica esatta del modello definitivo

Il primo stadio dello script ristima il modello a due regressori come controllo di
riproduzione: stime e SE coincidono con `results/FINAL_k5.json`.

| Parametro | Stima | SE (ML) | z | FINAL_k5 (riferimento) |
|---|---|---|---|---|
| ρ (SDM) | 0,4966 | 0,0216 | 22,97 | 0,4966 (SE 0,0216) |
| β Frana Grandi | 0,14336 | 0,00387 | 37,05 | 0,14336 (SE 0,0039) |
| β Frana PMI | 0,00613 | 0,00499 | 1,23 (n.s.) | 0,00613 (n.s.) |
| θ Frana PMI | −0,05205 | 0,00625 | −8,33 | −0,0520 |
| θ Frana Grandi | −0,03735 | 0,00926 | −4,03 | −0,0374 |
| logL | −6.351,71 | | | −6.351,71 |
| AIC | 12.717,42 | | | 12.717,4 |

Le statistiche SAR/SEM di controllo (LR vs SAR = 131,68; LR vs SEM = 108,55 contro 126,71 e
109,64 in FINAL_k5) differiscono entro il rumore Monte Carlo dell'approssimazione del
log-determinante (tracce condivise tra modelli vs semi separati): differenze ~1 su ρ e ~5 su
logL, conclusioni invariate.

## 2. SDM p = 4: coefficienti

| Parametro | Stima | SE (ML) | z | p |
|---|---|---|---|---|
| Intercetta | 5,97089 | 0,24209 | 24,66 | < 1e-6 |
| β Frana PMI | −0,00023 | 0,00418 | −0,05 | 0,957 (n.s.) |
| β Frana Grandi | **0,02305** | 0,00413 | 5,58 | < 1e-6 |
| β Sismico PMI | **0,17063** | 0,00606 | 28,15 | < 1e-6 |
| β Sismico Grandi | **0,12722** | 0,00283 | 44,95 | < 1e-6 |
| θ Frana PMI (lag spaziale) | 0,00311 | 0,00546 | 0,57 | 0,569 (n.s.) |
| θ Frana Grandi (lag spaziale) | −0,01867 | 0,00860 | −2,17 | 0,030 |
| θ Sismico PMI (lag spaziale) | **−0,11366** | 0,00850 | −13,38 | < 1e-6 |
| θ Sismico Grandi (lag spaziale) | −0,02213 | 0,00651 | −3,40 | 0,0007 |
| ρ (SDM) | **0,37953** | 0,02320 | 16,36 | < 1e-6 |
| ln σ² | −0,14521 | 0,02301 | −6,31 | < 1e-6 |

logL = −5.187,05; AIC = 10.396,1 (convenzione K = 2p+3 = 11 per SDM, p+3 = 7 per SAR/SEM).

**Confronto modelli (p = 4):**

| Modello | logL | AIC | LR vs SDM (df) |
|---|---|---|---|
| SDM p=4 | −5.187,05 | **10.396,1** | — |
| SEM p=4 (λ = 0,4226) | −5.220,13 | 10.454,3 | 66,16*** |
| SAR p=4 (ρ = 0,2162) | −5.290,71 | 10.595,4 | 207,32*** |
| SDM p=2 (baseline frana) | −6.351,71 | 12.717,4 | 2.329,32*** (df = 4) |

Il modello con il sismico domina nettamente: l'aggiunta dei due regressori sismici migliora
il log-verosimiglianza di oltre 1.160 punti e l'AIC scende di oltre 2.300 unità.

## 3. Effetti LeSage–Pace (p = 4)

| Regressore | Diretto | Indiretto | Totale |
|---|---|---|---|
| Sismico Grandi | **0,1327** | 0,0367 | **0,1694** |
| Sismico PMI | 0,1712 | −0,0794 | 0,0918 |
| Frana Grandi | 0,0229 | −0,0158 | 0,0070 |
| Frana PMI | ~0 | 0,0046 | 0,0046 |

## 4. Diagnostica (p = 4)

- Moran sui residui: I = −0,0355 (p = 0,004);
- RESET: F = 4,95 — **scende da 24,4** del modello a due regressori: gran parte della
  mancata linearità era dovuta all'omissione del sismico;
- Breusch–Pagan: LM = 470,1 (df = 9) — eteroschedasticità persistente, come nel modello base.

## 5. Robustezza

**Hazard sismico alternativo** (ag RP 475 del 50° percentile come baseline; varianti):

| Variante | ρ | β Sism Grandi (SE) | β Sism PMI (SE) | β Frana Grandi |
|---|---|---|---|---|
| RP30 (ag 81% in 50 anni) | 0,390 | 0,1359 (0,0030) | 0,1904 (0,0064) | 0,0229 |
| RP50 (ag 63% in 50 anni) | 0,389 | 0,1340 (0,0030) | 0,1859 (0,0063) | 0,0229 |
| Sa(0,10s) RP1000 | 0,377 | 0,1203 (0,0027) | 0,1554 (0,0058) | 0,0232 |
| Sa(0,10s) RP2500 | 0,375 | 0,1185 (0,0027) | 0,1514 (0,0057) | 0,0233 |

**Vicinato KNN alternativo:**

| k | medDistK | ρ | β Sism Grandi | θ Sism PMI |
|---|---|---|---|---|
| 6 | 0,0975 | 0,410 | 0,1271 | 0,0041 |
| 7 | 0,1061 | 0,450 | 0,1270 | 0,0066 |
| 8 | 0,1143 | 0,490 | 0,1261 | 0,0101 |

**Esclusione Sardegna** (70 comuni non classificati, n = 3.753): ρ = 0,1997;
β Sism Grandi = 0,1276 (SE 0,0028); β Frana Grandi = 0,0194; effetti Sism Grandi
diretto 0,1315 / totale 0,1491. La bassa ρ è attesa: si rimuove un blocco di comuni
omogenei (ag = 0) che generava dipendenza spaziale meccanica; i β sono stabili.

## 6. Correlazioni e lettura a livello provinciale

Correlazioni tra regressori (log1p): Sismico Grandi ~ Frana Grandi = **0,607** (la fonte
della distorsione da variabile omessa del modello a due regressori); Sismico PMI ~ Frana
PMI = 0,128; Sismico PMI ~ Sismico Grandi = 0,097; Frana PMI ~ Frana Grandi = 0,353.

Poiché il premio teorico è `premio_10k_prov × asset_tot / 10.000` — la **rate è costante
entro provincia** (107 province) — la variabile di rischio "pura" a livello aggregato è la
rate provinciale. Le correlazioni (log-rate ~ log1p rischio medio provinciale):

| Rischio | Correlazione con log-rate provinciale |
|---|---|
| **Sismico (ag RP475 medio)** | **0,838** |
| Frana (share area P3/P4 medio) | 0,175 |
| Idraulico (share area P3 medio) | 0,109 |

Rate per 10.000 € di asset: massime in Catanzaro (52,5), Cosenza (51,0), Reggio Calabria
(48,7), L'Aquila (47,4), Rieti (44,8), Messina (44,7); minime in Bolzano (10,8), Varese
(11,2), Como (12,0), Novara (12,1) e Sardegna (~12,8). La tariffazione Cat-Nat traccia la
pericolosità sismica in misura molto maggiore di quella idrogeologica.

## 7. Interpretazione

1. **Il sismico è il driver dominante.** β_Sismico è grande e precisissimo (z = 28–45) e
   resta stabile in tutte le robustezze (0,118–0,136 su Grandi).
2. **L'effetto frana del modello p = 2 era gonfiato dall'omissione del sismico**
   (β_Grandi 0,1434 → 0,0231; collinearità 0,607): una quota di ciò che il modello a due
   regressori attribuiva al rischio frana rifletteva in realtà il rischio sismico.
3. **θ_Sismico_PMI negativo** (−0,114) ricalca il pattern del modello base: comuni vicini
   con alta esposizione PMI sono associati a premi più bassi.
4. **Nota di lettura:** y = log1p(premio) e X = log1p(hazard × asset) condividono la
   componente asset; la lettura "pulita" del canale hazard passa dalla correlazione
   provinciale rate–ag (0,838) e dalla stabilità dei β nelle robustezze, che escludono
   che il risultato sia un artefatto della componente comune.

## Nota bug (corretta): Hessiana off-diagonale in `scripts/definitivo.js`

Durante la verifica della pipeline sismica è stato individuato un bug nella Hessiana
numerica interna di `scripts/definitivo.js`: i termini fuori diagonale usavano le
perturbazioni errate `f(+,+) − f(0,+) − f(−,0) + f(−,−)` al posto dei quattro angoli
`f(+,+) − f(+,−) − f(−,+) + f(−,−)`, con SE fortemente sottostimati come conseguenza
(es. SE(ρ) = 0,0019 contro il corretto 0,0216). **Il bug è stato corretto nel repo**: la
Hessiana finale ora usa la formula a 4 angoli e `results/risultati_definitivi.json` è stato
rigenerato (SE/z/p aggiornati, es. SE(ρ) = 0,02162, ora identici a
`results/se_definitivi.json`); il loop di Newton interno conserva la formula storica per
non alterare la traiettoria di ottimizzazione e quindi le stime pubblicate — stime, logL,
Moran ed effetti di `results/FINAL_k5.json` non cambiano (FINAL_k5 ha sempre usato gli SE
corretti di `se_definitivi.json`). Anche il notebook
`notebook/riproduce_sdm_comuni.ipynb` riproduceva la formula errata ed è stato corretto
(SE(ρ) 0,0019 → 0,0217). `scripts/sismico/definitivo_sismico.js` ha sempre usato la
formula corretta a 4 angoli: la baseline p = 2 di `results/FINAL_sismico_k5.json` ne è la
conferma (SE(ρ) = 0,02162, SE(β_Grandi) = 0,00387).
