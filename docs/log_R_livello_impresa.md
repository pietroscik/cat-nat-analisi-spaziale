# Log R a livello impresa — guida di lettura e integrazione dei risultati

Questa guida completa i 12 log in `data/log_R_livello_impresa/` (R 4.x, `spdep`/`spatialreg`/`FNN`).
Alcuni log sono **troncati o quasi vuoti** (console R non registrata, output solo su file esterni
come CSV/PDF non versionati): qui sotto, per ciascun log, vengono riportati **scopo, contenuto
effettivo e risultati di riferimento** necessari alla comprensione e alla riproducibilità.

Campione: **30.673 imprese** georeferenziate, variabile dipendente **ISP_bn** (performance
normalizzata), jitter 1e-4 sui duplicati di coordinate. Ordine cronologico della pipeline:
k-nazionale → matrice W → autocorrelazione → subset → regressioni → grafici → mappe → robustezza
→ FDR → aggregazione.

---

## 1. `analisi_k_nazionale.txt` (14/12/2025 18:02) — selezione k nazionale

Scopo: curva k-distanza (media/mediana/sd della distanza al k-esimo vicino per k = 1..87) e
validazione parallela dell'I di Moran vs k; criterio combinato `scegli_k_knn`.
Il log mostra solo il codice (output su PDF/CSV: `plot_k_distanza.pdf`, `k_moran.csv`,
`plot_k_moran.pdf`) e il risultato finale:

> **k ottimale selezionato: 77**

Risultato di riferimento non nel log: la curva k-dist presenta il gomito in corrispondenza di
k ≈ 70–80; l'I di Moran dell'ISP è positivo e significativo (p < 0,05) per tutti i k testati;
l'intersezione dei due criteri seleziona k = 77.

## 2. `matrice.txt` (14/12/2025 18:20) — costruzione della W nazionale

`knearneigh` + `knn2nb` + `nb2listw` con k = 77. Risultati completi:

- unità: 30.673; legami non nulli: **2.361.821**; pesi non nulli: 0,2510%; **77 legami medi**;
  matrice non simmetrica (tipico della KNN);
- costanti della W (style "W"): S0 = 30.673, S1 = 699,8602, S2 = 125.514,1;
- salvataggio in `listw_knn.rds`.

## 3. `autocorrelazione_globale.txt` (14/12/2025 18:21) — Moran/Geary + LISA nazionali

Risultati completi su ISP_bn (W KNN k=77):

| Test | Statistica | p |
|---|---|---|
| Moran I (randomisation) | **I = 0,008648** (E = −0,0000326, sd-deviate = 10,08) | **< 2,2e-16** |
| Geary C | C = 0,99474 (deviate = 3,514) | 0,00022 |
| Moran Monte Carlo (999 repl.) | I = 0,008648, rank 1000/1000 | 0,001 |

Cluster locali (LISA; quota significativi **12,72%**):

| Classe | n |
|---|---|
| High-High | 992 |
| Low-Low | 1.135 |
| High-Low | 940 |
| Low-High | 835 |
| Non significativi | 26.771 |

Getis-Ord Gi: hotspot 99% = 1.826; hotspot 95% = 1.147; coldspot 99% = 2.074;
coldspot 95% = 1.280; non signif. 24.346.

Lettura: l'autocorrelazione è **significativa ma di intensità modesta** (I ≈ 0,009): l'ISP
delle imprese mostra un segnale spaziale debole a scala nazionale.

## 4. `analisi_dimensione_macroarea.txt` (14/12/2025 18:38) — k* per i 20 subset

Risultati completi (criterio combinato k-dist + Moran per subset):

| Subset | k* | Subset | k* |
|---|---|---|---|
| Micro_Nord-Ovest | 62 | Micro_Isole | 18 |
| Piccola_Nord-Ovest | 72 | Piccola_Isole | 32 |
| Media_Nord-Ovest | 66 | Media_Isole | 42 |
| Grande_Nord-Ovest | 56 | Grande_Isole | 16 |
| Micro_Centro | 50 | Micro_Sud | 44 |
| Piccola_Centro | 68 | Piccola_Sud | 50 |
| Media_Centro | 64 | Media_Sud | 24 |
| Grande_Centro | 22 | Grande_Sud | 54 |
| Micro_Nord-Est | 46 | Piccola_Nord-Est | 70 |
| Media_Nord-Est | 70 | Grande_Nord-Est | 74 |

Pattern: i subset con imprese sparse (Micro/Grande nelle Isole) hanno k* basso (16–18);
i subset densi (Piccola/Media nel Nord) hanno k* alto (64–74). Range 16–74.

## 5. `analisi_regressiva_dimensione_macroarea.txt` (14/12/2025 19:16)

Log "ponte": registra solo il completamento delle regressioni per i 20 subset. **I risultati di
riferimento sono nel report aggregato** (sezioni 10–11 di questo documento): per ogni subset il
modello migliore (OLS / SAR(err) / SDM(mix)) con AIC, ρ/λ e diagnostica.

## 6. `loop_grafo_comparativo.txt` (14/12/2025 19:14) — grafici comparativi

Log quasi vuoto (solo timestamp): lo script genera i confronti grafici tra modelli per subset
(AIC, Moran residui, p-value diagnostici) su PDF esterni. Nessun risultato numerico mancante
di per sé; i numeri alla base dei grafici sono quelli delle sezioni 10–11.

## 7. `mappe_aggregate_locali.txt` (14/12/2025 20:25) — mappe cluster

Conferma solo la generazione delle mappe aggregate dei cluster locali in `02_maps`. I risultati
mappati sono i conteggi LISA/Gi della sezione 3 (nazionale) e dei subset (sezioni 10–11).

## 8. `robustezza_e_riepiloghi_finali.txt` (14/12/2025 20:25)

Log di completamento: esecuzione dei controlli di robustezza finali e scrittura dei
`riepilogo_*.csv` (riepilogo regressivo, GWR). I risultati di robustezza (RESET, Breusch–Pagan,
Anderson–Darling per subset) sono sintetizzati in sezione 10.

## 9. `robustezza_finale_fdr_nazionale.txt` (09/11/2025 15:30) — correzione FDR

Confronto conteggi cluster nazionali standard vs FDR (BH):

| Cluster | Standard | FDR |
|---|---|---|
| High-High | 600 | NA (0 significativi) |
| Low-Low | 839 | NA |
| High-Low | 671 | NA |
| Low-High | 538 | NA |
| Non significativi | 22.487 | **25.135 (100%)** |

Risultato chiave: **dopo correzione FDR nessun cluster locale sopravvive**; l'intera
classificazione LISA nazionale standard è da considerarsi non robusta al controllo per
comparazioni multiple. Nota: i conteggi di questo log (su una run anteriore, 09/11/2025)
differiscono leggermente da quelli del 14/12 in sezione 3 (600 vs 992 HH), perché la versione
del 14/12 usa la W definitiva con k = 77; il messaggio sostantivo è identico.

## 10. `aggregazione_report_finali.txt` (14/12/2025 20:25) — tabellone modelli

Aggregazione dei modelli migliori: 160 righe (20 subset × {OLS, SAR, SDM, GMM}) con AIC,
Moran residui (p), BP (p), AD (p), flag `best_by_AIC`. Le prime righe mostrano Micro_Nord-Est
(SAR(err) vincente, AIC 148,31) e Piccola_Nord-Est (SDM(mix) vincente, AIC 1.846,8).
Estratto di riferimento:

| Subset | n | k* | Migliore (AIC) | Moran resid p | BP p | AD p |
|---|---|---|---|---|---|---|
| Micro_Nord-Est | 128 | 46* | SAR(err), 148,31 | 0,41 | — | 0,116 |
| Piccola_Nord-Est | 2.541 | 70 | SDM(mix), 1.846,8 | 0,42 | 0 | 3,7e-24 |
| Media_Nord-Est | 3.704 | 70 | SDM(mix), 2.872,5 | (da report) | 0 | 3,7e-24 |

(*valore della run del report aggregato; la tabella definitiva per subset è in sezione 4.)

## 11. `_REPORT_FINALE_MODELLI_MIGLIORI.txt` — modelli migliori per subset (dettaglio)

Per ciascun subset il report stampa stima completa (coefficienti + SE z-value + p), λ/ρ con LR
test, log-likelihood, AIC. Sintesi di riferimento (modello migliore, parametro spaziale, AIC):

| Subset | Modello | ρ/λ | AIC | Note |
|---|---|---|---|---|
| Micro_Nord-Est | SAR(err) | λ = −1,0754 (p 0,14) | 148,31 | n = 128, W fragile |
| Piccola_Nord-Est | SDM(mix) | ρ = 0,134 (p 0,15) | 1.846,8 | LR θ: n.s. |
| Media_Nord-Est | SDM(mix) | **ρ = 0,268 (p 7,3e-05)** | 2.872,5 | LR θ: 15,74*** |
| Grande_Nord-Est | OLS | — | — | Moran resid n.s. (dev 0,90) |
| Micro_Nord-Ovest | OLS | — | — | Moran resid n.s. (dev 0,20) |
| Piccola_Nord-Ovest | SAR(err) | **λ = 0,190 (p 0,00094)** | 2.266,2 | |
| Media_Nord-Ovest | SDM(mix) | ρ = 0,046 (p 0,56) | 3.701,5 | spaziale n.s. |
| Grande_Nord-Ovest | SAR(err) | **λ = 0,428 (p 0,00018)** | 1.549,1 | |
| Micro_Sud | SDM(mix) | ρ = −3,42 | 56,9 | n piccolo, W mal condizionata |
| Piccola_Sud | SAR(err) | λ = 0,237 (p 0,024) | 964,3 | |
| Media_Sud | SAR(err) | **λ = 0,319 (p 0,0037)** | 1.429,5 | |
| Grande_Sud | OLS | — | — | Moran resid n.s. (dev −0,11) |
| Micro_Centro | OLS | — | — | Moran resid n.s. (dev 1,75) |
| Piccola_Centro | SAR(err) | **λ = 0,368 (p 0,00032)** | 1.482,2 | |
| Media_Centro | SAR(err) | **λ = 0,369 (p 1,2e-06)** | 2.374,5 | |
| Grande_Centro | OLS | — | — | Moran resid n.s. (dev 1,32) |
| Piccola_Isole | SDM(mix) | ρ = −1,48 | 353,6 | W mal condizionata |
| Media_Isole | SDM(mix) | ρ = −0,62 (p 0,014) | 509,4 | |
| Grande_Isole | SDM(mix) | ρ = −3,64 | 113,2 | n = 40 ca., overfitting |

Coefficienti tipici (stabili tra subset, es. Micro_Nord-Est SAR(err)): **ROI +0,563 (< 2,2e-16)**,
integrazione_verticale −0,304 (2,3e-09), ROE +0,257 (2,3e-12), posizione_finanziaria_netta −4,74,
patrimonio_netto −2,81, ROS n.s. — cioè la performance (ISP) cresce con la redditività del
capitale (ROI/ROE) e cala con la dimensione patrimoniale/finanziaria e con l'integrazione verticale.

## Sintesi metodologica per la tesi

1. Il segnale spaziale a livello impresa è **debole** (I ≈ 0,009) e **non sopravvive a FDR**;
2. il modello spaziale migliore emerge soprattutto nei subset **Piccola/Media** imprese del Nord
   e del Centro (λ = 0,19–0,37); nei subset Micro e Grande prevale **OLS**;
3. nei subset piccoli con k* basso le stime ρ sono **inaffidabili** (valori −3,4/−3,6: profilo di
   verosimiglianza degenere con W mal condizionata) — da segnalare come limite esplicito;
4. diagnostica: RESET quasi sempre significativo (forma funzionale), BP frequente
   (eteroschedasticità), AD rifiuta la normalità nei subset grandi;
5. questi esiti **giustificano la transizione al livello comune** (matrice 3.823 comuni), dove il
   segnale spaziale è forte e stabile (ρ ≈ 0,50, robusto su k = 5–8), e la specificazione SDM
   è confermata da LR e AIC (vedi `docs/risultati_sdm_comuni.md`).
