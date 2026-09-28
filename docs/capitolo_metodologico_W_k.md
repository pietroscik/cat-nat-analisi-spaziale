# Capitolo metodologico — Matrice W, scelta di k e specificazione dei modelli

## 1. Perché una matrice di pesi spaziali

Ogni modello di econometria spaziale richiede la pre-specificazione di una **matrice di pesi W**
che formalizza la struttura di vicinato tra unità. La specificazione di W è una scelta
metodologica, non un risultato di stima: condiziona il valore dei parametri (ρ, θ) e la
classificazione dei modelli (SAR vs SEM vs SDM). Per questo il capitolo documenta in modo
esplicito come W è stata costruita ai due livelli di analisi e come il parametro di densità k
è stato scelto con criteri oggettivi e replicabili.

Si adotta la struttura **KNN (k-nearest neighbours)**: ogni unità i è collegata ai suoi k unità più
vicine in termini di distanza euclidea tra coordinate (centroidi comunali a livello comune;
coordinate georeferenziate delle imprese a livello impresa). Rispetto alle matrici di contiguità
(o alle matrici a soglia di distanza), la KNN garantisce:

- **nessuna unità isolata** (problema rilevante a livello di impresa, dove i punti sono distribuiti
  in modo molto disomogeneo sul territorio e i duplicati di coordinate sono frequenti);
- **densità costante** per riga (ogni unità ha esattamente k legami), il che rende confrontabili
  le stime tra subset;
- standardizzazione per riga (style "W"), così che il ritardo spaziale Wy è una media locale.

## 2. Livello impresa (30.673 unità): k = 77

### 2.1 Il criterio di selezione

La scelta del k nazionale segue un **criterio combinato** su due strumenti:

1. **Curva k-distanza**: per k = 1..K si tracciano media, mediana e deviazione standard della
   distanza al k-esimo vicino. Il "gomito" della curva indica il punto oltre il quale aggiungere
   vicini significa includere unità sempre meno effettivamente vicine (perdite di omogeneità
   del vicinato). Con n = 30.673 il range esplorato arriva a K = floor(sqrt(n)·0,5) = 87.
2. **Stabilità di Moran**: per ciascun k si calcola l'I di Moran della variabile target (ISP) e
   si osserva la fascia di k in cui l'autocorrelazione è stabile e significativa.

L'intersezione dei due criteri seleziona **k = 77** (log `analisi_k_nazionale.txt`). La matrice
risultante (`matrice.txt`): 2.361.821 legami non nulli, 0,25% di pesi non nulli, 77 legami medi
per unità, costanti S0 = 30.673, S1 = 699,86, S2 = 125.514,1.

### 2.2 k per subset (Dimensione × Macroarea)

Ogni subset (20 combinazioni) ha ricevuto una propria matrice W KNN con k scelto con lo stesso
criterio, perché densità ottimale e intensità dell'autocorrelazione variano fortemente con la
geometria locale del campione. Range osservato: da k* = 16 (Grande_Isole) a k* = 74
(Grande_Nord-Est). Tabella completa in `docs/log_R_livello_impresa.md`.

**Lezione metodologica** (approfondita lì): nei subset piccoli e dispersi (Micro e Grande
imprese periferiche) il k* basso produce W mal condizionate, con stime ρ inaffidabili (valori
negativi estremi, segno di overfitting del profilo di verosimiglianza) e prevalenza dell'OLS;
nei subset grandi e densi (Piccola/Media imprese) il k* alto produce stime spaziali stabili.

## 3. Livello comune (3.823 unità): k = 5

### 3.1 Correzione delle coordinate

La qualità di una KNN dipende interamente dalla correttezza dei punti. Prima della costruzione
di W sono stati corretti i centroidi di quattro comuni con coordinate duplicate o anomale
(riportate a quelle ufficiali ISTAT): **Gazzo (VI), Lucignano (AR), Olgiate Olona (VA),
Telese Terme (BN)**. Senza la correzione questi comuni finivano nel vicinato di unità distantissime,
distorcendo i lag spaziali locali. Tutti i risultati definitivi (`results/FINAL_k5.json`,
`results/risultati_definitivi.json`) si riferiscono alla matrice con coordinate corrette.

### 3.2 Il criterio di selezione

Anche a livello comune si applica il criterio combinato k-dist + Moran
(`scripts/step1_kdist_moran.py`, griglia in `results/grid_results.json`):

- la curva k-distanza mostra il gomito in corrispondenza di k piccoli (4–6): la distanza mediana
  al k-esimo vicino cresce rapidamente oltre k = 5 (mediana ~0,088 gradi a k = 5, ~0,114 a k = 8);
- la stabilità dei parametri del modello (non solo la significatività di Moran) è stata valutata
  su una griglia k = 5, 6, 7, 8.

**k = 5** è la scelta definitiva per tre ragioni: (i) densità coerente con la scala comunale
italiana (un comune ha mediamente 5–6 comuni strettamente confinanti entro la soglia di
distanza del gomito); (ii) massima parsimonia (AIC cresce con k: 12.717 → 12.630, ma ρ è
monotono in k e assorbe la densità di W, quindi il guadagno di AIC non riflette struttura
spaziale aggiuntiva); (iii) stabilità dei coefficienti sostantivi su tutta la griglia
(β_Grandi = 0,142–0,144; θ_PMI = −0,042/−0,050 in ogni k). Il confronto con k = 77 dell'analisi
a livello impresa è naturale: il rapporto di densità riflette il rapporto tra n delle due
analisi (3.823 ≈ n/8 rispetto a 30.673) e la maggiore dispersione puntuale delle imprese.

## 4. Specificazione e stima a livello comune

Modello: **Spatial Durbin (SDM)**

```
y = ρWy + Xβ + WXθ + ε,   ε ~ N(0, σ²I)
y = log1p(Premio_Teorico_Comunale_EUR)
X = [log1p(Risk_Frana_Asset_PMI), log1p(Risk_Frana_Asset_Grandi)]
W = KNN k=5, row-standardized (coordinate corrette)
```

- **Stima**: Massima Verosimiglianza (script di riferimento `docs/script_spreg_sdm_catnat.py`
  con `spreg`; implementazione equivalente in ML pura in `scripts/definitivo.js`; errori standard
  dalla Hessiana numerica, `scripts/se_definitivi.js`).
- **Confronto tra modelli**: SAR (solo Wy), SEM (errore autorisgressivo), SDM; selezione con AIC
  e test Likelihood Ratio annidati (SDM vs SAR: θ = 0; SDM vs SEM: comune radice θ = −ρβ).
- **Effetti**: decomposizione LeSage–Pace con (I − ρW)^{-1}, riportando diretti, indiretti e totali.
- **Diagnostica**: Moran sui residui, RESET (forma funzionale), Breusch–Pagan (eteroschedasticità).

## 5. Limiti dichiarati

- l'I di Moran residuo, seppur piccolo in valore (−0,049), è significativo (p = 0,004): residui
  con leggera autocorrelazione negativa di corto raggio, tipica di W KNN con k basso;
- RESET significativo (F = 24,4): possibile forma funzionale non pienamente lineare nel
  log1p dei regressori; i segni e la gerarchia degli effetti restano stabili, quindi il
  contenuto sostantivo non è compromesso;
- eteroschedasticità (BP LM = 124,4): gli errori standard ML sono da considerare
  conservativamente minimi; per inferenza finale si raccomanda di affiancare stime robuste
  (es. SDM con errore standard di tipo White o GMM spaziale) — passaggio non ancora eseguito;
- la variabile dipendente a livello comune deriva da tariffe provinciali IVASS ripartite: per
  70 comuni sardi la ripartizione è approssimata (flag in matrice).
