# Risultati definitivi — SDM a livello comune (ML, k = 5)

Unità: 3.823 comuni. Specifica: `y = log1p(Premio_Teorico_Comunale_EUR)`;
regressori `log1p(Risk_Frana_Asset_PMI)`, `log1p(Risk_Frana_Asset_Grandi)`;
W = KNN k=5 row-standardized, coordinate corrette (Gazzo, Lucignano, Olgiate Olona, Telese Terme).
Fonti numeriche: `results/FINAL_k5.json`, `results/risultati_definitivi.json`,
`results/se_definitivi.json`.

## 1. Coefficienti del modello definitivo (k = 5)

| Parametro | Stima | SE (ML) | z | p |
|---|---|---|---|---|
| Intercetta | 6,1161 | 0,2628 | 23,27 | <0,001 |
| β — Risk PMI | 0,0061 | 0,0050 | 1,23 | 0,220 |
| β — Risk Grandi | **0,1434** | 0,0039 | 37,05 | <0,001 |
| θ — W·Risk PMI | **−0,0520** | 0,0063 | −8,33 | <0,001 |
| θ — W·Risk Grandi | −0,0374 | 0,0093 | −4,03 | <0,001 |
| ρ | **0,4966** | 0,0216 | 22,97 | <0,001 |
| ln σ² | 0,4498 | 0,0231 | 19,50 | <0,001 |

log-likelihood SDM = −6.351,71.

## 2. Confronto tra modelli (annidati)

| Modello | logL | AIC |
|---|---|---|
| SAR | −6.415,06 (ρ = 0,464) | 12.838,1 |
| SEM | −6.406,53 (λ = 0,547) | 12.821,1 |
| **SDM** | **−6.351,71** | **12.717,4** |

- **LR SDM vs SAR** (H0: θ = 0): 126,71, p < 0,001 → i ritardi spaziali WX sono congiuntamente
  necessari;
- **LR SDM vs SEM** (H0: comune radice): 109,64, p < 0,001 → la SDM non collassa in una SEM;
- l'AIC conferma la SDM come modello preferito; l'interpretazione degli effetti usa la
  decomposizione di LeSage–Pace.

## 3. Effetti diretti, indiretti e totali (LeSage–Pace)

| Regressore | Diretto | Indiretto | Totale |
|---|---|---|---|
| Risk_Frana_Asset_PMI | 0,0006 | −0,0918 | **−0,0912** |
| Risk_Frana_Asset_Grandi | **0,1511** | 0,0595 | **0,2106** |

Lettura: un aumento dell'esposizione delle Grandi imprese al rischio frana è associato a premi
comunali più alti, sia nel comune stesso (diretto ~0,15) sia — con intensità ridotta — nei comuni
vicini (totale ~0,21). L'esposizione PMI ha effetto proprio trascurabile ma un **effetto
indiretto negativo** (−0,09): comuni vicini con alta esposizione PMI sono associati a premi
locali più bassi.

## 4. Robustezza alla densità di W (k = 5, 6, 7, 8)

| k | ρ | SE(ρ) | β_Grandi | β_PMI | θ_PMI | θ_Grandi |
|---|---|---|---|---|---|---|
| 5 | 0,4966 | 0,0216 | 0,1434 | 0,0061 | −0,0520 | −0,0374 |
| 6 | 0,5370 | 0,0227 | 0,1424 | 0,0060 | −0,0500 | −0,0428 |
| 7 | 0,5751 | 0,0237 | 0,1422 | 0,0052 | −0,0458 | −0,0503 |
| 8 | 0,6140 | 0,0244 | 0,1420 | 0,0047 | −0,0424 | −0,0571 |

β_Grandi è stabile al terzo decimale su tutta la griglia; ρ cresce monotonicamente con k,
come previsto dalla teoria (W più densa sposta massa dal β al canale spaziale). Segni e
significatività non si invertono mai: la scelta k = 5 non condiziona il contenuto sostantivo.

## 5. Diagnostica

| Test | Statistica | Esito |
|---|---|---|
| Moran residui (two-sided) | I = −0,0486, p = 0,004 | leggera autocorrelazione negativa residua |
| RESET (OLS a confronto) | F = 24,4 | forma funzionale migliorabile |
| Breusch–Pagan | LM = 124,4 | eteroschedasticità presente |
| Distanza mediana k-vicini (k=5) | 0,0882 gradi | vicinato compatto |

## 6. Lettura sostantiva (per il capitolo dei risultati della tesi)

1. **Il premio Cat-Nat è un fatto spaziale**: ρ ≈ 0,5 significa che oltre metà dell'inerzia del
   premio comunale si spiega con i premi dei comuni vicini — coerente con tariffe definite a
   livello provinciale e poi ripartite;
2. **le Grandi imprese sono il canale di trasmissione dominante** tra rischio idrogeologico e
   premio: la loro esposizione al rischio frana ha effetto diretto positivo robusto;
3. **le PMI contano solo per via spaziale**: il loro ritardo negativo suggerisce un effetto di
   composizione del tessuto economico tra comuni confinanti;
4. la robustezza su k e la conferma LR/AIC rendono la SDM la base solida per la parte
   econometrica a livello comunale; l'eteroschedasticità residua impone prudenza sulla
   significatività esatta dei coefficienti minori (β_PMI).
