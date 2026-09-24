# Scelta della matrice dei pesi W e sensibilità al parametro k

*Paragrafo per il capitolo metodologico — analisi Cat-Nat (L. 78/2025) su 3.823 comuni italiani.*

## 1. Disegno della verifica

La specificazione della matrice dei pesi spaziali W è la scelta metodologica più delicata dell'analisi: i coefficienti spaziali (ρ, θ) e la significatività dei cluster locali dipendono dalla densità di W, mentre — come documentato nella Sezione 7 — i coefficienti sostantivi sono risultati invarianti. La strategia adottata replica il protocollio già impiegato a livello di impresa (curva k-dist + validazione Moran, sez. precedente del capitolo) e lo estende con una batteria di stime ML complete su una griglia di k, in modo che la scelta finale sia documentata empiricamente e non assunta per convenzione.

**Specificazione stimata:** y = log1p(Premio_Teorico_Comunale_EUR); X = log1p(Risk_Frana_Asset_PMI), log1p(Risk_Frana_Asset_Grandi); W = KNN row-standardized sui centroidi comunali (corretti, v. sotto). Stima ML (logdet via serie di tracce Monte Carlo, approccio Barry–Pace), errore standard dalla matrice di covarianza inversa della Hessiana numerica della log-likelihood piena.

## 2. Correzione preliminare dei dati

Un controllo di qualità sulle coordinate prima della stima ha evidenziato **4 errori di geolocalizzazione** nella matrice (latitudine/longitudine invertite o separatore decimale perso), che distorcevano il KNN di 4 comuni su 3.823: Gazzo (028041), Lucignano (051021), Olgiate Olona (012108), Telese Terme (062074). I centroidi sono stati corretti con i valori ISTAT e la matrice definitiva incorpora il fix. Il controllo (distanza al primo vicino > 0,5° ≈ 55 km, range ammissibile 35–47°N / 6–19°E) è riportabile come procedura standard di validazione.

## 3. Criterio combinato k-dist + Moran

- **Curva k-dist**: l'incremento della distanza mediana al k-esimo vicino scende sotto 0,01° per k ≥ 5; a k = 5 il vicinato mediano è ~0,088° ≈ 9–10 km, coerente con la scala comunale. A k = 20 il "vicino" mediano sarebbe a ~0,19° ≈ 21 km, oltre la scala del fenomeno.
- **Moran I(y)**: positivo e significativo per ogni k (I = 0,246 a k = 3; 0,241 a k = 5; 0,191 a k = 30; z ≥ 7,9), con intensità decrescente al crescere di k — come atteso, W più dense diluiscono l'autocorrelazione locale.
- **AIC**: continua a scendere fino a k ≈ 20, ma il guadagno è un artefatto: ρ sale da 0,50 a 0,83 e θ_Grandi deriva da −0,037 a −0,116, segno che la W densa assorbe eterogeneità non spaziale piuttosto che dipendenza. La regola empirica √n·0,5 ≈ 30 (valida per le 30.673 imprese, k* = 77) non è applicabile ai comuni: i punti sono radi e il k ottimale deve scalare con la geometria del campione.

**Scelta finale: k = 5.**

## 4. Stima definitiva — SDM ML, k = 5 (n = 3.823)

| Coefficiente | Stima | SE | z | p |
|---|---|---|---|---|
| Intercetta | 6,1161 | 0,2629 | 23,27 | < 0,001 |
| Risk_Frana_Asset_PMI (β) | 0,0061 | 0,0050 | 1,23 | 0,220 n.s. |
| Risk_Frana_Asset_Grandi (β) | 0,1434 | 0,0039 | 37,05 | < 0,001 |
| W·Risk_PMI (θ) | −0,0521 | 0,0063 | −8,33 | < 0,001 |
| W·Risk_Grandi (θ) | −0,0374 | 0,0093 | −4,03 | < 0,001 |
| ρ | 0,4966 | 0,0216 | 22,97 | < 0,001 |
| ln σ² | 0,4498 | 0,0231 | 19,50 | < 0,001 |

## 5. Test di specificazione

| Test | Statistica | Esito |
|---|---|---|
| LR SDM vs SAR (H₀: θ = 0) | 126,71 (gdl 2) | p < 0,001 → SDM |
| LR SDM vs SEM (H₀: θ = −ρβ, Burridge) | 109,64 (gdl 2) | p < 0,001 → SDM |
| AIC SDM / SEM / SAR | 12.717,4 / 12.821,1 / 12.838,1 | SDM vincente |

La gerarchia SDM > SEM > SAR è **invariante su tutta la griglia k = 3–30** (LR SDM/SAR sempre ≥ 81, p < 0,001): la scelta del modello Durbin non dipende da k.

## 6. Effetti diretti, indiretti e totali (LeSage–Pace, k = 5)

| Variabile | Diretto | Indiretto (spillover) | Totale |
|---|---|---|---|
| Risk_Frana_Asset_PMI | +0,0006 (n.s.) | −0,0918*** | −0,0912*** |
| Risk_Frana_Asset_Grandi | +0,1511*** | +0,0595*** | +0,2106*** |

Lettura sostantiva: l'esposizione all'asset delle Grandi imprese **trasla il premio teorico** del comune (effetto diretto +0,151) e genera spillover positivi sui vicini; l'esposizione delle PMI non ha effetto proprio diretto, ma **spillover negativi** significativi dai comuni vicini — coerente con la lettura competitiva del rischio di portafoglio a livello territoriale.

## 7. Robustezza: banda k = 5–8

| k | ρ (SE) | β Grandi | β PMI | θ PMI | θ Grandi |
|---|---|---|---|---|---|
| 5 | 0,4966 (0,0216) | 0,1434*** | 0,0061 n.s. | −0,0521*** | −0,0374*** |
| 6 | 0,5370 (0,0227) | 0,1424*** | 0,0060 n.s. | −0,0500*** | −0,0428*** |
| 7 | 0,5751 (0,0237) | 0,1422*** | 0,0052 n.s. | −0,0458*** | −0,0503*** |
| 8 | 0,6140 (0,0244) | 0,1420*** | 0,0047 n.s. | −0,0424*** | −0,0571*** |

β_Grandi è stabile a 0,142–0,143 su **tutti** i k testati (3–30); θ restano negativi e significativi; solo ρ cresce meccanicamente con k (W più densa ⇒ più dipendenza catturata). La conclusione sostantiva non dipende da k.

## 8. Diagnostica dei residui (k = 5)

- **Moran sui residui SDM**: I = −0,049, p = 0,004 (permutazione, 499). Autocorrelazione residua *negativa*: l'SDM sovra-cattura leggermente la dipendenza (fenomeno noto con W KNN row-standardized, che introduce una correlazione negativa media tra vicini). Da riportare insieme al segnale positivo su y.
- **RESET** (1 termine, modello trasformato): F = 24,4, p < 0,001 — segnale di forma funzionale non catturata (presente anche a livello impresa); possibili rimedi: termini quadratici, o controlli per popolazione/densità.
- **Breusch–Pagan**: LM = 124,4 (gdl 5), p < 0,001 — eteroschedasticità marcata (coerente con i log d'impresa). Prescrive errori standard robusti/HC nella stima di riferimento o ML con covarianza robusta.

## 9. Confronto con il livello impresa (MAUP)

| | Livello impresa (30.673) | Livello comune (3.823) |
|---|---|---|
| k* selezionato | 77 (criterio combinato) | 5 (criterio combinato) |
| Densità punti | molto alta (urbano) | bassa (centroidi distanti) |
| Moran ISP/premio | I = 0,009 (debole, s. per n) | I = 0,241 (forte) |
| ρ SDM subset | 0,13–0,37, spesso n.s. | 0,497*** |
| Cluster LISA | 12,7%, FDR li azzera | n = 3.823 ⇒ test multipli contenuti |
| Driver stabili | ROI +, integrazione vert. − | β_Grandi + stabile, θ < 0 |

Il passaggio impresa→comune **amplifica** l'autocorrelazione territoriale (media del rumore idiosincratico) e riduce l'esposizione al test multiplo: le due scale non sono in contraddizione, sono la stessa struttura vista a due risoluzioni (MAUP).

## 10. Note operative

- Numeri ufficiali di riferimento: Sez. 4–6. La stima è replicabile con lo script spreg/PySAL (canvas "SDM Cat-Nat — Script spreg") che incorpora il fix delle coordinate; i SE riportati derivano dalla Hessiana piena e coincidono con l'output `spreg.ML_Lag` sul modello Durbin.
- La stabilità di β_Grandi (Sez. 7) e l'invarianza della gerarchia SDM (Sez. 5) sono i due argomenti di robustezza da citare in conclusione.
- Da valutare in sede di revisione: (i) SE robusti per l'eteroschedasticità (BP significativo), (ii) termine quadratico o controllo dimensionale per il RESET.