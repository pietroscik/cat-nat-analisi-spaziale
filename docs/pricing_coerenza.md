# Modello di pricing (benchmark EAL) e coerenza asset

Estensione attuariale dell'analisi: un **benchmark di loss annuale attesa (EAL, expected
annual loss)** a tre hazard per comune, confrontato con la tariffa IVASS per misurare la
**coerenza della tariffazione Cat-Nat rispetto al rischio modellato**, più un'analisi di
**coerenza della componente asset** che separa il rischio "puro" (la rate provinciale)
dalla parte meccanica del premio (asset × rate).

**Avvertenza di metodo.** Il benchmark usa parametri fisici *illustrativi* (ordini di
grandezza della letteratura tecnica, dichiarati in §2 e perturbati in §8): **non è un
modello di pricing operativo** e non sostituisce un cat model. Ciò che misura è la
coerenza **relativa** della tariffazione: dove la tariffa segue il rischio modellato e
dove no. Il loss ratio va letto come indicatore di adeguatezza relativa, non come
economic loss ratio di portafoglio.

Output macchina: `results/pricing_benchmark.json` (parametri, calibrazione, tabella
province, regressione), `results/eal_comuni.csv` (3.823 righe), `docs/mappa_loss.svg`
(mappa di loss a due pannelli), `results/esposizione_tessuto.json` (peso della loss
sul tessuto produttivo, §5), `results/ep_curve.json` (AAL numerico e curva EP, §2.1),
`results/chi_paga.json` (ripartizione del premio, §5.4),
`results/spazializzazione_tariffa.json` + `docs/mappa_lisa_tariffa.svg` (Moran e
LISA della coerenza tariffaria, §6), `results/robustezza_tessuto.json` (verifiche
di robustezza, §5.5), `results/validazione_assunzioni.json` (validazione formale delle
assunzioni dichiarate, §8). **Quadro grafico**: `docs/grafici_pricing.svg` — sei pannelli
che visualizzano tutti i risultati precedenti (curva EP, quartili ISP, chi paga,
province estreme, robustezza, Moran), generato da `scripts/pricing/grafici.py`
leggendo i JSON: una vista dei dati, non un'origine. Tutto si rigenera con la
sequenza di comandi della §10 (solo
stdlib, deterministici).

## 1. Struttura del premio comunale

Il premio teorico comunale è `premio_10k_prov × asset_tot_EUR / 10.000`: la **rate è
costante entro provincia** (107 province), quindi tutta la variabilità comunale del
premio passa dagli asset. La decomposizione esatta del log-premio
(`log(premio) = log(asset) + log(rate) − log(10.000)`) sui 3.823 comuni:

| Componente | Quota di Var[log premio] |
|---|---|
| asset (dimensione esposta) | **95,9%** |
| rate (rischio provinciale) | 4,1% |
| correlazione asset–rate | −0,06 |

Il premio comunale è quindi quasi interamente un fatto di **esposizione**; la variabile
di pricing è la rate provinciale, ed è lì che va valutata la coerenza col rischio.

## 2. Benchmark EAL a tre hazard

Per ogni comune, rate di EAL per 10.000 € di asset (formalmente identiche alle tariffe
IVASS, così confrontabili):

- **sismico**: `EAL_rate(ag) = CURVE × min(1, MDR(ag)) / T`, con `T = 475` anni (design
  point 10% in 50), `MDR(ag) = min(1, 6·ag²)` (6,7% di danno al design point alla media
  nazionale pesata per asset, ag = 0,106; 37% a ag = 0,25 in Calabria) e `CURVE = 4`
  (punto di design + coda oltre il design point: con λ(ag) ∝ ag^−k, k = 3, il
  moltiplicatore è 1 + k/(k−2) = 4; la banda RP30–475 non è nel benchmark — §2.1 e §8);
- **frana**: `EAL_rate = 0,4% · share P3/P4` (damage-rate annuo atteso su area esposta);
- **idraulico**: `EAL_rate = 0,6% · share P3`.

Mix nazionale del benchmark (pesato per asset): **sismico 60,8%, idraulico 28,1%,
frana 11,2%** — l'ordine atteso per l'Italia, con il sismico dominante come nel modello
econometrico (§3.3 del README: l'omissione del sismico gonfiava β_Frana).

### 2.1 AAL numerico dalla curva di hazard e curva di eccedenza

Il fattore `CURVE = 4` vale quanto "punto di design + coda rara con k = 3"
(1 + k/(k−2) = 4): il benchmark integra gli eventi oltre il design point e lascia
fuori la banda RP30–475. La matrice estesa porta tre punti della curva MPS04 per
comune (ag a RP30, RP72, RP475), così l'AAL si stima per integrazione numerica
(`scripts/pricing/ep_curve.py`, `results/ep_curve.json`; correzione documentata in §8.1):

- **la banda RP30–RP475 vale da sola 4,53× il design point** (mediana per comune,
  p10–p90: 3,98–5,36; trapezoid sui punti noti, nessuna assunzione aggiuntiva): il 4
  assunto è numericamente vicino al moltiplicatore di banda, ma il moltiplicatore
  totale identificabile è **6,88** (banda + coda rara k = 3 con cap, p10–p90: 6,07–8,12);
- **AAL numerico sismico (RP≥30) = 1,71 mld €/anno** (trapezoid 1,14 + coda rara 0,58)
  contro 1,03 del sismico di benchmark: **1,66×** (sensibilità k = 2,5/3,5: 1,80/1,65
  mld, ratio 1,74/1,60). Il livello è assorbito dalla calibrazione c = 1,917 (§3) e la
  geografia del loss ratio è robusta a CURVE ×0,5/×2 (§9.1): è la misura di quanto il
  design point unico lascia fuori la banda RP30–475;
- **la coda frequente (RP<30) non è stimabile dai tre punti disponibili**: con qualunque
  legge di potenza λ(ag) ∝ ag^−k, k ≥ 2 l'integrale diverge, e la pendenza osservata del
  tratto RP30–RP72 è 3,70 (mediana) — ben oltre il confine di convergenza (§8.2). L'AAL
  numerico è quindi un **limite inferiore dichiarato** della parte RP≥30. La v1 di
  questa sezione riportava 2,30 mld € (2,2×) da una forma chiusa che era in realtà la
  coda rara da RP30 senza cap, col doppio conteggio della banda: valore ritirato,
  diagnosi completa in §8.1;
- **curva EP nazionale** (sismico puro, non calibrato): loss a scenario
  RP30 = 11,7 mld €, RP72 = 19,1, RP475 = **122,9 mld €** (banda epistemica 16/84:
  82–154 mld €): l'evento 1-in-475 vale **72 anni di AAL numerico** — il numero che
  dimensiona il rischio di coda rispetto al costo annuo atteso.

Aggregato nazionale: **EAL benchmark = 1,70 mld €/anno** su 1,43 mld € di asset
(rate media pesata 11,9 per 10.000 €), contro un premio teorico di 3,26 mld €/anno
(22,8 per 10.000 €): un fattore ~1,9 che assorbe insieme caricamenti, perils non
modellati (il perimetro DL 78/2025 va oltre sismico/frana/alluvione) e margine.

## 3. Calibrazione e loss ratio

Un **unico scalare di calibrazione** `c = 1,917` allinea la rate media pesata del
benchmark a quella IVASS (22,8 per 10.000 €): il benchmark calibrato e la tariffa
hanno lo stesso aggregato, e il confronto diventa puramente **relativo**.

`loss_ratio = rate_IVASS / rate_benchmark_calibrato` per comune
(3.820 comuni con ratio definito; 3 comuni sardi senza rischio modellato sono n/d):

| Statistica | Valore |
|---|---|
| mediana | 1,01 |
| p10 – p90 | 0,48 – 4,13 |
| comuni con ratio < 1 (tariffa sotto il benchmark) | 1.892 |
| comuni con ratio > 1 (tariffa sopra il benchmark) | 1.928 |
| comuni con ratio ≥ 1,5 | 1.222 (32%) |

La distribuzione è right-skewed per costruzione: nei comuni a basso hazard modellato
il benchmark tende a 0 mentre la tariffa ha un pavimento (~11 per 10.000 €), quindi
il ratio esplode — non è inefficienza, è il perimetro tariffario più ampio
(documentato in §9).

**Province estreme** (aggregato asset-weighted; tabella completa nel JSON):

| Sovraprezzo relativo (ratio alto) | | Sottoprezzo relativo (ratio basso) | |
|---|---|---|---|
| Agrigento | 4,59 | Treviso | 0,52 |
| Lecce | 3,53 | Udine | 0,51 |
| Enna | 3,51 | Grosseto | 0,51 |
| Monza e della Brianza | 3,41 | Forlì-Cesena | 0,49 |
| Rovigo | 3,29 | Rimini | 0,49 |
| Trapani | 3,14 | Pordenone | 0,47 |
| Biella | 3,12 | Valle d'Aosta | 0,28 |

Lettura: le province a ratio alto sono quelle a **basso sismico** (tariffa ≈ pavimento
di caricamento); quelle a ratio basso sono le province con forte **idraulico/frana non
riflesso nella tariffa** (prealpi venete e romagne, valle dell'Aosta) oppure sismico
alto ma tariffa contenuta (Vibo Valentia 0,45).

**EAL assoluta per provincia** (dove si concentra la loss attesa): Roma 144 M€/anno,
Bologna 120, Milano 88, Brescia 69, Verona 61, Napoli 61, Vicenza 60, Treviso 59
(l'esposizione pesa quanto l'hazard: la mappa di loss è dominata dalle aree
metropolitane).

## 4. Coerenza asset: la rate provinciale come variabile di pricing

Regressione OLS a livello provincia (n = 107) della log-rate IVASS sugli hazard pesati
per asset (log1p), con errori standard classici e **robusti (White HC1)**:

| Regressore | β | t | t robusto |
|---|---|---|---|
| log1p(ag medio pesato) | **+5,385** | +15,5 | **+16,7** |
| log1p(idraulico pesato) | +0,987 | +2,4 | +2,6 |
| log1p(frana pesato) | +0,203 | +0,7 (n.s.) | +0,7 (n.s.) |

R² = 0,708. **La tariffazione traccia la pericolosità sismica con elasticità ~5,4** —
più che quadratica rispetto all'ag, cioè più forte di quanto un danno ∝ ag²
comporterebbe — risponde debolmente all'idraulico (+1,0) e **non risponde al rischio
frana** (β ≈ 0,2, non significativo). È la conferma asset-free del risultato del SDM:
il premio è un fatto sismico, e il canale frana/idraulico della tariffazione è
sottodimensionato rispetto al benchmark (§3: ratio < 0,5 nelle province a forte
esposizione idrogeologica).

Confronto col modello econometrico: il SDM comunale (§3.2 del README) ha y e X che
condividono la componente asset; qui la rate provinciale la esclude per costruzione
(ed è costante entro provincia, quindi un SDM comunale sulla rate sarebbe degenere).
I due risultati convergono sul contenuto sostantivo (sismico dominante, frana
sottopagato), da angoli indipendenti.

## 5. Esposizione delle imprese performanti e peso della loss sul tessuto produttivo

Risposta attuariale alla domanda: **quanto pesa la loss Cat-Nat attesa sul valore
prodotto, e quanto sono esposte a questo peso le imprese performanti?** Lo script
`scripts/pricing/tessuto_produttivo.py` collega la EAL calibrata (§3) agli indicatori
del tessuto produttivo comunale della matrice estesa: `ISP_std_medio` (l'indicatore
composito di performance della tesi), EBITDA, dipendenti, numero imprese, struttura
dimensionale degli asset. Output: `results/esposizione_tessuto.json`.

### 5.1 Peso della loss sul tessuto

| Indicatore nazionale | Valore |
|---|---|
| EAL calibrata / EBITDA | **2,54%** (3,26 mld € su 128,5 mld €) |
| EAL per addetto | 895 €/anno |
| EAL per impresa | 106.408 €/anno |
| peso comunale (EAL/EBITDA): mediana / p90 / max | 2,44% / 7,50% / 70,7% |

Il peso comunale è fortemente right-skewed: la coda (max 70,7%, comune in Appennino
sismico con EBITDA piccolo) è dove un evento singolo più impattante può deprimere
per anni la capacità di ricostruzione del tessuto locale.

### 5.2 Esposizione per classe di performance (quartili ISP)

| Classe ISP | n comuni | rate bmk pesata | peso EAL/EBITDA | quota EAL / quota EBITDA |
|---|---|---|---|---|
| Q1 (meno performanti) | 949 | 30,74 | **3,91%** | 1,54 |
| Q2 | 950 | 27,67 | 3,07% | 1,21 |
| Q3 | 950 | 19,28 | 2,21% | 0,87 |
| Q4 (più performanti) | 950 | 25,02 | **2,06%** | **0,81** |

L'**intensità di esposizione** (quota di EAL nazionale / quota di EBITDA nazionale)
dice quanto pesa la loss di una classe rispetto al valore che quella classe produce:
i comuni con imprese *meno* performanti assorbono una quota di loss quasi **2 volte**
la loro quota di EBITDA (1,54 vs 0,81) — il gradiente è monotono su tutti e quattro i
quartili. I 77 comuni del cluster High-High dell'ISP espongono asset a rate medie
leggermente sopra il resto (24,53 vs 22,72 per 10.000 €: sono aree produttive
dell'Italia centro-meridionale con sismicità significativa), ma il loro peso-EBITDA (2,61%) resta vicino alla mediana
nazionale: l'alta produzione diluisce la loss.

Le correlazioni dicono che **la performance è geograficamente ortogonale al rischio
fisico** (ISP vs ag: +0,04; vs share frana/idraulico: ≈ 0): il gradiente di peso non
è "le imprese brave stanno dove non c'è rischio", è un fatto di composizione del
valore prodotto (EBITDA per asset più alto nei comuni performanti, stessa geografia
dell'hazard).

### 5.3 Modello predittivo del peso della loss

Regressione OLS del log-peso-EBITDA (n = 3.791 comuni; esclusi 24 senza ISP, 5 con
EBITDA ≤ 0 e 3 sardi con EAL calibrata = 0), SE robusti White HC1,
R² = 0,762:

| Regressore | β | t robusto |
|---|---|---|
| ISP (performance) | **−0,375** | **−21,8** |
| log1p(ag) | +14,68 | +71,9 |
| share frana | +2,76 | +29,1 |
| share idraulico | +4,30 | +25,1 |
| log(asset totale) | −0,038 | −4,7 |
| quota Grandi | +0,080 | +2,2 |

Lettura attuariale: **a parità di rischio fisico, dimensione e struttura dimensionale,
un incremento di 1 deviazione standard dell'ISP riduce il peso della loss del ~31%**
(e^−0,375 ≈ 0,69). La performance predice il peso della loss perché determina il
denominatore (valore prodotto) su una geografia dell'hazard che resta la stessa; gli
hazard confermano i segni attesi (sismico dominante, idraulico forte, frana presente),
la dimensione media degli asset ha effetto diluizione (−0,04 log-point per raddoppio),
e la quota di imprese Grandi ha un piccolo premio positivo (+0,08, t = 2,2): i comuni
a struttura dimensionale maggiore espongono asset più concentrati, e la loss pesa
leggermente di più — un segnale di **concentration risk** a livello comunale, marginale
ma statisticmente distinguibile.

A livello provinciale le posizioni estreme del peso-EBITDA sono **Vibo Valentia
(8,3%), Isernia (8,2%), Avellino (8,1%), Cosenza (7,8%)** — sismico calabro-lucano e
irpino su EBITDA bassi — contro **Monza e della Brianza (0,38%), Lecce (0,40%),
Brindisi (0,41%)** (tabella completa delle 107 province nel JSON, con la EAL per
addetto). Lo Spearman ISP–peso a livello comunale è −0,18: la relazione negativa
esiste ed è robusta, ma è un gradiente, non una legge.

### 5.4 Chi paga il premio: PMI vs Grandi imprese

La rate è costante entro provincia e si applica all'asset esposto, quindi la
ripartizione del premio segue gli asset (`scripts/pricing/chi_paga.py`,
`results/chi_paga.json`):

| Indicatore nazionale | Valore |
|---|---|
| Grandi imprese: quota delle imprese / quota del premio | 10,5% (3.648) / **63,7%** (2,08 mld €) |
| PMI: quota delle imprese / quota del premio | 89,5% (27.026) / 36,3% (1,18 mld €) |
| premio medio per impresa PMI | 43.884 €/anno |
| premio medio per impresa Grande | 569.617 €/anno (**×13** la PMI) |
| EAL calibrata su asset PMI / Grandi | 39,9% / 60,1% |
| incidenza comunale mediana: premio/EBITDA · premio/ricavi | 2,33% · 0,22% |

**Il 67% dei comuni (2.564 su 3.823) non ha Grandi imprese**: nei terzili di
quota_grandi dei comuni che le hanno, l'incidenza mediana non condizionata scende
da 2,41% (senza Grandi) a 2,15% (terzile alto) — i comuni dominati da Grandi
producono più EBITDA per asset e diluiscono l'incidenza. Nota di metodo: è una
media non condizionata, complementare al coefficiente **condizionato** +0,08 della
regressione del §5.3 (a parità di hazard, asset e dimensione la quota Grandi ha un
piccolo premio di esposizione): i due fatti convivono perché descrivono cose
diverse — il livello medio dell'incidenza e il suo gradiente marginale. L'incidenza
per classe non è calcolabile (EBITDA e ricavi in matrice sono totali di comune):
limite dichiarato.

### 5.5 Robustezza del risultato: il β della performance sopravvive a tutto

Cinque verifiche della specifica base, senza selezioni ex post
(`scripts/pricing/robustezza_tessuto.py`, `results/robustezza_tessuto.json`):

| Specifica | n | R² | β performance | t robusto |
|---|---|---|---|---|
| base (§5.3, ISP) | 3.791 | 0,762 | **−0,375** | −21,8 |
| winsorizzato log-peso 1%/99% | 3.791 | 0,777 | −0,361 | −23,1 |
| trim top 1% del peso | 3.754 | 0,763 | −0,359 | −22,4 |
| performance = ROA comunale (EBITDA/asset) | 3.791 | 0,844 | −10,51 (−0,41 per SD) | −44,7 |
| SLX (base + regressori col ritardo spaziale W) | 3.791 | 0,763 | −0,376 | −21,9 |

- il gradiente **non è portato dalla coda** del peso-EBITDA (winsorizzato e trim
  danno lo stesso β del base);
- **non è un fatto dell'ISP in sé**: con il ROA comunale — una misura di
  performance che non dipende dalla pipeline della tesi (correlazione con ISP:
  0,54) — segno e ordine di grandezza per deviazione standard reggono (−0,41 vs
  −0,25 per SD dell'ISP), con R² ancora più alto;
- **l'effetto è tutto locale**: nel SLX il ritardo spaziale della performance
  (W×ISP) è nullo (t = +0,1) — niente effetto di contesto, il peso della loss di un
  comune dipende dalla performance delle *sue* imprese, non da quelle dei vicini;
- i residui della specifica base sono autocorrelati (Moran I = 0,231, z = +24,6):
  l'OLS trascura una componente spaziale dell'hazard (geografia liscia), che però
  non sposta il coefficiente della performance (SLX: β invariato). Diagnostica
  documentata, non risolta con un SAR: la lettura sostantiva non cambia.

## 6. Spazializzazione della coerenza tariffaria: Moran e LISA del loss ratio

Il ponte tra i due mondi del repo: la geografia dell'adeguatezza tariffaria letta
con gli strumenti spaziali del SDM (`scripts/pricing/spazializzazione.py`,
`results/spazializzazione_tariffa.json`, mappa `docs/mappa_lisa_tariffa.svg`).
W: KNN k = 5 sui centroidi (come il SDM comunale), pesi di riga standardizzati;
Moran globale con statistica analitica Cliff–Ord **e** permutazioni (999, seed
fisso: deterministico); LISA con randomizzazione condizionata (499, p < 0,05).

| Variabile | n | Moran I | z | p (perm) |
|---|---|---|---|---|
| log loss ratio | 3.820 | **+0,649** | +69,2 | 0,001 |
| loss ratio grezzo | 3.820 | +0,0012 | +6,6 | 0,007 |
| log peso-EBITDA | 3.815 | +0,669 | +71,3 | 0,001 |
| rate benchmark calibrata | 3.823 | +0,763 | +81,3 | 0,001 |

- **l'adeguatezza tariffaria è fortemente clusterizzata** (I = 0,649 sul log): la
  distanza dal rischio modellato non è rumore comunale indipendente — è una
  struttura spaziale, coerente con il fatto che la rate è decisa a livello
  provinciale: comuni vicini condividono lo stesso errore di tariffazione. Il
  Moran del ratio grezzo (+0,001) è invece pilotato dalla coda (max 13.429,
  Soleminis): la trasformazione log è dichiarata ed è la lettura corretta;
- i cluster LISA significativi confermano la lettura della §3: **630 comuni HH**
  (sovraprezzo in cluster di sovraprezzo: Sardegna — dove il benchmark ≈ 0 e il
  ratio esplode — e Nord a basso sismico, effetto pavimento tariffario) contro
  **426 comuni LL** (sottoprezzo in cluster di sottoprezzo: Valle d'Aosta,
  Udine, Brescia, Treviso — l'idrogeologico non riflette in tariffa);
- la rate benchmark è più autocorrelata del loss ratio (0,76 vs 0,65): il rischio
  fisico è più liscio nello spazio dell'errore tariffario — la differenza è
  esattamente il disallineamento che la §4 misura con l'elasticità.

La mappa `docs/mappa_lisa_tariffa.svg` mostra i 3.820 comuni classificati per
cluster (HH rosso, LL blu, outliers spaziali arancio/celeste, non significativi
grigio).

## 7. Mappa di loss

`docs/mappa_loss.svg` — due pannelli sugli stessi 3.823 centroidi comunali:

- **Loss attesa (EAL benchmark, EUR/anno)**: 8 classi a quantili, dal beige al rosso
  scuro; mostra dove si concentra la loss attesa (metropoli + aree ad alta
  pericolosità);
- **Loss ratio (coerenza tariffa/rischio)**: classi fisse centrare su 1 (blu = tariffa
  sotto il benchmark, rosso = sopra); mostra la geografia dell'adeguatezza relativa.

Generata dallo stesso script (nessuna dipendenza); visualizzabile direttamente su
GitHub.

## 8. Validazione delle assunzioni

`scripts/pricing/validazione_assunzioni.py` (`results/validazione_assunzioni.json`)
sottopone ogni assunzione dichiarata a verifica empirica o di coerenza interna:
nessuna nuova ipotesi, verdetti riportati come vengono (solo stdlib, deterministico).

### 8.1 Correzione della coda: diagnosi di un errore di forma chiusa

La v1 della §2.1 stimava la "coda frequente" (RP<30) con la forma chiusa
`(k/(k−2))·MDR(ag30)/30`. La quadratura numerica dimostra (err 4·10⁻¹⁰) che
quell'integrale **non è la coda frequente: è la coda rara da RP30 a infinito, senza
cap di MDR** — sommata al trapezoid ricontava la banda RP30–475 già integrata. Il
totale v1 (2,30 mld €/anno, "2,2×") è **ritirato**. La forma chiusa corretta (coda
rara da RP475, cap attivo: `k/(k−2)·m475/475 − (2/(k−2))·λ_cap`) è verificata contro
quadratura su comuni rappresentativi (p10/p50/p90/max di ag): errore relativo ≤ 10⁻⁶.

### 8.2 Identificabilità della coda frequente (RP<30)

Con `λ(ag) ∝ ag^−k` l'integrale della coda frequente converge solo per **k < 2**. Le
pendenze locali osservate della curva MPS04 (n = 3.753 comuni con hazard):

| Tratto | k mediana | p10–p90 | comuni con k ≤ 2 |
|---|---|---|---|
| RP30–RP72 | 3,70 | 2,97–4,59 | 0,0% |
| RP72–RP475 | 2,11 | 1,77–2,69 | 38,9% |
| RP30–RP475 (OLS su 3 punti) | 2,36 | 1,97–2,94 | 14,7% |

La quadratura della coda frequente col k = 3 dichiarato cresce senza limite al
crescere del cutoff di λ (×6,3 a λ = 1/anno, ×40 a λ = 100, ×929 a λ = 10⁶, e oltre):
l'estrapolazione è dominata dal punto — non osservato — in cui la curva reale si
appiattisce. **Verdetto: non stimabile dai tre punti; dichiarata. L'AAL numerico è
un limite inferiore della parte RP≥30.**

### 8.3 Coda rara k=3

Assunzione dichiarata, non risolta dai dati: la pendenza del tratto RP72–475
osservato è più piatta di 3 (2,11 mediana), quindi la direzione del bias è al ribasso
(k più basso → coda più grande), ma sotto k=2 la chiusa non esiste. La sensibilità
k = 2,5/3,5 sul totale resta contenuta (**1,80/1,65 mld**, ±5%) perché il cap di MDR
taglia l'estrapolazione oltre RP ~19.600 anni (mediana di saturazione).

### 8.4 MDR quadratico

L'elasticità implementata è esattamente **2,000** (R² = 1: l'implementazione realizza
l'assunzione dichiarata, nessun comune saturo a RP475); l'elasticità IVASS osservata è
**+5,4** (§4): la quadratica è conservativa rispetto a come la tariffa traccia l'ag,
ma la vulnerabilità reale resta non osservata (parametro illustrativo, §9.1).

### 8.5 Diagnostica OLS del tessuto (replica della §5.3)

| Test | Risultato | Lettura |
|---|---|---|
| Jarque-Bera | p ≈ 0 | residui non normali (coda pesante): dichiarato |
| Breusch-Pagan | p ≈ 0 | eteroschedasticità: SE robusti HC1 già in uso |
| VIF | max 1,8 | nessuna collinearità rilevante |
| Cook (D > 4/n) | 234 osservazioni | β_ISP −0,375 → −0,371 senza le influenti |
| RESET (ISP²+ISP³) | p = 0,002 | non-linearità lieve statisticamente rilevabile |
| Dummies quartili ISP | Q2 −0,31, Q3 −0,06, Q4 −0,24 | tutti sotto Q1, non perfettamente monotoni |

Il β della performance sopravvive a tutta la diagnostica; la forma non è
perfettamente lineare (RESET e quartili), quindi il β va letto come riassunto medio —
coerente con la sua presentazione in §5.3 e con la robustezza della §5.5.

### 8.6 Trasformazione log, matrice W e seed

- **log del loss ratio**: skewness 61,8 (grezzo) → 1,07 (log): la trasformazione
  dichiarata della §6 è quantificata formalmente;
- **W KNN**: Moran del log loss ratio stabile su k = 3/5/7/10 (I = 0,68/0,65/0,63/0,60,
  z = 57/69/79/90): la clusterizzazione non è un artefatto della specifica della W.
  La replica esatta della §6 (k=5, 999 permutazioni, seed 42) dà delta_I = 0;
- **seed**: p = 0,001 per i seed 42/123/2024 (floor di risoluzione con 999
  permutazioni): l'inferenza permutativa non dipende dal seed.

### 8.7 Sintesi dei verdetti

| Assunzione | Esito |
|---|---|
| `CURVE = 4` | RICALIBRATA: moltiplicatore totale identificabile 6,88 (banda 4,53 + coda rara); livello assorbito da c = 1,917, geografia robusta (§9.1) |
| coda rara k = 3 | DICHIARATA: chiusa verificata contro quadratura; bias al ribasso documentato; sensibilità ±5% |
| coda frequente k = 3 (v1) | RITIRATA: errore di forma chiusa (doppio conteggio della banda); non stimabile → AAL = limite inferiore |
| MDR = 6·ag² | COERENZA INTERNA VERIFICATA (elasticità 2,000, R² = 1); illustrativa |
| OLS §5.3 | VALIDATA CON DIAGNOSTICA: β non portato da coda/influenti; non-linearità lieve dichiarata |
| log del ratio | VALIDATA (skewness 61,8 → 1,07) |
| W KNN k = 5 | VALIDATA (I stabile su k 3–10; replica esatta di §6) |
| seed permutazioni | VALIDATA (p = 0,001 su tre seed) |

## 9. Limiti e sensibilità

1. **Parametri illustrativi**: MDR, fattore curva e frequenze-danno sono ordini di
   grandezza dichiarati, non stime. La sensibilità (ogni parametro ×0,5 e ×2, più tutti
   insieme ×0,5) mostra Spearman del ranking dei loss ratio comuni **≥ 0,96 in ogni
   variante** (0,96–0,98): la geografia dell'adeguatezza relativa è robusta ai
   parametri, entro il loro ordine di grandezza.
2. **Perimetro**: la tariffa DL 78/2025 copre perils oltre i tre modellati; il
   pavimento tariffario (~11 per 10.000 €) spiega i ratio alti del Nord a basso
   sismico. Il loss ratio non è una misura di efficienza della tariffa ma di distanza
   dal rischio modellato.
3. **Niente vulnerabilità reale**: MDR uniforme per ag (niente microzonazione — vedi
   `docs/sismico_metodologia.md` §7 — né vulnerabilità per tipologia edilizia);
   Sardegna: benchmark ≈ 0, ratio non interpretabile (mediana 2,5 solo come ordine).
4. **Coda e moltiplicatore**: CURVE=4 copre design point + coda rara k=3 e lascia
   fuori la banda RP30–475 (4,53×), ora stimata separatamente (§2.1); la coda frequente
   RP<30 resta non identificabile dai tre punti (§8.2): l'AAL numerico è un limite
   inferiore; un cat model userebbe le percentili 16/84 (già in matrice) per l'incertezza.

## 10. Riproducibilità

```bash
# dalla root del repo (solo stdlib; ~2 + ~1 + ~1 + ~1 + ~1 + ~2 minuti)
python3 scripts/pricing/pricing_model.py
# output: results/pricing_benchmark.json, results/eal_comuni.csv, docs/mappa_loss.svg
python3 scripts/pricing/tessuto_produttivo.py
# output: results/esposizione_tessuto.json (richiede results/eal_comuni.csv)
python3 scripts/pricing/ep_curve.py
# output: results/ep_curve.json (AAL numerico e curva EP §2.1)
python3 scripts/pricing/chi_paga.py
# output: results/chi_paga.json (ripartizione del premio §5.4)
python3 scripts/pricing/spazializzazione.py
# output: results/spazializzazione_tariffa.json + docs/mappa_lisa_tariffa.svg (§6)
python3 scripts/pricing/robustezza_tessuto.py
# output: results/robustezza_tessuto.json (verifiche §5.5; richiede i due script sopra)
python3 scripts/pricing/validazione_assunzioni.py
# output: results/validazione_assunzioni.json (validazione delle assunzioni §8; richiede
# ep_curve.py, chi_paga.py e spazializzazione.py sopra)
python3 scripts/pricing/grafici.py
# output: docs/grafici_pricing.svg (quadro grafico; richiede tutti i JSON sopra)
```

Script deterministici (nessun campionamento; le permutazioni di Moran/LISA usano
un seed fisso, §6), percorsi configurabili via env (`MATRICE`, `OUT_CSV`,
`OUT_JSON`, `OUT_SVG` per pricing_model.py). Fonti: IVASS (tariffe provinciali),
INGV MPS04 (ag RP475/72/30 + percentili, Sa(0,10 s)), ISPRA (aree PAI P3/P4 e P3),
AIDA (asset per classe dimensione, aggregati a comune).
