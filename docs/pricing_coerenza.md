# Modello di pricing (benchmark EAL) e coerenza asset

Estensione attuariale dell'analisi: un **benchmark di loss annuale attesa (EAL, expected
annual loss)** a tre hazard per comune, confrontato con la tariffa IVASS per misurare la
**coerenza della tariffazione Cat-Nat rispetto al rischio modellato**, più un'analisi di
**coerenza della componente asset** che separa il rischio "puro" (la rate provinciale)
dalla parte meccanica del premio (asset × rate).

**Avvertenza di metodo.** Il benchmark usa parametri fisici *illustrativi* (ordini di
grandezza della letteratura tecnica, dichiarati in §2 e perturbati in §7): **non è un
modello di pricing operativo** e non sostituisce un cat model. Ciò che misura è la
coerenza **relativa** della tariffazione: dove la tariffa segue il rischio modellato e
dove no. Il loss ratio va letto come indicatore di adeguatezza relativa, non come
economic loss ratio di portafoglio.

Output macchina: `results/pricing_benchmark.json` (parametri, calibrazione, tabella
province, regressione), `results/eal_comuni.csv` (3.823 righe), `docs/mappa_loss.svg`
(mappa di loss a due pannelli), `results/esposizione_tessuto.json` (peso della loss
sul tessuto produttivo, §5). Tutto si rigenera con
`python3 scripts/pricing/pricing_model.py` seguito da
`python3 scripts/pricing/tessuto_produttivo.py` (solo stdlib, deterministici).

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
  (l'integrazione della curva di hazard sotto il design point: con danno ∝ ag², il solo
  punto RP475 sottostima l'AAL di ~4 volte);
- **frana**: `EAL_rate = 0,4% · share P3/P4` (damage-rate annuo atteso su area esposta);
- **idraulico**: `EAL_rate = 0,6% · share P3`.

Mix nazionale del benchmark (pesato per asset): **sismico 60,8%, idraulico 28,1%,
frana 11,2%** — l'ordine atteso per l'Italia, con il sismico dominante come nel modello
econometrico (§3.3 del README: l'omissione del sismico gonfiava β_Frana).

Aggregato nazionale: **EAL benchmark = 1,70 mld €/anno** su 1,43 mld € di asset
(rate media pesata 11,9 per 10.000 €), contro un premio teorico di 3,26 mld €/anno
(22,8 per 10.000 €): un fattore ~1,9 che assorbe insieme caricamenti, perils non
modellati (il perimetro DL 78/2025 va oltre sismico/frana/alluvione) e margine.

## 3. Calibrazione e loss ratio

Un **unico scalare di calibrazione** `c = 1,917` allinea la rate media pesata del
benchmark a quella IVASS (22,8 per 10.000 €): il benchmark calibrato e la tariffa
hanno lo stesso aggregato, e il confronto diventa puramente **relativo**.

`loss_ratio = rate_IVASS / rate_benchmark_calibrato` per comune
(1.920 comuni con ratio definito; 3 comuni sardi senza rischio modellato sono n/d):

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
(documentato in §7).

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

## 6. Mappa di loss

`docs/mappa_loss.svg` — due pannelli sugli stessi 3.823 centroidi comunali:

- **Loss attesa (EAL benchmark, EUR/anno)**: 8 classi a quantili, dal beige al rosso
  scuro; mostra dove si concentra la loss attesa (metropoli + aree ad alta
  pericolosità);
- **Loss ratio (coerenza tariffa/rischio)**: classi fisse centrare su 1 (blu = tariffa
  sotto il benchmark, rosso = sopra); mostra la geografia dell'adeguatezza relativa.

Generata dallo stesso script (nessuna dipendenza); visualizzabile direttamente su
GitHub.

## 7. Limiti e sensibilità

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
4. **Curva a un punto**: il fattore CURVE approssima l'integrazione della curva di
   hazard con un moltiplicatore costante; un cat model userebbe le percentili 16/84
   (già in matrice) per l'incertezza.

## 8. Riproducibilità

```bash
# dalla root del repo (solo stdlib; ~2 + ~1 secondi)
python3 scripts/pricing/pricing_model.py
# output: results/pricing_benchmark.json, results/eal_comuni.csv, docs/mappa_loss.svg
python3 scripts/pricing/tessuto_produttivo.py
# output: results/esposizione_tessuto.json (richiede results/eal_comuni.csv)
```

Script deterministici (nessun campionamento), percorsi configurabili via env
(`MATRICE`, `OUT_CSV`, `OUT_JSON`, `OUT_SVG` per pricing_model.py). Fonti: IVASS
(tariffe provinciali), INGV MPS04 (ag RP475), ISPRA (aree PAI P3/P4 e P3), AIDA
(asset per classe dimensione, aggregati a comune).
