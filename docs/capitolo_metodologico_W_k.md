# Capitolo metodologico — Dalla matrice dei dati alla domanda di ricerca: W, k e specificazione dei modelli

Questo capitolo esplicita il processo analitico nella sua sequenza logica — costruzione della
matrice dei dati, definizione della matrice di pesi spaziali W, scelta di k, specificazione dei
modelli — mantenendo il nesso tra i passaggi: ogni scelta successiva è condizionata dalle
precedenti e la sua correttezza è misurabile solo rispetto alle fonti e alle correlazioni
accademiche richiamate. L'oggetto finale dello studio **non è la performance d'impresa**, ma la
relazione tra **premi assicurativi Cat-Nat** e **copertura del rischio idrogeologico** a livello
territoriale (§5); la performance (ISP) entra solo come variabile esplorativa dell'analisi di
evoluzione a livello impresa (§3.2).

## 1. La matrice dei dati: specifica e costruzione

### 1.1 Unità di analisi e dataset

Le due analisi del repository poggiano su **due matrici distinte**, generate da un'unica pipeline
di allineamento delle fonti:

| Matrice | Unità | Dimensione | File |
|---|---|---|---|
| Livello impresa | impresa attiva georeferenziata | 30.673 × variabili di bilancio | (ambiente R della tesi) |
| Livello comune | comune italiano | **3.823 × 37** | `data/Matrice_Modello_Savelli_Final.csv` |

Fonti e ruolo di ciascuna (pipeline di costruzione in ordine di integrazione):

1. **AIDA (Bureau van Dijk)** — bilanci d'impresa: ROI, ROE, ROS, integrazione verticale,
   rendimento dei dipendenti, ricavi (log), posizione finanziaria netta, patrimonio netto,
   utile netto, indici di liquidità e di indipendenza finanziaria, durate medie crediti/debiti.
   Sono i regressori dell'analisi a livello impresa (§3) e la base per gli asset esposti.
2. **ISPRA — Carta della pericolosità da frana** — classi di pericolosità per comune, incrociate
   con gli **asset esposti** (PMI vs Grandi imprese) per generare le due variabili chiave
   `Risk_Frana_Asset_PMI` e `Risk_Frana_Asset_Grandi` a livello comune.
3. **IVASS — elaborazione Cat-Nat** — tariffe premiali teoriche per le 110 province
   (Decreto-Legge 78/2025), ripartite sui comuni per generare `Premio_Teorico_Comunale_EUR`.
4. **ISTAT** — georeferenziazione (centroidi comunali) e ricodifica territoriale.

### 1.2 Il processo di allineamento (con correzioni documentate)

Il passaggio da fonti eterogenee a matrice unica richiede decisioni che vengono qui rese
esplicite perché **condizionano tutti i risultati a valle**:

- **ricodifica delle fusioni comunali** (fonti con codici ISTAT di epoche diverse):
  024128→024103 (Sovizzo), 075098→075062 (Presicce-Acquarica), 048054→048052
  (Figline e Incisa Valdarno), 081025→081021 (Trapani/Misiliscemi);
- **Sardegna**: 70 comuni con tariffa provinciale ripartita secondo l'assetto provinciale 2025
  (flag `premio_provincia_appross` in matrice: segnala osservazioni con dipendente approssimata);
- **correzione dei centroidi** anomali/duplicati: **Gazzo (VI), Lucignano (AR), Olgiate Olona
  (VA), Telese Terme (BN)**. La correzione è *sensibile*: senza di essa questi comuni entrano
  nel vicinato KNN di unità distanti, inquinando i ritardi spaziali; tutti i risultati definitivi
  si riferiscono alla matrice corretta (verifica di sensibilità: commit "Matrice definitiva
  comuni: 3.823 × 37, coordinate corrette").

Le 37 colonne della matrice comune coprono: identificativi e geografia (codice ISTAT, provincia,
macroarea), il premio teorico comunale, le due esposizioni al rischio frana per classe di impresa,
e le variabili di controllo territoriali utilizzate nelle specificazioni della tesi.

## 2. La matrice di pesi spaziali W: definizione e giustificazione

Ogni modello di econometria spaziale (Anselin, 1988; LeSage & Pace, 2009; Elhorst, 2014) richiede
la pre-specificazione di una **matrice di pesi W** che formalizza il vicinato tra unità. W è una
scelta metodologica, non un risultato di stima: condiziona ρ e θ e la classificazione dei modelli
(SAR/SEM/SDM). Si adotta la struttura **KNN (k-nearest neighbours)** su distanza euclidea tra
coordinate, con standardizzazione per riga (style "W"), per tre ragioni:

- **nessuna unità isolata**, rilevante a livello impresa dove i punti sono disomogenei e i
  duplicati di coordinate frequenti (risolti con jitter 1e-4);
- **densità costante per riga** (ogni unità ha esattamente k legami): rende confrontabili le
  stime tra subset;
- il ritardo spaziale Wy è una **media locale**, interpretabile come contesto territoriale
  dell'unità i (LeSage & Pace, 2009, cap. 2).

## 3. Livello impresa: k = 77 e il ruolo dell'ISP

### 3.1 Scelta di k

Criterio combinato su due strumenti: **curva k-distanza** (media/mediana/sd della distanza al
k-esimo vicino; il "gomito" separa vicini omogenei da vicini arbitrari — criterio di tipo
k-nearest di Cover & Hart applicato alla selezione della banda) e **stabilità dell'I di Moran**
(Moran, 1950) della variabile target per k crescente. L'intersezione seleziona **k = 77**
(log `analisi_k_nazionale.txt`); matrice risultante: 2.361.821 legami, 0,251% pesi non nulli,
S0 = 30.673, S1 = 699,86, S2 = 125.514,1. Per i 20 subset Dimensione × Macroarea si ripete il
criterio con k* specifici (16–74; tabella in `docs/log_R_livello_impresa.md`).

### 3.2 Il contesto interpretativo dell'ISP

La variabile esplorata a livello impresa è l'**ISP — Indicatore Sintetico di Performance**,
indicatore composito sviluppato nella tesi. Costruzione e validazione sono documentate e
replicabili nel repository della tesi:
**[github.com/pietroscik/tesi-magistrale](https://github.com/pietroscik/tesi-magistrale)** —
in particolare `suddivisione_script/divisione e creazione ISP.R` (FASE 2: feature engineering e
calcolo), `suddivisione_script/ISP validazione e inferenza.R` (FASE 3: validazione e
finalizzazione) e gli output `03_validation/` (pesi e statistiche delle versioni ISP).

**Costruzione di base.** Le variabili di bilancio AIDA, organizzate in categorie (redditività,
solidità, produttività, liquidità, capitale circolante, rischio finanziario), sono normalizzate
min-max su scala 0–1000. L'ISP aggrega due sottogruppi con media ponderata:

| Sottogruppo | Componenti (peso) |
|---|---|
| **A — Redditività/Performance** | ROE (0,2727), EBITDA_su_vendite (0,3560), ROI (0,2386), rotazione_cap_investito (0,1327) |
| **B — Patrimoniale/Finanziario** | debt_equity_ratio (0,3162), debt_EBITDA_ratio (0,2703), totale_attivita (0,1583), PFN_EBITDA (0,2552) |

con pesi di aggregazione tra gruppi 0,4311 (A) + 0,5689 (B).

**Validazione.** La FASE 3 confronta strategie di pesatura per sottogruppo — pesi PCA (prima
componente), pesi da regressione (LM) e **LASSO** (`cv.glmnet`, lambda selezionato con
cross-validation), con i pesi tra gruppi stimati come quota di R² della regressione
(peso_A = R²_A/(R²_A+R²_B)). La versione definitiva è l'**ISP settoriale LASSO**
(`ISP_sett_lasso_norm`), con pesi specifici per sezione ATECO salvati in
`03_validation/risultati_settoriali_pesi.csv` (es. C-Manifatturiero: peso A = 0,7781;
K-Attività finanziarie: peso A = 0,5969): la performance è quindi pesata in modo diverso a
seconda del settore di attività, coerentemente con la letteratura degli **indicatori
compositi** (OECD/JRC, *Handbook on Constructing Composite Indicators*, 2008).

**Finalizzazione.** La versione usata in questo repo è **`ISP_bn`**: trasformazione robusta di
normalizzazione `bestNormalize`/OrderNorm (Peterson & Cavanaugh, 2019) dell'ISP settoriale
finale, scelta per avvicinare la distribuzione alla normalità richiesto dai test spaziali;
`ISP_std` e `lag_ISP_std` sono la standardizzazione e il ritardo spaziale usati nello
scatterplot di Moran.

**Conseguenze interpretative.** L'ISP è un **indice ordinale e relativo**: (i) i suoi valori
dipendono dalla normalizzazione min-max sul campione osservato e dai pesi stimati, quindi non
hanno unità di misura né lettura assoluta; (ii) confronti tra imprese di settori diversi
ereditano i pesi settoriali; (iii) essendo funzione monotona dei componenti di bilancio, le sue
correlazioni non aggiungono informazione causale oltre a quelle. Perciò l'analisi esplorativa
spaziale su ISP_bn (I di Moran = 0,00865, p < 2,2e-16; cluster LISA, Anselin 1995) è
**descrittiva**: dice *dove e quanto* la performance relativa è spazialmente associata, non
*perché*. Due cautele operative:

1. l'attribuzione causale richiede le regressioni di §3.3 (dove l'ISP è la dipendente e i
   fondamentali di bilancio i regressori) e il contesto della tesi;
2. i cluster LISA su un indice composito vanno letti con la prudenza imposta dal confronto
   multiplo: dopo correzione **FDR (Benjamini–Hochberg, 1995)** nessun cluster locale
   sopravvive (`robustezza_finale_fdr_nazionale.txt`) — esito che limita le conclusioni
   esplorative a livello di singola impresa e motivava la transizione al livello comunale
   (§4–5), dove la domanda di ricerca sui premi e le coperture trova risposta robusta.

### 3.3 Regressioni per subset

Per ogni subset Dimensione × Macroarea si stimano OLS, SAR(err), SDM(mix), GMM(err)
(`spatialreg`), con selezione per AIC e diagnostica (Moran residui, Breusch–Pagan 1979,
Anderson–Darling, Ramsey RESET 1969). Esiti di sintesi in `docs/log_R_livello_impresa.md`:
il modello spaziale è preferito nei subset Piccola/Media (λ = 0,19–0,37), l'OLS nei subset
Micro/Grande; le stime ρ nei subset piccoli con W mal condizionata sono inaffidabili
(ρ ≈ −3,4). Coefficienti stabili: **ROI +**, integrazione verticale −, ROE +.

## 4. Livello comune: k = 5

La qualità della KNN dipende dalle coordinate: prima di costruire W si applicano le correzioni
del §1.2. Il criterio combinato k-dist + Moran (`scripts/step1_kdist_moran.py`,
`results/grid_results.json`) mostra il gomito a k piccoli (distanza mediana 0,088° a k = 5,
0,114° a k = 8) e parametri sostantivi stabili sulla griglia k = 5–8 (β_Grandi = 0,142–0,144;
θ_PMI = −0,042/−0,050). **k = 5** per: densità coerente con la scala comunale (5–6 confinanti
stricti), parsimonia, e perché ρ cresce monotono con k assorbendo la densità di W (il guadagno
di AIC non riflette struttura aggiuntiva). Il confronto con k = 77 dell'analisi impresa riflette
il rapporto tra le n e la dispersione puntuale.

## 5. La domanda di ricerca: premi e copertura del rischio

### 5.1 Enunciato

L'interrogativo genuino del repository è: **la tariffazione Cat-Nat cattura l'esposizione al
rischio idrogeologico degli asset produttivi, e con quale geografia?** In formula: il premio
teorico comunale (proxy dell'offerta di copertura, DL 78/2025) risponde al rischio frana che
pesa sugli asset PMI e Grandi imprese, e questa risposta ha componenti locali e di vicinato?

### 5.2 Perché un modello spaziale è la specifica corretta

1. i premi nascono da tariffe **provinciali** ripartite sui comuni: per costruzione i comuni
   della stessa provincia condividono componente di premio → dipendenza spaziale attesa, da
   verificare con ρ (l'ipotesi è testabile, non assunta: Moran residui e LR la verificano);
2. il rischio idrogeologico è **spazialmente correlato per natura fisica**: un evento frana
   non rispetta i confini comunali, quindi l'esposizione dei comuni vicini contiene informazione
   rilevante per il premio locale (termini WX);
3. l'assenza di queste componenti produce stime biased e inconsistenti (omissione di
   dipendenza spaziale, Anselin 1988), con effetti diretti e indiretti non separabili.

Da qui la **SDM** come specifica di partenza, perché nests SAR e SEM (LeSage & Pace, 2009):
se i ritardi WX non servono, la SDM collassa in SAR (test LR su θ = 0); se vale la comune
radice θ = −ρβ, collassa in SEM. L'evidenza definitiva (§5.3) conferma la SDM: **ρ = 0,4966***
(z = 22,97), **β_Grandi = 0,1434*** (z = 37,05), θ_PMI = −0,052***, LR vs SAR = 126,71***,
LR vs SEM = 109,64***, AIC: SDM < SEM < SAR; effetti LeSage–Pace: Grandi diretto +0,151 /
totale +0,211; PMI totale −0,091 (indiretto dominante). Dettaglio: `docs/risultati_sdm_comuni.md`.

### 5.3 Lettura assicurativa (per il capitolo dei risultati)

- il premio è **fortemente spaziale** (ρ ≈ 0,5): tariffare comune per comune ignorando il
  vicinato sottostima la co-movimento dei premi; per il regolatore è argomento a favore di
  tariffe coordinate a scala sovracomunale;
- **le Grandi imprese sono il canale rischio→prezzo**: l'esposizione delle Grandi imprese al
  rischio frana incrementa il premio (elasticità totale ~0,21), suggerendo che il mercato
  Cat-Nat prezza soprattutto i grandi asset concentrati;
- **le PMI entrano solo per via indiretta e con segno negativo** (−0,09): i comuni con alta
  esposizione PMI nei dintorni sono associati a premi locali minori; possibile lettura
  assicurativa: l'esposizione diffusa e policentrica delle PMI è meno captata dalla tariffazione,
  un tema di *insurability* delle PMI rispetto al rischio idrogeologico;
- implicazione per le **coperture**: la differenza di sensibilità del premio alle due classi
  di asset indica un potenziale divario tra rischio effettivo portato dalle PMI e prezzo della
  copertura disponibile per esse.

## 6. Limiti dichiarati

- Moran residuo −0,049 (p = 0,004): lieve autocorrelazione negativa di corto raggio, tipica di
  KNN con k basso;
- RESET F = 24,4: forma funzionale migliorabile nel log1p; segni e gerarchia stabili;
- Breusch–Pagan LM = 124,4: eteroschedasticità → SE ML da considerare conservativamente;
  raccomandata stima robusta (White/GMM spaziale) come sviluppo;
- 70 comuni sardi con premio approssimato (flag in matrice);
- i risultati del livello impresa su ISP non sono confrontabili direttamente con quelli del
  livello comune (unità, variabile dipendente e scopi diversi): il primo descrive la
  performance, il secondo risponde alla domanda sui premi e le coperture.
