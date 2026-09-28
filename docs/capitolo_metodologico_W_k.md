# Capitolo metodologico — Dalla matrice dei dati alla domanda di ricerca: W, k e specificazione dei modelli

Questo capitolo esplicita il processo analitico nella sua sequenza logica — costruzione della
matrice dei dati, **transizioni di scala** (impresa→comune e provincia→comune), definizione
della matrice di pesi spaziali W, scelta di k, specificazione dei modelli — mantenendo il nesso
tra i passaggi: ogni scelta successiva è condizionata dalle precedenti e la sua correttezza è
misurabile solo rispetto alle fonti e alle correlazioni accademiche richiamate. L'oggetto
finale dello studio **non è la performance d'impresa**, ma la relazione tra **premi assicurativi
Cat-Nat** e **copertura del rischio idrogeologico** a livello territoriale (§5); la performance
(ISP) entra solo come variabile aggregata di contesto (§1.3) e come variabile dell'analisi di
evoluzione a livello impresa, il cui studio è già stato completato nella tesi (§3.2).

## 1. La matrice dei dati: specifica e costruzione

### 1.1 Unità di analisi e dataset

Le due analisi del repository poggiano su **due matrici distinte**, generate da un'unica pipeline
di allineamento delle fonti:

| Matrice | Unità | Dimensione | File |
|---|---|---|---|
| Livello impresa | impresa attiva georeferenziata | 30.673 × variabili di bilancio | (repo tesi-magistrale) |
| Livello comune | comune italiano | **3.823 × 37** | `data/Matrice_Modello_Savelli_Final.csv` |

La matrice a livello comune **riassume** quella a livello impresa e **arricchisce** con i dati
territoriali (ISTAT/ISPRA) e assicurativi (IVASS): senza esplicitare *come* si passa di scala la
matrice resterebbe una specificazione mancante — le variabili comunali apparirebbero come dati
primitivi senza regola di costruzione, e ogni risultato a valle (W, k, SDM) sarebbe
non riproducibile. Le due transizioni sono documentate al §1.3.

Fonti e ruolo di ciascuna (pipeline di costruzione in ordine di integrazione):

1. **AIDA (Bureau van Dijk)** — bilanci d'impresa: attivo totale, EBITDA, ricavi, dipendenti,
   integrazione verticale, redditività. Ricondotti ai comuni via §1.3a.
2. **ISPRA — Piano di Assetto Idrogeologico (PAI)** e aree a pericolosità idraulica — classi
   P3/P4 (frana) e P3 (idraulico) per comune: superfici e popolazione esposte.
3. **IVASS — elaborazione Cat-Nat** — tariffe premiali teoriche per le 110 province
   (Decreto-Legge 78/2025), espresse **per 10.000 € di asset esposto**: ripartite sui comuni
   via §1.3b.
4. **ISTAT** — georeferenziazione (centroidi comunali, confini `Com01012025_WGS84`),
   popolazione 2018, superfici, ricodifica territoriale.

### 1.2 Il processo di allineamento (con correzioni documentate)

Il passaggio da fonti eterogenee a matrice unica richiede decisioni che vengono qui rese
esplicite perché **condizionano tutti i risultati a valle**:

- **ricodifica delle fusioni comunali** (fonti con codici ISTAT di epoche diverse):
  024128→024103 (Sovizzo), 075098→075062 (Presicce-Acquarica), 048054→048052
  (Figline e Incisa Valdarno), 081025→081021 (Trapani/Misiliscemi);
- **Sardegna**: 70 comuni con tariffa provinciale ripartita secondo l'assetto provinciale 2025
  (flag `premio_provincia_appross`);
- **correzione dei centroidi** anomali/duplicati: **Gazzo (VI), Lucignano (AR), Olgiate Olona
  (VA), Telese Terme (BN)**. La correzione è *sensibile*: senza di essa questi comuni entrano
  nel vicinato KNN di unità distanti, inquinando i ritardi spaziali; tutti i risultati definitivi
  si riferiscono alla matrice corretta.

### 1.3 Le transizioni di scala: come si arriva al comune

La matrice comunale non è un rilievo primario: ogni colonna deriva da una regola di
aggregazione esplicita da un livello di osservazione superiore o inferiore. Senza queste regole
la matrice soffre di **specificazione mancante**; con esse, ogni colonna è tracciabile alla
fonte e riproducibile.

#### a) Transizione impresa → comune (aggregazione bottom-up)

*Assegnazione*: ogni impresa AIDA georeferenziata è assegnata al comune di appartenenza tramite
**spatial join** dei punti-coordinata con i poligoni comunali ISTAT
(`Com01012025_WGS84.shp`, `st_join` con predicato `st_within` — stessa macchina della FASE 2 del
repo [tesi-magistrale](https://github.com/pietroscik/tesi-magistrale), usata lì per la
correzione delle macroaree). Le imprese sono classificate per dimensione secondo i criteri UE
2023 (Micro/Piccola/Media/Grande, su dipendenti, ricavi e attivo — vedi `divisione e creazione
ISP.R`); **PMI = Micro + Piccola + Media**, contrapposte alle **Grandi**.

*Regole di aggregazione per variabile* (dalla matrice impresa a quella comune):

| Colonna comunale | Regola | Significato |
|---|---|---|
| `n_imprese`, `n_micro`…`n_grande` | conteggio per classe | struttura imprenditoriale del comune |
| `dipendenti` | somma | occupazione |
| `asset_PMI_EUR`, `asset_grandi_EUR`, `asset_tot_EUR` | somma (totale = PMI + Grandi, identità verificata su tutti i comuni) | base di esposizione al rischio |
| `EBITDA_migl_EUR`, `ricavi_migl_EUR` | somma | dimensione economica |
| `integrazione_verticale_media` | media delle imprese | struttura produttiva |
| `ISP_std_medio` | media dell'ISP_std delle imprese del comune | contesto di performance (dall'ISP della tesi, §3.2) |
| `LISA_cluster_moda`, `Gi_bin_moda` | **moda** del cluster LISA / bin Gi* delle imprese del comune | trasposizione del risultato spaziale della tesi al comune |

#### b) Transizione provincia → comune (ripartizione top-down delle tariffe IVASS)

Le tariffe Cat-Nat IVASS sono definite a livello **provinciale** (110 province riconciliate)
come **tasso premiale per 10.000 € di asset esposto** (`premio_10k_prov`). La ripartizione sul
comune usa come chiave l'asset aggregato dalla transizione (a):

```
Premio_Teorico_Comunale_EUR = premio_10k_prov × asset_tot_EUR / 10.000
```

Formula verificata su tutti i 3.823 comuni della matrice definitiva. Per i **70 comuni sardi**
la tariffa provinciale è ripartita secondo l'assetto provinciale 2025 (flag
`premio_provincia_appross` = "yes"): il premio comunale resta una *quantità approssimata*.

#### c) Il rischio a livello comune (incrocio delle due fonti)

L'esposizione al rischio idrogeologico è costruita incrociando l'hazard territoriale (ISPRA)
con gli asset aggregati (transizione a):

```
hazard_frana_share    = PAI_area_P3P4_kmq / SUP_kmq        (quota di superficie a pericolosità frana P3/P4)
Risk_Frana_Asset_X    = hazard_frana_share × asset_X_EUR   (X = PMI, Grandi; idem per l'idraulico con IDR_area_P3)
```

*Esempio verificato* (Milano): `PAI_area_P3P4_kmq = 4,801891` su `SUP_kmq = 181,6727` →
`hazard_idraulico_share = 0,026432`; `Risk_Idraulico_Asset_PMI = 0,026432 × 43.788.803.521 =
1.157.405.937 €`; `premio_10k_prov = 15,05` su `asset_tot = 234.669.494.797` →
`Premio_Teorico_Comunale_EUR = 353.177.637,63`.

#### d) Conseguenza metodologica: MAUP e tracciabilità

Aggregare impresa→comune e ripartire provincia→comune espone l'analisi al **problema
dell'unità areale modificabile** (MAUP, Openshaw 1984) e alla **fallacia ecologica**
(Robinson 1950): i risultati valgono a scala comunale e non sono riportabili alle singole
imprese (e viceversa). Le regole §1.3a–c rispondono in modo trasparente: (i) ogni colonna ha
una formula esplicita; (ii) la dipendente `Premio_Teorico_Comunale_EUR` **incorpora per
costruzione** il tasso provinciale e gli asset — quindi la SDM del §5 stima *come la tariffa
risponde alla composizione del rischio sul territorio*, non un prezzo di mercato indipendente;
questo va dichiarato nell'interpretazione (§5.2); (iii) l'identità `asset_tot = asset_PMI +
asset_grandi` e la verifica della formula del premio su tutti i comuni sono controlli di
coerenza interni della matrice.

Le 37 colonne coprono quindi: geografia e identificativi (1–7), aggregati d'impresa dalla
transizione a (8–19), dati territoriali ISTAT/ISPRA (20–28), rischio incrociato dalla
transizione c (29–31), tariffazione IVASS dalla transizione b (32–34), e i riassuntori
dell'ISP dalla tesi (35–37).

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
k-esimo vicino; il "gomito" separa vicini omogenei da vicini arbitrari) e **stabilità dell'I di
Moran** (Moran, 1950) della variabile target per k crescente. L'intersezione seleziona **k = 77**
(log `analisi_k_nazionale.txt`); matrice risultante: 2.361.821 legami, 0,251% pesi non nulli,
S0 = 30.673, S1 = 699,86, S2 = 125.514,1. Per i 20 subset Dimensione × Macroarea si ripete il
criterio con k* specifici (16–74; tabella in `docs/log_R_livello_impresa.md`).

### 3.2 Il contesto interpretativo dell'ISP

La variabile esplorata a livello impresa è l'**ISP — Indicatore Sintetico di Performance**,
indicatore composito sviluppato nella tesi. **Nomenclatura**: nel README del repo
[tesi-magistrale](https://github.com/pietroscik/tesi-magistrale) compare come "Indice di
Sostenibilità Potenziale" — un errore di autocorrezione da correggere in lettura; nel commento
dello script compare come "Indicatore Sintetico Posizione", denominazione che richiama la
**fase di studio della distribuzione finale**: l'indice viene analizzato tramite la
trasformazione **OrderNorm** (`bestNormalize`), che ordina i ranghi e li ridistribuisce
valutando la **differenza interquartile e la distanza dalla media centrale** — di qui il
riferimento alla "posizione" dell'impresa nella distribuzione. Il nome corretto dell'indicatore
è **Indicatore Sintetico di Performance**.

Costruzione e validazione sono documentate e replicabili nel repository della tesi:
`suddivisione_script/divisione e creazione ISP.R` (FASE 2: feature engineering e calcolo),
`ISP validazione e inferenza.R` (FASE 3: validazione e finalizzazione) e gli output
`03_validation/` (pesi e statistiche delle versioni ISP).

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

**Finalizzazione.** La versione usata è **`ISP_bn`**: trasformazione robusta
`bestNormalize`/OrderNorm (Peterson & Cavanaugh, 2019) dell'ISP settoriale finale, scelta per
avvicinare la distribuzione alla normalità richiesto dai test spaziali; `ISP_std` e
`lag_ISP_std` sono la standardizzazione e il ritardo spaziale usati nello scatterplot di Moran.

**Conseguenze interpretative.** L'ISP è un **indice ordinale e relativo**: (i) i suoi valori
dipendono dalla normalizzazione min-max sul campione osservato e dai pesi stimati, quindi non
hanno unità di misura né lettura assoluta; (ii) confronti tra imprese di settori diversi
ereditano i pesi settoriali; (iii) essendo funzione monotona dei componenti di bilancio, le sue
correlazioni non aggiungono informazione causale oltre a quelle. Perciò l'analisi esplorativa
spaziale su ISP_bn (I di Moran = 0,00865, p < 2,2e-16; cluster LISA, Anselin 1995) è
**descrittiva**: dice *dove e quanto* la performance relativa è spazialmente associata, non
*perché*.

**Ruolo in questo repo.** Lo studio dell'ISP è già stato completato nella tesi: qui non viene
ri-stimato, ma **riutilizzato** in due forme derivate dalla transizione §1.3a — come media
comunale (`ISP_std_medio`, colonna di contesto della matrice) e come modà dei cluster
(`LISA_cluster_moda`, `Gi_bin_moda`) — e come variabile dipendente dell'analisi di evoluzione
a livello impresa i cui log sono documentati in `docs/log_R_livello_impresa.md`. Dopo
correzione **FDR (Benjamini–Hochberg, 1995)** nessun cluster locale dell'ISP sopravvive
(`robustezza_finale_fdr_nazionale.txt`): esito che limita le conclusioni esplorative a livello
di singola impresa e motivava la transizione al livello comunale (§4–5), dove la domanda di
ricerca sui premi e le coperture trova risposta robusta.

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

1. i premi nascono da tariffe **provinciali** ripartite sui comuni (§1.3b): per costruzione i
   comuni della stessa provincia condividono il tasso → dipendenza spaziale attesa, da
   verificare con ρ (l'ipotesi è testabile, non assunta: Moran residui e LR la verificano);
2. il rischio idrogeologico è **spazialmente correlato per natura fisica**: un evento frana
   non rispetta i confini comunali, quindi l'esposizione dei comuni vicini contiene informazione
   rilevante per il premio locale (termini WX);
3. l'assenza di queste componenti produce stime biased e inconsistenti (omissione di
   dipendenza spaziale, Anselin 1988), con effetti diretti e indiretti non separabili;
4. la dipendente **incorpora per costruzione** il tasso provinciale e gli asset (§1.3d):
   la lettura corretta dei parametri è *quanto la tariffazione provincialmente uniforme
   lascia spazio a una risposta al rischio locale* — ed è esattamente ciò che θ e gli effetti
   indiretti misurano.

Da qui la **SDM** come specifica di partenza, perché nests SAR e SEM (LeSage & Pace, 2009):
se i ritardi WX non servono, la SDM collassa in SAR (test LR su θ = 0); se vale la comune
radice θ = −ρβ, collassa in SEM. L'evidenza definitiva conferma la SDM: **ρ = 0,4966***
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
- 70 comuni sardi con premio approssimato (flag in matrice, §1.3b);
- MAUP/fallacia ecologica (§1.3d): i risultati valgono a scala comunale;
- i risultati del livello impresa su ISP non sono confrontabili direttamente con quelli del
  livello comune (unità, variabile dipendente e scopi diversi): il primo descrive la
  performance, il secondo risponde alla domanda sui premi e le coperture.
