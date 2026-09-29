# Terzo hazard sismico: fonti, estrazione dei dati e matching

Metodologia completa dell'estensione con la pericolosità sismica (MPS04, Istituto Nazionale
di Geofisica e Vulcanologia) come terzo hazard del modello a livello comune, accanto a
frana e idraulico. I risultati sono in `docs/risultati_sdm_sismico.md`; il quadro macchina
in `results/FINAL_sismico_k5.json`. Tutti gli script citati sono in `scripts/sismico/`
(Python stdlib puro + un solo script Node).

## 1. Fonti INGV (MPS04)

Le mappe di pericolosità sismica di riferimento per l'Italia sono le elaborazioni MPS04
(Meletti et al., 2004; Gruppo di lavoro MPS04, 2004) pubblicate dall'INGV.

| Griglia | Contenuto | Formato | Punti |
|---|---|---|---|
| `italia_ag_002` | **ag** (accelerazione massima del suolo, g), RP 475 anni (10% in 50), con 16°/84° percentile | testo, passo 0,02° | 104.564 |
| `SA_0475` | **Sa(T=0,10 s)**, RP 475 (10% in 50), fogli al 16°/50°/84° percentile | `.xls` (BIFF8) | 16.852 |
| `SA_1000` | Sa(T=0,10 s), RP 1.000 anni (5% in 100) | `.xls` (BIFF8) | 16.852 |
| `SA_2500` | Sa(T=0,10 s), RP 2.500 anni (2% in 50) | `.xls` (BIFF8) | 16.852 |
| ag "81%" | ag RP ~30 anni, con 16°/84° percentile | `.xlsx` fornito in elaborazione | 16.852 |
| ag "63%" | ag RP ~72 anni, con 16°/84° percentile | `.xlsx` fornito in elaborazione | 16.852 |

URL di download (griglie non archiviate nel repo, ~12 MB l'una):

- ag RP475: `https://zonesismiche.mi.ingv.it/elaborazioni/dati/italia_ag_002_txt.zip`
- Sa RP475/RP1000/RP2500: `https://esse1.mi.ingv.it/data/SA_0475.zip`,
  `https://esse1.mi.ingv.it/data/SA_1000.zip`, `https://esse1.mi.ingv.it/data/SA_2500.zip`

Verifiche di coerenza eseguite sulle griglie: L'Aquila Sa(0,10 s) RP475 (50° perc) = 0,533,
Milano = 0,116; rapporto Sa(0,10 s)/ag ≈ 2 coerente tra i punti; griglie 16.852 punti con
ID concordi tra i diversi tempi di ritorno; join con la griglia ag completa senza anomalie.

## 2. Estrazione: perché un parser dedicato

I file `.xls` SA_* dell'INGV hanno una struttura OLE2 **difettosa** (la continuazione DIFAT
è sfasata di un settore e i blocchi FAT di coda non sono raggiungibili): i lettori standard
(LibreOffice, pandas/xlrd, Excel) li aprono in modo parziale o non li aprono. La pipeline
usa quindi due strumenti stdlib-only:

- **`parse_biff.py`** — parser minimale OLE2 + BIFF8 (NUMBER, RK, MULRK, LABELSST con SST e
  CONTINUE, LABEL, MULBLANK) che estrae tutte le celle di tutti i fogli;
- **`recover_sa.py`** — recovery specifico: ricostruisce il FAT dai settori contigui a
  partire dal primo settore DIFAT valido, rilegge la directory direttamente da settore e
  ricuce lo stream del Workbook (scartando il settore-header OLE duplicato). Estrae i 9
  fogli (3 file × 3 percentili) con **0 anomalie**;
- **`xlsx_to_csv.py`** — conversione dei `.xlsx` delle griglie ag "81%"/"63%" (zip + XML,
  stdlib).

Uso: `python3 recover_sa.py SA_0475.xls prefisso` produce `prefisso_SA_475_{16,50,84}percentile.csv`.
Il 50° percentile alimenta il matching; 16°/84° restano disponibili per l'analisi
d'incertezza. Una riga incompleta per foglio (percentile 84) è trascurabile e scartata dai
controlli di completezza.

## 3. Matching comune → griglia (`match_sismico.py`)

- **Hash spaziale** con celle di 0,2° per restringere la ricerca, poi punto di griglia più
  vicino a ciascun comune (raggio massimo ammesso **5 km**, controllo di qualità);
- griglia ag RP475 a passo 0,02°: distanza massima osservata **1,91 km**; griglie
  Sa/ag 16.852 punti (passo ~0,05°): distanza massima ~4,2 km;
- esito: **3.823/3.823 comuni abbinati, 0 respinti**;
- output `data/sismico/matrice_sismica_ingv.csv`: 22 colonne (ag e Sa ai diversi tempi di
  ritorno con percentili, distanze di match, flag `sismico_non_classificato`).

**Regola Sardegna.** La Sardegna non è classificata sismicamente (OPCM 3274/2003 e
successive: zona "non classificata") e le griglie INGV non la coprono. I 70 comuni con
`COD_REG = 20` ricevono `ag = 0` (e analogamente Sa = 0) e `sismico_non_classificato = 1`,
senza tentare il match — evitando abbinamenti spuri a centinaia di km. La robustezza
"esclusione Sardegna" (n = 3.753) verifica che le stime non dipendano da questo blocco.

## 4. Matrice estesa (`build_matrice_v2.py`)

`data/Matrice_Modello_Savelli_Final_sismico.csv`: 3.823 comuni × **53 colonne** (37 della
matrice definitiva + 16 nuove), join su `PRO_COM`:

- `ag_RP475` (+16°/84° perc), `ag_RP30`, `ag_RP72`, `Sa01_RP475`, `Sa01_RP1000`, `Sa01_RP2500`;
- `sismico_non_classificato`;
- `hazard_sismico = ag_RP475` (in g);
- `Risk_Sismico_Asset_PMI = ag_RP475 × asset_PMI_EUR` e
  `Risk_Sismico_Asset_Grandi = ag_RP475 × asset_grandi_EUR`.

**Unità e confronto dei coefficienti.** Per frana e idraulico l'hazard è una **quota di
area comunale** (share, adimensionale 0–1); per il sismico è un'**accelerazione** (g).
Le due famiglie di regressori hanno quindi scale diverse: i β non sono confrontabili in
magnitudine tra hazard, mentre segni, significatività e stabilità nelle robustezze sono
confrontabili. Le robustezze sui tempi di ritorno (RP30/RP72/Sa RP1000/Sa RP2500) testano
la sensibilità alla scelta dell'intensità.

## 5. Stima (`definitivo_sismico.js`)

Stesso motore ML del modello definitivo (tracce Monte Carlo T=45, M=120, seed fisso;
sezione aurea per ρ; Newton numerico; convenzioni AIC: K = 2p+3 per SDM, K = p+3 per
SAR/SEM; effetti LeSage–Pace via matrice (I−ρW)^{-1}), con:

- **baseline p = 2** come controllo di riproduzione (replica esatta di
  `results/FINAL_k5.json`, §1 di `docs/risultati_sdm_sismico.md`);
- **SDM p = 4** con i quattro regressori Frana/Sismico × PMI/Grandi e confronto SAR/SEM;
- robustezze: hazard alternativi (RP30, RP72, Sa RP1000, Sa RP2500), k = 6/7/8,
  esclusione Sardegna;
- errori standard via Hessiana numerica finale con la **formula corretta a 4 angoli**
  (vedi la nota bug in `docs/risultati_sdm_sismico.md`: la Hessiana interna di
  `scripts/definitivo.js` usava perturbazioni errate fuori diagonale).

## 6. Pipeline di rigenerazione completa

```bash
# 1. scaricare e scompattare le griglie (URL del §1) in data/sismico/griglie/
cd data/sismico/griglie

# 2. estrarre i fogli SA (recovery) e convertire le griglie ag .xlsx
python3 scripts/sismico/recover_sa.py SA_0475.xls sa0475    # idem SA_1000, SA_2500
python3 scripts/sismico/xlsx_to_csv.py ag_81percento.xlsx ag_81_RP30.csv

# 3. matching comune → griglia
python3 scripts/sismico/match_sismico.py

# 4. matrice estesa 53 colonne
python3 scripts/sismico/build_matrice_v2.py

# 5. stima SDM p=4 (baseline p=2 inclusa come controllo) -> results/FINAL_sismico_k5.json
node scripts/sismico/definitivo_sismico.js
```

Tutti i percorsi sono configurabili via variabili d'ambiente (`MATRICE`, `AG10_TXT`,
`AG81_CSV`, `AG63_CSV`, `SA475`, `SA1000`, `SA2500`, `OUT`, `SISMICA`, `MATRIX`, `OUTJSON`).
Il log-determinante è approssimato per tracce Monte Carlo: tolleranza attesa ~1e-3 su ρ e
logL tra esecuzioni con PRNG diversi.
