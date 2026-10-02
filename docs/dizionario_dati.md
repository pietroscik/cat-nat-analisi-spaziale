# Dizionario dei dati — matrici a livello comune

Descrizione delle colonne di `data/Matrice_Modello_Savelli_Final.csv` (3.823 comuni × 37
colonne) e delle 16 colonne sismiche aggiunte in
`data/Matrice_Modello_Savelli_Final_sismico.csv` (3.823 × 53). Le formule di costruzione
e le transizioni di scala sono documentate nel capitolo metodologico (§1.3) e nel
`README.md` (§2); qui il significato colonna per colonna.

Unità di misura: gli importi monetari sono in **EUR** (migliaia di EUR dove indicato);
le quote di area (hazard frana/idraulico) sono **adimensionali** (0–1); l'accelerazione
sismica (ag, Sa) è in **g**.

## 1. Identificativi e geografia (1–7)

| Colonna | Descrizione |
|---|---|
| `PRO_COM` | codice ISTAT numerico del comune (chiave di join; 3.823 valori unici) |
| `COMUNE` | denominazione del comune |
| `COD_PROV` | codice ISTAT della provincia |
| `Provincia` | denominazione della provincia (107 province in matrice su 110 riconciliate IVASS) |
| `COD_REG` | codice ISTAT della regione (20 = Sardegna, non classificata sismicamente) |
| `lat`, `long` | coordinate del centroide comunale; corrette per Gazzo (VI), Lucignano (AR), Olgiate Olona (VA), Telese Terme (BN); fusioni comunali ricodificate (Sovizzo, Presicce-Acquarica, Figline e Incisa Valdarno, Trapani/Misiliscemi) |

## 2. Imprese — transizione impresa → comune (8–19)

Aggregazione bottom-up delle 30.673 imprese AIDA georeferenziate via spatial join sui
poligoni ISTAT. PMI = Micro + Piccola + Media.

| Colonna | Descrizione |
|---|---|
| `n_imprese` | numero di imprese geolocalizzate nel comune |
| `n_micro`, `n_piccola`, `n_media`, `n_grande` | conteggi per classe dimensione |
| `dipendenti` | somma dei dipendenti |
| `asset_PMI_EUR` | somma degli asset delle PMI (EUR) |
| `asset_grandi_EUR` | somma degli asset delle Grandi imprese (EUR) |
| `asset_tot_EUR` | somma complessiva; identità verificata `asset_tot = asset_PMI + asset_grandi` (scarto massimo 0 EUR) |
| `EBITDA_migl_EUR` | somma dell'EBITDA (migliaia di EUR) |
| `ricavi_migl_EUR` | somma dei ricavi (migliaia di EUR) |
| `integrazione_verticale_media` | media comunale dell'indice di integrazione verticale |

## 3. Territorio e popolazione — ISTAT/ISPRA (20–26)

| Colonna | Descrizione |
|---|---|
| `SUP_kmq` | superficie comunale (kmq) |
| `POP_2018` | popolazione residente 2018 |
| `DENSPOP` | densità di popolazione = `POP_2018 / SUP_kmq` (ab/kmq) |
| `PAI_area_P3P4_kmq` | area in pericolosità da frana P3/P4 (kmq, ISPRA) |
| `IDR_area_P3_kmq` | area in pericolosità idraulica P3 (kmq, ISPRA) |
| `PAI_Pop_P3P4` | popolazione esposta a frana P3/P4 |
| `IDR_Pop_P3` | popolazione esposta a idraulico P3 |

## 4. Hazard — quote di area (27–28)

| Colonna | Descrizione |
|---|---|
| `hazard_frana_share` | quota di area comunale in pericolosità frana = `PAI_area_P3P4_kmq / SUP_kmq` (0–1) |
| `hazard_idraulico_share` | quota di area comunale in pericolosità idraulica = `IDR_area_P3_kmq / SUP_kmq` (0–1) |

## 5. Rischio incrociato — hazard × asset, EUR (29–31)

| Colonna | Descrizione |
|---|---|
| `Risk_Frana_Asset_PMI` | `hazard_frana_share × asset_PMI_EUR` (esposizione effettiva delle PMI) |
| `Risk_Frana_Asset_Grandi` | `hazard_frana_share × asset_grandi_EUR`; **= 0 in 3.216 comuni su 3.823** (nessuna Grande impresa esposta e/o quota di area nulla): con `log1p` questi comuni formano il gruppo di riferimento del modello |
| `Risk_Idraulico_Asset_PMI` | `hazard_idraulico_share × asset_PMI_EUR` |

## 6. Premio teorico — transizione provincia → comune, IVASS (32–34)

| Colonna | Descrizione |
|---|---|
| `premio_10k_prov` | tariffa provinciale IVASS per 10.000 EUR di asset esposto |
| `premio_provincia_appross` | flag `"si"` per i 70 comuni sardi con tariffa ripartita per l'assetto provinciale 2025 (`"no"` negli altri) |
| `Premio_Teorico_Comunale_EUR` | `premio_10k_prov × asset_tot_EUR / 10.000` (formula verificata su tutti i 3.823 comuni, errore relativo < 1e-5); variabile dipendente del modello (via `log1p`) |

## 7. Indicatori d'impresa aggregati — dalla tesi (35–37)

| Colonna | Descrizione |
|---|---|
| `ISP_std_medio` | media comunale dell'ISP standardizzato (Indicatore Sintetico di Performance, indicatore composito della tesi) |
| `LISA_cluster_moda` | moda del cluster LISA/Gi\* delle imprese del comune (`High-High`, `High-Low`, `Low-High`, `Low-Low`, `Non signif`, vuoto se senza imprese) |
| `Gi_bin_moda` | moda della classificazione binaria Gi\* (`Hotspot 95%`/`99%`, `Coldspot 95%`/`99%`, `Non signif`; vuoto = comune senza imprese classificabili) |

## 8. Colonne sismiche — MPS04 INGV (38–53)

Aggiunte da `scripts/sismico/build_matrice_v2.py` (join su `PRO_COM` con
`data/sismico/matrice_sismica_ingv.csv`; matching al punto di griglia più vicino, vedi
`docs/sismico_metodologia.md`). La **Sardegna** non è classificata sismicamente (OPCM
3274/2003): `ag`/`Sa` = 0 e `sismico_non_classificato = 1` senza tentare il match.

| Colonna | Descrizione |
|---|---|
| `ag_RP475` | accelerazione massima del suolo (g), RP 475 anni — 10% in 50 anni; griglia `italia_ag_002` a passo 0,02° (104.565 punti) |
| `ag_RP475_16perc`, `ag_RP475_84perc` | ag RP475 al 16° e 84° percentile |
| `ag_RP30` | ag 81% in 50 anni (RP ~30 anni), con `ag_RP30_16perc`/`ag_RP30_84perc` |
| `ag_RP72` | ag 63% in 50 anni (RP ~50 anni; la denominazione RP72 è d'archivio), con `ag_RP72_16perc`/`ag_RP72_84perc` |
| `Sa01_RP475`, `Sa01_RP1000`, `Sa01_RP2500` | spettro di risposta Sa(T = 0,10 s) ai RP 475 (10% in 50), 1.000 (5% in 50) e 2.500 anni (2% in 50) |
| `sismico_non_classificato` | flag 1 = Sardegna / match respinto (ag = 0) |
| `hazard_sismico` | = `ag_RP475` (g); accelerazione, non quota di area: i β delle due famiglie di hazard non sono confrontabili in magnitudine |
| `Risk_Sismico_Asset_PMI` | `ag_RP475 × asset_PMI_EUR` |
| `Risk_Sismico_Asset_Grandi` | `ag_RP475 × asset_grandi_EUR` |

## Nota su dati terzi

Tutte le colonne derivano da elaborazioni aggregate di dati di terzi parti (AIDA/Bureau
van Dijk, ISTAT, ISPRA, INGV, IVASS) archiviate nel solo scopo di riproducibilità
scientifica: vedi `LICENSE` per la nota sui diritti dei dati originali.
