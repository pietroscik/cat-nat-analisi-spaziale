# Recupero e validazione dati per la chiusura della coda frequente RP<30 (ramo dev)

Fase di validazione del ramo `dev`: si recuperano le fonti primarie INGV per
rispondere alla domanda aperta della **coda frequente RP<30** (dichiarata non
stimabile dalla curva EP in `docs/pricing_coerenza.md` §2.1 e §8.2) e si valuta
se, e con quali dati, la banda può essere chiusa. Questo documento **recensisce
il tutto**: cosa esiste pubblicamente, cosa è stato scaricato, come è stato
validato, e quale protocollo ne deriva. Script: `scripts/esplorazione/validazione_cpti15.py`
(sole stdlib, deterministico); output: `results/esplorazione_cpti15.json`.

## 1. Recupero hazard: verdetto MPS04 — sotto RP30 non esiste nulla di pubblico

La §"fonti" di `docs/sismico_metodologia.md` documenta le griglie MPS04 usate:
ag a RP475 (`italia_ag_002_txt.zip`, ancora online: HTTP 200, ~1,5 MB,
Last-Modified 2006) e le griglie ag "81% in 50 anni" (≈ RP30) e "63%" (≈ RP50–72)
**"fornite in elaborazione"**, cioè senza URL pubblico. Verifiche fatte in questa
fase:

- catalogo CKAN INGV (`https://data.ingv.it`, `package_show` id `70` e `193`):
  i metadati esistono ma `num_resources = 0` — **nessun file scaricabile**;
- probe diretti sotto `https://zonesismiche.mi.ingv.it/elaborazioni/dati/`
  con i nomi plausibili (`italia_ag_81_002_txt.zip`, varianti 81/63, zip/txt):
  tutti **404**; l'unico URL vivo resta la griglia RP475.

**Verdetto**: il punto di pericolo più frequente pubblicato in MPS04 è l'81% in
50 anni, che equivale a **RP≈30**. **Non esistono punti di hazard pubblici con
periodo di ritorno inferiore a 30 anni**: la pendenza della coda frequente
(`k` locale nella banda RP<30) **non è misurabile dalle mappe di pericolo**.
Le pendenze MPS04 misurate (§8.2 di `pricing_coerenza.md`: k=3,70 in RP30–72,
k=2,11 in RP72–475) si fermano al bordo RP30 e lì resta il muro informativo.
L'unica strada per la banda RP<30 è **dal lato eventi**: frequenze empiriche
del catalogo (CPTI15) × fragilità × esposizione.

## 2. Recupero eventi: catalogo CPTI15 v4.0

Scaricato il catalogo parametrico dei terremoti italiani di riferimento:

- URL: `https://emidius.mi.ingv.it/CPTI15-DBMI15/data/CPTI15_v4.0.xlsx`
- md5: `c7e9a48ec7b348fec669ecf21906e3bc`, 1.363.375 byte;
- citazione: Rovida A., Locati M., Camassi R., Lolli B., Gasperini P. (a cura
  di), *CPTI15 v4.0*, INGV — uso scientifico con citazione obbligatoria;
- **il file non è archiviato nel repo** (convenzione del repo: i dati non si
  committano, si documenta l'URL); lo script lo cerca in `data/sismico/` o nel
  percorso indicato dalla variabile d'ambiente `CPTI15_XLSX`.

Nota di parsing: il catalogo è nel **foglio 3** del workbook (il parser di
`scripts/sismico/xlsx_to_csv.py` legge solo il primo foglio, che qui è "about"):
`validazione_cpti15.py` individua il foglio giusto cercando l'intestazione
contenente `MwDef` — robusto a riordinamenti dei fogli.

## 3. Validazione del catalogo (`results/esplorazione_cpti15.json`)

**Copertura**: 4.894 righe-evento (1005–2020), di cui **4.737 con Mw definito**
(157 senza, quasi tutte medievali); Mw 2,22–7,32. Aggiornamento 2020: contiene
l'intera sequenza 2016–2017 e la serie 2012 Emilia.

**Tassi empirici** (eventi/anno, finestre di completezza crescente; quella di
riferimento è 1980–2020, strumentale, completa per M≥4,5):

| Magnitudo | 1980–2020 (rif.) | 1960–2020 | 1900–2020 | Periodo di ritorno (1980–2020) |
|---|---|---|---|---|
| Mw ≥ 4,5 | 8,56/anno | 8,72/anno | 7,79/anno | — |
| Mw ≥ 5,0 | 2,24/anno | 2,34/anno | 2,46/anno | ~0,45 anni |
| Mw ≥ 5,5 | 0,71/anno | 0,71/anno | 0,62/anno | **~1,4 anni** |
| Mw ≥ 6,0 | 0,15/anno | 0,16/anno | 0,15/anno | **~6,8 anni** |
| Mw ≥ 6,5 | 0,049/anno | 0,033/anno | 0,058/anno | ~20,5 anni |

I tassi di M≥5,5–6,0 sono **stabili** al variare della finestra (1960/1980/2000):
la stima non dipende dalla scelta della soglia di completezza. I sei M≥6
1980–2020 sono tutti eventi noti e identificati: Irpinia-Basilicata 1980
(Mw 6,81), L'Aquila 2009 (6,29), Pianura emiliana 2012 (6,09) e la sequenza
2016 (Monti della Laga 6,18; Visso/Valnerina 6,07; Norcia/Valnerina 6,61).
L'Umbria-Marche 1997 resta appena sotto soglia (Mw 5,97).

**Eventi noti** (controllo puntuale, 5/5 superati): L'Aquila 2009 Mw 6,29;
Amatrice 2016 Mw 6,18; Emilia 2012 Mw 6,09; Irpinia 1980 Mw 6,81; Molise 2002
Mw 5,74 (atteso 5,72). Aree epicentrali corrette ("Aquilano", "Monti della
Laga", "Pianura emiliana", "Irpinia-Basilicata", "Molise").

**Coerenza interna Mw→Io**: mediana di Io crescente monotona per classe di Mw
(Io 5 per M∈[4,0;4,5) fino a Io 10 per M≥6; Io max 11). Il catalogo è
internamente coerente.

**Nota di parsing critica**: il campo `IoDef` può essere un intervallo (es.
`6-7`): va convertito al punto medio, non con `float()` — altrimenti ~1.310
eventi vengono scartati silenziosamente (bug trovato e corretto in questa fase).

## 4. Lettura attuariale: gli eventi dannosi cadono DENTRO la banda RP<30

Il fatto chiave che cambia la lettura: gli eventi **dannosi** M≥6 ricorrono
storicamente **ogni ~6,8 anni** — un periodo di ritorno ben dentro la banda
RP<30 che la curva EP non copre. Cioè:

1. la coda frequente **non è un'astrazione da estrapolare**: è fatta di eventi
   osservati e ricorrenti (6 eventi M≥6 in 41 anni);
2. la chiusura corretta è **event-based**: tassi CPTI15 × fragilità ×
   esposizione, con l'esposizione già disponibile nella matrice del modello
   (`asset_tot_EUR` per comune);
3. il bound max-ent di `coda_frequente_prior.py` (AAL frequente mediana
   **3,18 mld €/anno**, totale 4,7× il benchmark CURVE=4) va trattato come
   **limite superiore di convenienza**, non come stima: estrapolare la curva EP
   con k=3,70 dalla banda RP30–72 giù fino a RP<30 attribuisce alla banda
   frequente più AAL di quanto l'intero benchmark attribuisca al sismico, e
   non è ancorato ad alcun evento osservato. Il confronto con i tassi empirici
   di CPTI15 è esattamente il collaudo che manca.

## 5. Protocollo della fase di validazione (prossimi passi dev)

1. **Stima event-based della coda frequente** (`esplorazione_coda_eventi.py`):
   per bin di magnitudo (5,0–5,5–6,0–6,5+), tasso empirico CPTI15 × intensità
   attesa (Io mediano per classe, già validato) × fragilità. Fonte fragilità:
   funzioni di fragilità per edifici italiani (Rosti et al. 2021, URM e RC),
   coerenti con l'MDR del benchmark; sensibilità su Due Ombrature.
2. **Gate di decisione**: confrontare l'AAL frequente event-based con
   (a) il limite inferiore RP≥30 = **1,713 mld €/anno** (totale 1,66× benchmark),
   (b) il bound max-ent = 3,18 mld, (c) il benchmark CURVE=4 = 1,035 mld.
   Se event-based << max-ent, il max-ent si archivia come bound lasco e il
   totale RP≥0 si stima come 1,71 + AAL_f(event-based); se invece convergono,
   la curva EP va ricalibrata.
3. **DBMI15 come passo successivo**: il database macrosismico (danni osservati
   per comune, stesso sito INGV) permette di sostituire la fragilità analitica
   con intensità osservate → danno, validando la catena eventi→danno sui dati
   storici (es. Irpinia 1980, Umbria-Marche 1997).
4. **Cautela interazioni** (principio stabilito nel ramo dev): ogni termine
   aggiunto va centrato, con VIF confrontato alla base 1,8 e stabilità LOPO
   verificata (barra β_ISP [−0,377;−0,373]); su questo lato eventi non si
   tocca il modello tariffario pubblicato su `main`.

## 6. Riproducibilità

Lo script è deterministico (doppia esecuzione: JSON identico, md5
`c1e516fe301be3d4390949cd12172694`) e usa solo stdlib. Il catalogo non è nel
repo: si scarica dall'URL di §2 (oppure si imposta `CPTI15_XLSX`), poi

```bash
python3 scripts/esplorazione/validazione_cpti15.py
```

genera `results/esplorazione_cpti15.json` (tassi su 5 finestre, eventi noti,
coerenza Mw→Io, verdetto). Nessun output di `main` cambia: l'esplorazione sta
sul ramo `dev` e confluis su `main` solo dopo validazione via PR.
