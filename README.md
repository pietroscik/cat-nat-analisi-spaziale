# cat-nat-analisi-spaziale

Pipeline di analisi spaziale per la tesi magistrale — **Cat-Nat (Legge 78/2025)**:
dal livello impresa (30.673 unità AIDA) al livello comune (3.823 comuni),
con stima SDM (Spatial Durbin Model) ML del premio teorico comunale.

## Struttura

```
data/
  Matrice_Modello_Savelli_Final.csv     # matrice definitiva 3.823 comuni x 37 col
                                        # (coordinate corrette: Gazzo, Lucignano,
                                        #  Olgiate Olona, Telese Terme)
  log_R_livello_impresa/               # 12 log R (spatialreg) dell'analisi
                                        # precedente a livello impresa (evoluzione)
docs/
  capitolo_metodologico_W_k.md          # paragrafo tesi: scelta W, sensibilità k,
                                        # stima definitiva, robustezza, MAUP
  script_spreg_sdm_catnat.py           # script spreg/PySAL definitivo
                                        # (da eseguire in locale per i numeri ufficiali)
scripts/
  step1_kdist_moran.py                 # curva k-dist + Moran vs k (comuni)
  sdm_grid.js                           # griglia SDM/SAR/SEM k=3..30 (pure JS)
  definitivo.js                         # stima definitiva SDM ML k=5..8,
                                        # Hessiana, effetti, Moran perm, RESET, BP
  se_definitivi.js                      # SE via inversa Hessiana piena (validato)
results/
  grid_results.json                     # esiti griglia k=3..30
  risultati_definitivi.json             # stime k=5..8 complete
  se_definitivi.json                    # coefficienti + SE validati (2 strade)
  FINAL_k5.json                         # numeri ufficiali di riferimento (k=5)
```

## Risultato definitivo (SDM ML, k=5, coord. corrette, n=3.823)

- ρ = 0,4966*** (SE 0,0216) | β_Grandi = 0,1434*** | β_PMI = 0,0061 n.s.
- θ_PMI = −0,0521*** | θ_Grandi = −0,0374***
- LR vs SAR = 126,71*** | LR vs SEM (Burridge) = 109,64*** → SDM vincente
- Effetti LeSage–Pace: Grandi dir +0,151*** / tot +0,211***; PMI dir n.s. /
  spillover −0,092***
- Robustezza: β_Grandi stabile 0,142–0,145 su k=3–30; gerarchia SDM invariante
- Diagnostica: Moran residui −0,049 (p=0,004, sovra-cattura lieve);
  RESET F=24,4***; BP LM=124,4*** → usare SE robusti/HC

## Riproducibilità

1. `scripts/step1_kdist_moran.py` → criterio k-dist + Moran (Python stdlib)
2. `scripts/sdm_grid.js` → griglia k (Node, nessuna dipendenza)
3. `scripts/definitivo.js` → stima definitiva (Node, nessuna dipendenza)
4. `docs/script_spreg_sdm_catnat.py` → replica con spreg/PySAL
   (`pip install pandas numpy libpysal spreg esda`)

Fonti dati: AIDA (valori grezzi, integrazione_verticale), ISTAT/ISPRA
Dissesto 2018, IVASS CatNat (110 province), centroidi ISTAT
(opendatasicilia/comuni-italiani).
