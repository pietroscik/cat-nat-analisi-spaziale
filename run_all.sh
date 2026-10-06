#!/bin/bash
# =================================================================
# Cat-Nat Analisi Spaziale - pipeline di riproduzione completa
# =================================================================
# Esegue la catena deterministica del repository (gli stessi passi
# della CI, .github/workflows/ci.yml):
#   Python (stdlib)  pricing, tessuto, EP, chi paga, spazializzazione,
#                    robustezza, validazione, grafici SVG
#   Node             SE robusti HC1 del SDM (scripts/se_robusti.js,
#                    riusa scripts/definitivo.js)
#
# Utilizzo:
#   bash run_all.sh                solo pipeline deterministica
#   bash run_all.sh --verify       + verifica git diff sugli output (come la CI)
#   bash run_all.sh --notebook     + notebook end-to-end (nbconvert, requirements.txt)
#   bash run_all.sh --riferimento  + selezione k e stima spreg/PySAL di riferimento
#                                   (richiede spreg>=1.3, libpysal>=4.7)
#
# Requisiti: Python 3.9+ (la catena pricing è solo stdlib), Node 18+.
# =================================================================
set -euo pipefail
cd "$(dirname "$0")"

VERIFY=0; NOTEBOOK=0; RIF=0
for arg in "$@"; do
  case "$arg" in
    --verify)     VERIFY=1 ;;
    --notebook)   NOTEBOOK=1 ;;
    --riferimento) RIF=1 ;;
    *) echo "flag sconosciuto: $arg (usi: --verify, --notebook, --riferimento)"; exit 1 ;;
  esac
done

command -v python3 >/dev/null 2>&1 || { echo "ERRORE: python3 non trovato"; exit 1; }
command -v node    >/dev/null 2>&1 || { echo "ERRORE: node non trovato (serve per scripts/se_robusti.js)"; exit 1; }

echo "== Cat-Nat: pipeline deterministica (pricing -> validazione -> SE robusti) =="
echo "-- [1/3] benchmark pricing e loss ratio (pricing_model)"
python3 scripts/pricing/pricing_model.py
echo "-- [2/3] esposizione del tessuto e approfondimenti attuariali"
python3 scripts/pricing/tessuto_produttivo.py
python3 scripts/pricing/ep_curve.py
python3 scripts/pricing/chi_paga.py
python3 scripts/pricing/spazializzazione.py
python3 scripts/pricing/robustezza_tessuto.py
python3 scripts/pricing/validazione_assunzioni.py
echo "-- [3/3] grafici SVG e SE robusti HC1 del SDM (Node)"
python3 scripts/pricing/grafici.py
python3 scripts/grafici_sdm.py
node scripts/se_robusti.js

if [ "$RIF" = 1 ]; then
  echo "-- opzionale: selezione k (curva k-dist + Moran) e stima di riferimento spreg/PySAL"
  python3 scripts/step1_kdist_moran.py
  python3 scripts/script_spreg_sdm_catnat.py
fi
if [ "$NOTEBOOK" = 1 ]; then
  echo "-- opzionale: notebook end-to-end (nbconvert --execute, assert vs FINAL_k5)"
  jupyter nbconvert --to notebook --execute \
    --output /tmp/riproduce_sdm_comuni_eseguito.ipynb \
    notebook/riproduce_sdm_comuni.ipynb
fi
if [ "$VERIFY" = 1 ]; then
  echo "-- verifica di riproducibilita' (git diff, come la CI)"
  git diff --quiet -- results/ docs/ || {
    echo "ERRORE: output non riproducibili (git diff)"; git diff --stat -- results/ docs/; exit 1; }
fi

echo "== pipeline completata: risultati in results/, grafici in docs/ =="
