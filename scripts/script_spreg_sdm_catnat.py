# =====================================================================
# Spatial Durbin Model — Cat-Nat (Legge 78/2025) — STIMA DEFINITIVA
# Input: Matrice_Modello_Savelli_Final.csv (3.823 comuni, 37 colonne)
#        NB: il workbook canvas "matrice-modello-savelli" ha già le
#        coordinate corrette per i 4 comuni; il fix è ridondante qui sotto.
# Requisiti: pip install pandas numpy libpysal spreg esda
#
# Risultati ML definitivi (replicati in pure JS: logdet Barry-Pace,
# tracce MC, Newton sulla log-likelihood piena, SE da Hessiana inversa):
#   SDM k=5: rho=0.4966*** (SE 0.0216, z=22.97)  logL=-6351.71
#   beta_PMI=0.0061 n.s. (SE 0.0050) | beta_Grandi=0.1434*** (SE 0.0039)
#   theta_PMI=-0.0521*** (SE 0.0063) | theta_Grandi=-0.0374*** (SE 0.0093)
#   SAR: rho=0.4641 logL=-6415.06 | SEM: lam=0.5471 logL=-6406.53
#   LR SDM vs SAR = 126.71*** | LR SDM vs SEM (Burridge) = 109.64***
#   AIC: SDM 12717.4 < SEM 12821.1 < SAR 12838.1 -> SDM vincente
#   Effetti: Grandi dir=+0.1511 ind=+0.0595 tot=+0.2106
#            PMI   dir=+0.0006 (n.s.) ind=-0.0918 tot=-0.0912
#   Moran residui SDM: I=-0.0486, p(perm 499, 2 code)=0.004
#   RESET F=24.4*** (forma funzionale) | BP LM=124.4*** (eterosched.)
#   -> usare SE robusti/HC; vedere canvas "metodologia-w-sensibilita-k"
# =====================================================================

import numpy as np
import pandas as pd
from libpysal.weights import KNN, DistanceBand
from esda.moran import Moran, Moran_Local
import spreg

# ---------------------------------------------------------------
# 1. DATI
# ---------------------------------------------------------------
CSV_PATH = "data/Matrice_Modello_Savelli_Final.csv"
df = pd.read_csv(CSV_PATH, dtype={"PRO_COM": str})

# fix coordinate errate nel CSV originale (già corretto nel workbook canvas)
_COORD_FIX = {
    "28041": (45.581, 11.706548),    # Gazzo (VI)
    "51021": (43.274463, 11.746),    # Lucignano (AR)
    "12108": (45.631, 8.889841),     # Olgiate Olona (VA)
    "62074": (41.216414, 14.527),    # Telese Terme (BN)
}
for pc, (la, lo) in _COORD_FIX.items():
    m = df["PRO_COM"].str.strip() == pc
    df.loc[m, "lat"] = la
    df.loc[m, "long"] = lo

Y_NAME  = "Premio_Teorico_Comunale_EUR"
X_NAMES = ["Risk_Frana_Asset_PMI", "Risk_Frana_Asset_Grandi"]

df = df.dropna(subset=[Y_NAME] + X_NAMES).reset_index(drop=True)

# trasformazione log (molte osservazioni = 0 -> log1p)
y = np.log1p(df[Y_NAME]).values.reshape(-1, 1)
X = np.log1p(df[X_NAMES]).values

# ---------------------------------------------------------------
# 2. MATRICE DEI PESI W — KNN k=5 sui centroidi (lat/long in matrice)
#    (stessa W usata nella stima ML replicata: k=5, row-standardized)
# ---------------------------------------------------------------
coords = df[["long", "lat"]].values
w = KNN.from_array(coords, k=5)
w.transform = "r"
print(f"W: {w.n} unità, isolati: {w.islands}")

# ---------------------------------------------------------------
# 3. MORAN'S I (Y e X) + LISA
# ---------------------------------------------------------------
for var in [Y_NAME] + X_NAMES:
    mi = Moran(df[var].values, w, permutations=999)
    print(f"\n--- Moran's I: {var} ---")
    print(f"I = {mi.I:.4f} | z = {mi.z_norm:.3f} | p = {mi.p_norm:.4f}")

# ---------------------------------------------------------------
# 4. STIMA SDM (ML)  y = rho*Wy + X*beta + WX*theta + a + eps
#    spreg non ha ML-SDM dedicato: si usa ML_Lag su [X, WX]
# ---------------------------------------------------------------
WX = w.sparse @ X
X_durbin = np.hstack([X, WX])
names_durbin = X_NAMES + [f"W_{v}" for v in X_NAMES]

model_durbin = spreg.ML_Lag(
    y=y, X=X_durbin, w=w, method="full", epsilon=1e-7,
    name_y="ln_premio", name_x=names_durbin, name_w="knn5",
)
print("\n================ SDM (ML) ================")
print(model_durbin.summary)

# confronti per i test di specificazione
model_sar = spreg.ML_Lag(y=y, X=X, w=w, method="full", name_w="knn5")
model_sem = spreg.ML_Error(y=y, X=X, w=w, method="full", name_w="knn5")

ll_sdm = model_durbin.ll if hasattr(model_durbin, "ll") else model_durbin.logll
ll_sar = model_sar.ll     if hasattr(model_sar, "ll")     else model_sar.logll
ll_sem = model_sem.ll     if hasattr(model_sem, "ll")     else model_sem.logll

from scipy import stats as sstats
LR_sar = -2 * (ll_sar - ll_sdm)   # H0: theta = 0          (2 gdl)
LR_sem = -2 * (ll_sem - ll_sdm)   # H0 Burridge: theta=-rho*beta (2 gdl)
print(f"LR SDM vs SAR = {LR_sar:.3f}, p = {1 - sstats.chi2.cdf(LR_sar, 2):.6f}")
print(f"LR SDM vs SEM = {LR_sem:.3f}, p = {1 - sstats.chi2.cdf(LR_sem, 2):.6f}")

# ---------------------------------------------------------------
# 5. EFFETTI DIRETTI/INDIRETTI/TOTALI
#    (valori definitivi attesi: PMI dir=+0.0006 n.s., ind=-0.0918; Grandi dir=+0.1511, ind=+0.0595)
# ---------------------------------------------------------------
def impacts_sdm(beta_k, theta_k, rho, w_sparse, n_iter=200):
    """Effetti via traccia MC della serie (I-rho*W)^-1 (beta I + theta W)."""
    n = w_sparse.shape[0]
    rng = np.random.default_rng(777)
    tot = 0.0
    for _ in range(20):
        u = rng.choice([-1.0, 1.0], size=n)
        c = beta_k * u + theta_k * (w_sparse @ u)
        z = c.copy()
        wk = c.copy()
        for _ in range(n_iter):
            wk = rho * (w_sparse @ wk)
            z = z + wk
        tot += u @ z
    direct = tot / (20 * n)
    total = (beta_k + theta_k) / (1 - rho)   # esatto per W row-standardized
    return direct, total - direct, total

betas = model_durbin.betas.flatten()
for j, name in enumerate(X_NAMES):
    d, i_, t_ = impacts_sdm(betas[1 + j], betas[3 + j], betas[-1], w.sparse)
    print(f"{name}: direct={d:.4f} indirect={i_:.4f} total={t_:.4f}")

# ---------------------------------------------------------------
# 6. MORAN SUI RESIDUI SDM
#    (valore definitivo atteso: I = -0.0486, p(perm, 2 code) = 0.004 — sovra-cattura lieve, v. nota nel testo)
# ---------------------------------------------------------------
resid = model_durbin.u.flatten() if hasattr(model_durbin, "u") else None
if resid is not None:
    mi_r = Moran(resid, w, permutations=999)
    print(f"Moran residui SDM: I = {mi_r.I:.4f}, p = {mi_r.p_sim:.4f}")

# ---------------------------------------------------------------
# 7. EXPORT
# ---------------------------------------------------------------
df["resid_sdm"] = resid if resid is not None else np.nan
df.to_csv("Matrice_SDM_risultati.csv", index=False)
print("\nSalvato: Matrice_SDM_risultati.csv")
