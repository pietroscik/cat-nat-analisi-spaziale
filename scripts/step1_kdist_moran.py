import csv, math

# --- carica dati ---
rows = list(csv.DictReader(open('tool-results/data-analysis/scelta-k/matrice.csv')))
def f(s):
    try: return float(s)
    except: return None
lat = [f(r['lat']) for r in rows]; lon = [f(r['long']) for r in rows]
y   = [math.log1p(f(r['Premio_Teorico_Comunale_EUR'])) for r in rows]
x1  = [math.log1p(f(r['Risk_Frana_Asset_PMI'])) for r in rows]
x2  = [math.log1p(f(r['Risk_Frana_Asset_Grandi'])) for r in rows]
n = len(rows); print('n =', n)

KMAX = 40
pts = list(zip(lon, lat))
# k-vicini ordinati per ogni punto (brute force)
nn_idx = [None]*n
for i in range(n):
    xi, yi = pts[i]
    d = [ ( (xi-pts[j][0])**2 + (yi-pts[j][1])**2, j) for j in range(i)] + \
        [ ((xi-pts[j][0])**2 + (yi-pts[j][1])**2, j) for j in range(i+1, n)]
    d.sort()
    nn_idx[i] = [j for _, j in d[:KMAX]]

# --- curva k-dist (distanza al k-esimo vicino) ---
print('\nk  d_mean(k-esimo)  d_median  d_sd   dist_media_knn_cum')
dk = []
for k in range(1, KMAX+1):
    ds = [ math.sqrt((pts[i][0]-pts[nn_idx[i][k-1]][0])**2 + (pts[i][1]-pts[nn_idx[i][k-1]][1])**2) for i in range(n)]
    m = sum(ds)/n
    med = sorted(ds)[n//2]
    sd = (sum((d-m)**2 for d in ds)/n)**0.5
    dk.append((k, m, med, sd))
for k, m, med, sd in dk[:40]:
    print(f'{k:3d} {m:.5f} {med:.5f} {sd:.5f}')

# --- Moran I vs k ---
def moran(v, k):
    num = 0.0; S0 = 0.0
    mv = sum(v)/n
    denom = sum((a-mv)**2 for a in v)
    for i in range(n):
        wsum = 1.0/k
        for j in nn_idx[i][:k]:
            num += wsum*(v[i]-mv)*(v[j]-mv)
            S0 += wsum
    I = (n/S0)*num/denom
    EI = -1/(n-1)
    # varianza sotto randomizzazione (semplificata: normal approx)
    s1 = sum((1.0/k + 1.0/k)**2 for _ in range(0)) # placeholder
    # var normal assumption: n * S1 / ((n-1)*S0) etc. uso approssimazione normale standard
    A = n*(n*n-3*n+3)*1.0 - 0  # trascuro: uso z con formula classica sotto normalita'
    b2 = n*sum((a-mv)**4 for a in v)/denom**2
    s1v = 0.0
    # S1 = sum_i sum_j (wij+wji)^2 ; con RS-KNN wij=1/k per i->j
    # approssimo con 2*S0 (simmetrica) per z
    S1 = 2*S0
    S2 = sum((sum(1.0/k for _ in range(k)))**2 for _ in range(n))
    varI = ( (n*((n*n-3*n+3)*S1 - n*S2 + 3*S0*S0)) -
             b2*(( (n*n-n)*S1 - 2*n*S2 + 6*S0*S0 )) ) / ((n-1)*(n-2)*(n-3)*S0*S0) - EI*EI
    z = (I-EI)/varI**0.5 if varI>0 else float('nan')
    return I, z

print('\nk | I(y) z(y) | I(x1) z(x1) | I(x2) z(x2)')
for k in [3,4,5,6,7,8,9,10,12,15,18,20,25,30,35,40]:
    Iy, zy = moran(y,k); I1,z1 = moran(x1,k); I2,z2 = moran(x2,k)
    print(f'{k:3d} | {Iy:+.4f} {zy:+7.3f} | {I1:+.4f} {z1:+7.3f} | {I2:+.4f} {z2:+7.3f}')

# salva indici per step2
import pickle
pickle.dump({'nn_idx': nn_idx, 'y': y, 'x1': x1, 'x2': x2, 'n': n},
            open('tool-results/data-analysis/scelta-k/data.pkl','wb'))
print('\nOK salvato')
