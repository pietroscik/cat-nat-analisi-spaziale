#!/usr/bin/env python3
"""Esplorazione (ramo dev): banda a priori della coda frequente RP<30.

La §8.2 di docs/pricing_coerenza.md dichiara la coda RP<30 non stimabile dai
tre punti della matrice (l'integrale diverge per k>=2). Questo esperimento
formalizza l'ignoranza invece di lasciarla qualitativa: distribuzione
predittiva a massima entropia sul pendine locale k della coda frequente,
momenti stimati con Monte Carlo deterministico (LCG a seed fisso).

Modello: nella banda RP in [RP_min, 30) la curva EP e' una potenza
    L(RP) = L30 * (RP/30)^(2/k)
con L30 = loss nazionale a scenario RP30 (results/ep_curve.json) e 2/k
l'esponente ereditato dalla forma λ(ag) ∝ ag^(-k) (stessa convenzione di
ep_curve.py). Il contributo frequente all'AAL:
    AAL_f(k) = (L30 / 30^(2/k)) * [30^(2/k - 1) - RP_min^(2/k - 1)] / (2/k - 1)
(finale, chiusa; per 2/k = 1 si riduce a (L30/30)*log(30/RP_min)). Il cap di
MDR non vincola mai sotto RP30 (la saturazione mediana e' a RP ~19.600 anni),
quindi la chiusa e' esatta.

Prior di massima entropia, tre varianti:
  P1 "ancorata": lognormale su ln k con mediana = pendine mediano osservato
     sul tratto adiacente RP30-72 (3,695) e sigma dalla dispersione empirica
     p10-p90 dello stesso tratto (validazione_assunzioni.json);
  P2 "appiattita": come P1 ma con mediana = pendine RP72-475 (2,111): il prior
     piu' prudente fra quelli compatibili con i dati osservati;
  P3 "scala-free": log-uniforme su k in [1, 20], massima entropia senza ancore.

Inoltre la griglia deterministica AAL_f(k) e il calcolo inverso: quale k
renderebbe il contributo frequente compatibile con l'esperienza storica.

Output: results/esplorazione_coda_frequente.json. Ramo dev: esplorazione.
"""
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
EP_JSON = os.path.join(ROOT, 'results', 'ep_curve.json')
VAL_JSON = os.path.join(ROOT, 'results', 'validazione_assunzioni.json')
OUT_JSON = os.path.join(ROOT, 'results', 'esplorazione_coda_frequente.json')

N_MC = 200_000
SEED = 42


# ------------------------------------------------------------- LCG deterministico
class LCG:
    """PCG-style a 64 bit con seed fisso: riproducibile al byte."""

    def __init__(self, seed):
        self.s = (seed * 6364136223846793005 + 1442695040888963407) % (1 << 64)

    def u01(self):
        self.s = (self.s * 6364136223846793005 + 1442695040888963407) % (1 << 64)
        return (self.s >> 11) / float(1 << 53)


def norm_ppf(p):
    """ inversa della normale standard (Acklam, asintotica alle code)."""
    a = [-3.969683028665376e+01, 2.209471984487870e+02, -2.759285104694865e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798265878e+02,
         6.680131188771972e+01, -1.328068912466706e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.373661411651208e+00, 2.938159225199150e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]
    plow, phigh = 0.02425, 1 - 0.02425
    if p < plow:
        q = math.sqrt(-2 * math.log(p))
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
    if p > phigh:
        q = math.sqrt(-2 * math.log(1 - p))
        return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
    q = p - 0.5
    r = q * q
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / \
           (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1)


# ------------------------------------------------------------- chiusa del contributo
def aal_f(k, l30, rp_min=1.0):
    """Contributo AAL della banda [rp_min, 30) con esponente 2/k (chiusa esatta)."""
    e = 2.0 / k
    if abs(e - 1.0) < 1e-12:
        return (l30 / 30.0) * math.log(30.0 / rp_min)
    num = 30.0 ** (e - 1.0) - rp_min ** (e - 1.0)
    return (l30 / 30.0 ** e) * num / (e - 1.0)


def aal_f_quad(k, l30, rp_min=1.0, nstep=100000):
    """Verifica per trapezio in RP lineare della stessa chiusa (err atteso ~1e-9)."""
    e = 2.0 / k
    s = 0.0
    for i in range(nstep):
        ra = rp_min + (30.0 - rp_min) * i / nstep
        rb = rp_min + (30.0 - rp_min) * (i + 1) / nstep
        s += ((ra / 30.0) ** e / (ra * ra) + (rb / 30.0) ** e / (rb * rb)) / 2.0 * (rb - ra)
    return l30 * s


def quantili(v, ps=(0.1, 0.25, 0.5, 0.75, 0.9)):
    v = sorted(v)
    n = len(v)
    return {('p%d' % round(p * 100)): round(v[min(n - 1, int(p * n))], 3) for p in ps}


def main():
    with open(EP_JSON, encoding='utf-8') as f:
        ep = json.load(f)
    with open(VAL_JSON, encoding='utf-8') as f:
        val = json.load(f)

    l30 = ep['ep_nazionale_sismica']['RP30']
    tot_identificato = ep['aal_numerico_sismico']['totale_k3_EUR']
    bmk = ep['aal_numerico_sismico']['benchmark_CURVE4_EUR']
    pend = val['3_coda_frequente_identificabilita']['pendenze_locali']
    k_adj = pend['RP30_72']          # tratto adiacente alla banda mancante
    k_rar = pend['RP72_475']

    # verifica chiusa vs quadratura su tre k
    check = {}
    for k in (1.0, 2.0, 3.695):
        a, q = aal_f(k, l30), aal_f_quad(k, l30)
        check['k=%g' % k] = {'chiusa': round(a, 2), 'quadratura': round(q, 2),
                            'err_rel': round(abs(a - q) / a, 12)}

    # sigma dei prior dalla dispersione empirica p10-p90 (tratto RP30-72)
    sig1 = (math.log(k_adj['p90']) - math.log(k_adj['p10'])) / (2 * 1.2816)
    sig2 = (math.log(k_rar['p90']) - math.log(k_rar['p10'])) / (2 * 1.2816)

    priors = {
        'P1_ancorata_RP30_72': ('lognormale', math.log(k_adj['mediana']), sig1),
        'P2_appiattita_RP72_475': ('lognormale', math.log(k_rar['mediana']), sig2),
        'P3_loguniforme_1_20': ('logunif', 1.0, 20.0),
    }

    risultati = {}
    for nome, (tipo, p1, p2) in priors.items():
        rng = LCG(SEED)
        vals = []
        k_sum = 0.0
        for _ in range(N_MC):
            if tipo == 'lognormale':
                k = math.exp(p1 + p2 * norm_ppf(rng.u01()))
            else:
                k = math.exp(math.log(p1) + rng.u01() * (math.log(p2) - math.log(p1)))
            k_sum += k
            vals.append(aal_f(k, l30))
        tots = [tot_identificato + v for v in vals]
        risultati[nome] = {
            'tipo': ('lognormale mediana %.3f sigma %.3f' % (math.exp(p1), p2)) if tipo == 'lognormale'
                    else ('log-uniforme [%.0f, %.0f]' % (p1, p2)),
            'AAL_frequente_EUR': quantili(vals),
            'AAL_totale_EUR': quantili(tots),
            'rapporto_totale_su_benchmark': quantili([t / bmk for t in tots]),
            'k_medio_MC': round(k_sum / N_MC, 4),
        }
        vs, ts = sorted(vals), sorted(tots)
        print('%-22s AAL_f p50=%.0f p10=%.0f p90=%.0f  tot p50=%.0f (%.1fx bmk)' % (
            nome, vs[N_MC // 2], vs[N_MC // 10], vs[9 * N_MC // 10],
            ts[N_MC // 2], ts[N_MC // 2] / bmk))

    # griglia deterministica: la "curva dell'ignoranza"
    griglia = []
    for k in (1.0, 1.5, 2.0, 2.5, 3.0, 3.695, 4.459, 6.0, 10.0, 20.0):
        griglia.append({'k': k, 'AAL_f_EUR': round(aal_f(k, l30), 0),
                        'AAL_tot_EUR': round(tot_identificato + aal_f(k, l30), 0),
                        'x_benchmark': round((tot_identificato + aal_f(k, l30)) / bmk, 2)})

    # calcolo inverso: quale k rende AAL_f = frazione dell'identificato
    bersagli = {}
    for fr in (0.1, 0.25, 0.5, 1.0):
        target = fr * tot_identificato
        lo, hi = 0.05, 200.0
        for _ in range(200):
            mid = math.sqrt(lo * hi)
            if aal_f(mid, l30) < target:
                lo = mid
            else:
                hi = mid
        bersagli['AAL_f=%d%%_identificato' % round(fr * 100)] = {'k_richiesto': round(math.sqrt(lo * hi), 3),
                                                                'rapporto_su_pendine_osservato': round(math.sqrt(lo * hi) / k_adj['mediana'], 3)}

    # dominio di integrazione: sensibilita' a RP_min
    dom = {}
    for rpmin in (1.0, 2.0, 5.0):
        dom['RP_min=%g' % rpmin] = round(aal_f(k_adj['mediana'], l30, rpmin), 0)

    # divergenza: integrando fino a RP->0 (λ->infinito) il contributo e' finito
    # solo per k<2: quota di massa del prior nella regione divergente
    def frac_ge2(tipo, p1, p2, n=100000):
        rng = LCG(SEED)
        c = 0
        for _ in range(n):
            k = (math.exp(p1 + p2 * norm_ppf(rng.u01())) if tipo == 'lognormale'
                 else math.exp(math.log(p1) + rng.u01() * (math.log(p2) - math.log(p1))))
            if k >= 2.0:
                c += 1
        return c / n

    diverg = {nome: round(frac_ge2(*v[:3]), 4) for nome, v in priors.items()}

    out = {
        'esperimento': "banda a priori della coda frequente RP<30: momento entropico previsionato + Monte Carlo (ramo dev)",
        'ancoraggi': {'L30_nazionale_EUR': l30,
                      'AAL_identificato_RP30_inf_EUR': tot_identificato,
                      'benchmark_CURVE4_EUR': bmk,
                      'pendine_RP30_72': k_adj, 'pendine_RP72_475': k_rar},
        'chiusa': {'formula': 'AAL_f(k) = (L30/30^(2/k)) * [30^(2/k-1) - RP_min^(2/k-1)]/(2/k-1), L(RP)=L30*(RP/30)^(2/k)',
                   'cap_MDR_non_vincola': 'saturazione mediana del cap a RP ~19.600 anni: sotto RP30 MDR<1 sempre',
                   'verifica_quadratura': check},
        'prior': {'P1': 'max-ent lognormale, mediana = pendine osservato RP30-72 (3,695), sigma dalla dispersione p10-p90 empirica',
                  'P2': 'max-ent lognormale, mediana = pendine RP72-475 (2,111), il prior finito piu prudente compatibile coi dati',
                  'P3': 'max-ent log-uniforme su k in [1,20], senza ancore (scala-free)',
                  'MC': '%d estrazioni, LCG 64 bit seed %d, deterministico e riproducibile' % (N_MC, SEED)},
        'predittiva': risultati,
        'griglia_deterministica': griglia,
        'calcolo_inverso': bersagli,
        'sensibilita_dominio': dom,
        'massa_prior_k_oltre_2': diverg,
        'verdetto': '',
    }

    tot_p50_p1 = risultati['P1_ancorata_RP30_72']['AAL_totale_EUR']['p50']
    k_inv25 = bersagli['AAL_f=25%_identificato']['k_richiesto']
    out['verdetto'] = (
        "il momento entropico previsionato formalizza la non-stimabilita' invece di risolverla: anche col "
        "prior piu' informato (P1, ancorato al pendine osservato RP30-72 = 3,70) il contributo RP<30 predetto "
        "e' enorme (mediana %.2f mld) e il totale sale a %.1fx il benchmark; per riportare il contributo "
        "frequente al 25%% dell'AAL identificato servirebbe k=%.2f, circa 1/%d del pendine osservato: "
        "una curva che sotto RP30 diventi tanto piu' ripida e' in direzione OPPOSTA al gradiente osservato "
        "(k cresce verso RP piu' frequenti: 2,11 su RP72-475 -> 3,70 su RP30-72). L'extrapolazione della "
        "curva EP contraddice l'esperienza e la coda RP<30 va chiusa con dati evento (catalogo storico + "
        "fragilita'), non con la statistica della curva: l'1,71 mld resta il limite inferiore identificato."
        % (risultati['P1_ancorata_RP30_72']['AAL_frequente_EUR']['p50'] / 1e9,
           tot_p50_p1 / bmk, k_inv25, round(k_adj['mediana'] / k_inv25)))

    with open(OUT_JSON, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print('scritto:', OUT_JSON)


if __name__ == '__main__':
    main()
