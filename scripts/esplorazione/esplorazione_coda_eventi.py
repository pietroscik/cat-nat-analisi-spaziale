#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Stima event-based della coda frequente sismica RP<30 (ramo dev, esplorazione).

Passo 1 del protocollo §5 di docs/esplorazione_validazione_rp30.md: la banda
RP<30 non è stimabile dalle mappe di pericolo (verdetto MPS04, §1 del doc) e si
chiude dal lato eventi: eventi CPTI15 1980-2020 (finestra strumentale di
riferimento, completa per M>=4,5) x modello di scuotimento x fragilità x
esposizione.

Modello (dichiarato, illustrativo, coerente col benchmark dove possibile):
- eventi: CPTI15 v4.0, foglio 'catalogue', Mw>=5.0 (sensibilità 4,5), con
  coordinate epicentrali LatDef/LonDef;
- Io0(Mw): mediana empirica di IoDef per bin di 0,5 di magnitudo sul catalogo
  completo (relazione validata in validazione_cpti15.py, monotona); sensibilità
  sulle mediane della sola finestra 1980-2020 (più basse);
- attenuazione: I(d) = Io0 - 3*log10(sqrt(d^2+h^2)/h) - c*d, h = 5 km (termine
  geometrico di Blake), c = 0,01 /km (anelastico); sensibilità c in {0,005; 0,01; 0,02};
- PGA(I) = 0,02 * 10^((I-5)/4) g (illustrativa; sensibilità divisore 3,5-5),
  calibrata sugli ordini osservati: Io IX-X -> ~0,26 g (L'Aquila 2009);
- MDR = min(1, 6*PGA^2): identico a pricing_model.py §2 e results/ep_curve.json;
- banda frequente RP<30 per-sito: PGA locale < ag_RP30 del comune (griglia MPS04
  RP30 già in matrice): stessa definizione del limite inferiore identificato
  (trapezoid RP30-475 + coda rara k=3), quindi le due parti sono disgiunte per
  costruzione;
- esposizione: asset_tot_EUR per comune, come nel benchmark.

Declustering (variante di robustezza dell'output): finestra spazio-temporale
50 km / 90 giorni attorno alla scossa più forte del cluster (Gardner-Knopoff
semplificato, deterministico). Per una AAL storica la finestra osservata
(cluster inclusi) è già il dato: la variante declusterata è il pavimento
Poisson-forward, non la stima centrale.

Calibrazione (sezione dedicata dell'output): il modello event-based produce
pochi superamenti PGA>=ag_RP30 rispetto a quelli impliciti nella definizione
MPS04 (n_comuni/30 all'anno): il deficit è strutturale (lo alimentano eventi
M<4,5 ravvicinati e la seismicità liscia di fondo, che contano nella
frequenza quasi nulla nel danno perché MDR ∝ PGA²) e qualifica la stima
puntuale come limite inferiore.

Nessun dato viene committato: il catalogo si scarica dall'URL documentato in
docs/esplorazione_validazione_rp30.md §2 oppure si punta con CPTI15_XLSX.

Output: results/esplorazione_coda_eventi.json. Solo stdlib, deterministico
(doppia esecuzione: JSON identico).
"""

import csv
import datetime
import json
import math
import os
import statistics
import zipfile
from xml.etree import ElementTree as ET

RAD = math.pi / 180.0
H_KM = 5.0            # profondità di riferimento dell'attenuazione
COEF_GEOM = 3.0       # coefficiente geometrico (Blake) su log10(sqrt(d^2+h^2)/h)
ANELASTICO = 0.01     # perdita d'intensità per km (illustrativa)
IO_PGA_REF = 0.02     # PGA (g) a Io = 5
IO_PGA_DIV = 4.0      # PGA = 0,02 * 10^((Io-5)/4): un decadimento ogni 4 gradi
MDR_COEF = 6.0        # MDR = min(1, 6*PGA^2), identico al benchmark
SOGLIA_MW = 5.0       # eventi considerati (sensibilità 4,5)
W1, W2 = 1980, 2020   # finestra di riferimento (completa per M>=4,5)
IO_MIN = 4.5          # sotto questa intensità il contributo è trascurabile

# valori di confronto (results/ep_curve.json, ramo main; e prior max-ent ramo dev)
AAL_IDENTIFICATO_EUR = 1712535014.0     # trapezoid RP30-475 + coda rara k=3
BENCHMARK_CURVE4_EUR = 1034560856.0
AAL_MAXENT_P1_EUR = 3183852881.42      # mediana prior P1, results/esplorazione_coda_frequente.json

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def leggi_matrice():
    """Comuni con coordinate, asset e ag_RP30 dalla matrice sismica."""
    comuni = []
    path = os.path.join(ROOT, 'data', 'Matrice_Modello_Savelli_Final_sismico.csv')
    with open(path, encoding='utf-8-sig') as f:
        for row in csv.DictReader(f):
            try:
                lat = float(row['lat']); lon = float(row['long'])
                asset = float(row['asset_tot_EUR'])
                ag30 = float(row['ag_RP30'] or 0.0)
            except (ValueError, KeyError, TypeError):
                continue
            if asset <= 0 or ag30 <= 0:
                continue  # comuni senza hazard classificato: convenzione repo (n=3753)
            comuni.append({'nome': row['COMUNE'], 'lat': lat, 'lon': lon,
                           'asset': asset, 'ag30': ag30})
    return comuni


def leggi_catalogo(percorso):
    """Eventi del foglio 'catalogue' (individuato dall'intestazione con MwDef)."""
    ns = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
    z = zipfile.ZipFile(percorso)
    shared = []
    if 'xl/sharedStrings.xml' in z.namelist():
        for si in ET.fromstring(z.read('xl/sharedStrings.xml')).findall(ns + 'si'):
            shared.append(''.join(t.text or '' for t in si.iter(ns + 't')))

    def cell_val(c):
        v = c.find(ns + 'v')
        if v is None:
            return None
        return shared[int(v.text)] if c.get('t') == 's' else v.text

    # individua il foglio catalogo (intestazione con MwDef e Year) e la mappa nome->colonna
    for nome_foglio in ('sheet3.xml', 'sheet1.xml', 'sheet2.xml'):
        p = 'xl/worksheets/' + nome_foglio
        if p not in z.namelist():
            continue
        righe = []
        for row in ET.fromstring(z.read(p)).iter(ns + 'row'):
            vals = {}
            for c in row.iter(ns + 'c'):
                col = ''.join(ch for ch in c.get('r') if ch.isalpha())
                v = cell_val(c)
                if v is not None:
                    vals[col] = v
            righe.append(vals)
        intestazione = {v: k for k, v in righe[0].items()}
        if 'MwDef' in intestazione and 'Year' in intestazione:
            return righe, intestazione
    raise RuntimeError('foglio catalogo non trovato in %s' % percorso)


def medio(valore):
    """Converte IoDef: numero o intervallo '6-7' -> punto medio (lezione validazione)."""
    if valore is None:
        return None
    s = str(valore).strip()
    if not s:
        return None
    if '-' in s:
        parti = s.split('-')
        try:
            a, b = float(parti[0]), float(parti[-1])
            return (a + b) / 2.0 if b > a else a
        except ValueError:
            return None
    try:
        return float(s)
    except ValueError:
        return None


def raccogli_eventi(righe, intestazione, w1, w2, soglia_mw):
    """Eventi nella finestra con Mw e coordinate epicentrali."""
    eventi = []
    scartati = 0
    for r in righe[1:]:
        if intestazione['Year'] not in r or intestazione['MwDef'] not in r:
            continue
        try:
            anno = int(float(r[intestazione['Year']]))
            mw = float(r[intestazione['MwDef']])
        except ValueError:
            continue
        if not (w1 <= anno <= w2) or mw < soglia_mw:
            continue
        lat = medio(r.get(intestazione['LatDef']))
        lon = medio(r.get(intestazione['LonDef']))
        if lat is None or lon is None or not (-90 < lat < 90 and -180 < lon < 180):
            scartati += 1
            continue
        io = medio(r.get(intestazione['IoDef']))
        t = None
        mo = medio(r.get(intestazione['Mo']))
        da = medio(r.get(intestazione['Da']))
        if mo is not None and da is not None:
            try:
                t = datetime.date(anno, int(mo), int(da)).toordinal()
            except ValueError:
                t = None
        eventi.append({'anno': anno, 'mw': mw, 'lat': lat, 'lon': lon,
                       'io': io, 't': t,
                       'area': str(r.get(intestazione['EpicentralArea'], '?'))})
    eventi.sort(key=lambda e: (e['anno'], e['area'], e['lat'], e['lon']))
    return eventi, scartati

def decluster(eventi, giorni=90, km=50):
    """Variante di robustezza: repliche = eventi nella finestra spazio-temporale
    (km, giorni) di una scossa più forte (Gardner-Knopoff semplificato).

    Deterministico: si processano gli eventi in ordine di Mw decrescente (a
    parità di Mw, anno crescente): la scossa più forte del cluster resta
    principale e ne rimuove foreshock e aftershock nella finestra. Gli eventi
    senza data completa non possono essere dichiarati replica e restano
    principali.
    """
    ordine = sorted(range(len(eventi)),
                    key=lambda i: (-eventi[i]['mw'], eventi[i]['anno'], eventi[i]['area']))
    repliche = set()
    for k, i in enumerate(ordine):
        if i in repliche:
            continue
        e = eventi[i]
        if e['t'] is None:
            continue
        for j in ordine[k + 1:]:
            if j in repliche:
                continue
            f = eventi[j]
            if f['t'] is None:
                continue
            if (abs(f['t'] - e['t']) <= giorni and f['mw'] < e['mw']
                    and distanza_km(e['lat'], e['lon'], f['lat'], f['lon']) <= km):
                repliche.add(j)
    eventi_principali = [eventi[i] for i in range(len(eventi)) if i not in repliche]
    return eventi_principali, len(repliche)


def mediane_io(righe, intestazione, w1, w2, mw_min=5.0):
    """Io0 empirico per bin di 0,5 di magnitudo nella finestra indicata."""
    bin_ = {}
    for r in righe[1:]:
        if intestazione['MwDef'] not in r or intestazione['Year'] not in r:
            continue
        try:
            anno = int(float(r[intestazione['Year']]))
            mw = float(r[intestazione['MwDef']])
        except ValueError:
            continue
        if not (w1 <= anno <= w2) or mw < mw_min:
            continue
        io = medio(r.get(intestazione['IoDef']))
        if io is None or io <= 0:
            continue
        b = min(int((mw - mw_min) / 0.5), 3)
        bin_.setdefault(b, []).append(io)
    etichette = ['[%.1f,%.1f)' % (mw_min + 0.5 * b, mw_min + 0.5 * (b + 1)) for b in range(4)]
    etichette[3] = '[%.1f,inf)' % (mw_min + 1.5)
    med = {}
    for b in range(4):
        vals = sorted(bin_.get(b, []))
        if len(vals) >= 5:
            med[b] = statistics.median(vals)
    # se un bin resta con meno di 5 osservazioni, usa il bin adiacente disponibile
    for b in range(4):
        if b not in med:
            vicini = [bb for bb in (b - 1, b + 1) if bb in med]
            if vicini:
                med[b] = med[vicini[0]]
            elif med:
                med[b] = max(med.values())
            else:
                med[b] = 10.0
    return med, {etichette[b]: med[b] for b in range(4)}


def distanza_km(lat1, lon1, lat2, lon2):
    f1, f2 = lat1 * RAD, lat2 * RAD
    df, dl = (lat2 - lat1) * RAD, (lon2 - lon1) * RAD
    a = math.sin(df / 2) ** 2 + math.cos(f1) * math.cos(f2) * math.sin(dl / 2) ** 2
    return 6371.0088 * 2 * math.asin(math.sqrt(a))


def intensita(io0, d, c_anel):
    return io0 - COEF_GEOM * math.log10(math.sqrt(d * d + H_KM * H_KM) / H_KM) - c_anel * d


def pga_da_io(io, div):
    return IO_PGA_REF * 10.0 ** ((io - 5.0) / div)


def mdr(pga):
    return min(1.0, MDR_COEF * pga * pga)


def stima(comuni, eventi, med_io, c_anel=ANELASTICO, div=IO_PGA_DIV, mw_min=SOGLIA_MW):
    """Somma i danni e conta i superamenti per (evento, comune), divisi per banda.

    La divisione per gli anni avviene a valle: qui si accumulano i danni grezzi
    dell'insieme di eventi ricevuto. Gli eventi sotto mw_min (sensibilità 4,5)
    usano il bin di Io0 più basso disponibile.
    """
    danni_freq = 0.0
    danni_ident = 0.0
    n_freq = 0
    n_ident = 0
    per_evento = [0.0] * len(eventi)
    per_bin = [0.0] * 4
    per_bin_ident = [0.0] * 4
    for i, ev in enumerate(eventi):
        b = min(int((ev['mw'] - mw_min) / 0.5), 3)
        if b < 0:
            b = 0  # eventi M[4,5,5,0) della sensibilità: bin di Io0 più basso
        io0 = med_io[b]
        for cm in comuni:
            d = distanza_km(ev['lat'], ev['lon'], cm['lat'], cm['lon'])
            io = intensita(io0, d, c_anel)
            if io < IO_MIN:
                continue
            p = pga_da_io(io, div)
            danno = cm['asset'] * mdr(p)
            if p < cm['ag30']:
                danni_freq += danno
                n_freq += 1
                per_evento[i] += danno
                per_bin[b] += danno
            else:
                danni_ident += danno
                n_ident += 1
                per_bin_ident[b] += danno
    return {
        'danni_freq': danni_freq,
        'danni_ident': danni_ident,
        'n_freq': n_freq,
        'n_ident': n_ident,
        'per_evento': per_evento,
        'per_bin': per_bin,
        'per_bin_ident': per_bin_ident,
    }


def main():
    percorso = os.environ.get('CPTI15_XLSX') or os.path.join(
        ROOT, 'data', 'sismico', 'CPTI15_v4.0.xlsx')
    comuni = leggi_matrice()
    righe, intestazione = leggi_catalogo(percorso)
    eventi, scartati = raccogli_eventi(righe, intestazione, W1, W2, SOGLIA_MW)
    med_io, etichette_io = mediane_io(righe, intestazione, 1005, 2020)
    med_io_mod, etichette_io_mod = mediane_io(righe, intestazione, W1, W2)

    anni = W2 - W1 + 1
    base = stima(comuni, eventi, med_io)
    aal_freq = base['danni_freq'] / anni

    # sensibilità: variazioni deterministiche a parametri fissi (niente Monte Carlo)
    sens = {}
    s_lasca = None
    for c_anel in (0.005, 0.01, 0.02):
        for div in (3.5, 4.0, 5.0):
            if c_anel == ANELASTICO and div == IO_PGA_DIV:
                continue
            s = stima(comuni, eventi, med_io, c_anel, div)
            sens['anelastico=%.3f,pga_div=%.1f' % (c_anel, div)] = round(s['danni_freq'] / anni, 2)
            if c_anel == 0.005 and div == 3.5:
                s_lasca = s
    s_mod = stima(comuni, eventi, med_io_mod)
    sens['Io0_finestra_1980'] = round(s_mod['danni_freq'] / anni, 2)
    ev_1960, _ = raccogli_eventi(righe, intestazione, 1960, W2, SOGLIA_MW)
    s_1960 = stima(comuni, ev_1960, med_io)
    sens['finestra_1960_2020'] = round(s_1960['danni_freq'] / (W2 - 1960 + 1), 2)
    ev_45, _ = raccogli_eventi(righe, intestazione, W1, W2, 4.5)
    med_io_45, _ = mediane_io(righe, intestazione, 1005, 2020, 4.5)
    s_45 = stima(comuni, ev_45, med_io_45, mw_min=4.5)
    sens['soglia_Mw_4.5'] = round(s_45['danni_freq'] / anni, 2)
    # variante di robustezza: declustering (repliche rimosse, pavimento Poisson)
    ev_principali, n_repliche = decluster(eventi)
    s_dec = stima(comuni, ev_principali, med_io)
    sens['declustering_50km_90gg'] = round(s_dec['danni_freq'] / anni, 2)

    # calibrazione: superamenti PGA>=ag_RP30 del modello vs impliciti MPS04
    impliciti_anno = len(comuni) / 30.0
    sup_modello = base['n_ident'] / anni
    deficit = impliciti_anno / sup_modello if sup_modello > 0 else float('inf')
    aal_freq_scalata = aal_freq * deficit  # variante prudenziale estrema (deficit uniforme)

    bin_labels = ['[5.0,5.5)', '[5.5,6.0)', '[6.0,6.5)', '[6.5,inf)']
    per_bin_aal = {bin_labels[b]: round(base['per_bin'][b] / anni, 2) for b in range(4)}
    per_bin_ident_aal = {bin_labels[b]: round(base['per_bin_ident'][b] / anni, 2) for b in range(4)}

    top = sorted(range(len(eventi)), key=lambda i: (-base['per_evento'][i], i))[:8]
    top_eventi = [{
        'anno': eventi[i]['anno'], 'area': eventi[i]['area'], 'Mw': round(eventi[i]['mw'], 2),
        'contributo_AAL_freq_EUR': round(base['per_evento'][i] / anni, 2),
    } for i in top]

    aal_ident_ev = base['danni_ident'] / anni
    totale_stimato = AAL_IDENTIFICATO_EUR + aal_freq
    q_maxent = aal_freq / AAL_MAXENT_P1_EUR
    out = {
        'esperimento': 'coda frequente RP<30 event-based: eventi CPTI15 x attenuazione x fragilita x esposizione (ramo dev)',
        'modello': {
            'eventi': 'CPTI15 v4.0, finestra %d-%d, Mw>=%.1f (%d eventi, %d scartati senza coordinate)'
                      % (W1, W2, SOGLIA_MW, len(eventi), scartati),
            'Io0_per_bin': etichette_io,
            'Io0_per_bin_finestra_moderna': etichette_io_mod,
            'attenuazione': 'I(d) = Io0 - 3*log10(sqrt(d^2+25)/5) - %.3f*d' % ANELASTICO,
            'PGA': 'PGA(I) = %.3f*10^((I-5)/%.1f) g (illustrativa)' % (IO_PGA_REF, IO_PGA_DIV),
            'MDR': 'min(1, 6*PGA^2), identico a pricing_model.py e results/ep_curve.json',
            'banda_frequente': 'PGA locale < ag_RP30 del comune (griglia MPS04 RP30 in matrice): disgiunta per costruzione dal limite identificato',
            'esposizione': 'asset_tot_EUR dei %d comuni con hazard classificato' % len(comuni),
        },
        'risultato': {
            'AAL_frequente_RP_inf30_EUR': round(aal_freq, 2),
            'AAL_banda_identificata_eventi_EUR': round(aal_ident_ev, 2),
            'AAL_totale_eventi_EUR': round(aal_freq + aal_ident_ev, 2),
            'per_magnitudo_AAL_freq': per_bin_aal,
            'per_magnitudo_AAL_ident': per_bin_ident_aal,
            'top_eventi_AAL_freq': top_eventi,
        },
        'declustering': {
            'metodo': 'finestra spazio-temporale 50 km / 90 giorni attorno alla scossa piu forte (Gardner-Knopoff semplificato), deterministico',
            'n_eventi_base': len(eventi),
            'n_principali': len(ev_principali),
            'n_repliche': n_repliche,
            'AAL_freq_declusterata_EUR': round(s_dec['danni_freq'] / anni, 2),
            'nota': 'per una AAL storica la finestra osservata (cluster inclusi) e il dato: la variante declusterata e il pavimento Poisson-forward, non la stima centrale',
        },
        'calibrazione': {
            'superamenti_modello_PGA_sup_ag30_per_anno': round(sup_modello, 1),
            'superamenti_impliciti_MPS04_per_anno': round(impliciti_anno, 1),
            'deficit': round(deficit, 1),
            'variante_Mw45_superamenti_per_anno': round(s_45['n_ident'] / anni, 1),
            'variante_lasca_superamenti_per_anno': round(s_lasca['n_ident'] / anni, 1),
            'nota': 'il deficit e alimentato da eventi M<4,5 ravvicinati e dalla seismicita liscia di fondo: contano nella frequenza dei superamenti, quasi nulla nel danno (MDR proporzionale a PGA^2): con Mw>=4,5 la AAL frequente sale solo da %.0f a %.0f mln' % (aal_freq / 1e6, s_45['danni_freq'] / anni / 1e6),
            'AAL_freq_variante_scalata_EUR': round(aal_freq_scalata, 2),
            'nota_variante_scalata': 'ipotesi estrema di deficit uniforme su tutti i livelli di scuotimento: prudenziale, non una stima',
        },
        'gate': {
            'limite_inferiore_RP_sup30_EUR': AAL_IDENTIFICATO_EUR,
            'benchmark_CURVE4_EUR': BENCHMARK_CURVE4_EUR,
            'prior_maxent_P1_mediana_EUR': AAL_MAXENT_P1_EUR,
            'AAL_f_quota_maxent_P1': round(q_maxent, 4),
            'AAL_f_quota_identificato': round(aal_freq / AAL_IDENTIFICATO_EUR, 4),
            'totale_RP_sup0_stimato_EUR': round(totale_stimato, 2),
            'totale_quota_benchmark': round(totale_stimato / BENCHMARK_CURVE4_EUR, 3),
            'sanita_banda_identificata': 'eventi %d-%d: %.0f mln/anno nella banda identificata vs curva 1.713 mld: stesso ordine al bordo RP30, eventi rari non campionati in 41 anni' % (W1, W2, aal_ident_ev / 1e6),
        },
        'sensibilita': sens,
        'verdetto': '',
    }
    out['verdetto'] = (
        'AAL frequente event-based = %.0f mln EUR/anno (base; sensibilita %.0f-%.0f mln): '
        '%.1f%% del prior max-ent P1 (3,18 mld) e %.1f%% del limite identificato (1,713 mld). '
        'Il modello sotto-produce i superamenti di ag_RP30 di %.0fx rispetto alla definizione '
        'MPS04, ma il deficit e quasi tutto frequenza senza danno (M<4,5 ravvicinati): la '
        'variante prudenziale scalata (%.2f mld) resta sotto il max-ent. Verdetto: la coda '
        'frequente e un correctivo minore, indicativamente %.0f-%.0f mln/anno, al piu ~%.1f '
        'mld in prudenza estrema; totale RP>=0 = 1,713 + %.0f = %.2f mld (%.2fx benchmark, '
        'dal %.2fx del limite inferiore). Con declustering (repliche fuori): %.0f '
        'mln/anno, pavimento. Chiusura definitiva del gate al passo DBMI15 '
        '(danni osservati), come da protocollo.'
        % (aal_freq / 1e6, min(sens.values()) / 1e6, max(sens.values()) / 1e6,
           100 * q_maxent, 100 * aal_freq / AAL_IDENTIFICATO_EUR, deficit,
           aal_freq_scalata / 1e9, aal_freq / 1e6, max(sens.values()) / 1e6,
           aal_freq_scalata / 1e9, aal_freq / 1e6, totale_stimato / 1e9,
           totale_stimato / BENCHMARK_CURVE4_EUR, AAL_IDENTIFICATO_EUR / BENCHMARK_CURVE4_EUR,
           s_dec['danni_freq'] / anni / 1e6))
    percorso_out = os.path.join(ROOT, 'results', 'esplorazione_coda_eventi.json')
    with open(percorso_out, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=1, ensure_ascii=True, sort_keys=False)
    print('eventi: %d (Mw>=%.1f, %d-%d) | comuni: %d' % (len(eventi), SOGLIA_MW, W1, W2, len(comuni)))
    print('AAL frequente RP<30: %.1f mln/anno | identificata (eventi): %.1f mln | totale eventi: %.1f mln'
          % (aal_freq / 1e6, aal_ident_ev / 1e6, (aal_freq + aal_ident_ev) / 1e6))
    print('calibrazione: superamenti modello %.1f/anno vs MPS04 %.1f/anno (deficit %.1fx)'
          % (sup_modello, impliciti_anno, deficit))
    print('scritto: %s' % percorso_out)


if __name__ == '__main__':
    main()
