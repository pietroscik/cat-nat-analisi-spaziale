// =============================================================
// SDM ESTESO Cat-Nat — terzo hazard SISMICO (MPS04 INGV)
// Baseline p=2 (frana PMI/Grandi, replica FINAL_k5) vs p=4
// (+ Risk_Sismico_Asset_PMI/Grandi = ag_RP475 * asset).
// ML full-likelihood, tracce MC T=45 M=120 (come definitivo.js),
// Hessian numerico, effetti LeSage-Pace (MC mrep=80, nser=120),
// Moran residui (permutazione), RESET, BP.
// Robustezza: RP30/RP50/Sa01_RP1000/Sa01_RP2500, esclusione Sardegna, k=6,7,8.
// Convenzioni AIC: SDM K=2p+3, SAR/SEM K=p+3.
// =============================================================
const fs = require('fs');

const MATRIX = process.env.MATRIX || 'data/Matrice_Modello_Savelli_Final_sismico.csv';
const OUTJSON = process.env.OUTJSON || 'results/FINAL_sismico_k5.json';

// ---- lettura CSV (con gestione dei quoted fields) ----
const raw = fs.readFileSync(MATRIX, 'utf8');
const recs = []; let cur = '', inq = false;
for (const ch of raw) {
  if (inq) { cur += ch; if (ch === '"') inq = false; }
  else if (ch === '"') { inq = true; }
  else if (ch === '\n') { recs.push(cur); cur = ''; }
  else if (ch !== '\r') cur += ch;
}
recs.push(cur);
const hdr = recs[0].split(','); const col = {}; hdr.forEach((h, i) => col[h.trim()] = i);
const rows = recs.slice(1).filter(r => r.length).map(l => l.split(','));

// ---- coordinate-fix inline (come definitivo.js) ----
const fix = { '28041': [45.581, 11.706548], '51021': [43.274463, 11.746], '12108': [45.631, 8.889841], '62074': [41.216414, 14.527] };

// ---- vettori dati (tutti i candidati) ----
let n;
let lon, lat, y, xi = {};
function loadAll(subset) {
  // subset: array di indici di comuni (default tutti)
  const idx = subset || rows.map((_, i) => i);
  n = idx.length;
  lon = new Float64Array(n); lat = new Float64Array(n); y = new Float64Array(n);
  const mk = () => new Float64Array(n);
  xi = {
    franaPMI: mk(), franaGra: mk(), sismPMI: mk(), sismGra: mk(),
    sismPMI_RP30: mk(), sismGra_RP30: mk(),
    sismPMI_RP50: mk(), sismGra_RP50: mk(),
    sismPMI_Sa1000: mk(), sismGra_Sa1000: mk(),
    sismPMI_Sa2500: mk(), sismGra_Sa2500: mk(),
  };
  for (let a = 0; a < n; a++) {
    const i = idx[a];
    let lo = parseFloat(rows[i][col['long']]), la = parseFloat(rows[i][col['lat']]);
    const f = fix[(rows[i][col['PRO_COM']] || '').trim()]; if (f) { la = f[0]; lo = f[1]; }
    lon[a] = lo; lat[a] = la;
    y[a] = Math.log1p(parseFloat(rows[i][col['Premio_Teorico_Comunale_EUR']]));
    const aPMI = parseFloat(rows[i][col['asset_PMI_EUR']]);
    const aGra = parseFloat(rows[i][col['asset_grandi_EUR']]);
    xi.franaPMI[a] = Math.log1p(parseFloat(rows[i][col['Risk_Frana_Asset_PMI']]));
    xi.franaGra[a] = Math.log1p(parseFloat(rows[i][col['Risk_Frana_Asset_Grandi']]));
    const ag = parseFloat(rows[i][col['ag_RP475']]);
    xi.sismPMI[a] = Math.log1p(ag * aPMI);
    xi.sismGra[a] = Math.log1p(ag * aGra);
    xi.sismPMI_RP30[a] = Math.log1p(parseFloat(rows[i][col['ag_RP30']]) * aPMI);
    xi.sismGra_RP30[a] = Math.log1p(parseFloat(rows[i][col['ag_RP30']]) * aGra);
    xi.sismPMI_RP50[a] = Math.log1p(parseFloat(rows[i][col['ag_RP72']]) * aPMI);
    xi.sismGra_RP50[a] = Math.log1p(parseFloat(rows[i][col['ag_RP72']]) * aGra);
    xi.sismPMI_Sa1000[a] = Math.log1p(parseFloat(rows[i][col['Sa01_RP1000']]) * aPMI);
    xi.sismGra_Sa1000[a] = Math.log1p(parseFloat(rows[i][col['Sa01_RP1000']]) * aGra);
    xi.sismPMI_Sa2500[a] = Math.log1p(parseFloat(rows[i][col['Sa01_RP2500']]) * aPMI);
    xi.sismGra_Sa2500[a] = Math.log1p(parseFloat(rows[i][col['Sa01_RP2500']]) * aGra);
  }
}

// ---- KNN KMAX=8 ----
const KMAX = 8;
let nn;
function buildNN() {
  nn = new Int32Array(n * KMAX);
  for (let i = 0; i < n; i++) {
    const ds = [];
    for (let j = 0; j < n; j++) { if (j === i) continue;
      const dx = lon[i] - lon[j], dy = lat[i] - lat[j]; ds.push([dx * dx + dy * dy, j]); }
    ds.sort((a, b) => a[0] - b[0]);
    for (let t = 0; t < KMAX; t++) nn[i * KMAX + t] = ds[t][1];
  }
}
function kdistMed(k) { const a = []; for (let i = 0; i < n; i++) { const j = nn[i * KMAX + k - 1];
  a.push(Math.sqrt((lon[i] - lon[j]) ** 2 + (lat[i] - lat[j]) ** 2)); }
  a.sort((p, q) => p - q); return a[n >> 1]; }

// ---- W ----
function buildW(k) { const wi = new Int32Array(n * k), ww = new Float64Array(n * k).fill(1 / k);
  for (let i = 0; i < n; i++) for (let t = 0; t < k; t++) wi[i * k + t] = nn[i * KMAX + t];
  return { wi, ww, k }; }
function Wm(w, v, out) { out.fill(0); const { wi, ww, k } = w;
  for (let i = 0; i < n; i++) { let s = 0; const b = i * k; for (let t = 0; t < k; t++) s += ww[b + t] * v[wi[b + t]]; out[i] = s; } }

// ---- tracce MC (identico a definitivo.js) ----
let seed = 987654321;
function rnd() { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff; }
function traces(w, T, m) { const tr = new Float64Array(T + 1);
  const u = new Float64Array(n), z = new Float64Array(n), tmp = new Float64Array(n);
  for (let r = 0; r < m; r++) {
    for (let i = 0; i < n; i++) { const a = rnd(), b = rnd() || 1e-12;
      u[i] = Math.sqrt(-2 * Math.log(b)) * Math.cos(2 * Math.PI * a); }
    z.set(u);
    for (let t = 1; t <= T; t++) { Wm(w, z, tmp); z.set(tmp);
      let s = 0; for (let i = 0; i < n; i++) s += u[i] * z[i]; tr[t] += s; } }
  for (let t = 1; t <= T; t++) tr[t] /= m; return tr; }
function logdet(tr, rho, T) { let s = 0;
  for (let t = 1; t <= T; t++) s += (t % 2 ? 1 : -1) * Math.pow(rho, t) * tr[t] / t; return s; }

// ---- OLS ----
function ols(yv, Z, nc) {
  const A = []; for (let a = 0; a < nc; a++) A.push(new Float64Array(nc)); const b = new Float64Array(nc);
  for (let a = 0; a < nc; a++) { for (let c2 = a; c2 < nc; c2++) { let s = 0; for (let i = 0; i < n; i++) s += Z[a][i] * Z[c2][i]; A[a][c2] = s; A[c2][a] = s; }
    let s = 0; for (let i = 0; i < n; i++) s += Z[a][i] * yv[i]; b[a] = s; }
  const M = A.map((r, ii) => Float64Array.from([...Array.from(r), b[ii]]));
  for (let c2 = 0; c2 < nc; c2++) { let p = c2; for (let r = c2 + 1; r < nc; r++) if (Math.abs(M[r][c2]) > Math.abs(M[p][c2])) p = r;
    [M[c2], M[p]] = [M[p], M[c2]];
    for (let r = c2 + 1; r < nc; r++) { const f = M[r][c2] / M[c2][c2]; for (let cc = c2; cc <= nc; cc++) M[r][cc] -= f * M[c2][cc]; } }
  const beta = new Float64Array(nc);
  for (let r = nc - 1; r >= 0; r--) { let s = M[r][nc]; for (let cc = r + 1; cc < nc; cc++) s -= M[r][cc] * beta[cc]; beta[r] = s / M[r][r]; }
  let sse = 0; for (let i = 0; i < n; i++) { let f = 0; for (let a = 0; a < nc; a++) f += Z[a][i] * beta[a]; const r = yv[i] - f; sse += r * r; }
  return { beta, sse }; }

// ---- SDM generico: Xs = [x1..xp], nomi in names ----
function sdmFull(w, Xs, tr) {
  const p = Xs.length;
  const T = 45;
  const Wy = new Float64Array(n); Wm(w, y, Wy);
  const Wxs = Xs.map(x => { const o = new Float64Array(n); Wm(w, x, o); return o; });
  const ones = new Float64Array(n).fill(1);
  const Z = [ones, ...Xs, ...Wxs];       // nc = 2p+1
  const nc = 2 * p + 1;
  const nll = (par) => {
    const rho = par[nc], lns2 = par[nc + 1];
    let sse = 0;
    for (let i = 0; i < n; i++) {
      let f = par[0];
      for (let a = 0; a < p; a++) f += par[1 + a] * Xs[a][i] + par[1 + p + a] * Wxs[a][i];
      const e = y[i] - rho * Wy[i] - f; sse += e * e; }
    const s2 = Math.exp(lns2);
    return n / 2 * (Math.log(2 * Math.PI) + lns2) + sse / (2 * s2) - logdet(tr, rho, T);
  };
  // concentrate su rho (golden section)
  const conc = (rho) => { const yst = new Float64Array(n);
    for (let i = 0; i < n; i++) yst[i] = y[i] - rho * Wy[i];
    const { sse } = ols(yst, Z, nc);
    return -n / 2 * (Math.log(2 * Math.PI * sse / n) + 1) + logdet(tr, rho, T); };
  let a = -0.95, b = 0.95; const gr = 0.6180339887498949;
  let c = b - gr * (b - a), d = a + gr * (b - a);
  for (let it = 0; it < 80; it++) { if (conc(c) > conc(d)) b = d; else a = c; c = b - gr * (b - a); d = a + gr * (b - a); }
  const rho0 = (a + b) / 2;
  const yst = new Float64Array(n);
  for (let i = 0; i < n; i++) yst[i] = y[i] - rho0 * Wy[i];
  const { beta, sse } = ols(yst, Z, nc);
  const s2 = sse / n;
  const np = nc + 2;
  let par = Float64Array.from([...Array.from(beta), rho0, Math.log(s2)]);
  // Newton con Hessian numerico
  // NOTA: l'Hessiana usa la formula corretta a 4 angoli f(+,+)-f(+,-)-f(-,+)+f(-,-)
  // (il vecchio definitivo.js usava perturbazioni errate fuori diagonale: bug
  //  documentato, per questo il repo adotta gli SE di se_definitivi.js)
  const h = 1e-5;
  const hess = (par, hstep) => {
    const np_ = par.length;
    const Hh = []; for (let r = 0; r < np_; r++) Hh.push(new Float64Array(np_));
    const fbase = nll(par);
    for (let r = 0; r < np_; r++) {
      const qpr = Float64Array.from(par); qpr[r] += hstep;
      const qmr = Float64Array.from(par); qmr[r] -= hstep;
      Hh[r][r] = (nll(qpr) - 2 * fbase + nll(qmr)) / (hstep * hstep);
      for (let c2 = r + 1; c2 < np_; c2++) {
        const qpp = Float64Array.from(par); qpp[r] += hstep; qpp[c2] += hstep;
        const qmm = Float64Array.from(par); qmm[r] -= hstep; qmm[c2] -= hstep;
        const qpm = Float64Array.from(par); qpm[r] += hstep; qpm[c2] -= hstep;
        const qmp = Float64Array.from(par); qmp[r] -= hstep; qmp[c2] += hstep;
        Hh[r][c2] = (nll(qpp) - nll(qpm) - nll(qmp) + nll(qmm)) / (4 * hstep * hstep);
        Hh[c2][r] = Hh[r][c2];
      }
    }
    return Hh;
  };
  const grad = new Float64Array(np);
  for (let it = 0; it < 200; it++) {
    const f0 = nll(par);
    for (let j = 0; j < np; j++) { const pj = par[j]; par[j] += h; grad[j] = (nll(par) - f0) / h; par[j] = pj; }
    const H = hess(par, h);
    const M = H.map((r, ii) => Float64Array.from([...Array.from(r), -grad[ii]]));
    for (let c2 = 0; c2 < np; c2++) { let piv = c2;
      for (let r = c2 + 1; r < np; r++) if (Math.abs(M[r][c2]) > Math.abs(M[piv][c2])) piv = r;
      [M[c2], M[piv]] = [M[piv], M[c2]];
      for (let r = c2 + 1; r < np; r++) { const f = M[r][c2] / M[c2][c2]; for (let cc = c2; cc <= np; cc++) M[r][cc] -= f * M[c2][cc]; } }
    const step = new Float64Array(np);
    for (let r = np - 1; r >= 0; r--) { let s = M[r][np]; for (let cc = r + 1; cc < np; cc++) s -= M[r][cc] * step[cc]; step[r] = s / M[r][r]; }
    let lam = 1.0, ok = false;
    for (let t2 = 0; t2 < 20; t2++) { const p2 = new Float64Array(np);
      for (let j = 0; j < np; j++) p2[j] = par[j] - lam * step[j];
      if (nll(p2) < f0) { par = p2; ok = true; break; } lam /= 2; }
    if (!ok || Math.max(...Array.from(grad).map(Math.abs)) < 1e-7) break;
  }
  // Hessian finale -> cov (hf=1e-3; formula corretta a 4 angoli, come se_definitivi.js)
  const f0 = nll(par);
  const hf = 1e-3;
  const H = hess(par, hf);
  // inv(H) via Gauss-Jordan
  const A2 = H.map(r => Float64Array.from(Array.from(r)));
  const I = Array.from({ length: np }, (_, i2) => Float64Array.from({ length: np }, (_, j) => i2 === j ? 1 : 0));
  for (let c2 = 0; c2 < np; c2++) { let piv = c2;
    for (let r = c2 + 1; r < np; r++) if (Math.abs(A2[r][c2]) > Math.abs(A2[piv][c2])) piv = r;
    [A2[c2], A2[piv]] = [A2[piv], A2[c2]]; [I[c2], I[piv]] = [I[piv], I[c2]];
    const f = A2[c2][c2];
    for (let cc = 0; cc < np; cc++) { A2[c2][cc] /= f; I[c2][cc] /= f; }
    for (let r = 0; r < np; r++) { if (r === c2) continue; const g = A2[r][c2];
      if (g !== 0) { for (let cc = 0; cc < np; cc++) { A2[r][cc] -= g * A2[c2][cc]; I[r][cc] -= g * I[c2][cc]; } } } }
  const se = new Float64Array(np);
  for (let j = 0; j < np; j++) se[j] = Math.sqrt(Math.max(I[j][j], 0));
  const logL = -f0;
  const rho = par[nc];
  const resid = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    let f2 = par[0];
    for (let a = 0; a < p; a++) f2 += par[1 + a] * Xs[a][i] + par[1 + p + a] * Wxs[a][i];
    resid[i] = y[i] - rho * Wy[i] - f2; }
  // effetti LeSage-Pace (MC come definitivo.js)
  const effects = {};
  for (let a = 0; a < p; a++) {
    const betaK = par[1 + a], thetaK = par[1 + p + a];
    const total = (betaK + thetaK) / (1 - rho);
    const u = new Float64Array(n), tmp = new Float64Array(n), tmp2 = new Float64Array(n);
    let acc = 0; const mrep = 80, nser = 120;
    for (let r = 0; r < mrep; r++) {
      for (let i = 0; i < n; i++) { const a2 = rnd(), b2 = rnd() || 1e-12; u[i] = Math.sqrt(-2 * Math.log(b2)) * Math.cos(2 * Math.PI * a2); }
      Wm(w, u, tmp);
      const c = new Float64Array(n);
      for (let i = 0; i < n; i++) c[i] = betaK * u[i] + thetaK * tmp[i];
      const z = Float64Array.from(c);
      const wgt = Float64Array.from(c);
      for (let t = 1; t <= nser; t++) { Wm(w, wgt, tmp2); for (let i = 0; i < n; i++) wgt[i] = rho * tmp2[i];
        for (let i = 0; i < n; i++) z[i] += wgt[i]; }
      let val = 0; for (let i = 0; i < n; i++) val += u[i] * z[i];
      acc += val / n; }
    const direct = acc / mrep;
    effects[a] = { direct, indirect: total - direct, total };
  }
  return { par: Array.from(par), se: Array.from(se), logL, resid, effects, rho, Wy, Wxs, Z, p };
}

// ---- SAR/SEM generici ----
function fitSARSEM(w, Xs, tr, model) {
  const p = Xs.length; const T = 45;
  const Wy = new Float64Array(n); Wm(w, y, Wy);
  const Wxs = Xs.map(x => { const o = new Float64Array(n); Wm(w, x, o); return o; });
  const ones = new Float64Array(n).fill(1);
  const logL = (rho) => {
    const yst = new Float64Array(n);
    if (model === 'sar') { for (let i = 0; i < n; i++) yst[i] = y[i] - rho * Wy[i];
      const { sse } = ols(yst, [ones, ...Xs], p + 1);
      return -n / 2 * (Math.log(2 * Math.PI * sse / n) + 1) + logdet(tr, rho, T); }
    const Axs = Xs.map((x, a) => { const o = new Float64Array(n); for (let i = 0; i < n; i++) o[i] = x[i] - rho * Wxs[a][i]; return o; });
    for (let i = 0; i < n; i++) yst[i] = y[i] - rho * Wy[i];
    const { sse } = ols(yst, [ones, ...Axs], p + 1);
    return -n / 2 * (Math.log(2 * Math.PI * sse / n) + 1) + logdet(tr, rho, T); };
  let a = -0.95, b = 0.95; const gr = 0.6180339887498949;
  let c = b - gr * (b - a), d = a + gr * (b - a);
  for (let it = 0; it < 80; it++) { if (logL(c) > logL(d)) b = d; else a = c; c = b - gr * (b - a); d = a + gr * (b - a); }
  const rho = (a + b) / 2;
  let fit;
  if (model === 'sar') { const yst = new Float64Array(n); for (let i = 0; i < n; i++) yst[i] = y[i] - rho * Wy[i];
    fit = ols(yst, [ones, ...Xs], p + 1); }
  else { const yst = new Float64Array(n); for (let i = 0; i < n; i++) yst[i] = y[i] - rho * Wy[i];
    const Axs = Xs.map((x, a2) => { const o = new Float64Array(n); for (let i = 0; i < n; i++) o[i] = x[i] - rho * Wxs[a2][i]; return o; });
    fit = ols(yst, [ones, ...Axs], p + 1); }
  return { rho, logL: logL(rho), beta: Array.from(fit.beta), sse: fit.sse };
}

// ---- Moran con permutazione (come definitivo.js) ----
function moranPerm(v, w, nperm) {
  const { ww, k } = w;
  const m = Array.from(v).reduce((a, b) => a + b, 0) / n;
  let denom = 0;
  for (let i = 0; i < n; i++) { const d = v[i] - m; denom += d * d; }
  const numv = new Float64Array(n);
  for (let i = 0; i < n; i++) { const d = v[i] - m; let s = 0; const b = i * k;
    for (let t = 0; t < k; t++) s += ww[b + t] * v[w.wi[b + t]];
    numv[i] = d * (s - m); }
  const norm = n / n; // S0 = n (righe standardizzate)
  const I = norm * numv.reduce((a, b) => a + b, 0) / denom;
  const idx = Array.from({ length: n }, (_, i) => i);
  let count = 0, countLo = 0;
  const vv = new Float64Array(n);
  for (let p2 = 0; p2 < nperm; p2++) {
    for (let i = n - 1; i > 0; i--) { const j = Math.floor(rnd() * (i + 1)); [idx[i], idx[j]] = [idx[j], idx[i]]; }
    for (let i = 0; i < n; i++) vv[i] = v[idx[i]];
    let num2 = 0;
    for (let i = 0; i < n; i++) { const d = vv[i] - m; let s = 0; const b = i * k;
      for (let t = 0; t < k; t++) s += ww[b + t] * vv[w.wi[b + t]];
      num2 += d * (s - m); }
    const Ip = norm * num2 / denom;
    if (Ip >= I) count++;
    if (Ip <= I) countLo++;
  }
  const pUp = (count + 1) / (nperm + 1), pLo = (countLo + 1) / (nperm + 1);
  return { I, pUp, pLo, p: Math.min(1, 2 * Math.min(pUp, pLo)) };
}

// ---- p-value normali / chi2 (come definitivo.js) ----
function pval(z) { return 2 * (1 - (() => { const t = Math.abs(z); const tt = 1 / (1 + 0.2316419 * t);
  const d = 0.3989423 * Math.exp(-t * t / 2);
  return 1 - d * tt * (0.3193815 + tt * (-0.3565638 + tt * (1.781478 + tt * (-1.821256 + tt * 1.330274)))); })()); }
function chi2cdf(x, df) {
  if (df === 2) return 1 - Math.exp(-x / 2);
  if (df === 1) { const z = Math.sqrt(x);
    const Phi = (t) => { const t2 = 1 / (1 + 0.2316419 * Math.abs(t));
      const d = 0.3989423 * Math.exp(-t * t / 2);
      let p = d * t2 * (0.3193815 + t2 * (-0.3565638 + t2 * (1.781478 + t2 * (-1.821256 + t2 * 1.330274))));
      return t > 0 ? 1 - p : p; };
    return 2 * Phi(z) - 1; }
  // Wilson-Hilferty
  const zwh = (Math.pow(x / df, 1 / 3) - (1 - 2 / (9 * df))) / Math.sqrt(2 / (9 * df));
  const Phi = (t) => { const tt = 1 / (1 + 0.2316419 * Math.abs(t));
    const d = 0.3989423 * Math.exp(-t * t / 2);
    let pp = d * tt * (0.3193815 + tt * (-0.3565638 + tt * (1.781478 + tt * (-1.821256 + tt * 1.330274))));
    return t > 0 ? 1 - pp : pp; };
  return Phi(zwh);
}

// ---- correlazioni tra regressori ----
function corr(a, b) { let ma = 0, mb = 0; for (let i = 0; i < n; i++) { ma += a[i]; mb += b[i]; } ma /= n; mb /= n;
  let num = 0, da = 0, db = 0;
  for (let i = 0; i < n; i++) { num += (a[i] - ma) * (b[i] - mb); da += (a[i] - ma) ** 2; db += (b[i] - mb) ** 2; }
  return num / Math.sqrt(da * db); }

// ================================ RUN ================================
const out = { note: 'SDM esteso con terzo hazard sismico (MPS04 INGV); y=log1p(premio EUR); X=log1p(risk)' };

function reportCoefs(sdm, names, label) {
  const p = sdm.p; const nc = 2 * p + 1;
  const all = ['interc', ...names.map(x => x), ...names.map(x => 'W_' + x), 'rho', 'ln_sigma2'];
  const coefs = [];
  for (let j = 0; j < nc + 2; j++) { const est = sdm.par[j], se = sdm.se[j];
    coefs.push({ name: all[j], est: +est.toFixed(5), se: +se.toFixed(5), z: +(est / se).toFixed(3), p: +pval(est / se).toFixed(6) }); }
  return coefs;
}

// ---------- A. setup completo ----------
loadAll();
buildNN();
out.n = n;
console.log('n =', n);

const w5 = buildW(5);
const tr5 = traces(w5, 45, 120);
out.medDistK = { k5: kdistMed(5) };

// correlazioni log1p(risk) (multicollinearita' attesa)
out.corr = {
  sismPMI_franaPMI: +corr(xi.sismPMI, xi.franaPMI).toFixed(3),
  sismGra_franaGra: +corr(xi.sismGra, xi.franaGra).toFixed(3),
  sismPMI_sismGra: +corr(xi.sismPMI, xi.sismGra).toFixed(3),
  franaPMI_franaGra: +corr(xi.franaPMI, xi.franaGra).toFixed(3),
};

// ---------- B. baseline p=2 (verifica vs FINAL_k5) ----------
const X2 = [xi.franaPMI, xi.franaGra];
const sdm2 = sdmFull(w5, X2, tr5);
const sar2 = fitSARSEM(w5, X2, tr5, 'sar'), sem2 = fitSARSEM(w5, X2, tr5, 'sem');
out.baseline_p2 = {
  coefs: reportCoefs(sdm2, ['Risk_Frana_Asset_PMI', 'Risk_Frana_Asset_Grandi']),
  logL: sdm2.logL, aic: -2 * sdm2.logL + 2 * 7,
  sar: { rho: sar2.rho, logL: sar2.logL, aic: -2 * sar2.logL + 2 * 4 },
  sem: { rho: sem2.rho, logL: sem2.logL, aic: -2 * sem2.logL + 2 * 4 },
  lrSAR: 2 * (sdm2.logL - sar2.logL), lrSEM: 2 * (sdm2.logL - sem2.logL),
  effects: { franaPMI: sdm2.effects[0], franaGra: sdm2.effects[1] },
};
console.log('\n== BASELINE p=2 (atteso: rho~0.4966, b_Grandi~0.1434, logL~-6351.7, AIC~12717.4) ==');
console.log('rho=%.4f  b_Grandi=%.5f  logL=%.2f  AIC=%.1f  LR_SAR=%.2f  LR_SEM=%.2f',
  sdm2.rho, sdm2.par[2], sdm2.logL, out.baseline_p2.aic, out.baseline_p2.lrSAR, out.baseline_p2.lrSEM);

// ---------- C. modello esteso p=4 ----------
const X4 = [xi.franaPMI, xi.franaGra, xi.sismPMI, xi.sismGra];
const names4 = ['Risk_Frana_Asset_PMI', 'Risk_Frana_Asset_Grandi', 'Risk_Sismico_Asset_PMI', 'Risk_Sismico_Asset_Grandi'];
const sdm4 = sdmFull(w5, X4, tr5);
const sar4 = fitSARSEM(w5, X4, tr5, 'sar'), sem4 = fitSARSEM(w5, X4, tr5, 'sem');
const K4 = { sdm: 11, sarsem: 7 };
out.sdm_p4 = {
  coefs: reportCoefs(sdm4, names4),
  logL: sdm4.logL, aic: -2 * sdm4.logL + 2 * K4.sdm,
  sar: { rho: sar4.rho, logL: sar4.logL, aic: -2 * sar4.logL + 2 * K4.sarsem },
  sem: { rho: sem4.rho, logL: sem4.logL, aic: -2 * sem4.logL + 2 * K4.sarsem },
  lr_vs_SAR: 2 * (sdm4.logL - sar4.logL), lr_vs_SEM: 2 * (sdm4.logL - sem4.logL),
  effects: { franaPMI: sdm4.effects[0], franaGra: sdm4.effects[1], sismPMI: sdm4.effects[2], sismGra: sdm4.effects[3] },
};
// nested: SDM4 vs SDM2 (df = 4: 2 beta + 2 theta sismici)
const lrNested = 2 * (sdm4.logL - sdm2.logL);
out.sdm_p4.lr_vs_baseline_p2 = lrNested;
out.sdm_p2_logL = sdm2.logL; out.sdm_p2_aic = out.baseline_p2.aic;
console.log('\n== SDM p=4 (frana + sismico) ==');
console.log('rho=%.4f  logL=%.2f  AIC=%.1f', sdm4.rho, sdm4.logL, out.sdm_p4.aic);
for (const c of out.sdm_p4.coefs) console.log(`${c.name}\t${c.est}\t(${c.se})\tz=${c.z}\tp=${c.p}`);
console.log('LR SDM4 vs SDM2 = %.2f (df=4, p=%.3e)', lrNested, 1 - chi2cdf(lrNested, 4));
console.log('LR SDM4 vs SAR4 = %.2f (df=4) | LR SDM4 vs SEM4 = %.2f (df=4)', out.sdm_p4.lr_vs_SAR, out.sdm_p4.lr_vs_SEM);
console.log('AIC: SDM4=%.1f SDM2=%.1f SAR4=%.1f SEM4=%.1f', out.sdm_p4.aic, out.baseline_p2.aic, out.sdm_p4.sar.aic, out.sdm_p4.sem.aic);
for (const [nm, e] of Object.entries(out.sdm_p4.effects))
  console.log(`effetti ${nm}: direct=${e.direct.toFixed(4)} indirect=${e.indirect.toFixed(4)} total=${e.total.toFixed(4)}`);

// diagnostica k=5 sul modello esteso
const mr4 = moranPerm(sdm4.resid, w5, 499);
out.sdm_p4.moranResid = mr4;
// RESET
{
  const rho = sdm4.rho;
  const yst = new Float64Array(n);
  for (let i = 0; i < n; i++) yst[i] = y[i] - rho * sdm4.Wy[i];
  const Z = sdm4.Z; const nc = 9;
  const { beta: b1, sse: sse1 } = ols(yst, Z, nc);
  const yhat = new Float64Array(n);
  for (let i = 0; i < n; i++) { let f = 0; for (let a = 0; a < nc; a++) f += Z[a][i] * b1[a]; yhat[i] = f; }
  const yh2 = new Float64Array(n); for (let i = 0; i < n; i++) yh2[i] = yhat[i] * yhat[i];
  const { sse: sse2 } = ols(yst, [...Z, yh2], nc + 1);
  const F = ((sse1 - sse2) / 1) / (sse2 / (n - (nc + 2)));
  out.sdm_p4.reset = { F };
  console.log('RESET F=%.3f | Moran residui I=%.4f p=%.4f', F, mr4.I, mr4.p);
}
// BP
{
  const s2 = Math.exp(sdm4.par[10]);
  const e2 = new Float64Array(n); for (let i = 0; i < n; i++) e2[i] = sdm4.resid[i] ** 2 / s2;
  const Z = sdm4.Z; const nc = 9;
  const { sse: sBP } = ols(e2, Z, nc);
  const mE = e2.reduce((a, b) => a + b, 0) / n;
  let sstBP = 0; for (let i = 0; i < n; i++) sstBP += (e2[i] - mE) ** 2;
  const LM = n * (1 - sBP / sstBP);
  out.sdm_p4.bp = { LM, df: nc };
  console.log('BP LM=%.2f (df=%d)', LM, nc);
}

// ---------- D. robustezza: hazard sismico alternativo (stessa W) ----------
out.robust_hazard = {};
const variants = [
  ['RP30', xi.sismPMI_RP30, xi.sismGra_RP30],
  ['RP50', xi.sismPMI_RP50, xi.sismGra_RP50],
  ['Sa01_RP1000', xi.sismPMI_Sa1000, xi.sismGra_Sa1000],
  ['Sa01_RP2500', xi.sismPMI_Sa2500, xi.sismGra_Sa2500],
];
for (const [tag, xp, xg] of variants) {
  const s = sdmFull(w5, [xi.franaPMI, xi.franaGra, xp, xg], tr5);
  out.robust_hazard[tag] = {
    rho: s.rho, logL: s.logL, aic: -2 * s.logL + 2 * 11,
    b_franaGra: s.par[2], b_sismGra: s.par[4], b_sismPMI: s.par[3],
    th_sismGra: s.par[6], th_sismPMI: s.par[5],
    se_sismGra: s.se[4], se_sismPMI: s.se[3],
  };
  console.log('robust %s: rho=%.4f b_sismGra=%.5f (se %.5f) b_sismPMI=%.5f logL=%.2f',
    tag, s.rho, s.par[4], s.se[4], s.par[3], s.logL);
}

// ---------- E. robustezza: k=6,7,8 ----------
out.robust_k = {};
for (const k of [6, 7, 8]) {
  const wk = buildW(k);
  const trk = traces(wk, 45, 120);
  const s = sdmFull(wk, X4, trk);
  out.robust_k[k] = {
    rho: s.rho, logL: s.logL, aic: -2 * s.logL + 2 * 11,
    b_franaGra: s.par[2], b_sismGra: s.par[4], b_sismPMI: s.par[3],
    th_sismGra: s.par[6], th_sismPMI: s.par[5],
    se_sismGra: s.se[4], se_sismPMI: s.se[3],
    medDistK: kdistMed(k),
  };
  console.log('robust k=%d: rho=%.4f b_franaGra=%.5f b_sismGra=%.5f (se %.5f)', k, s.rho, s.par[2], s.par[4], s.se[4]);
}

// ---------- F. robustezza: esclusione Sardegna ----------
{
  const idx = [];
  for (let i = 0; i < rows.length; i++) {
    if (parseInt(rows[i][col['COD_REG']]) !== 20) idx.push(i);
  }
  loadAll(idx);
  buildNN();
  const wns = buildW(5);
  const trns = traces(wns, 45, 120);
  const s = sdmFull(wns, [xi.franaPMI, xi.franaGra, xi.sismPMI, xi.sismGra], trns);
  out.robust_noSardegna = {
    n: n, rho: s.rho, logL: s.logL, aic: -2 * s.logL + 2 * 11,
    b_franaGra: s.par[2], b_sismGra: s.par[4], b_sismPMI: s.par[3],
    th_sismGra: s.par[6], th_sismPMI: s.par[5],
    se_sismGra: s.se[4], se_sismPMI: s.se[3],
    effects: { franaGra: s.effects[1], sismGra: s.effects[3] },
  };
  console.log('robust noSardegna: n=%d rho=%.4f b_franaGra=%.5f b_sismGra=%.5f (se %.5f)',
    n, s.rho, s.par[2], s.par[4], s.se[4]);
}

fs.writeFileSync(OUTJSON, JSON.stringify(out, null, 1));
console.log('\nOK salvato', OUTJSON);
