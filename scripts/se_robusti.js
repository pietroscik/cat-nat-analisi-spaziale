// =============================================================
// SE ROBUSTI (SANDWICH HUBER–WHITE, HC1) PER IL SDM Cat-Nat
// -----------------------------------------------------------------
// Estende scripts/se_definitivi.js: stessa stima ML (modulo di
// scripts/definitivo.js, stesse tracce Monte Carlo -> i SE ML
// ricalcolati replicano esattamente results/se_definitivi.json)
// piu' covarianza sandwich robusta a eteroschedasticita' (il
// Breusch-Pagan del modello k=5 dava LM = 124,37, vedi
// results/FINAL_k5.json):
//
//   V_HC1 = (n/(n-7)) * C * S * C,   C = H_nll^{-1} (7x7),
//   S = sum_i s_i s_i' con score per osservazione della logL:
//     s_i(gamma_j) = z_ij e_i / sigma^2
//         (gamma = interc, b_PMI, b_Grandi, th_PMI, th_Grandi)
//     s_i(rho)     = e_i (Wy)_i / sigma^2 + (1/n) d/drho log|I-rhoW|
//     s_i(lns2)    = (e_i^2 / sigma^2 - 1) / 2
//   Il termine jacobiano d/drho log|I-rhoW| (non separabile per
//   osservazione) e' ripartito uniformemente: contributo O(1/n).
//   Con W assente la formula si riduce a White HC1 per OLS.
//
// Verifica interna: le stime e i SE ML ricalcolati devono coincidere
// con results/se_definitivi.json (errore relativo < 1e-9, atteso ~0).
// Output: results/se_robusti_hc1.json (deterministico, k = 5..8).
// =============================================================
const fs = require('fs');
const src = fs.readFileSync(process.env.DEFINITIVO || 'scripts/definitivo.js','utf8');
const mod = src.split('// ============================ RUN ============================')[0];
eval(mod + `
const REF = JSON.parse(fs.readFileSync('results/se_definitivi.json','utf8'));
const NAMES = ['interc','b_PMI','b_Grandi','th_PMI','th_Grandi','rho','lns2'];

// Phi normale standard (Zelen-Severo, stessa formula di chi2cdf in definitivo.js)
function PhiN(t){ const t2=1/(1+0.2316419*Math.abs(t));
  const d=0.3989423*Math.exp(-t*t/2);
  const p=d*t2*(0.3193815+t2*(-0.3565638+t2*(1.781478+t2*(-1.821256+t2*1.330274))));
  return t>0?1-p:p; }
function pTwo(z){ return 2*(1-PhiN(Math.abs(z))); }  // z grande -> 0 (underflow, < 1e-300)

// inversa 7x7 (Gauss-Jordan con pivoting, come se_definitivi.js)
function inv7(H){
  const A=H.map(r=>Float64Array.from(Array.from(r)));
  const I=Array.from({length:7},(_,i)=>Float64Array.from({length:7},(_,j)=>i===j?1:0));
  for(let c2=0;c2<7;c2++){ let piv=c2;
    for(let r=c2+1;r<7;r++)if(Math.abs(A[r][c2])>Math.abs(A[piv][c2]))piv=r;
    [A[c2],A[piv]]=[A[piv],A[c2]];[I[c2],I[piv]]=[I[piv],I[c2]];
    const f=A[c2][c2];
    for(let cc=0;cc<7;cc++){A[c2][cc]/=f;I[c2][cc]/=f;}
    for(let r=0;r<7;r++){ if(r===c2)continue; const g=A[r][c2];
      if(g!==0)for(let cc=0;cc<7;cc++){A[r][cc]-=g*A[c2][cc];I[r][cc]-=g*I[c2][cc];} } }
  return I; }
function mat7mul(A,B){ const R=[]; for(let a=0;a<7;a++){R.push(new Float64Array(7));
  for(let b2=0;b2<7;b2++){let s=0;for(let j=0;j<7;j++)s+=A[a][j]*B[j][b2];R[a][b2]=s;} } return R; }

const RES={}; let maxPDiff=0, maxSeDiff=0;
for(const kk of [5,6,7,8]){
  const w=buildW(kk);
  const sdm=sdmFull(w);   // invariato: stessa sequenza PRNG di se_definitivi.js (include gli effetti MC)
  const p0=Float64Array.from(sdm.p);
  const T=45; const tr=sdm.tr;
  const Wy=sdm.Wy,Wx1=sdm.Wx1,Wx2=sdm.Wx2;
  const nll=(q)=>{ const rho=q[5],lns2=q[6]; let sse=0;
    for(let i=0;i<n;i++){ let f=q[0]+q[1]*x1[i]+q[2]*x2[i]+q[3]*Wx1[i]+q[4]*Wx2[i];
      const e=y[i]-rho*Wy[i]-f; sse+=e*e; }
    const s2=Math.exp(lns2);
    return n/2*(Math.log(2*Math.PI)+lns2)+sse/(2*s2)-logdet(tr,rho,T); };
  // Hessiana numerica di nll (identica a se_definitivi.js: h=1e-3, 4 angoli)
  const h=1e-3;
  const f0=nll(p0);
  const H=[]; for(let r=0;r<7;r++)H.push(new Float64Array(7));
  for(let r=0;r<7;r++){ for(let c2=r;c2<7;c2++){
    if(r===c2){ const qp=Float64Array.from(p0);qp[r]+=h;
      const qm=Float64Array.from(p0);qm[r]-=h;
      H[r][c2]=(nll(qp)-2*f0+nll(qm))/(h*h); }
    else { const qpp=Float64Array.from(p0);qpp[r]+=h;qpp[c2]+=h;
      const qmm=Float64Array.from(p0);qmm[r]-=h;qmm[c2]-=h;
      const qpm=Float64Array.from(p0);qpm[r]+=h;qpm[c2]-=h;
      const qmp=Float64Array.from(p0);qmp[r]-=h;qmp[c2]+=h;
      H[r][c2]=(nll(qpp)-nll(qpm)-nll(qmp)+nll(qmm))/(4*h*h); H[c2][r]=H[r][c2]; } } }
  const C=inv7(H);
  // sanity: H*inv ~ I
  let herr=0; for(let r=0;r<7;r++)for(let c2=0;c2<7;c2++){ let s=0;
    for(let j=0;j<7;j++)s+=H[r][j]*C[j][c2]; if(r===c2)s-=1; if(Math.abs(s)>herr)herr=Math.abs(s); }
  // SE ML ricalcolati + verifica vs se_definitivi.json
  const seML=new Float64Array(7);
  for(let j=0;j<7;j++){ seML[j]=Math.sqrt(Math.max(C[j][j],0));
    maxSeDiff=Math.max(maxSeDiff,Math.abs(seML[j]-REF[kk].se[j])/REF[kk].se[j]);
    maxPDiff=Math.max(maxPDiff,Math.abs(p0[j]-REF[kk].p[j])/Math.max(1e-3,Math.abs(REF[kk].p[j]))); }
  console.log('k='+kk+' (H*inv-I)max='+herr.toExponential(2));
  // ---- sandwich: score per osservazione della log-likelihood ----
  const rho=p0[5], s2=Math.exp(p0[6]);
  const e=sdm.resid;               // y - rho*Wy - Z*gamma a p0
  let dld=0; for(let t=1;t<=T;t++) dld+=(t%2?1:-1)*Math.pow(rho,t-1)*tr[t]; // d/drho log|I-rhoW| (serie, stesse tracce)
  const dldn=dld/n;
  const S=[]; for(let a=0;a<7;a++)S.push(new Float64Array(7));
  for(let i=0;i<n;i++){
    const ei=e[i], g0=ei/s2, grho=ei*Wy[i]/s2+dldn, gl=(ei*ei/s2-1)/2;
    const g=[g0, x1[i]*g0, x2[i]*g0, Wx1[i]*g0, Wx2[i]*g0, grho, gl];
    for(let a=0;a<7;a++){ const ga=g[a];
      for(let b2=a;b2<7;b2++) S[a][b2]+=ga*g[b2]; } }
  for(let a=1;a<7;a++)for(let b2=0;b2<a;b2++)S[a][b2]=S[b2][a];
  const CS=mat7mul(C,S); const V=mat7mul(CS,C);
  const hc1=n/(n-7);
  const seHC=[],zHC=[],pHC=[],ratio=[];
  for(let j=0;j<7;j++){
    const v=Math.max(hc1*V[j][j],0); const sj=Math.sqrt(v);
    seHC.push(sj); zHC.push(p0[j]/sj); pHC.push(pTwo(p0[j]/sj));
    ratio.push(sj/seML[j]); }
  RES[kk]={p:Array.from(p0),se_ml:Array.from(seML),se_hc1:seHC,z_hc1:zHC,p_hc1:pHC,rapporto_se_hc1_su_ml:ratio,hc1_factor:hc1,names:NAMES.slice()};
  console.log('  '+'parametro'.padEnd(11)+'stima'.padStart(10)+'SE ML'.padStart(10)+'SE HC1'.padStart(10)+'z HC1'.padStart(9)+'p HC1'.padStart(10));
  for(let j=0;j<7;j++) console.log('  '+NAMES[j].padEnd(11)+p0[j].toFixed(5).padStart(10)+seML[j].toFixed(5).padStart(10)+seHC[j].toFixed(5).padStart(10)+zHC[j].toFixed(2).padStart(9)+(pHC[j]===0?'<1e-300':pHC[j].toFixed(4)).padStart(10));
}
if(maxSeDiff>1e-9||maxPDiff>1e-9) throw new Error('MISMATCH vs se_definitivi.json: se='+maxSeDiff+' p='+maxPDiff);
const FINALQ = JSON.parse(fs.readFileSync('results/FINAL_k5.json','utf8'));
const META={ n_comuni:n,
  metodo:"covarianza sandwich Huber-White robusta a eteroschedasticita' con correzione a campione finito HC1 (n/(n-7))",
  basato_su:"scripts/definitivo.js via scripts/se_robusti.js (stessa stima ML e stesse tracce Monte Carlo: stime e SE ML replicano results/se_definitivi.json)",
  score_per_osservazione:"s_i(gamma)=z_ij e_i/sigma^2; s_i(rho)=e_i(Wy)_i/sigma^2+(1/n)d/drho log|I-rhoW| (jacobiano ripartito uniformemente); s_i(lns2)=(e_i^2/sigma^2-1)/2",
  breusch_pagan_LM_k5:FINALQ.k5.BP_LM,
  breusch_pagan_fonte:"results/FINAL_k5.json: eteroschedasticita' rilevata (LM=124,37, df=5) che il sandwich corregge",
  max_rel_diff_vs_se_definitivi:Math.max(maxSeDiff,maxPDiff) };
fs.writeFileSync(process.env.OUTJSON || 'results/se_robusti_hc1.json', JSON.stringify(Object.assign({meta:META},RES), null, 1));
console.log('\\nmax rel diff vs se_definitivi.json = '+Math.max(maxSeDiff,maxPDiff).toExponential(2));
console.log('scritto: '+(process.env.OUTJSON || 'results/se_robusti_hc1.json'));
`);
