// Griglia k per SDM ML — Cat-Nat comuni (n=3823), pure JS
const fs = require('fs');

// ---- CSV parse ----
const raw = fs.readFileSync('tool-results/data-analysis/scelta-k/matrice.csv','utf8');
const recs = []; let cur = '', inq = false;
for (const ch of raw) {
  if (inq) { cur += ch; if (ch === '"') inq = false; }
  else if (ch === '"') { inq = true; }
  else if (ch === '\n') { recs.push(cur); cur = ''; }
  else if (ch !== '\r') cur += ch;
}
recs.push(cur);
const hdr = recs[0].split(',');
const col = {}; hdr.forEach((h,i)=>col[h.trim()]=i);
const rows = recs.slice(1).filter(r=>r.length).map(l=>l.split(','));
const n = rows.length;
const num = (r,h)=>parseFloat(r[col[h]]);

// fix coordinate errate (lat/lon scambiati o virgola persa)
const fix = {'28041':[45.581,11.706548],'51021':[43.274463,11.746],'12108':[45.631,8.889841],'62074':[41.216414,14.527]};
let lon = new Float64Array(n), lat = new Float64Array(n), y = new Float64Array(n),
    x1 = new Float64Array(n), x2 = new Float64Array(n);
for (let i=0;i<n;i++){
  let lo = num(rows[i],'long'), la = num(rows[i],'lat');
  const f = fix[rows[i][col['PRO_COM']].trim()];
  if (f){ la=f[0]; lo=f[1]; }
  lon[i]=lo; lat[i]=la;
  y[i]=Math.log1p(num(rows[i],'Premio_Teorico_Comunale_EUR'));
  x1[i]=Math.log1p(num(rows[i],'Risk_Frana_Asset_PMI'));
  x2[i]=Math.log1p(num(rows[i],'Risk_Frana_Asset_Grandi'));
}

// ---- KNN (kmax=30) ----
const KMAX=30;
const nn = new Int32Array(n*KMAX);
const nd = new Float64Array(n*KMAX);
for (let i=0;i<n;i++){
  const ds=[];
  for (let j=0;j<n;j++){
    if (j===i) continue;
    const dx=lon[i]-lon[j], dy=lat[i]-lat[j];
    ds.push([dx*dx+dy*dy, j]);
  }
  // partial selection: sort first KMAX (full sort ok in JS)
  ds.sort((a,b)=>a[0]-b[0]);
  for (let k=0;k<KMAX;k++){ nn[i*KMAX+k]=ds[k][1]; nd[i*KMAX+k]=Math.sqrt(ds[k][0]); }
}

// ---- statistica: distanza mediana al k-esimo vicino ----
function kdist(k){
  const arr=Array.from(nd.filter((_,idx)=>(idx%KMAX)===k-1));
  arr.sort((a,b)=>a-b);
  return {median:arr[n>>1], mean:arr.reduce((a,b)=>a+b,0)/n};
}

// ---- Moran I ----
function moran(v,k){
  const m=Array.from(v).reduce((a,b)=>a+b,0)/n;
  let denom=0,num=0;
  for(let i=0;i<n;i++){const d=v[i]-m; denom+=d*d;}
  for(let i=0;i<n;i++){const d=v[i]-m;
    for(let t=0;t<k;t++) num+=d*(v[nn[i*KMAX+t]]-m);}
  return (n/k)*num/denom;
}

// ---- W sparse row-standardized ----
function buildW(k){
  const w=new Float64Array(n*k), wi=new Int32Array(n*k), wv=1.0/k;
  for(let i=0;i<n;i++) for(let t=0;t<k;t++){ wi[i*k+t]=nn[i*KMAX+t]; w[i*k+t]=wv; }
  return {wi,w,k};
}
function Wm(w,v,out){ // out = W v
  const {wi,w:ww,k}=w; out.fill(0);
  for(let i=0;i<n;i++){ let s=0; const b=i*k;
    for(let t=0;t<k;t++) s+=ww[b+t]*v[wi[b+t]];
    out[i]=s; }
}

// ---- tracce Monte Carlo tr(W^t), t=1..T ----
let seed=1234567;
function rnd(){ seed=(seed*1103515245+12345)&0x7fffffff; return seed/0x7fffffff; }
function traces(w,T,m){
  const tr=new Float64Array(T+1);
  const u=new Float64Array(n), z=new Float64Array(n), tmp=new Float64Array(n);
  const ww=w.w;
  for(let r=0;r<m;r++){
    let nrm=0;
    for(let i=0;i<n;i++){ // gauss da due uniformi
      const a=rnd(), b=rnd()||1e-12;
      u[i]=Math.sqrt(-2*Math.log(b||1e-12))*Math.cos(2*Math.PI*a); nrm+=u[i]*u[i];
    }
    z.set(u);
    for(let t=1;t<=T;t++){
      Wm(w,z,tmp); z.set(tmp);
      let s=0; for(let i=0;i<n;i++) s+=u[i]*z[i];
      tr[t]+= s;
    }
  }
  for(let t=1;t<=T;t++) tr[t]/=m;
  return tr;
}
function logdet(tr,rho,T){ // log|I-rho W| = sum_{t>=1} (-1)^{t+1} rho^t tr(W^t)/t
  let s=0;
  for(let t=1;t<=T;t++) s+= (t%2?1:-1)*Math.pow(rho,t)*tr[t]/t;
  return s;
}

// ---- OLS via normal equations ----
function ols(yv,Z,nc){ // Z: array di colonne Float64Array
  const A=[]; for(let a=0;a<nc;a++) A.push(new Float64Array(nc));
  const b=new Float64Array(nc);
  for(let a=0;a<nc;a++){
    for(let c2=a;c2<nc;c2++){ let s=0; for(let i=0;i<n;i++) s+=Z[a][i]*Z[c2][i]; A[a][c2]=s; A[c2][a]=s; }
    let s=0; for(let i=0;i<n;i++) s+=Z[a][i]*yv[i]; b[a]=s;
  }
  // gauss elim
  const M=A.map((r,ii)=>{const rr=Float64Array.from([...Array.from(r),b[ii]]); return rr;});
  for(let c2=0;c2<nc;c2++){
    let p=c2; for(let r=c2+1;r<nc;r++) if(Math.abs(M[r][c2])>Math.abs(M[p][c2])) p=r;
    [M[c2],M[p]]=[M[p],M[c2]];
    for(let r=c2+1;r<nc;r++){ const f=M[r][c2]/M[c2][c2];
      for(let cc=c2;cc<=nc;cc++) M[r][cc]-=f*M[c2][cc]; }
  }
  const beta=new Float64Array(nc);
  for(let r=nc-1;r>=0;r--){ let s=M[r][nc];
    for(let cc=r+1;cc<nc;cc++) s-=M[r][cc]*beta[cc];
    beta[r]=s/M[r][r]; }
  let sse=0;
  for(let i=0;i<n;i++){ let f=0; for(let a=0;a<nc;a++) f+=Z[a][i]*beta[a];
    const r=yv[i]-f; sse+=r*r; }
  return {beta,sse};
}

function fitModel(k, model){ // model: 'sdm' | 'sar' | 'sem'
  const w=buildW(k);
  const Wy=new Float64Array(n); Wm(w,y,Wy);
  const Wx1=new Float64Array(n); Wm(w,x1,Wx1);
  const Wx2=new Float64Array(n); Wm(w,x2,Wx2);
  const T=45, M=100;
  const tr=traces(w,T,M);
  const ones=new Float64Array(n).fill(1);
  let Z, nc;
  if(model==='sdm'){ Z=[ones,x1,x2,Wx1,Wx2]; nc=5; }
  else if(model==='sar'){ Z=[ones,x1,x2]; nc=3; }
  else { Z=null; }
  const logL=(rho)=>{
    if(model==='sem'){
      const yst=new Float64Array(n), x1st=new Float64Array(n), x2st=new Float64Array(n);
      for(let i=0;i<n;i++){ yst[i]=y[i]-rho*Wy[i]; x1st[i]=x1[i]-rho*Wx1[i]; x2st[i]=x2[i]-rho*Wx2[i]; }
      const {sse}=ols(yst,[ones,x1st,x2st],3);
      return -n/2*(Math.log(2*Math.PI*sse/n)+1)+logdet(tr,rho,T);
    }
    const yst=new Float64Array(n);
    for(let i=0;i<n;i++) yst[i]=y[i]-rho*Wy[i];
    const {sse}=ols(yst,Z,nc);
    return -n/2*(Math.log(2*Math.PI*sse/n)+1)+logdet(tr,rho,T);
  };
  // golden section
  let a=-0.95,b=0.95;
  const gr=0.6180339887498949;
  let c2=b-gr*(b-a), d=a+gr*(b-a);
  for(let it=0;it<60;it++){
    if(logL(c2)>logL(d)) b=d; else a=c2;
    c2=b-gr*(b-a); d=a+gr*(b-a);
  }
  const rho=(a+b)/2;
  const L=logL(rho);
  // SE numerico da hessiana
  const h=1e-4;
  const Lm=logL(rho-h), Lp=logL(rho+h);
  const d2=(Lp-2*L+Lm)/(h*h);
  const seRho=Math.sqrt(-1/d2);
  // beta a rho ottimale
  let beta=null, sse=null;
  if(model==='sem'){
    const yst=new Float64Array(n),x1st=new Float64Array(n),x2st=new Float64Array(n);
    for(let i=0;i<n;i++){ yst[i]=y[i]-rho*Wy[i]; x1st[i]=x1[i]-rho*Wx1[i]; x2st[i]=x2[i]-rho*Wx2[i]; }
    ({beta,sse}=ols(yst,[ones,x1st,x2st],3));
  } else {
    const yst=new Float64Array(n);
    for(let i=0;i<n;i++) yst[i]=y[i]-rho*Wy[i];
    ({beta,sse}=ols(yst,Z,nc));
  }
  const p=(model==='sdm')?6:4;
  const aic=-2*L+2*p;
  return {k,model,rho,seRho,logL:L,aic,beta:Array.from(beta),sse};
}

// ---- selezione ----
const grid=[3,4,5,6,7,8,9,10,12,15,20,25,30];
const results=[];
console.log('k | med_dk | I(y) | SDM: rho(se) logL AIC | SAR: rho logL AIC | SEM: lam logL AIC');
for(const k of grid){
  const kd=kdist(k);
  const Iy=moran(y,k);
  const sdm=fitModel(k,'sdm'), sar=fitModel(k,'sar'), sem=fitModel(k,'sem');
  const LR=2*(sdm.logL-sar.logL); // SDM vs SAR
  results.push({k,med:kd.median,Iy,sdm,sar,sem,LR});
  console.log(`${k}\t${kd.median.toFixed(5)}\t${Iy.toFixed(4)}\t`+
    `SDM ${sdm.rho.toFixed(4)}(${sdm.seRho.toFixed(4)}) ${sdm.logL.toFixed(2)} ${sdm.aic.toFixed(1)}\t`+
    `SAR ${sar.rho.toFixed(4)} ${sar.logL.toFixed(2)} ${sar.aic.toFixed(1)}\t`+
    `SEM ${sem.rho.toFixed(4)} ${sem.logL.toFixed(2)} ${sem.aic.toFixed(1)}\tLR(SDM/SAR)=${LR.toFixed(2)}`);
}
fs.writeFileSync('tool-results/data-analysis/scelta-k/grid_results.json',JSON.stringify(results,null,1));
console.log('OK');
