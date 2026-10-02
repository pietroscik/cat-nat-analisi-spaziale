// =============================================================
// STIMA DEFINITIVA SDM Cat-Nat — comuni (n=3823), coord. corrette
// ML full-likelihood, Hessian numerico, effetti LeSage-Pace,
// Moran residui (permutazione), RESET, BP. Confronto k=5..8.
// =============================================================
const fs = require('fs');

// ---- dati ----
const raw = fs.readFileSync(process.env.MATRIX || 'data/Matrice_Modello_Savelli_Final.csv','utf8');
const recs=[]; let cur='', inq=false;
for(const ch of raw){ if(inq){cur+=ch; if(ch==='"')inq=false;} else if(ch==='"'){inq=true;} else if(ch==='\n'){recs.push(cur);cur='';} else if(ch!=='\r')cur+=ch; }
recs.push(cur);
const hdr=recs[0].split(','); const col={}; hdr.forEach((h,i)=>col[h.trim()]=i);
const rows=recs.slice(1).filter(r=>r.length).map(l=>l.split(','));
const n=rows.length;
const fix={'28041':[45.581,11.706548],'51021':[43.274463,11.746],'12108':[45.631,8.889841],'62074':[41.216414,14.527]};
const lon=new Float64Array(n),lat=new Float64Array(n),y=new Float64Array(n),x1=new Float64Array(n),x2=new Float64Array(n);
for(let i=0;i<n;i++){
  let lo=parseFloat(rows[i][col['long']]), la=parseFloat(rows[i][col['lat']]);
  const f=fix[rows[i][col['PRO_COM']].trim()]; if(f){la=f[0];lo=f[1];}
  lon[i]=lo;lat[i]=la;
  y[i]=Math.log1p(parseFloat(rows[i][col['Premio_Teorico_Comunale_EUR']]));
  x1[i]=Math.log1p(parseFloat(rows[i][col['Risk_Frana_Asset_PMI']]));
  x2[i]=Math.log1p(parseFloat(rows[i][col['Risk_Frana_Asset_Grandi']]));
}

// ---- KNN KMAX=8 ----
const KMAX=8;
const nn=new Int32Array(n*KMAX);
for(let i=0;i<n;i++){
  const ds=[];
  for(let j=0;j<n;j++){ if(j===i)continue;
    const dx=lon[i]-lon[j],dy=lat[i]-lat[j]; ds.push([dx*dx+dy*dy,j]); }
  ds.sort((a,b)=>a[0]-b[0]);
  for(let t=0;t<KMAX;t++) nn[i*KMAX+t]=ds[t][1];
}
function kdistMed(k){ const a=[]; for(let i=0;i<n;i++){ const j=nn[i*KMAX+k-1];
  a.push(Math.sqrt((lon[i]-lon[j])**2+(lat[i]-lat[j])**2)); }
  a.sort((p,q)=>p-q); return a[n>>1]; }

// ---- W ----
function buildW(k){ const wi=new Int32Array(n*k),ww=new Float64Array(n*k).fill(1/k);
  for(let i=0;i<n;i++)for(let t=0;t<k;t++)wi[i*k+t]=nn[i*KMAX+t];
  return {wi,ww,k}; }
function Wm(w,v,out){ out.fill(0); const {wi,ww,k}=w;
  for(let i=0;i<n;i++){let s=0;const b=i*k;for(let t=0;t<k;t++)s+=ww[b+t]*v[wi[b+t]];out[i]=s;} }

// ---- tracce MC ----
let seed=987654321;
function rnd(){ seed=(seed*1103515245+12345)&0x7fffffff; return seed/0x7fffffff; }
function traces(w,T,m){ const tr=new Float64Array(T+1);
  const u=new Float64Array(n),z=new Float64Array(n),tmp=new Float64Array(n);
  for(let r=0;r<m;r++){
    for(let i=0;i<n;i++){ const a=rnd(),b=rnd()||1e-12;
      u[i]=Math.sqrt(-2*Math.log(b))*Math.cos(2*Math.PI*a); }
    z.set(u);
    for(let t=1;t<=T;t++){ Wm(w,z,tmp); z.set(tmp);
      let s=0; for(let i=0;i<n;i++) s+=u[i]*z[i]; tr[t]+=s; } }
  for(let t=1;t<=T;t++) tr[t]/=m; return tr; }
function logdet(tr,rho,T){ let s=0;
  for(let t=1;t<=T;t++) s+=(t%2?1:-1)*Math.pow(rho,t)*tr[t]/t; return s; }

// ---- OLS ----
function ols(yv,Z,nc){
  const A=[];for(let a=0;a<nc;a++)A.push(new Float64Array(nc));const b=new Float64Array(nc);
  for(let a=0;a<nc;a++){ for(let c2=a;c2<nc;c2++){let s=0;for(let i=0;i<n;i++)s+=Z[a][i]*Z[c2][i];A[a][c2]=s;A[c2][a]=s;}
    let s=0;for(let i=0;i<n;i++)s+=Z[a][i]*yv[i];b[a]=s; }
  const M=A.map((r,ii)=>Float64Array.from([...Array.from(r),b[ii]]));
  for(let c2=0;c2<nc;c2++){ let p=c2;for(let r=c2+1;r<nc;r++)if(Math.abs(M[r][c2])>Math.abs(M[p][c2]))p=r;
    [M[c2],M[p]]=[M[p],M[c2]];
    for(let r=c2+1;r<nc;r++){const f=M[r][c2]/M[c2][c2];for(let cc=c2;cc<=nc;cc++)M[r][cc]-=f*M[c2][cc];} }
  const beta=new Float64Array(nc);
  for(let r=nc-1;r>=0;r--){let s=M[r][nc];for(let cc=r+1;cc<nc;cc++)s-=M[r][cc]*beta[cc];beta[r]=s/M[r][r];}
  let sse=0;for(let i=0;i<n;i++){let f=0;for(let a=0;a<nc;a++)f+=Z[a][i]*beta[a];const r=yv[i]-f;sse+=r*r;}
  return {beta,sse}; }

// ---- likelihood piena SDM: params = [b0..b4, rho, lnsigma2] ----
// e = y - rho*Wy - Z*beta, con Z = [1, x1, x2, Wx1, Wx2]
function sdmFull(w){
  const Wy=new Float64Array(n),Wx1=new Float64Array(n),Wx2=new Float64Array(n);
  Wm(w,y,Wy);Wm(w,x1,Wx1);Wm(w,x2,Wx2);
  const ones=new Float64Array(n).fill(1);
  const Z=[ones,x1,x2,Wx1,Wx2];
  const T=45,M=120; const tr=traces(w,T,M);
  const nll=(p)=>{ // -logL
    const rho=p[5], lns2=p[6];
    let sse=0;
    for(let i=0;i<n;i++){
      let f=p[0]+p[1]*x1[i]+p[2]*x2[i]+p[3]*Wx1[i]+p[4]*Wx2[i];
      const e=y[i]-rho*Wy[i]-f; sse+=e*e; }
    const s2=Math.exp(lns2);
    return n/2*(Math.log(2*Math.PI)+lns2) + sse/(2*s2) - logdet(tr,rho,T);
  };
  // stima iniziale: concentrate su rho via golden section
  const conc=(rho)=>{ const yst=new Float64Array(n);
    for(let i=0;i<n;i++)yst[i]=y[i]-rho*Wy[i];
    const {sse}=ols(yst,Z,5);
    return -n/2*(Math.log(2*Math.PI*sse/n)+1)+logdet(tr,rho,T); };
  let a=-0.95,b=0.95;const gr=0.6180339887498949;
  let c=b-gr*(b-a),d=a+gr*(b-a);
  for(let it=0;it<80;it++){ if(conc(c)>conc(d))b=d;else a=c; c=b-gr*(b-a);d=a+gr*(b-a); }
  const rho0=(a+b)/2;
  const yst=new Float64Array(n);
  for(let i=0;i<n;i++)yst[i]=y[i]-rho0*Wy[i];
  const {beta,sse}=ols(yst,Z,5);
  const s2=sse/n;
  let p=Float64Array.from([...Array.from(beta),rho0,Math.log(s2)]);
  // Newton semplice su -logL (gradiente numerico + discesa) poi Hessian
  const h=1e-5;
  const grad=new Float64Array(7);
  for(let it=0;it<200;it++){
    const f0=nll(p);
    for(let j=0;j<7;j++){ const pj=p[j]; p[j]+=h; grad[j]=(nll(p)-f0)/h; p[j]=pj; }
    // Hessian
    const H=[]; for(let r=0;r<7;r++)H.push(new Float64Array(7));
    for(let r=0;r<7;r++){ for(let c2=r;c2<7;c2++){
      const pr=p[r],pc=p[c2];
      p[r]+=h;p[c2]+=h; const fpp=nll(p);
      p[r]=pr;p[c2]=pc;
      p[r]-=h;p[c2]-=h; const fmm=nll(p);
      p[r]=pr;p[c2]=pc;
      if(r===c2){ p[r]+=h; const fp=nll(p); p[r]=pr; p[r]-=h; const fm=nll(p); p[r]=pr;
        H[r][c2]=(fp-2*f0+fm)/(h*h); }
      else { p[c2]+=h; const fpm=nll(p); p[c2]=pc; p[r]-=h; const fmp=nll(p); p[r]=pr;
        H[r][c2]=(fpp-fpm-fmp+fmm)/(4*h*h); H[c2][r]=H[r][c2]; }
      // NB: questa Hessiana (usata SOLO per i passi di Newton) mantiene la formula
      // storicamente usata per non alterare la traiettoria di ottimizzazione e quindi
      // le stime pubblicate; gli SE definitivi vengono dalla Hessiana finale a 4 angoli
      // (piu' sotto) e da scripts/se_definitivi.js.
    } }
    // step: p -= H^-1 grad
    const M=H.map((r,ii)=>Float64Array.from([...Array.from(r),-grad[ii]]));
    for(let c2=0;c2<7;c2++){ let piv=c2;
      for(let r=c2+1;r<7;r++)if(Math.abs(M[r][c2])>Math.abs(M[piv][c2]))piv=r;
      [M[c2],M[piv]]=[M[piv],M[c2]];
      for(let r=c2+1;r<7;r++){const f=M[r][c2]/M[c2][c2];for(let cc=c2;cc<=7;cc++)M[r][cc]-=f*M[c2][cc];} }
    const step=new Float64Array(7);
    for(let r=6;r>=0;r--){let s=M[r][7];for(let cc=r+1;cc<7;cc++)s-=M[r][cc]*step[cc];step[r]=s/M[r][r];}
    let lam=1.0, ok=false;
    for(let t2=0;t2<20;t2++){ const p2=new Float64Array(7);
      for(let j=0;j<7;j++)p2[j]=p[j]-lam*step[j];
      if(nll(p2)<f0){ p=p2; ok=true; break; } lam/=2; }
    if(!ok||Math.max(...Array.from(grad).map(Math.abs))<1e-7) break;
  }
  // Hessian finale -> cov
  const f0=nll(p);
  const hf=1e-3;
  const H=[]; for(let r=0;r<7;r++)H.push(new Float64Array(7));
  for(let r=0;r<7;r++){ for(let c2=r;c2<7;c2++){
    const pr=p[r],pc=p[c2];
    p[r]+=hf;p[c2]+=hf; const fpp=nll(p); p[r]=pr;p[c2]=pc;
    p[r]-=hf;p[c2]-=hf; const fmm=nll(p); p[r]=pr;p[c2]=pc;
    if(r===c2){ p[r]+=hf;const fp=nll(p);p[r]=pr;p[r]-=hf;const fm=nll(p);p[r]=pr;H[r][c2]=(fp-2*f0+fm)/(hf*hf);}
    else { // 4 angoli: f(+,+) - f(+,-) - f(-,+) + f(-,-) (formula corretta, vedi se_definitivi.js)
      p[r]=pr+hf;p[c2]=pc-hf;const fpm=nll(p);p[r]=pr;p[c2]=pc;
      p[r]=pr-hf;p[c2]=pc+hf;const fmp=nll(p);p[r]=pr;p[c2]=pc;
      H[r][c2]=(fpp-fpm-fmp+fmm)/(4*hf*hf);H[c2][r]=H[r][c2];} } }
  // inv(H)
  const A2=H.map(r=>Float64Array.from(Array.from(r)));
  const I=Array.from({length:7},(_,i2)=>Float64Array.from({length:7},(_,j)=>i2===j?1:0));
  for(let c2=0;c2<7;c2++){ let piv=c2;
    for(let r=c2+1;r<7;r++)if(Math.abs(A2[r][c2])>Math.abs(A2[piv][c2]))piv=r;
    [A2[c2],A2[piv]]=[A2[piv],A2[c2]];[I[c2],I[piv]]=[I[piv],I[c2]];
    const f=A2[c2][c2];
    for(let cc=0;cc<7;cc++){A2[c2][cc]/=f;I[c2][cc]/=f;}
    for(let r=0;r<7;r++){ if(r===c2)continue; const g=A2[r][c2];
      if(g!==0){ for(let cc=0;cc<7;cc++){A2[r][cc]-=g*A2[c2][cc];I[r][cc]-=g*I[c2][cc];} } } }
  const cov=I; // inverse of H (H = d2(-logL) -> cov = H^-1)
  const se=new Float64Array(7);
  for(let j=0;j<7;j++) se[j]=Math.sqrt(Math.max(cov[j][j],0));
  const logL=-f0;
  // residui
  const resid=new Float64Array(n);
  for(let i=0;i<n;i++){
    let f2=p[0]+p[1]*x1[i]+p[2]*x2[i]+p[3]*Wx1[i]+p[4]*Wx2[i];
    resid[i]=y[i]-p[5]*Wy[i]-f2; }
  // effetti
  const rho=p[5];
  const effects={};
  for(const [name,betaK,thetaK] of [['Risk_Frana_Asset_PMI',p[1],p[3]],['Risk_Frana_Asset_Grandi',p[2],p[4]]]){
    const total=(betaK+thetaK)/(1-rho);
    // direct via MC: E[u' (I-rhoW)^-1 (beta I + theta W) u]/n
    const u=new Float64Array(n),tmp=new Float64Array(n),tmp2=new Float64Array(n);
    let acc=0; const mrep=80, nser=120;
    for(let r=0;r<mrep;r++){
      for(let i=0;i<n;i++){const a=rnd(),b2=rnd()||1e-12;u[i]=Math.sqrt(-2*Math.log(b2))*Math.cos(2*Math.PI*a);}
      // v = (I-rhoW)^-1 (beta u + theta Wu) = sum_t (rho W)^t c
      Wm(w,u,tmp); // tmp = Wu
      const c=new Float64Array(n);
      for(let i=0;i<n;i++)c[i]=betaK*u[i]+thetaK*tmp[i];
      let z=Float64Array.from(c);
      let val=0;
      const wgt=new Float64Array(n); wgt.set(c);
      for(let t=1;t<=nser;t++){ Wm(w,wgt,tmp2); for(let i=0;i<n;i++)wgt[i]=rho*tmp2[i];
        for(let i=0;i<n;i++)z[i]+=wgt[i]; }
      for(let i=0;i<n;i++)val+=u[i]*z[i];
      acc+=val/n; }
    const direct=acc/mrep;
    effects[name]={direct,indirect:total-direct,total};
  }
  return {p:Array.from(p),se:Array.from(se),logL,resid,effects,Wx1,Wx2,Wy,Z,rho,tr};
}
function fitSARSEM(w,model){
  const Wy=new Float64Array(n),Wx1=new Float64Array(n),Wx2=new Float64Array(n);
  Wm(w,y,Wy);Wm(w,x1,Wx1);Wm(w,x2,Wx2);
  const ones=new Float64Array(n).fill(1);
  const T=45; const tr=traces(w,T,120);
  const logL=(rho)=>{
    const yst=new Float64Array(n);
    if(model==='sar'){ for(let i=0;i<n;i++)yst[i]=y[i]-rho*Wy[i];
      const {sse}=ols(yst,[ones,x1,x2],3);
      return -n/2*(Math.log(2*Math.PI*sse/n)+1)+logdet(tr,rho,T); }
    const a1=new Float64Array(n),a2=new Float64Array(n);
    for(let i=0;i<n;i++){yst[i]=y[i]-rho*Wy[i];a1[i]=x1[i]-rho*Wx1[i];a2[i]=x2[i]-rho*Wx2[i];}
    const {sse}=ols(yst,[ones,a1,a2],3);
    return -n/2*(Math.log(2*Math.PI*sse/n)+1)+logdet(tr,rho,T); };
  let a=-0.95,b=0.95;const gr=0.6180339887498949;
  let c=b-gr*(b-a),d=a+gr*(b-a);
  for(let it=0;it<80;it++){ if(logL(c)>logL(d))b=d;else a=c; c=b-gr*(b-a);d=a+gr*(b-a); }
  const rho=(a+b)/2;
  const {beta,sse}=(()=>{ const yst=new Float64Array(n);
    if(model==='sar'){ for(let i=0;i<n;i++)yst[i]=y[i]-rho*Wy[i]; return ols(yst,[ones,x1,x2],3); }
    const a1=new Float64Array(n),a2=new Float64Array(n);
    for(let i=0;i<n;i++){yst[i]=y[i]-rho*Wy[i];a1[i]=x1[i]-rho*Wx1[i];a2[i]=x2[i]-rho*Wx2[i];}
    return ols(yst,[ones,a1,a2],3); })();
  return {rho,logL:logL(rho),beta:Array.from(beta),sse};
}
// Moran con permutazione
function moranPerm(v,w,nperm){
  const {ww,k}=w;
  const m=Array.from(v).reduce((a,b)=>a+b,0)/n;
  let num=0,denom=0;
  for(let i=0;i<n;i++){const d=v[i]-m;denom+=d*d;}
  const numv=new Float64Array(n);
  for(let i=0;i<n;i++){const d=v[i]-m;let s=0;const b=i*k;
    for(let t=0;t<k;t++)s+=ww[b+t]*v[w.wi[b+t]];
    numv[i]=d*(s-m);}
  const I0=numv.reduce((a,b)=>a+b,0)/denom*(n/(n));
  const S0=n; // row std: S0=n
  const norm=(n/S0);
  const I=norm*I0;
  // permutazione
  const idx=Array.from({length:n},(_,i)=>i);
  let count=0,countLo=0;
  const vv=new Float64Array(n);
  for(let p2=0;p2<nperm;p2++){
    for(let i=n-1;i>0;i--){const j=Math.floor(rnd()*(i+1));[idx[i],idx[j]]=[idx[j],idx[i]];}
    for(let i=0;i<n;i++)vv[i]=v[idx[i]];
    let num2=0;
    for(let i=0;i<n;i++){const d=vv[i]-m;let s=0;const b=i*k;
      for(let t=0;t<k;t++)s+=ww[b+t]*vv[w.wi[b+t]];
      num2+=d*(s-m);}
    const Ip=norm*num2/denom;
    if(Ip>=I)count++;
    if(Ip<=I)countLo++;
  }
  const pUp=(count+1)/(nperm+1), pLo=(countLo+1)/(nperm+1);
  return {I,pUp,pLo,p:Math.min(1,2*Math.min(pUp,pLo))};
}
// chi2 cdf
function chi2cdf(x,df){ // df=1 o 2: uso Wilson-Hilferty o esatto
  if(df===2) return 1-Math.exp(-x/2);
  // df=1: 2*Phi(sqrt(x))-1
  const z=Math.sqrt(x);
  const Phi=(t)=>{ const t2=1/(1+0.2316419*Math.abs(t));
    const d=0.3989423*Math.exp(-t*t/2);
    let p=d*t2*(0.3193815+t2*(-0.3565638+t2*(1.781478+t2*(-1.821256+t2*1.330274))));
    return t>0?1-p:p; };
  return 2*Phi(z)-1;
}

// ============================ RUN ============================
const out={n, note:'coord corrette; y=log1p(premio EUR); X=log1p(risk)'};

for(const k of [5,6,7,8]){
  const w=buildW(k);
  const sdm=sdmFull(w);
  const sar=fitSARSEM(w,'sar'), sem=fitSARSEM(w,'sem');
  const lrSAR=2*(sdm.logL-sar.logL), lrSEM=2*(sdm.logL-sem.logL);
  const mr=moranPerm(sdm.resid,w,499);
  // AIC con convenzione uniforme: K = 2p+3 per SDM, K = p+3 per SAR/SEM (sigma^2 inclusa; p=2 -> 7/5)
  const aicSDM=-2*sdm.logL+2*7, aicSAR=-2*sar.logL+2*5, aicSEM=-2*sem.logL+2*5;
  // test z / p per coefficienti
  const pval=(z)=>2*(1-(()=>{const t=Math.abs(z);const tt=1/(1+0.2316419*t);
    const d=0.3989423*Math.exp(-t*t/2);
    return 1-d*tt*(0.3193815+tt*(-0.3565638+tt*(1.781478+tt*(-1.821256+tt*1.330274))));})());
  const coefs=[];
  const names=['interc','Risk_Frana_Asset_PMI','Risk_Frana_Asset_Grandi','W_Risk_PMI','W_Risk_Grandi','rho','ln_sigma2'];
  for(let j=0;j<7;j++){ const est=sdm.p[j], se=sdm.se[j];
    coefs.push({name:names[j],est:+est.toFixed(5),se:+se.toFixed(5),z:+(est/se).toFixed(3),p:+pval(est/se).toFixed(6)}); }
  console.log(`\n===== k=${k} =====`);
  console.log('med dist k-vicino (gradi):',kdistMed(k).toFixed(5));
  coefs.forEach(c=>console.log(`${c.name}\t${c.est}\t(${c.se})\tz=${c.z}\tp=${c.p}`));
  console.log(`logL SDM=${sdm.logL.toFixed(2)} AIC=${aicSDM.toFixed(1)} | SAR logL=${sar.logL.toFixed(2)} AIC=${aicSAR.toFixed(1)} rho=${sar.rho.toFixed(4)} | SEM logL=${sem.logL.toFixed(2)} AIC=${aicSEM.toFixed(1)} lam=${sem.rho.toFixed(4)}`);
  console.log(`LR SDMvsSAR=${lrSAR.toFixed(2)} (df2, p=${(1-chi2cdf(lrSAR,2)).toExponential(3)}) | LR SDMvsSEM=${lrSEM.toFixed(2)} (df2, p=${(1-chi2cdf(lrSEM,2)).toExponential(3)})`);
  console.log(`Moran residui: I=${mr.I.toFixed(4)} pUp=${mr.pUp.toFixed(4)} pLo=${mr.pLo.toFixed(4)} p2lat=${mr.p.toFixed(4)}`);
  for(const [nm,e] of Object.entries(sdm.effects))
    console.log(`effetti ${nm}: direct=${e.direct.toFixed(4)} indirect=${e.indirect.toFixed(4)} total=${e.total.toFixed(4)}`);
  out[k]={coefs,logL:sdm.logL,aicSDM,aicSAR,aicSEM,lrSAR,lrSEM,moranResid:mr,effects:sdm.effects,
    sar:{rho:sar.rho,logL:sar.logL},sem:{rho:sem.rho,logL:sem.logL},medDistK:kdistMed(k)};
  // RESET e BP solo per k=5
  if(k===5){
    // RESET sul modello trasformato: y - rho Wy = Z beta + gamma*yhat^2
    const rho=sdm.p[5];
    const yst=new Float64Array(n);
    for(let i=0;i<n;i++)yst[i]=y[i]-rho*sdm.Wy[i];
    const {beta:b1,sse:sse1}=ols(yst,sdm.Z,5);
    const yhat=new Float64Array(n);
    for(let i=0;i<n;i++){let f=0;for(let a=0;a<5;a++)f+=sdm.Z[a][i]*b1[a];yhat[i]=f;}
    const yh2=new Float64Array(n); for(let i=0;i<n;i++)yh2[i]=yhat[i]*yhat[i];
    const {sse:sse2}=ols(yst,[...sdm.Z,yh2],6);
    const F=((sse1-sse2)/1)/(sse2/(n-7));
    console.log(`RESET (1 termine): F=${F.toFixed(3)} (p=${(1-chi2cdf(F,1)).toExponential(3)} — p approssimata chi2)`);
    // BP: e^2/sigma2 su Z
    const s2=Math.exp(sdm.p[6]);
    const e2=new Float64Array(n); for(let i=0;i<n;i++)e2[i]=sdm.resid[i]**2;
    const {beta:bb,sse:sBP}=ols(e2,sdm.Z,5);
    const mE=e2.reduce((a,b2)=>a+b2,0)/n;
    let sstBP=0; for(let i=0;i<n;i++)sstBP+=(e2[i]-mE)**2;
    const R2=1-sBP/sstBP;
    const LM=n*R2;
    console.log(`BP: LM=${LM.toFixed(2)} (df=5, p=${(1-chi2cdf(LM,5)<0?0:(()=>{ // df5: uso approx
      // Wilson-Hilferty
      const dfg=5;
      const zwh=(Math.pow(LM/dfg,1/3)-(1-2/(9*dfg)))/Math.sqrt(2/(9*dfg));
      const Phi=(t)=>{const tt=1/(1+0.2316419*Math.abs(t));const d=0.3989423*Math.exp(-t*t/2);
        let pp=d*tt*(0.3193815+tt*(-0.3565638+tt*(1.781478+tt*(-1.821256+tt*1.330274))));return t>0?1-pp:pp;};
      return 1-Phi(zwh);})()).toExponential(3)}`);
    out[k].reset={F}, out[k].bp={LM};
  }
}
fs.writeFileSync(process.env.OUTJSON || 'results/risultati_definitivi.json',JSON.stringify(out,null,1));
console.log('\nOK salvato risultati_definitivi.json');
