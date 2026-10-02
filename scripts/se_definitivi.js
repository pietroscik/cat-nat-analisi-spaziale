const fs=require('fs');
const src=fs.readFileSync(process.env.DEFINITIVO || 'scripts/definitivo.js','utf8');
const mod=src.split('// ============================ RUN ============================')[0];
eval(mod+`
const res={};
for(const k of [5,6,7,8]){
  const w=buildW(k);
  const sdm=sdmFull(w);
  const p0=Float64Array.from(sdm.p);
  const T=45; const tr=sdm.tr;
  const Wy=sdm.Wy,Wx1=sdm.Wx1,Wx2=sdm.Wx2;
  const nll=(q)=>{ const rho=q[5],lns2=q[6]; let sse=0;
    for(let i=0;i<n;i++){ let f=q[0]+q[1]*x1[i]+q[2]*x2[i]+q[3]*Wx1[i]+q[4]*Wx2[i];
      const e=y[i]-rho*Wy[i]-f; sse+=e*e; }
    const s2=Math.exp(lns2);
    return n/2*(Math.log(2*Math.PI)+lns2)+sse/(2*s2)-logdet(tr,rho,T); };
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
  const A=H.map(r=>Float64Array.from(Array.from(r)));
  const I=Array.from({length:7},(_,i2)=>Float64Array.from({length:7},(_,j)=>i2===j?1:0));
  for(let c2=0;c2<7;c2++){ let piv=c2;
    for(let r=c2+1;r<7;r++)if(Math.abs(A[r][c2])>Math.abs(A[piv][c2]))piv=r;
    [A[c2],A[piv]]=[A[piv],A[c2]];[I[c2],I[piv]]=[I[piv],I[c2]];
    const f=A[c2][c2];
    for(let cc=0;cc<7;cc++){A[c2][cc]/=f;I[c2][cc]/=f;}
    for(let r=0;r<7;r++){ if(r===c2)continue; const g=A[r][c2];
      if(g!==0)for(let cc=0;cc<7;cc++){A[r][cc]-=g*A[c2][cc];I[r][cc]-=g*I[c2][cc];} } }
  // sanity: H*inv ~ I
  let err=0; for(let r=0;r<7;r++)for(let c2=0;c2<7;c2++){ let s=0;
    for(let j=0;j<7;j++)s+=H[r][j]*I[j][c2]; if(r===c2)s-=1; err=Math.max(err,Math.abs(s)); }
  const se=Array.from({length:7},(_,i)=>Math.sqrt(Math.max(I[i][i],0)));
  const names=['interc','b_PMI','b_Grandi','th_PMI','th_Grandi','rho','lns2'];
  console.log('k='+k+' (H*inv-I)max='+err.toExponential(2));
  names.forEach((nm,j)=>console.log('  '+nm+'  est='+p0[j].toFixed(5)+'  se='+se[j].toFixed(5)+'  z='+(p0[j]/se[j]).toFixed(3)));
  res[k]={p:Array.from(p0),se,names};
}
fs.writeFileSync(process.env.OUTJSON || 'results/se_definitivi.json',JSON.stringify(res,null,1));
`);
