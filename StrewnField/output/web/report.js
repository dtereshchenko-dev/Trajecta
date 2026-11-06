(function(root) {
  'use strict';
  function quantile(values, weights, p) {
    const rows=values.map((v,i)=>[v,weights?.[i]??1]).filter(r=>Number.isFinite(r[0])&&r[1]>0).sort((a,b)=>a[0]-b[0]);
    if(!rows.length) return null;
    if(rows.every(r=>r[1]===rows[0][1])) {
      const index=(rows.length-1)*p, i=Math.floor(index);
      return rows[i][0]+(rows[Math.min(i+1,rows.length-1)][0]-rows[i][0])*(index-i);
    }
    const total=rows.reduce((s,r)=>s+r[1],0);
    let before=0, last=rows[0][0], lastP=0;
    for(const [value,weight] of rows) {
      const at=(before+weight/2)/total;
      if(p<=at) return last+(value-last)*Math.max(0,(p-lastP)/(at-lastP));
      before+=weight;last=value;lastP=at;
    }
    return last;
  }
  function describe(values,weights) {
    let sum=0,total=0,n=0;
    values.forEach((v,i)=>{const w=weights?.[i]??1;if(Number.isFinite(v)&&w>0){sum+=v*w;total+=w;n++;}});
    return {n,mean:total?sum/total:null,p05:quantile(values,weights,.05),median:quantile(values,weights,.5),p95:quantile(values,weights,.95)};
  }
  function histogram(values,weights,bins=28,log=false) {
    const rows=values.map((v,i)=>[v,weights?.[i]??1]).filter(r=>Number.isFinite(r[0])&&r[1]>0&&(!log||r[0]>0));
    if(!rows.length) return {x:[],y:[],width:[],ranges:[]};
    let lo=Infinity,hi=-Infinity;
    for(const row of rows) {row[0]=log?Math.log10(row[0]):row[0];lo=Math.min(lo,row[0]);hi=Math.max(hi,row[0]);}
    if(hi-lo<=32*Number.EPSILON*Math.max(1,Math.abs(lo),Math.abs(hi))) {
      const value=log?10**((lo+hi)/2):(lo+hi)/2;
      return {x:[value],y:[1],width:[0],ranges:[[value,value]],singleValue:value};
    }
    const step=(hi-lo)/bins,y=Array(bins).fill(0),x=[],width=[],ranges=[];
    for(const [value,w] of rows) y[Math.min(bins-1,Math.floor((value-lo)/step))]+=w;
    const total=y.reduce((a,b)=>a+b,0);
    for(let i=0;i<bins;i++) {
      const a=log?10**(lo+i*step):lo+i*step,b=log?10**(lo+(i+1)*step):lo+(i+1)*step;
      x.push((a+b)/2);width.push((b-a)*.96);ranges.push([a,b]);y[i]/=total;
    }
    return {x,y,width,ranges};
  }
  function floats(encoded,shuffled=false) {
    const bytes=Uint8Array.from(atob(encoded),c=>c.charCodeAt(0));
    return readFloats(bytes,shuffled);
  }
  function readFloats(bytes,shuffled=false) {
    const count=bytes.length/4;
    if(!Number.isInteger(count))throw new Error('Invalid packed numeric array');
    if(shuffled) {
      const linear=new Uint8Array(bytes.length);
      for(let i=0;i<count;i++)for(let b=0;b<4;b++)linear[i*4+b]=bytes[b*count+i];
      bytes=linear;
    }
    const view=new DataView(bytes.buffer),result=new Float32Array(count);
    for(let i=0;i<result.length;i++) result[i]=view.getFloat32(i*4,true);
    return result;
  }
  async function inflate(encoded) {
    if(typeof DecompressionStream==='undefined') throw new Error('This report requires a browser with DecompressionStream support. Use a current browser.');
    const bytes=Uint8Array.from(atob(encoded),c=>c.charCodeAt(0));
    const stream=new Blob([bytes]).stream().pipeThrough(new DecompressionStream('deflate'));
    return new Uint8Array(await new Response(stream).arrayBuffer());
  }
  async function unpack(encoded) {
    return readFloats(await inflate(encoded),true);
  }
  async function readFlight(text) {
    const payload=JSON.parse(text);
    if(!payload.encoding)return payload;
    if(payload.encoding!=='base64-deflate-json')throw new Error('Unknown trajectory payload encoding');
    return JSON.parse(new TextDecoder().decode(await inflate(payload.payload)));
  }
  function windProfile(rows) {
    return rows.filter(row=>row.slice(0,3).every(Number.isFinite)).map(([height,east,north])=>{
      const speed=Math.hypot(east,north);
      return {height,speed,direction:speed<1e-6?null:(Math.atan2(east,north)*180/Math.PI+360)%360};
    }).sort((a,b)=>a.height-b.height);
  }
  const metrics=[
    ['Total landed mass','kg',r=>r.mass_balance_kg.landed],
    ['Surviving fragments','count',r=>r.ground_count],
    ['Tracked representatives','count',r=>r.tracked_count],
    ['Largest surviving fragment','kg',r=>r.largest_mass_kg],
    ['First breakup height','km WGS84',r=>r.first_breakup_km],
    ['Breakup events','weighted count',r=>r.breakup_weighted_count??r.breakups.reduce((s,b)=>s+b.multiplicity,0)],
    ['Initial mass','kg',r=>r.mass_balance_kg.initial],
    ['Ablated mass','kg',r=>r.mass_balance_kg.ablated],
    ['Dust mass','kg',r=>r.mass_balance_kg.dust],
    ['Below-cutoff mass','kg',r=>r.mass_balance_kg.below_cutoff],
    ['Unpropagated remainder','kg',r=>r.mass_balance_kg.residual]
  ];
  root.StrewnReport={quantile,describe,histogram,floats,unpack,readFlight,metrics,windProfile};
  if(typeof module!=='undefined') module.exports=root.StrewnReport;
})(globalThis);
