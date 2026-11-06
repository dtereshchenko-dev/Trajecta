const assert=require('node:assert/strict');
const R=require('../output/web/report.js');
assert.equal(R.quantile([0,2,4,6],null,.5),3);
assert.equal(R.quantile([null,30,40],null,.5),35);
assert.equal(R.quantile([],null,.5),null);
assert.equal(R.quantile([0,100],[0,1],.5),100);
assert.deepEqual(R.describe([0,2,4],null),{n:3,mean:2,p05:.2,median:2,p95:3.8});
assert.equal(R.describe([1,3],[1,3]).mean,2.5);
const h=R.histogram([1,1,10],[3,2,1],2,true);
assert.equal(h.y[0],5/6);assert.equal(h.y[1],1/6);
assert.deepEqual(R.histogram([null,0,-1],null,20,true).x,[]);
for(const log of [false,true]) {
  const fixed=R.histogram([59,59+1e-14],null,28,log);
  assert.ok(Math.abs(fixed.singleValue-59)<1e-12);assert.deepEqual(fixed.y,[1]);
}
assert.equal(R.histogram([2.2048,2.2058],null).singleValue,undefined,'Real mass variation must not be collapsed');
const wind=R.windProfile([[4,-1,0],[0,0,1],[3,0,-1],[2,1,0],[1,0,0],[5,NaN,0]]);
assert.deepEqual(wind.map(p=>p.direction),[0,null,90,180,270]);
assert.deepEqual(wind.map(p=>p.speed),[1,0,1,1,1]);
assert.deepEqual(R.windProfile([]),[]);
const bytes=Buffer.alloc(8);bytes.writeFloatLE(1.5,0);bytes.writeFloatLE(-3.25,4);
assert.deepEqual([...R.floats(bytes.toString('base64'))],[1.5,-3.25]);
if(process.argv[2]) {
  const fs=require('node:fs'),html=fs.readFileSync(process.argv[2],'utf8');
  const match=html.match(/<script id="field-data" type="application\/json">([\s\S]*?)<\/script>/);
  assert.ok(match);const data=JSON.parse(match[1]);
  const impactIndex=new Map(data.impacts.map(p=>[`${p[5]}:${p[6]}`,p]));
  for(const run of data.runs) {
    const marker=`<script id="${run.flight.element}" type="application/json">`,start=html.indexOf(marker)+marker.length;
    let flight=JSON.parse(html.slice(start,html.indexOf('</script>',start)));
    if(flight.encoding==='base64-deflate-json')flight=JSON.parse(require('node:zlib').inflateSync(Buffer.from(flight.payload,'base64')));
    run.flight=flight;
  }
  assert.equal(data.runs.length,data.summary.completed);
  for(const run of data.runs) {
    assert.equal(run.flight.verified,true,`Run ${run.run_id}: ${run.flight.error}`);
    assert.equal(run.flight.trace_selection?.total??run.flight.traces.length,run.particles);
    assert.equal(run.flight.trace_selection?.shown??run.flight.traces.length,run.flight.traces.length);
    if(run.flight.trace_selection)assert.ok(run.flight.traces.length<=run.flight.trace_selection.limit);
    const ids=new Set(run.flight.traces.map(t=>t.id));
    for(const t of run.flight.traces)if(t.generation>0)assert.ok(ids.has(t.parent_id),'Preview must retain ancestors');
    const samples=R.floats(require('node:zlib').inflateSync(Buffer.from(run.flight.samples,'base64')).toString('base64'),true);
    const stride=data.trace_columns.length;
    for(const t of run.flight.traces) {
      assert.ok(t.count<=96&&t.count>=1);
      assert.ok(t.offset+t.count*stride<=samples.length);
      for(let i=1;i<t.count;i++)assert.ok(samples[t.offset+i*stride]>=samples[t.offset+(i-1)*stride]);
      assert.ok(Math.abs(samples[t.offset]-t.first[0])<1e-3);
      assert.ok(Math.abs(samples[t.offset+(t.count-1)*stride]-t.last[0])<1e-3);
      if(stride>=7&&t.end_event==='ground') {
        const impact=impactIndex.get(`${run.run_id}:${t.id}`);
        if(impact){assert.ok(Math.abs(t.last[5]-impact[0])<1e-7);assert.ok(Math.abs(t.last[6]-impact[1])<1e-7);}
        else assert.ok(data.impact_preview?.shown<data.impact_preview?.total);
      }
    }
    const impacts=data.impacts.filter(p=>p[5]===run.run_id);
    if(!data.impact_preview||data.impact_preview.shown===data.impact_preview.total) {
      assert.equal(impacts.length,run.tracked_count);
      assert.ok(Math.abs(impacts.reduce((s,p)=>s+p[3],0)-run.ground_count)<1e-8);
      assert.ok(Math.abs(impacts.reduce((s,p)=>s+p[2]*p[3],0)-run.mass_balance_kg.landed)<1e-8);
    }
  }
  console.log(JSON.stringify({runs:data.runs.length,verified:data.runs.filter(r=>r.flight.verified).length,
    medianLanded:R.describe(data.runs.map(r=>r.mass_balance_kg.landed)).median,
    medianCount:R.describe(data.runs.map(r=>r.ground_count)).median,
    medianFirstBreakupKm:R.describe(data.runs.map(r=>r.first_breakup_km)).median}));
}
console.log('Report statistics and packed-trace checks passed.');
