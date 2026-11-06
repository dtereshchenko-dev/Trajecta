const assert=require('node:assert/strict');
const {deflateSync}=require('node:zlib');
const R=require('../output/web/report.js');

(async()=>{
  const flight={verified:true,traces:[{id:1,first:[0,95,12000000],last:[620,0,.002]}],breakup_rows:[{multiplicity:3}],wind:[[0,1,2,0]]};
  const text=JSON.stringify(flight);
  assert.deepEqual(await R.readFlight(text),flight,'Legacy plain payload remains readable');
  const packed=JSON.stringify({encoding:'base64-deflate-json',payload:deflateSync(text).toString('base64')});
  assert.deepEqual(await R.readFlight(packed),flight,'Compressed JSON retains double-precision endpoints');
  await assert.rejects(()=>R.readFlight('{"encoding":"unknown"}'),/Unknown/);
  await assert.rejects(()=>R.readFlight('{"encoding":"base64-deflate-json","payload":"broken"}'));
  assert.equal(R.metrics[5][2]({breakup_weighted_count:123}),123);
  assert.equal(R.metrics[5][2]({breakups:[{multiplicity:3}]}),3);
  console.log('Lazy compressed payload, legacy compatibility and complete breakup metrics passed.');
})().catch(error=>{console.error(error);process.exitCode=1;});
