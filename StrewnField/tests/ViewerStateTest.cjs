const assert=require('node:assert/strict');
const fs=require('node:fs');
const vm=require('node:vm');
const path=require('node:path');
const report=fs.readFileSync(process.argv[2],'utf8');
let payload=report.match(/<script id="field-data" type="application\/json">([\s\S]*?)<\/script>/)[1];
const data=JSON.parse(payload),elements=new Map(),graphs=new Map(),events={},panes=new Map(),tileEvents={},calls={print:0,pixels:0,flightLoads:0,fit:0,resize:0,exports:0,onlineAdded:0,onlineCreated:0};
if(process.env.STREWN_BASEMAP_TEST==='offline')data.options.offline=true;
if(process.env.STREWN_BASEMAP_TEST==='none')data.options.basemap='none';
payload=JSON.stringify(data);
assert.ok(report.includes(fs.readFileSync(path.join(__dirname,'../output/web/viewer.js'),'utf8')),'Report must embed the tested viewer source');
assert.ok(report.includes(fs.readFileSync(path.join(__dirname,'../output/web/print.js'),'utf8')),'Report must embed the tested print source');
assert.ok(report.includes(fs.readFileSync(path.join(__dirname,'../output/web/graticule.js'),'utf8')),'Report must embed the tested coordinate grid');
assert.ok(report.includes(fs.readFileSync(path.join(__dirname,'../output/web/map-flight.js'),'utf8')),'Report must embed the tested animation controls');
class ClassList {
  constructor(){this.values=new Set();}
  add(value){this.values.add(value);}remove(value){this.values.delete(value);}contains(value){return this.values.has(value);}
  toggle(value,on){if(on??!this.contains(value))this.add(value);else this.remove(value);}
}
class Element {
  constructor(id='') {this.id=id;this.children=[];this.style={};this.dataset={};this.value='';this.checked=false;this.attrs={};this.textContent='';this.classList=new ClassList();}
  get visible(){const section=this.section&&elements.get(this.section);return !section||section.classList.contains('active');}
  get clientWidth(){return this.visible?650:0;}get clientHeight(){return this.visible?370:0;}get offsetWidth(){return this.clientWidth;}
  set id(id){this._id=id;if(id)elements.set(id,this);}get id(){return this._id;}
  set innerHTML(text){for(const match of text.matchAll(/id="([^"]+)"/g))new Element(match[1]);}
  append(...values){this.children.push(...values);}appendChild(v){this.append(v);}replaceChildren(...v){this.children=v;}add(v){this.children.push(v);}
  setAttribute(k,v){this.attrs[k]=v;}createCaption(){const n=new Element();this.append(n);return n;}createTHead(){return this.createCaption();}createTBody(){return this.createCaption();}insertRow(){return this.createCaption();}insertCell(){return this.createCaption();}
  removeAttribute(k){delete this.attrs[k];}querySelectorAll(){return [];}cloneNode(){const n=new Element();n.textContent=this.textContent;n.children=this.children.map(c=>c.cloneNode());return n;}
  removeAllListeners(){}on(){}
  toDataURL(){return 'data:image/png;base64,dGVzdA==';}
  getContext(){return {createImageData:(w,h)=>({data:new Uint8ClampedArray(w*h*4)}),putImageData:image=>{calls.pixels=image.data.filter(v=>v!==0).length;},
    fillRect(){},fillText(){},drawImage(){},strokeRect(){},beginPath(){},moveTo(){},lineTo(){},stroke(){},arc(){},fill(){},
    save(){},restore(){},rect(){},clip(){},setLineDash(){},scale(){},measureText:s=>({width:s.length*7})};}
}
const template=fs.readFileSync(path.join(__dirname,'../output/web/viewer.html'),'utf8');
const rates=html=>[...html.match(/<select id="map-rate"[^>]*>([\s\S]*?)<\/select>/)[1].matchAll(/<option value="([^"]+)"/g)].map(m=>m[1]);
assert.deepEqual(rates(template),['1','5','20','100']);
assert.deepEqual(rates(report),rates(template));
for(const match of template.matchAll(/id="([^"]+)"/g))new Element(match[1]);
for(const section of template.matchAll(/<section id="([^"]+)"[^>]*>([\s\S]*?)<\/section>/g)) {
  elements.get(section[1]).section=section[1];
  for(const match of section[2].matchAll(/id="([^"]+)"/g))elements.get(match[1]).section=section[1];
}
elements.get('overview').classList.add('active');
for(const run of data.runs) {
  const marker=`<script id="${run.flight.element}" type="application/json">`,start=report.indexOf(marker)+marker.length;
  new Element(run.flight.element).textContent=report.slice(start,report.indexOf('</script>',start));
}
const tabs=['overview','field','flight','wind'].map(name=>{const n=new Element();n.dataset.view=name;return n;});
const body=new Element();
const document={body,addEventListener:(name,handler)=>events['doc:'+name]=handler,getElementById:id=>elements.get(id),createElement:()=>new Element(),querySelectorAll:q=>q==='.view'?tabs.map(t=>elements.get(t.dataset.view)):tabs};
elements.get('field-data').textContent=payload;
for(const [id,value] of Object.entries({quantity:'count',opacity:'.65','point-scope':'all',fragment:'all',start:'breakup','print-scope':'current','map-rate':'1'}))elements.get(id).value=value;
for(const id of ['density','hpd','log-time','coordinates'])elements.get(id).checked=true;
const group=()=>({clearLayers(){},addTo(){return this;},bindPopup(){return this;},bindTooltip(){return this;},on(){return this;}});
const map={setView(){},createPane:name=>{const pane=new Element();panes.set(name,pane);return pane;},getPane:()=>new Element(),getSize:()=>({x:elements.get('map').clientWidth,y:elements.get('map').clientHeight}),containerPointToLayerPoint:p=>p,
  containerPointToLatLng:p=>({lng:data.bounds?data.bounds[0][1]+p[0]/650*(data.bounds[1][1]-data.bounds[0][1]):0,
    lat:data.bounds?data.bounds[1][0]-p[1]/370*(data.bounds[1][0]-data.bounds[0][0]):0}),
  latLngToContainerPoint:p=>({x:650*(p[1]-data.bounds[0][1])/(data.bounds[1][1]-data.bounds[0][1]),y:370*(data.bounds[1][0]-p[0])/(data.bounds[1][0]-data.bounds[0][0])}),
  on(names,handler){for(const name of names.split(' '))(events['map:'+name]??=[]).push(handler);},fitBounds(){calls.fit++;assert.ok(elements.get('map').clientWidth>36&&elements.get('map').clientHeight>36,'Leaflet fitBounds cannot fit a hidden map with padding');},invalidateSize(){calls.resize++;},attributionControl:{addAttribution(){}}};
const control=()=>({addTo(){this.onAdd();}});control.scale=()=>({addTo(){}});
const L={map:options=>map,control,imageOverlay:(url,bounds,options)=>{assert.equal(options.pane,'cached-imagery');return group();},
  tileLayer:(url,options)=>{calls.onlineCreated++;assert.ok(url.includes('/tile/{z}/{y}/{x}'));assert.equal(options.pane,'online-imagery');
    assert.ok(options.maxNativeZoom>data.tiles[0]?.z||!data.tiles.length);
    return {on(name,handler){tileEvents[name]=handler;return this;},addTo(){calls.onlineAdded++;return this;}};},
  layerGroup:group,geoJSON:group,circleMarker:group,DomUtil:{create:()=>new Element(),setPosition(){}}};
const Plotly={async react(node,traces,layout){graphs.set(node.id,{traces,layout});},async relayout(){throw new Error('Printing must not resize live charts');},Plots:{async resize(){}},
  async toImage(spec,options){if(calls.failExport){calls.failExport=false;throw new Error('Test export failure');}calls.exports++;
    assert.equal(options.format,'svg');assert.equal(spec.layout.width,520);assert.equal(spec.layout.autosize,false);
    assert.ok(spec.layout.height>=260);return 'data:image/svg+xml,'+encodeURIComponent('<svg xmlns="http://www.w3.org/2000/svg" width="520" height="260"/>');}};
let animateFrames=false,nextFrame=0;const animationQueue=new Map();
const context={document,window:{addEventListener:(name,handler)=>events[name]=handler,print:()=>{if(calls.failPrint)throw new Error('Test print failure');calls.print++;}},
  JSON:{stringify:JSON.stringify,parse:text=>{const value=JSON.parse(text);if(value.traces)calls.flightLoads++;return value;}},
  L,Plotly,StrewnReport:{...require('../output/web/report.js'),readFlight:async text=>{calls.flightLoads++;return require('../output/web/report.js').readFlight(text);}},Option:class extends Element{constructor(text,value){super();this.textContent=text;this.value=value;}},
  Image:class extends Element{set src(value){this._src=value;queueMicrotask(()=>this.onload());}},
  structuredClone,setTimeout,clearTimeout,requestAnimationFrame:handler=>{if(!animateFrames){handler(0);return 0;}const id=++nextFrame;animationQueue.set(id,handler);return id;},cancelAnimationFrame:id=>animationQueue.delete(id),console};
(async()=>{
  vm.runInNewContext(fs.readFileSync(path.join(__dirname,'../output/web/graticule.js'),'utf8'),context);
  vm.runInNewContext(fs.readFileSync(path.join(__dirname,'../output/web/trajectory.js'),'utf8'),context);
  vm.runInNewContext(fs.readFileSync(path.join(__dirname,'../output/web/map-flight.js'),'utf8'),context);
  let lastGround=null,lastProfile=null,profileFrames=0;
  const drawGround=context.StrewnTrajectory.draw,drawProfile=context.StrewnTrajectory.drawProfile;
  context.StrewnTrajectory.draw=(...args)=>{const states=drawGround(...args);lastGround={model:args[1],time:args[2],states};return states;};
  context.StrewnTrajectory.drawProfile=(...args)=>{
    assert.equal(args[1],lastGround.model);assert.equal(args[2],lastGround.time);assert.equal(args[6],lastGround.states);
    profileFrames++;lastProfile=drawProfile(...args);return lastProfile;
  };
  let grids=0;const drawGrid=context.StrewnGraticule.draw;
  context.StrewnGraticule.draw=(...args)=>{grids++;return drawGrid(...args);};
  vm.runInNewContext(fs.readFileSync(path.join(__dirname,'../output/web/print.js'),'utf8'),context);
  let printContext;
  const prepare=context.StrewnPrint.prepare;
  context.StrewnPrint.prepare=async args=>{printContext=args;return prepare(args);};
  await vm.runInNewContext(fs.readFileSync(path.join(__dirname,'../output/web/viewer.js'),'utf8'),context);
  assert.equal(body.dataset.ready,'true',elements.get('warnings').textContent);
  if(data.entry_summary) {
    const entryTable=elements.get('entry-statistics').children[0];
    assert.equal(entryTable.children[0].textContent,'Atmospheric entry');
    assert.equal(entryTable.children[2].children.length,10);
    assert.ok(elements.get('entry-note').textContent.includes(`${data.entry_summary.positive_weight_entries} positive-weight entries`));
    assert.ok(elements.get('entry-note').textContent.includes('not atmospheric MC resampling'));
  }
  assert.equal(calls.flightLoads,1,'Only selected realization is decoded');
  assert.equal(elements.get('run').children.length,data.runs.length);
  assert.equal(graphs.size,16);
  const online=data.options.basemap==='esri_imagery'&&!data.options.offline;
  assert.equal(calls.onlineCreated,online?1:0);assert.equal(calls.onlineAdded,0,'Hidden map must not download tiles around its initial world view');
  assert.equal(panes.get('cached-imagery').style.zIndex,200);assert.equal(panes.get('online-imagery').style.zIndex,201);
  await elements.get('print').onclick();assert.equal(calls.print,1);await events.afterprint();
  await tabs.find(t=>t.dataset.view==='field').onclick();
  assert.equal(calls.onlineAdded,online?1:0);
  if(online){tileEvents.tileerror();assert.ok(elements.get('map-meta').textContent.includes('unavailable'));tileEvents.loading();tileEvents.load();assert.ok(elements.get('map-meta').textContent.includes('Online imagery'));}
  assert.ok(calls.pixels>0);
  const gridCanvas=panes.get('coordinates').children[0];
  assert.equal(gridCanvas.style.display,'block');assert.ok(grids>0);
  let previousGrid=grids;elements.get('coordinates').checked=false;elements.get('coordinates').onchange();
  assert.equal(gridCanvas.style.display,'none');assert.equal(grids,previousGrid);
  elements.get('density').checked=false;elements.get('density').onchange();
  elements.get('coordinates').checked=true;elements.get('coordinates').onchange();
  assert.equal(grids,previousGrid+1,'Grid must not depend on density visibility');
  for(const event of ['move','zoomend','resize']) {
    previousGrid=grids;for(const handler of events['map:'+event])handler();
    assert.equal(grids,previousGrid+1,'Grid must update with '+event);
  }
  elements.get('density').checked=true;elements.get('density').onchange();
  const index=Math.min(124,data.runs.length-1),run=data.runs[index];
  const flight=await require('../output/web/report.js').readFlight(elements.get(run.flight.element).textContent);
  await tabs.find(t=>t.dataset.view==='flight').onclick();
  elements.get('run').value=String(index);await elements.get('run').onchange();
  assert.equal(calls.flightLoads,2);
  assert.ok(elements.get('run-meta').textContent.includes('Run '+run.run_id));
  assert.ok(elements.get('run-meta').textContent.includes('Landed '+run.mass_balance_kg.landed.toLocaleString('en-US',{maximumSignificantDigits:6})));
  assert.equal(elements.get('fragment').children.length,flight.traces.length+1);
  assert.equal(graphs.get('mass').layout.yaxis.type,'log');
  assert.equal(graphs.get('mass').layout.xaxis.type,'log');
  elements.get('log-time').checked=false;await elements.get('log-time').onchange();
  assert.equal(graphs.get('mass').layout.xaxis.type,'linear');
  const fragment=flight.traces.find(t=>t.end_event==='ground');
  if(fragment){elements.get('fragment').value=String(fragment.id);await elements.get('fragment').onchange();
    for(const trace of graphs.get('mass').traces)for(const text of trace.text??[])if(text)assert.ok(text.includes(`Fragment ${fragment.id},`));}
  await tabs.find(t=>t.dataset.view==='wind').onclick();
  assert.ok(elements.get('wind-meta').textContent.includes('Run '+run.run_id));
  assert.equal(graphs.get('wind-components').traces[0].x[0],flight.wind[0][1]);
  const polar=graphs.get('wind-direction');
  assert.equal(polar.layout.polar.angularaxis.rotation,90);assert.equal(polar.layout.polar.angularaxis.direction,'clockwise');
  assert.ok(polar.traces.every(t=>t.type==='scatterpolar'));
  const windMarkers=polar.traces.find(t=>t.marker?.showscale);
  assert.equal(windMarkers.r.length,flight.wind.length);
  for(let i=0;i<flight.wind.length;i++) {
    const row=flight.wind[i],speed=Math.hypot(row[1],row[2]);
    assert.equal(windMarkers.r[i],row[0]);assert.equal(windMarkers.marker.color[i],speed);
    assert.equal(windMarkers.theta[i],speed<1e-6?null:(Math.atan2(row[1],row[2])*180/Math.PI+360)%360);
  }
  assert.ok(polar.traces.some(t=>t.marker?.symbol==='square'));assert.ok(polar.traces.some(t=>t.marker?.symbol==='triangle-up'));
  elements.get('print-scope').value='all';await elements.get('print').onclick();assert.equal(calls.print,2);await events.afterprint();
  assert.equal(elements.get('print').disabled,false);
  assert.equal(body.classList.contains('print-ready'),false);
  for(const tab of tabs)for(const scope of ['current','all']) {
    await tab.onclick();elements.get('print-scope').value=scope;
    const before=calls.print,fit=calls.fit,resize=calls.resize;await elements.get('print').onclick();
    assert.equal(calls.print,before+1,`${tab.dataset.view}, ${scope}`);
    assert.equal(body.classList.contains('print-ready'),true);
    const entryPage=data.entry_summary?1:0;
    assert.equal(elements.get('print-pages').children.length,scope==='all'?8+entryPage:({overview:4+entryPage,field:1,flight:2,wind:1})[tab.dataset.view]);
    if(data.entry_summary&&(scope==='all'||tab.dataset.view==='overview')) {
      const entry=elements.get('print-pages').children[0];
      assert.ok(entry.children[0].textContent.endsWith('Atmospheric entry'));
      assert.equal(entry.children[2].children[0].children[2].children.length,10);
    }
    events.resize();assert.equal(calls.resize,resize,'Print preparation must preserve live map dimensions');
    await events.afterprint();assert.equal(calls.fit,fit,'Printing must not change map zoom');
    assert.equal(elements.get('print-pages').children.length,0);
    assert.equal(elements.get('print-error').textContent,'');
    assert.equal(body.classList.contains('print-ready'),false);
    assert.equal(elements.get(tab.dataset.view).classList.contains('active'),true);
    assert.equal(elements.get('print').disabled,false);
  }
  for(const failure of ['failExport','failPrint']) {
    calls[failure]=true;const before=calls.print;await elements.get('print').onclick();
    assert.equal(calls.print,before);
    assert.ok(elements.get('print-error').textContent.startsWith('Print failed: Test '));
    assert.equal(body.classList.contains('print-ready'),false);
    assert.equal(elements.get('print').disabled,false);
    calls[failure]=false;await elements.get('print').onclick();await events.afterprint();
    assert.equal(calls.print,before+1);assert.equal(elements.get('print-error').textContent,'');
  }
  assert.equal(elements.get('warnings').textContent,data.warnings.join('\n'));
  if(data.trace_columns.length>=7) {
    await tabs.find(t=>t.dataset.view==='field').onclick();elements.get('trajectory').checked=true;await elements.get('trajectory').onchange();
    assert.equal(elements.get('map-playback').hidden,false);assert.equal(elements.get('map-play').disabled,false);
    assert.ok(elements.get('map-flight-status').textContent.includes('Ground + side projections'));
    assert.equal(elements.get('profile-view').hidden,false);assert.equal(elements.get('map-ground-title').hidden,false);
    assert.ok(elements.get('map-views').classList.contains('two-plane'));assert.ok(profileFrames>0);
    const step=now=>{const jobs=[...animationQueue.values()];animationQueue.clear();for(const job of jobs)job(now);};
    animateFrames=true;elements.get('map-play').onclick();step(0);step(100);
    assert.ok(Number(elements.get('map-time').value)>0);assert.equal(elements.get('map-play').attrs['aria-pressed'],'true');
    elements.get('map-play').onclick();assert.equal(animationQueue.size,0);
    const playbackTimes=[lastGround.model.start,lastGround.model.firstBreakup+6].filter(t=>t<lastGround.model.end);
    for(const start of playbackTimes)for(const rate of rates(template)) {
      elements.get('map-time').value=String(start);elements.get('map-time').oninput();
      elements.get('map-rate').value=rate;
      const before=Number(elements.get('map-time').value);
      elements.get('map-play').onclick();step(0);step(100);
      assert.ok(Math.abs(Number(elements.get('map-time').value)-Math.min(lastGround.model.end,before+.1*Number(rate)))<1e-9,
        `Fixed ${rate}x playback at ${start}s must not change after breakup`);
      if(elements.get('map-play').attrs['aria-pressed']==='true')elements.get('map-play').onclick();
      assert.equal(animationQueue.size,0);
    }
    elements.get('map-rate').value='1';
    elements.get('map-time').value=elements.get('map-time').max;elements.get('map-time').oninput();
    assert.ok(elements.get('map-flight-status').textContent.includes('0 airborne'));
    if(lastGround.states.length) {
      const point=lastProfile.point(lastGround.states[0].row);
      elements.get('map-profile').onclick({offsetX:point.x,offsetY:point.y});
      assert.equal(elements.get('profile-details').hidden,false);
      assert.ok(elements.get('profile-details').children[0].children[0].textContent.includes('Fragment'));
    }
    elements.get('map-reset').onclick();assert.equal(Number(elements.get('map-time').value),0);
    assert.equal(elements.get('profile-details').hidden,true);
    const breakup=lastGround.model.markers.find(m=>m.kind==='breakup');
    if(breakup) {
      const point=lastProfile.point(breakup.row);
      elements.get('map-profile').onclick({offsetX:point.x,offsetY:point.y});
      assert.equal(lastGround.time,breakup.row[0]);assert.equal(elements.get('profile-details').hidden,false);
      assert.ok(elements.get('profile-details').children[0].children[0].textContent.includes('Breakup'));
    }
    const framesBeforeResize=profileFrames;for(const handler of events['map:resize'])handler();step(200);
    assert.ok(profileFrames>framesBeforeResize);
    elements.get('map-run').value='0';await elements.get('map-run').onchange();
    assert.ok(elements.get('map-flight-status').textContent.startsWith('Run '+data.runs[0].run_id));
    assert.equal(lastGround.time,Number(elements.get('map-time').min));
    step(300);
    elements.get('map-play').onclick();document.hidden=true;events['doc:visibilitychange']();document.hidden=false;assert.equal(animationQueue.size,0);
    elements.get('map-play').onclick();events.beforeprint();assert.equal(animationQueue.size,0);
    elements.get('map-play').onclick();await tabs.find(t=>t.dataset.view==='wind').onclick();assert.equal(animationQueue.size,0);
    animateFrames=false;elements.get('trajectory').checked=false;await elements.get('trajectory').onchange();assert.equal(elements.get('map-playback').hidden,true);
    assert.equal(elements.get('profile-view').hidden,true);assert.equal(elements.get('map-ground-title').hidden,true);
    assert.equal(elements.get('map-views').classList.contains('two-plane'),false);
  }
  await tabs.find(t=>t.dataset.view==='field').onclick();elements.get('print-scope').value='current';
  for(const enabled of [false,true]) {
    elements.get('coordinates').checked=enabled;previousGrid=grids;
    await elements.get('print').onclick();assert.equal(grids,previousGrid+(enabled?1:0),'Print must respect the grid checkbox');
    await events.afterprint();
  }
  if(process.env.STREWN_PLOTLY_CHECK)await require('./PlotlyRenderTest.cjs')(printContext.plots,process.env.STREWN_PLOTLY_CHECK);
  if(process.argv[3]) {
    const native=require('@napi-rs/canvas'),render={document:{createElement:()=>native.createCanvas(1,1)},Image:native.Image,setTimeout,clearTimeout,structuredClone};
    vm.runInNewContext(fs.readFileSync(path.join(__dirname,'../output/web/graticule.js'),'utf8'),render);
    vm.runInNewContext(fs.readFileSync(path.join(__dirname,'../output/web/print.js'),'utf8'),render);
    const uri=await render.StrewnPrint.mapImage(printContext,{quantity:'count',density:true,hpd:true,points:false,opacity:.65,pointScope:'all'});
    const png=Buffer.from(uri.split(',')[1],'base64'),img=await native.loadImage(png),c=native.createCanvas(img.width,img.height),ctx=c.getContext('2d');ctx.drawImage(img,0,0);
    assert.equal(img.width,1600);assert.equal(img.height,892);
    const pixels=ctx.getImageData(0,0,1600,820).data;
    let red=0,coloured=0;for(let i=0;i<pixels.length;i+=4){if(pixels[i]>220&&pixels[i+1]<120&&pixels[i+2]<140)red++;if(Math.max(...pixels.subarray(i,i+3))-Math.min(...pixels.subarray(i,i+3))>30)coloured++;}
    assert.ok(red>100,'HPD contours must be visible');assert.ok(coloured>1600*820*.2,'Embedded imagery and density must not be blank');
    fs.writeFileSync(process.argv[3],png);console.log('Native canvas print-map composition verified: '+process.argv[3]);
  }
  console.log(`Viewer state tests passed: ${data.runs.length} runs; selection, fragment filter, log axis, wind, map data and print hooks. DOM stubs, not a visual browser test.`);
})().catch(error=>{console.error(error);process.exitCode=1;});
