(function(){
  'use strict';
  const el=id=>document.getElementById(id),data=JSON.parse(el('mov-data').textContent),h=OrbitalReportMap;
  const finite=x=>typeof x==='number'&&Number.isFinite(x),fmt=x=>finite(x)?Number(x).toLocaleString('en-US',{maximumSignificantDigits:8}):'\u2014';
  const total=data.nodes.reduce((s,n)=>s+n.weight,0),hits=data.nodes.filter(n=>n.status===5&&n.entry);
  const probability=n=>total>0?n.weight/total:0;
  const color=n=>n.status===5?(n.entry?'#ed5551':'#07808a'):[2,4,6,7].includes(n.status)?'#bf8a28':'#a6b1b7';
  el('subtitle').textContent=`${data.object} | ${data.sampling} | ${data.target_name} | ${fmt(data.altitude_km)} km`;
  el('warnings').textContent=[...(data.cancelled?['Cancelled calculation. No final probability.']:[]),...data.warnings].join('\n');
  const p=finite(data.probability_lower)?data.probability_lower===data.probability_upper?fmt(data.probability_lower):`${fmt(data.probability_lower)} to ${fmt(data.probability_upper)}`:'Unavailable';
  [['Fitted nodes',data.nodes.filter(n=>finite(n.chi_square)).length],['Selected nodes',data.nodes.filter(n=>n.weight>0).length],['Crossing nodes',hits.length],['Weighted P',p]].forEach(([label,value])=>{
    const box=document.createElement('div');box.className='metric';const name=document.createElement('div'),v=document.createElement('div');name.textContent=label;v.className='value';v.textContent=value;box.append(name,v);el('metrics').append(box);
  });
  el('statistics').textContent=`${h.utc(data.start_utc)} to ${h.utc(data.end_utc)} UTC. RA15; ${data.workers} workers; chi <= ${data.sigma_max}; H <= ${data.H_max}. `+
    'Probability is a weighted MOV integral, not the fraction of crossing nodes. Conditional on the sampled domain and converged fits. '+
    `Outer-edge weight: ${fmt(data.outer_boundary_weight)}. A single grid does not establish convergence.`;
  const config={...data};delete config.nodes;delete config.preliminary;el('configuration').textContent=JSON.stringify(config,null,2);
  let page=0,plotted=[],map,markers;
  function draw(){
    const canvas=el('range-plot'),w=canvas.clientWidth,height=canvas.clientHeight,dpr=devicePixelRatio||1;
    canvas.width=w*dpr;canvas.height=height*dpr;const c=canvas.getContext('2d');c.scale(dpr,dpr);c.font='12px sans-serif';
    const nodes=el('preliminary').checked?data.preliminary:data.nodes,log=el('log-range').checked;
    const points=nodes.filter(n=>n.attributable[4]>0&&finite(n.attributable[5])).map(n=>({n,x:log?Math.log10(n.attributable[4]):n.attributable[4],y:n.attributable[5]}));
    plotted=[];if(!points.length){c.fillText('No grid nodes',60,60);return;}
    let x0=Infinity,x1=-Infinity,y0=Infinity,y1=-Infinity;
    for(const p of points){x0=Math.min(x0,p.x);x1=Math.max(x1,p.x);y0=Math.min(y0,p.y);y1=Math.max(y1,p.y);}
    const dx=Math.max(x1-x0,1e-12),dy=Math.max(y1-y0,1e-12);x0-=dx*.04;x1+=dx*.04;y0-=dy*.04;y1+=dy*.04;
    const px=x=>72+(x-x0)/(x1-x0)*(w-100),py=y=>height-52-(y-y0)/(y1-y0)*(height-80);
    c.strokeStyle='#ced9de';c.fillStyle='#526a75';
    for(let i=0;i<=4;i++){const x=x0+(x1-x0)*i/4,y=y0+(y1-y0)*i/4;c.beginPath();c.moveTo(px(x),28);c.lineTo(px(x),height-52);c.moveTo(72,py(y));c.lineTo(w-28,py(y));c.stroke();c.textAlign='center';c.fillText((log?10**x:x).toPrecision(3),px(x),height-32);c.textAlign='right';c.fillText(y.toPrecision(3),65,py(y)+4);}
    c.textAlign='center';c.fillText('Range, AU',w/2,height-8);c.save();c.translate(14,height/2);c.rotate(-Math.PI/2);c.fillText('Range rate, AU/day',0,0);c.restore();
    points.sort((a,b)=>a.n.weight-b.n.weight).forEach(p=>{const x=px(p.x),y=py(p.y);c.fillStyle=color(p.n);c.globalAlpha=p.n.status===1?.3:.85;c.beginPath();c.arc(x,y,p.n.weight>0?2.7:1.7,0,2*Math.PI);c.fill();plotted.push({x,y,n:p.n});});c.globalAlpha=1;
  }
  el('range-plot').onmousemove=e=>{const r=e.currentTarget.getBoundingClientRect(),x=e.clientX-r.left,y=e.clientY-r.top;let best=null,d=64;for(const p of plotted){const t=(p.x-x)**2+(p.y-y)**2;if(t<d){d=t;best=p.n;}}el('selected').textContent=best?`Node ${best.row}:${best.column} | rho ${fmt(best.attributable[4])} AU | rho-dot ${fmt(best.attributable[5])} AU/day | chi ${fmt(best.chi)} | weight ${fmt(probability(best))} | ${data.status_codes[best.status]}`:'';};
  el('log-range').onchange=draw;el('preliminary').onchange=draw;el('preliminary').disabled=!data.preliminary.length;
  el('log-range').checked=data.sampling==='grid';new ResizeObserver(draw).observe(el('range-plot'));
  function table(){
    const mode=el('filter').value,rows=data.nodes.filter(n=>mode==='all'||mode==='crossing'&&n.entry&&n.status===5||mode==='failed'&&[2,6].includes(n.status)||mode==='selected'&&n.weight>0);
    const pages=Math.max(1,Math.ceil(rows.length/100));page=Math.max(0,Math.min(page,pages-1));
    const t=document.createElement('table'),header=t.createTHead().insertRow();['Node','Range, AU','Rate, AU/day','Chi','H','Normalized weight','Status','Crossing UTC','Latitude N','Longitude E'].forEach(v=>{const th=document.createElement('th');th.textContent=v;header.append(th);});
    const body=t.createTBody();rows.slice(page*100,(page+1)*100).forEach(n=>{const row=body.insertRow();[`${n.row}:${n.column}`,fmt(n.attributable[4]),fmt(n.attributable[5]),fmt(n.chi),fmt(n.H),fmt(probability(n)),data.status_codes[n.status],n.entry?h.utc(n.entry.epoch_utc):'',n.entry?fmt(n.entry.latitude_deg):'',n.entry?fmt(n.entry.longitude_deg):''].forEach(v=>row.insertCell().textContent=v);row.title=n.error||'';});
    el('node-table').replaceChildren(t);el('count').textContent=`${rows.length} nodes`;el('page-count').textContent=`${page+1} / ${pages}`;el('previous').disabled=page===0;el('next').disabled=page+1===pages;
  }
  el('filter').onchange=()=>{page=0;table();};el('previous').onclick=()=>{page--;table();};el('next').onclick=()=>{page++;table();};
  function initMap(){
    if(map)return;
    const points=hits.map(n=>n.entry),profile=h.mapProfile(data.target_name,points),centre=h.longitudeFrame(points);
    const crs=profile.geographic?L.extend({},L.CRS.EPSG4326,{distance:(a,b)=>L.CRS.Earth.distance(a,b)*(profile.radius||6371000)/6371000}):L.CRS.EPSG3857;
    map=L.map('map',{preferCanvas:true,zoomAnimation:false,crs,minZoom:0,maxZoom:22,wheelDebounceTime:120});markers=L.layerGroup().addTo(map);
    hits.forEach(n=>{const e=n.entry,content=document.createElement('div');content.textContent=`Node ${n.row}:${n.column} | ${h.utc(e.epoch_utc)} UTC | weight ${fmt(probability(n))}`;L.circleMarker([e.latitude_deg,h.unwrap(e.longitude_deg,centre)],{radius:3.5,color:'#702925',fillColor:'#ff5b57',fillOpacity:.9,weight:.7}).bindPopup(content).addTo(markers);});
    el('fit').onclick=()=>points.length?map.fitBounds(points.map(p=>[p.latitude_deg,h.unwrap(p.longitude_deg,centre)]),{padding:[40,40],maxZoom:13,animate:false}):map.setView([0,0],2);el('fit').onclick();
    const imagery=h.imageryLayer(L,profile),overview=profile.kind==='wmts'?h.imageryLayer(L,profile,true):null;
    const note=`${data.body_fixed_frame}; ${data.body_fixed_model}. Red points are weighted MOV nodes, not equally probable clones. Coordinates work offline; optional imagery requires internet. ${data.altitude_km>0?'Altitude crossings, not ground-fall locations.':''}`;
    el('map-note').textContent=note;el('map-title').textContent=`${data.target_name} boundary at ${fmt(data.altitude_km)} km`;
    el('imagery').disabled=!imagery;if(!imagery)el('imagery').checked=false;if(imagery){imagery.on('tileerror',()=>el('map-note').textContent=note+' Some imagery is unavailable.');el('imagery').onchange=()=>{if(el('imagery').checked){overview?.addTo(map);imagery.addTo(map);}else{imagery.remove();overview?.remove();}};if(el('imagery').checked){overview?.addTo(map);imagery.addTo(map);}}
    if(profile.radius)L.control.scale({imperial:false}).addTo(map);
    const pane=map.createPane('coordinates');pane.style.zIndex=450;pane.style.pointerEvents='none';const canvas=L.DomUtil.create('canvas','coordinate-canvas');pane.append(canvas);
    function grid(){canvas.style.display=el('coordinates').checked?'block':'none';if(!el('coordinates').checked)return;const size=map.getSize(),dpr=devicePixelRatio||1;canvas.width=size.x*dpr;canvas.height=size.y*dpr;canvas.style.width=size.x+'px';canvas.style.height=size.y+'px';L.DomUtil.setPosition(canvas,map.containerPointToLayerPoint([0,0]));const c=canvas.getContext('2d');c.scale(dpr,dpr);StrewnGraticule.draw(c,size.x,size.y,(x,y)=>map.containerPointToLatLng([x,y]),(lat,lon)=>{const p=map.latLngToContainerPoint([lat,lon]);return[p.x,p.y];});}
    map.on('moveend zoomend resize',grid);el('coordinates').onchange=grid;grid();
  }
  document.querySelectorAll('[data-view]').forEach(b=>b.onclick=()=>{document.querySelectorAll('[data-view]').forEach(v=>v.setAttribute('aria-selected',String(v===b)));document.querySelectorAll('.view').forEach(v=>v.classList.toggle('active',v.id===b.dataset.view));if(b.dataset.view==='field'){initMap();map.invalidateSize(false);}if(b.dataset.view==='overview')draw();});
  el('print').onclick=()=>window.print();table();draw();
})();
