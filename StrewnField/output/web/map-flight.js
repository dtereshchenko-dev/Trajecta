(function(root) {
  'use strict';
  function create({map,L,palette,load,choose}) {
    const el=id=>document.getElementById(id),controls=['map-play','map-reset','map-time','map-rate','map-fit-flight'];
    el('map-flight-status').setAttribute('aria-live','off');
    const pane=map.createPane('flight');pane.style.zIndex=460;pane.style.pointerEvents='none';
    const canvas=L.DomUtil.create('canvas','flight-canvas');canvas.setAttribute('aria-hidden','true');pane.append(canvas);
    const markers=L.layerGroup().addTo(map);
    const profile=el('map-profile');let profileView=null;
    let model=null,run=null,time=0,playing=false,frameId=null,last=null,revision=0,drawn=[],visible=false,pending=false;
    function pause() {
      playing=false;if(frameId!==null)cancelAnimationFrame(frameId);frameId=null;last=null;
      el('map-play').setAttribute('aria-pressed','false');el('map-play').title='Play';el('map-play').setAttribute('aria-label','Play');
      el('map-play').classList.remove('playing');
    }
    function status() {
      if(!model)return;
      el('map-time').value=String(time);el('map-clock').textContent=`${time.toFixed(2)} / ${model.end.toFixed(2)} s`;
      const airborne=drawn.filter(s=>!s.landed).length,landed=drawn.length-airborne;
      el('map-flight-status').textContent=`Run ${run.run_id} | ${model.tracks.length} / ${model.total} tracks | ${airborne} airborne, ${landed} landed representatives | Ground + side projections | Colour: mass (log) | Interpolated preview`;
    }
    function draw() {
      canvas.style.display=visible&&model?'block':'none';
      if(!model){profile.width=profile.width;profileView=null;}
      if(!visible||!model||!el('map').clientWidth||!el('map').clientHeight)return;
      const size=map.getSize(),ratio=Math.min(2,window.devicePixelRatio||1);
      canvas.width=Math.round(size.x*ratio);canvas.height=Math.round(size.y*ratio);canvas.style.width=size.x+'px';canvas.style.height=size.y+'px';
      L.DomUtil.setPosition(canvas,map.containerPointToLayerPoint([0,0]));const ctx=canvas.getContext('2d');ctx.scale(ratio,ratio);
      drawn=StrewnTrajectory.draw(ctx,model,time,(lat,lon)=>map.latLngToContainerPoint([lat,lon]),palette);status();
      const width=profile.clientWidth,height=profile.clientHeight;
      if(width&&height) {
        profile.width=Math.round(width*ratio);profile.height=Math.round(height*ratio);
        const side=profile.getContext('2d');side.scale(ratio,ratio);
        profileView=StrewnTrajectory.drawProfile(side,model,time,width,height,palette,drawn);
      }
    }
    function clearDetails(){el('profile-details').hidden=true;el('profile-details').replaceChildren();}
    function seek(value){pause();clearDetails();if(model){time=Math.max(model.start,Math.min(model.end,value));draw();}}
    function tick(now) {
      if(!playing)return;
      const dt=last===null?0:Math.min(.25,Math.max(0,(now-last)/1000));last=now;
      const rate=Number(el('map-rate').value);
      time=Math.min(model.end,time+dt*rate);draw();
      if(time>=model.end)pause();else frameId=requestAnimationFrame(tick);
    }
    function fit(){if(model&&visible&&el('map').clientWidth)map.fitBounds(model.bounds,{padding:[45,45],animate:false});}
    function popup(marker) {
      const r=marker.row,node=document.createElement('div'),title=document.createElement('strong'),details=document.createElement('div');
      title.textContent=marker.label;details.style.whiteSpace='pre-line';
      details.textContent=`t = ${r[0].toFixed(3)} s\nWGS84: ${r[1].toFixed(2)} km\n${r[5].toFixed(5)}, ${r[6].toFixed(5)}\n${(r[2]*1000).toPrecision(5)} g | ECEF ${r[3].toFixed(3)} km/s`;
      if(r[7]!==undefined)details.textContent+=`\nAir-relative ${r[7].toFixed(3)} km/s`;
      if(Number.isFinite(marker.pressure))details.textContent+=`\nModel pressure q = ${(marker.pressure/1e6).toFixed(3)} MPa`;
      if(marker.kind==='dark')details.textContent+='\nSpeed threshold, not an observed end of light.';
      node.append(title,details);return node;
    }
    function annotate() {
      markers.clearLayers();if(!visible||!model)return;
      for(const m of model.markers) {
        const marker=L.circleMarker([m.row[5],m.row[6]],{radius:m.kind==='height'?3:6,color:'#fff',weight:1.5,
          fillColor:m.kind==='breakup'?'#e77c42':m.kind==='dark'?'#3cc9da':'#fff2bc',fillOpacity:1}).bindPopup(popup(m));
        if(m.kind==='height')marker.bindTooltip(m.label,{permanent:true,direction:'top',className:'flight-height-label'});
        marker.on('click',()=>seek(m.row[0]));marker.addTo(markers);
      }
    }
    async function refresh(fitAfter=false) {
      const token=++revision;pause();clearDetails();model=null;drawn=[];markers.clearLayers();draw();
      controls.forEach(id=>el(id).disabled=true);
      if(!visible)return;
      el('map-flight-status').textContent='Loading saved trajectory...';
      try {
        const result=await load();if(token!==revision)return;
        run=result.run;model=StrewnTrajectory.prepare(result.flight,result.values,result.stride,120,result.threshold,run.entry);
        time=model.start;el('map-time').min=String(model.start);el('map-time').max=String(model.end);
        controls.forEach(id=>el(id).disabled=false);annotate();if(fitAfter)fit();draw();
      }catch(error){if(token===revision)el('map-flight-status').textContent=error.message;}
    }
    el('trajectory').onchange=()=>{
      visible=el('trajectory').checked;el('map-playback').hidden=!visible;el('profile-view').hidden=!visible;el('map-ground-title').hidden=!visible;
      el('map-views').classList.toggle('two-plane',visible);map.invalidateSize({animate:false});return refresh(visible);
    };
    el('map-play').onclick=()=>{
      if(playing){pause();return;}if(!model||!visible)return;if(time>=model.end)time=model.start;
      clearDetails();
      playing=true;el('map-play').setAttribute('aria-pressed','true');el('map-play').setAttribute('aria-label','Pause');el('map-play').title='Pause';el('map-play').classList.add('playing');
      frameId=requestAnimationFrame(tick);
    };
    el('map-reset').onclick=()=>seek(model?.start??0);el('map-time').oninput=()=>seek(Number(el('map-time').value));
    el('map-fit-flight').onclick=fit;el('map-run').onchange=()=>choose(Number(el('map-run').value));
    map.on('move zoomend resize',()=>{if(!pending){pending=true;requestAnimationFrame(()=>{pending=false;draw();});}});
    map.on('click',event=>{
      if(!visible||!model||!event.containerPoint)return;
      let nearest=null,distance=13;
      for(const s of drawn){const p=map.latLngToContainerPoint([s.row[5],s.row[6]]),d=Math.hypot(p.x-event.containerPoint.x,p.y-event.containerPoint.y);if(d<distance){nearest=s;distance=d;}}
      if(nearest){pause();L.popup().setLatLng([nearest.row[5],nearest.row[6]]).setContent(popup({row:nearest.row,label:`Fragment #${nearest.track.id} | x${nearest.track.multiplicity}${nearest.landed?' | Landed':''}`})).openOn(map);}
    });
    profile.onclick=event=>{
      if(!visible||!model||!profileView)return;
      let nearest=null,distance=13;
      const candidates=[...drawn.map(s=>({row:s.row,label:`Fragment #${s.track.id} | x${s.track.multiplicity}${s.landed?' | Landed':''}`})),...model.markers.filter(m=>m.kind!=='height')];
      for(const candidate of candidates) {
        const p=profileView.point(candidate.row),d=Math.hypot(p.x-event.offsetX,p.y-event.offsetY);
        if(d<distance){nearest=candidate;distance=d;}
      }
      if(nearest){if(nearest.kind)seek(nearest.row[0]);else pause();el('profile-details').replaceChildren(popup(nearest));el('profile-details').hidden=false;}
    };
    document.addEventListener('visibilitychange',()=>{if(document.hidden)pause();});
    window.addEventListener('beforeprint',pause);
    return {pause,refresh,draw};
  }
  root.StrewnMapFlight={create};
})(globalThis);
