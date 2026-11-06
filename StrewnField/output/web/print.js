(function(root) {
  'use strict';
  const el=id=>document.getElementById(id);
  const create=(tag,text,cls)=>{const n=document.createElement(tag);if(text)n.textContent=text;if(cls)n.className=cls;return n;};
  const canvas=(w,h)=>{const c=create('canvas');c.width=w;c.height=h;return c;};
  const mercator=lat=>(1-Math.asinh(Math.tan(Math.max(-85.05112878,Math.min(85.05112878,lat))*Math.PI/180))/Math.PI)/2;
  const latitude=y=>Math.atan(Math.sinh(Math.PI*(1-2*y)))*180/Math.PI;

  function projection(bounds,width,height,padding=30) {
    const [[south,west],[north,east]]=bounds,x=(west+east+360)/720,y=(mercator(north)+mercator(south))/2;
    const scale=Math.min((width-2*padding)/Math.max(1e-12,(east-west)/360),(height-2*padding)/Math.max(1e-12,mercator(south)-mercator(north)));
    const left=x-width/(2*scale),top=y-height/(2*scale);
    return {scale,world:(wx,wy)=>[(wx-left)*scale,(wy-top)*scale],
      point:(lat,lon)=>[((lon+180)/360-left)*scale,(mercator(lat)-top)*scale],
      coordinate:(px,py)=>({lat:latitude(top+py/scale),lng:(left+px/scale)*360-180})};
  }
  function image(url,alt='') {
    return new Promise((resolve,reject)=>{
      const img=new Image();img.alt=alt;
      const timer=setTimeout(()=>reject(new Error('Timed out preparing print image')),30000);
      img.onload=()=>{clearTimeout(timer);resolve(img);};
      img.onerror=()=>{clearTimeout(timer);reject(new Error('Cannot decode print image'));};img.src=url;
    });
  }
  function contourLines(geometry) {
    if(geometry.type==='LineString')return [geometry.coordinates];
    if(geometry.type==='MultiLineString'||geometry.type==='Polygon')return geometry.coordinates;
    if(geometry.type==='MultiPolygon')return geometry.coordinates.flat();
    return [];
  }
  async function mapImage({data,densityImage,palette,peaks,runId},settings) {
    const width=1600,height=820,footer=72,c=canvas(width,height+footer),ctx=c.getContext('2d');
    ctx.fillStyle='#e5ebed';ctx.fillRect(0,0,width,height);
    ctx.font='18px Arial';ctx.fillStyle='#334a54';
    if(!data.bounds){ctx.fillText('No surviving impacts',40,50);return c.toDataURL('image/png');}
    const p=projection(data.bounds,width,height);
    for(const tile of data.tiles) {
      const img=await image(tile.url),n=2**tile.z;
      const [x,y]=p.world(tile.x/n,tile.y/n),[right,bottom]=p.world((tile.x+1)/n,(tile.y+1)/n);
      ctx.drawImage(img,Math.floor(x),Math.floor(y),Math.ceil(right)-Math.floor(x),Math.ceil(bottom)-Math.floor(y));
    }
    if(settings.density&&data.grid) {
      const overlay=canvas(width,height),oc=overlay.getContext('2d');
      oc.putImageData(densityImage(oc,width,height,p.coordinate,settings.quantity),0,0);
      ctx.globalAlpha=settings.opacity;ctx.drawImage(overlay,0,0);ctx.globalAlpha=1;
    }
    if(settings.hpd)for(const f of data.contours.features) {
      if(f.properties.quantity!==settings.quantity)continue;
      ctx.strokeStyle=f.properties.hpd>.98?'#f75464':f.properties.hpd>.7?'#ead53f':data.tiles.length?'#f4ffff':'#476670';ctx.lineWidth=2.5;
      for(const line of contourLines(f.geometry)) {
        ctx.beginPath();line.forEach(([lon,lat],i)=>{const [x,y]=p.point(lat,lon);if(i)ctx.lineTo(x,y);else ctx.moveTo(x,y);});ctx.stroke();
      }
    }
    if(settings.points)for(const point of data.impacts) {
      if(settings.pointScope==='selected'&&point[5]!==runId)continue;
      const [x,y]=p.point(point[0],point[1]);ctx.beginPath();ctx.arc(x,y,Math.min(7,2+Math.log10(1+point[2]*1000)),0,2*Math.PI);
      ctx.fillStyle='#ca893dcc';ctx.fill();ctx.strokeStyle='#fff';ctx.lineWidth=.7;ctx.stroke();
    }
    if(settings.coordinates!==false)root.StrewnGraticule.draw(ctx,width,height,p.coordinate,p.point,{fontSize:19,corner:115});
    ctx.fillStyle='#fff';ctx.fillRect(0,height,width,footer);
    ctx.strokeStyle='#abbcc4';ctx.lineWidth=1;ctx.strokeRect(.5,.5,width-1,height-1);
    const center=p.coordinate(width/2,height/2),metresPerPixel=40075016.68557849*Math.cos(center.lat*Math.PI/180)/p.scale;
    const target=160*metresPerPixel,unit=10**Math.floor(Math.log10(target));
    const metres=([5,2,1].find(v=>v*unit<=target)??1)*unit,barWidth=metres/metresPerPixel;
    ctx.strokeStyle='#26343a';ctx.lineWidth=3;ctx.beginPath();ctx.moveTo(24,height+40);ctx.lineTo(24+barWidth,height+40);ctx.stroke();
    ctx.fillStyle='#26343a';ctx.fillText(metres>=1000?`${metres/1000} km`:`${metres} m`,24,height+27);
    ctx.fillText(data.tiles.length?'Imagery (c) Esri | Embedded tiles':'No cached imagery',24,height+65);
    const legendX=width-470;
    ctx.fillText(`${settings.quantity==='mass'?'Mass':'Number'}-weighted probability / cell`,legendX,height+23);
    for(let x=0;x<240;x++){ctx.fillStyle=`rgb(${palette(x/239).join(',')})`;ctx.fillRect(legendX+x,height+32,1,12);}
    ctx.fillStyle='#26343a';ctx.fillText('0',legendX,height+65);ctx.fillText(peaks[settings.quantity].toExponential(3),legendX+148,height+65);
    for(const [i,label,color] of [[0,'99%','#f75464'],[1,'80%','#ba9d00'],[2,'50%','#476670']]) {
      ctx.fillStyle=color;ctx.fillRect(legendX+270+i*65,height+33,20,3);ctx.fillStyle='#26343a';ctx.fillText(label,legendX+270+i*65,height+65);
    }
    return c.toDataURL('image/png');
  }
  function plotSpec(model,height=260) {
    const spec=structuredClone(model),layout=spec.layout;
    const margin=layout.polar?{l:40,r:85,t:60,b:80}:{l:68,r:22,t:42,b:layout.showlegend?90:57};
    Object.assign(layout,{width:520,height,autosize:false,margin,
      font:{family:'Arial, sans-serif',size:11,color:'#26343a'},paper_bgcolor:'#fff',plot_bgcolor:'#fff'});
    if(layout.polar)delete layout.plot_bgcolor;
    layout.title={...layout.title,x:.13,font:{size:13}};
    for(const key of ['xaxis','yaxis'])if(layout[key]) {
      layout[key].automargin=true;
      if(layout[key].type==='log')Object.assign(layout[key],{dtick:1,tickformat:'.3~g'});
    }
    if(layout.showlegend)layout.legend={orientation:'h',x:0,y:-50/(height-margin.t-margin.b),yanchor:'top',font:{size:10}};
    return spec;
  }
  function copy(id) {
    const node=el(id).cloneNode(true);node.removeAttribute('id');
    node.querySelectorAll('[id]').forEach(n=>n.removeAttribute('id'));return node;
  }
  async function prepare(context) {
    const {data,scope,plots}=context,pages=create('div');
    function page(title,meta='') {
      const node=create('section',null,'print-page');node.append(create('h2',`${data.title} | ${title}`));
      if(meta)node.append(create('p',meta,'print-meta'));pages.append(node);return node;
    }
    async function charts(target,ids,height=260) {
      const row=create('div',null,'print-charts');target.append(row);
      for(const id of ids) {
        const model=plots.get(id);if(!model)continue;
        const url=await Plotly.toImage(plotSpec(model,height),{format:'svg',width:520,height});
        const img=await image(url,model.layout.title.text);img.width=520;img.height=height;row.append(img);
      }
    }
    if(scope.includes('overview')) {
      if(data.entry_summary)page('Atmospheric entry',el('subtitle').textContent).append(copy('entry-statistics'),create('p',el('entry-note').textContent,'print-meta'));
      const stats=page('Ensemble statistics',el('subtitle').textContent);
      if(el('warnings').textContent)stats.append(create('p',el('warnings').textContent,'print-meta'));
      stats.append(copy('metrics'),create('p',el('statistics-note').textContent,'print-meta'),copy('statistics'));
      const ids=['landed','counts','largest','breakup-height','mass-spectrum','convergence','mass-count','mass-budget','impact-speed','impact-time'];
      for(let i=0;i<ids.length;i+=4)await charts(page(`Ensemble distributions ${i/4+1}`,el('subtitle').textContent),ids.slice(i,i+4));
    }
    if(scope.includes('field')) {
      const settings={quantity:el('quantity').value,density:el('density').checked,hpd:el('hpd').checked,
        points:el('points').checked,coordinates:el('coordinates').checked,opacity:Number(el('opacity').value),pointScope:el('point-scope').value};
      const target=page('Impact field',`${settings.quantity==='mass'?'Mass':'Number'}-weighted ensemble | Density: ${settings.density?'on':'off'} | HPD: ${settings.hpd?'on':'off'} | Impacts: ${settings.points?settings.pointScope==='selected'?'run '+context.runId:'all runs':'off'} | Grid: ${settings.coordinates?'WGS84':'off'} | Full ensemble extent`);
      const img=await image(await mapImage(context,settings),'Impact field with imagery, density and HPD');img.className='print-map';target.append(img);
      target.append(create('p',el('map-meta').textContent,'print-meta'));
    }
    if(scope.includes('flight')) {
      const target=page('Selected realization',el('run-meta').textContent+' | '+el('run-print').textContent);
      if(el('run-error').textContent)target.append(create('p',el('run-error').textContent,'print-meta'));
      await charts(target,['altitude','mass','speed','range']);
      const key=create('div',null,'print-mass-key'),c=canvas(256,12),ctx=c.getContext('2d');
      for(let x=0;x<256;x++){ctx.fillStyle=`rgb(${context.palette(x/255).join(',')})`;ctx.fillRect(x,0,1,12);}
      key.append(create('span','Instantaneous fragment mass: '+el('mass-scale').textContent),await image(c.toDataURL(),'Mass colour scale'),create('span','x Breakup'));
      target.append(key);
      const events=page('Realization details',el('run-meta').textContent);
      events.append(create('p',el('trace-note').textContent,'print-meta'),copy('events'),copy('parameters'));
    }
    if(scope.includes('wind'))await charts(page('Wind at selected entry',el('wind-meta').textContent),['wind-components','wind-direction'],500);
    el('print-pages').replaceChildren(...pages.children);
  }
  const api={prepare,projection,plotSpec,mapImage};
  if(typeof module==='object'&&module.exports)module.exports=api;else root.StrewnPrint=api;
})(globalThis);
