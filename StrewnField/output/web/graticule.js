(function(root) {
  'use strict';
  function ticks(min,max,pixels,spacing) {
    if(![min,max,pixels,spacing].every(Number.isFinite)||max<=min||pixels<=0||spacing<=0)return {step:0,values:[]};
    const target=(max-min)/Math.max(1,Math.min(32,pixels/spacing)),unit=10**Math.floor(Math.log10(target));
    const step=([1,2,5,10].find(v=>v*unit>=target*(1-1e-12))??10)*unit;
    const first=Math.ceil(min/step-1e-9),last=Math.floor(max/step+1e-9),values=[];
    for(let i=first;i<=last&&values.length<64;i++)values.push(i*step);
    return {step,values};
  }
  function label(value,step,axis) {
    const decimals=Math.max(0,Math.min(9,Math.ceil(-Math.log10(step)-1e-9)));
    if(axis==='longitude')value=((value+180)%360+360)%360-180;
    const rounded=Number(value.toFixed(decimals)),absolute=Math.abs(rounded);
    const direction=rounded===0||(axis==='longitude'&&absolute===180)?'':axis==='longitude'?(rounded<0?' W':' E'):(rounded<0?' S':' N');
    return absolute.toFixed(decimals)+'\u00b0'+direction;
  }
  function layout(width,height,coordinate,point,fontSize=13) {
    const nw=coordinate(0,0),se=coordinate(width,height),center=coordinate(width/2,height/2);
    const longitude=ticks(nw.lng,se.lng,width,Math.max(120,fontSize*10));
    const latitude=ticks(Math.max(-85.05112878,se.lat),Math.min(85.05112878,nw.lat),height,fontSize*6);
    const project=(axis,source)=>source.values.map(value=>({value,
      position:axis==='longitude'?point(center.lat,value)[0]:point(value,center.lng)[1],
      label:label(value,source.step,axis)})).filter(t=>Number.isFinite(t.position)&&t.position>=0&&t.position<=(axis==='longitude'?width:height));
    return {longitudeStep:longitude.step,latitudeStep:latitude.step,
      longitudes:project('longitude',longitude),latitudes:project('latitude',latitude)};
  }
  function draw(ctx,width,height,coordinate,point,{fontSize=13,corner=90}={}) {
    const grid=layout(width,height,coordinate,point,fontSize),padding=fontSize*.35;
    ctx.save();ctx.beginPath();ctx.rect(0,0,width,height);ctx.clip();
    ctx.beginPath();
    for(const t of grid.longitudes){ctx.moveTo(t.position,0);ctx.lineTo(t.position,height);}
    for(const t of grid.latitudes){ctx.moveTo(0,t.position);ctx.lineTo(width,t.position);}
    ctx.setLineDash([4,5]);ctx.strokeStyle='rgba(255,255,255,.6)';ctx.lineWidth=2;ctx.stroke();
    ctx.strokeStyle='rgba(30,49,58,.5)';ctx.lineWidth=.7;ctx.stroke();ctx.setLineDash([]);
    ctx.font=`${fontSize}px Arial, sans-serif`;ctx.textBaseline='middle';ctx.textAlign='left';
    function text(t,x,y) {
      const w=ctx.measureText(t.label).width;
      ctx.fillStyle='rgba(255,255,255,.9)';ctx.fillRect(x-padding,y-fontSize/2-padding,w+2*padding,fontSize+2*padding);
      ctx.fillStyle='#263e49';ctx.fillText(t.label,x,y);
    }
    let right=corner;
    for(const t of grid.longitudes) {
      const w=ctx.measureText(t.label).width,x=t.position-w/2;
      if(x-padding<right||x+w+padding>width-6)continue;
      text(t,x,fontSize*.9+padding);right=x+w+padding+fontSize;
    }
    let bottom=corner;
    for(const t of [...grid.latitudes].sort((a,b)=>a.position-b.position)) {
      const top=t.position-fontSize/2-padding;
      if(top<bottom||t.position+fontSize/2+padding>height-36)continue;
      text(t,6+padding,t.position);bottom=t.position+fontSize/2+padding+fontSize;
    }
    ctx.restore();return grid;
  }
  const api={ticks,label,layout,draw};
  if(typeof module==='object'&&module.exports)module.exports=api;else root.StrewnGraticule=api;
})(globalThis);
