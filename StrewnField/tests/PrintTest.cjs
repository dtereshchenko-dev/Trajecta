const assert=require('node:assert/strict');
const P=require('../output/web/print.js');
const G=require('../output/web/graticule.js');
let grids=0;
global.StrewnGraticule={draw(...args){grids++;return G.draw(...args);}};
for(const bounds of [[[49.78,.72],[49.86,.9]],[[-21.4,22.5],[-21,24]],[[84,-20],[85,20]],[[0,179],[1,181]]]) {
  for(const [width,height] of [[1600,820],[400,600]]) {
    const p=P.projection(bounds,width,height);
    for(const lat of [bounds[0][0],bounds[1][0]])for(const lon of [bounds[0][1],bounds[1][1]]) {
      const [x,y]=p.point(lat,lon),back=p.coordinate(x,y);
      assert.ok(x>=29.99&&x<=width-29.99);assert.ok(y>=29.99&&y<=height-29.99);
      assert.ok(Math.abs(back.lat-lat)<1e-9);assert.ok(Math.abs(back.lng-lon)<1e-9);
    }
  }
}
const model={data:[{x:[1,1000],y:[100,.001]}],layout:{title:{text:'Mass'},xaxis:{type:'log'},yaxis:{type:'log'},showlegend:true}};
const before=JSON.stringify(model),result=P.plotSpec(model);
assert.equal(JSON.stringify(model),before,'Export must not mutate live data or axes');
assert.equal(result.layout.width,520);assert.equal(result.layout.height,260);assert.equal(result.layout.autosize,false);
assert.equal(result.layout.xaxis.dtick,1);assert.equal(result.layout.yaxis.dtick,1);
assert.ok(result.layout.height>result.layout.margin.t+result.layout.margin.b+100);
for(const height of [260,390,500]) {
  const l=P.plotSpec(model,height).layout,plotHeight=height-l.margin.t-l.margin.b;
  const legendTop=height-l.margin.b-l.legend.y*plotHeight;
  assert.ok(legendTop+20<height,'Legend must fit even on taller wind charts');
}
const rowHeightMm=260/520*((273-2)/2);
assert.ok(2*rowHeightMm+3+35<186,'Two chart rows plus headers must fit A4 landscape');
assert.ok(273*892/1600+30<186,'Map with legend and metadata must fit A4 landscape');
const polar=P.plotSpec({data:[{type:'scatterpolar',r:[0,50],theta:[359,1]}],layout:{title:{text:'Wind'},showlegend:true,polar:{angularaxis:{rotation:90,direction:'clockwise'}}}},500);
assert.equal(polar.layout.polar.angularaxis.direction,'clockwise');assert.ok(polar.layout.margin.r>=80);
assert.ok(500/520*((273-2)/2)+30<186,'Tall wind plots must fit the print page');
console.log('Print projection, fixed page dimensions and immutable chart export tests passed.');
const calls={tiles:0,points:0,density:0,paths:0};
global.document={createElement:()=>({getContext:()=>({fillRect(){},fillText(){},strokeRect(){},beginPath(){},moveTo(){},
  lineTo(){calls.paths++;},stroke(){},arc(){calls.points++;},fill(){},putImageData(){},drawImage(){calls.tiles++;},
  save(){},restore(){},rect(){},clip(){},setLineDash(){},measureText:s=>({width:s.length*10})}),toDataURL:()=> 'data:image/png;base64,test'})};
global.Image=class {set src(value){queueMicrotask(()=>value==='bad'?this.onerror():this.onload());}};
const context={data:{bounds:[[49,.5],[50,1]],grid:{},tiles:[{x:1,y:1,z:2,url:'embedded'}],
  contours:{features:[{properties:{quantity:'mass',hpd:.99},geometry:{type:'LineString',coordinates:[[.6,49.1],[.9,49.9]]}}]},
  impacts:[[49.2,.7,.1,1,1,0],[49.3,.8,.2,1,1,1]]},peaks:{mass:1},palette:()=>[0,128,100],runId:1,
  densityImage:(ctx,w,h,coordinate,quantity)=>{assert.equal(quantity,'mass');assert.equal(w,1600);assert.equal(h,820);calls.density++;return {};}};
(async()=>{
  const settings={quantity:'mass',density:true,hpd:true,points:true,coordinates:false,opacity:.5,pointScope:'selected'};
  await P.mapImage(context,settings);
  assert.deepEqual(calls,{tiles:2,points:1,density:1,paths:2},'Tiles, density, HPD and selected-run impacts must be composed');
  calls.points=0;calls.density=0;
  await P.mapImage(context,{...settings,pointScope:'all',density:false,hpd:false});
  assert.equal(calls.points,2);assert.equal(calls.density,0);
  assert.equal(grids,0);
  await P.mapImage(context,{...settings,coordinates:true,density:false,hpd:false,points:false});
  assert.equal(grids,1,'Print grid must work independently of the other map layers');
  context.data.tiles[0].url='bad';await assert.rejects(P.mapImage(context,settings),/Cannot decode print image/);
  context.data.bounds=null;assert.ok((await P.mapImage(context,settings)).startsWith('data:image/png'));
  console.log('Print layer toggles, selected-run markers, missing-data and image-failure tests passed.');
})().catch(error=>{console.error(error);process.exitCode=1;});
