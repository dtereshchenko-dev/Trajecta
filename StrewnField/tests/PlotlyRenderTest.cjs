const assert=require('node:assert/strict');
const fs=require('node:fs');
const path=require('node:path');

module.exports=async function(models,outputDirectory) {
  const {JSDOM}=require('jsdom'),native=require('@napi-rs/canvas');
  const dom=new JSDOM('<!doctype html><html><body></body></html>',{pretendToBeVisual:true,runScripts:'outside-only'}),w=dom.window;
  w.URL.createObjectURL=()=> 'blob:test';
  w.HTMLCanvasElement.prototype.getContext=function(type) {
    if(type!=='2d')return null;
    this.backing??=native.createCanvas(this.width||1,this.height||1);return this.backing.getContext('2d');
  };
  const measure=native.createCanvas(1,1).getContext('2d');
  w.SVGElement.prototype.getBBox=function() {
    const style=w.getComputedStyle(this),size=parseFloat(style.fontSize)||12;
    measure.font=`${size}px Arial`;return {x:0,y:-size*.8,width:measure.measureText(this.textContent||'').width,height:size};
  };
  w.SVGElement.prototype.getComputedTextLength=function(){return this.getBBox().width;};
  w.SVGElement.prototype.getBoundingClientRect=function() {
    const b=this.getBBox();return {...b,left:b.x,right:b.x+b.width,top:b.y,bottom:b.y+b.height};
  };
  w.eval(fs.readFileSync(path.join(__dirname,'../output/web/vendor/plotly/plotly.min.js'),'utf8'));
  assert.equal(w.Plotly.version,'3.1.0');assert.ok(w.Plotly.PlotSchema.get().traces.scatterpolar);
  const print=require('../output/web/print.js');
  fs.mkdirSync(outputDirectory,{recursive:true});
  try {
    for(const [id,height] of [['wind-components',500],['wind-direction',500],['counts',260]]) {
      const model=w.JSON.parse(JSON.stringify(print.plotSpec(models.get(id),height)));
      assert.equal(w.Plotly.validate(model.data,model.layout),undefined,id+' must use supported Plotly attributes');
      const uri=await w.Plotly.toImage(model,{format:'svg',width:520,height});
      const svg=decodeURIComponent(uri.slice(uri.indexOf(',')+1));
      assert.ok(svg.includes('width="520"'));assert.ok(!svg.includes('NaN'));
      if(id==='wind-direction'){assert.ok(svg.includes('polarlayer'));assert.ok(svg.includes('Wind speed'));}
      const img=await native.loadImage(Buffer.from(svg)),c=native.createCanvas(520,height),ctx=c.getContext('2d');ctx.drawImage(img,0,0);
      fs.writeFileSync(path.join(outputDirectory,id+'.png'),c.toBuffer('image/png'));
    }
  } finally {w.close();}
  console.log('Real Plotly schema and SVG-export checks passed under jsdom with approximate text metrics, not a browser layout test.');
};
