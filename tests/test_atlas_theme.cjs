// Theme resolution must survive navigation, system changes and unavailable storage.
const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict');
const source=fs.readFileSync(require('node:path').join(__dirname,'../analyses/atlas/design/atlas_theme.js'),'utf8');
function browser({saved=null,dark=false,blocked=false}={}){
 const handlers={},events={},changes=[],storage=new Map(saved===null?[]:[['hf-atlas-theme-v1',saved]]);
 const select={value:'',addEventListener:(name,fn)=>handlers[name]=fn};let loaded=false;
 const media={matches:dark,addEventListener:(name,fn)=>changes.push(fn)};
 const document={readyState:'loading',documentElement:{dataset:{},style:{}},getElementById:()=>loaded?select:null,addEventListener:(name,fn)=>events[name]=fn};
 const window={matchMedia:()=>media,addEventListener:(name,fn)=>events[name]=fn};
 const localStorage={getItem:k=>{if(blocked)throw Error('denied');return storage.get(k)||null;},setItem:(k,v)=>{if(blocked)throw Error('denied');storage.set(k,v);}};
 vm.runInNewContext(source,{window,document,localStorage});
 return {document,storage,select,ready:()=>{loaded=true;events.DOMContentLoaded();},change:value=>{select.value=value;handlers.change();},system:dark=>{media.matches=dark;changes.forEach(fn=>fn());},sync:value=>events.storage({key:'hf-atlas-theme-v1',newValue:value})};
}
const a=browser({dark:true});assert.equal(a.document.documentElement.dataset.theme,'dark');a.ready();assert.equal(a.select.value,'system');
a.system(false);assert.equal(a.document.documentElement.dataset.theme,'light');
a.change('dark');assert.equal(a.storage.get('hf-atlas-theme-v1'),'dark');a.system(false);assert.equal(a.document.documentElement.dataset.theme,'dark');
const next=browser({saved:a.storage.get('hf-atlas-theme-v1'),dark:false});assert.equal(next.document.documentElement.dataset.theme,'dark');next.ready();assert.equal(next.select.value,'dark');
a.sync('light');assert.equal(a.select.value,'light');assert.equal(a.document.documentElement.dataset.theme,'light');
a.change('system');a.system(true);assert.equal(a.document.documentElement.dataset.theme,'dark');a.system(false);assert.equal(a.document.documentElement.dataset.theme,'light');
a.change('invalid');assert.equal(a.storage.get('hf-atlas-theme-v1'),'system');
const denied=browser({dark:true,blocked:true});denied.ready();denied.change('light');assert.equal(denied.document.documentElement.dataset.theme,'light');assert.equal(denied.document.documentElement.style.colorScheme,'light');
const corrupt=browser({saved:'invalid',dark:true});assert.equal(corrupt.document.documentElement.dataset.theme,'dark');corrupt.ready();assert.equal(corrupt.select.value,'system');
console.log('PASS: first-paint system theme, explicit override, navigation persistence, cross-tab sync, system changes and blocked storage.');
