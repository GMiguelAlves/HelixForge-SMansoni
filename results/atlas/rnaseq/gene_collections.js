/* Shared, reference-bound collections. No assay matrices are loaded here. */
(function(root,factory){const api=factory();if(typeof module==='object'&&module.exports)module.exports=api;else root.GeneCollections=api;})(typeof window==='object'?window:globalThis,()=>{
  'use strict';
  const esc=v=>String(v??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  function parse(text,reference){
    text=text.trim();if(!text)return [];
    if(text.startsWith('{')||text.startsWith('[')){const obj=JSON.parse(text);if(obj.reference_id&&obj.reference_id!==reference)throw Error('This file uses a different reference assembly.');const ids=Array.isArray(obj)?obj:obj.ids;if(!Array.isArray(ids)||ids.some(x=>typeof x!=='string'))throw Error('Collection JSON must contain an array of gene IDs.');return ids;}
    const lines=text.split(/\r?\n/).filter(x=>x.trim()&&!x.trim().startsWith('#'));
    if(!lines.length)return [];
    const cells=line=>line.split(/[\t,;]/).map(x=>x.trim().replace(/^"|"$/g,''));
    const header=cells(lines[0]),column=header.findIndex(x=>/^(gene_?id|resolved_gene_id|id)$/i.test(x));
    return column>=0?lines.slice(1).map(line=>cells(line)[column]||'').filter(Boolean):text.split(/[\s,;]+/).filter(Boolean);
  }
  function resolve(tokens,catalog,choices={}){
    const current=new Map(catalog.genes.map(([id])=>[id.toLowerCase(),id])),seen=new Set(),ids=[],unknown=[],ambiguous=[],mapped=[];let duplicates=0;
    for(const raw of tokens){const token=raw.trim();if(!token)continue;const key=token.toLowerCase(),canonical=current.get(key),candidates=canonical?[canonical]:(catalog.aliases[key]||[]).filter(id=>current.has(id.toLowerCase()));
      let id=candidates.length===1?candidates[0]:choices[key];if(!candidates.includes(id))id=undefined;
      if(candidates.length===0){if(!unknown.includes(token))unknown.push(token);continue;}
      if(!id){if(!ambiguous.some(a=>a.token.toLowerCase()===key))ambiguous.push({token,candidates});continue;}
      if(seen.has(id)){duplicates++;continue;}seen.add(id);ids.push(id);if(id.toLowerCase()!==key)mapped.push({token,id});
    }
    return {ids,unknown,ambiguous,mapped,duplicates};
  }
  function regionMatches(region,scope,genes){
    if(scope===null)return true;
    if(region.kind==='consensus')return String(region.genes||'').split(';').some(id=>scope.has(id));
    return genes.some(g=>scope.has(g.id)&&g.chrom===region.chrom&&g.start<region.end&&region.start<g.end);
  }
  function download(name,text,type='text/plain'){const url=URL.createObjectURL(new Blob([text],{type})),a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);}
  async function mount(host,{onChange=()=>{}}={}){
    host.innerHTML='<p role="status">Loading gene collections…</p>';
    const response=await fetch('data/gene_collections.json');if(!response.ok)throw Error('Gene collection catalog could not be loaded.');
    const catalog=await response.json(),key='hf-atlas-collections-v1',valid=new Set(catalog.genes.map(r=>r[0])),rna=new Set(catalog.genes.filter(r=>r[1]).map(r=>r[0]));
    let saved=[],state={key:'all',name:'All genes',ids:null},draft=null,choices={},storageOK=true,loadWarning='';
    function read(){try{const raw=JSON.parse(localStorage.getItem(key)||'null');if(!raw)return;if(raw.reference_id!==catalog.reference_id){loadWarning='Saved collections use another reference; import a compatible file.';return;}saved=(raw.saved||[]).filter(x=>Array.isArray(x.ids)&&x.ids.every(id=>valid.has(id)));if(raw.active&&(raw.active.ids===null||Array.isArray(raw.active.ids)&&raw.active.ids.length>0&&raw.active.ids.every(id=>valid.has(id))))state=raw.active;else loadWarning='The saved selection is incompatible with this catalog; all genes are selected.';}catch{storageOK=false;}}
    read();
    const query=new URLSearchParams(location.search),preset=query.get('collection');
    function presetState(value){if(value==='all')return {key:'all',name:'All genes',ids:null};if(value==='machinery')return {key:value,...catalog.presets.machinery};if(value.startsWith('group:')&&catalog.groups[value.slice(6)])return {key:value,name:value.slice(6),ids:catalog.groups[value.slice(6)]};return null;}
    if(preset&&presetState(preset))state=presetState(preset);
    host.classList.add('gene-collection');
    host.innerHTML=`<div class="gc-heading"><div><h2>Gene collection</h2></div><details class="gc-info"><summary aria-label="About gene collections">i</summary><div>Use the same collection in RNA-seq and ChIP-seq. Each gene remains a separate row. The machinery preset contains 118 resolved genes; the seven excluded records are omitted. Functional groups are shortcuts, not evidence tiers. Lists and names are saved only in this browser. Download a file to move a collection to another device.</div></details></div><div class="gc-bar"><label>Collection<select data-gc="select"></select></label><button type="button" data-gc="edit">Edit gene list</button><button type="button" data-gc="export">Export collection</button><a class="gc-assay" href="${location.pathname.endsWith('chipseq.html')?'atlas.html#collection':'chipseq.html#collection'}">${location.pathname.endsWith('chipseq.html')?'View in RNA-seq':'View in ChIP-seq'} →</a></div><div data-gc="summary" class="gc-summary" role="status"></div><details data-gc="editor" class="gc-editor"><summary>List editor</summary><div class="gc-editor-grid"><div><label>Gene IDs<textarea data-gc="input" rows="5" placeholder="Smp_244880&#10;Smp_179650"></textarea></label><label>Or load a list<input data-gc="file" type="file" accept=".txt,.tsv,.csv,.json"></label></div><div><label>Collection name<input data-gc="name" maxlength="80" value="My collection"></label><button type="button" data-gc="check">Check IDs</button><div data-gc="preview" role="status"></div><div class="gc-actions"><button type="button" data-gc="apply" disabled>Use matched genes</button><button type="button" data-gc="save" disabled>Save collection</button></div></div></div></details><div data-gc="message" role="status"></div>`;
    const $=name=>host.querySelector(`[data-gc="${name}"]`);
    function menus(){$('select').innerHTML='<option value="all">All genes</option><option value="machinery">Epigenetic machinery · 118</option>'+Object.keys(catalog.groups).map(name=>`<option value="group:${esc(name)}">↳ ${esc(name)} · ${catalog.groups[name].length}</option>`).join('')+saved.map((x,i)=>`<option value="saved:${i}">${esc(x.name)} · ${x.ids.length}</option>`).join('')+(state.key==='custom'?`<option value="custom">${esc(state.name)} · ${state.ids.length}</option>`:'');$('select').value=state.key;}
    function persist(){try{localStorage.setItem(key,JSON.stringify({reference_id:catalog.reference_id,active:state,saved}));}catch{storageOK=false;}}
    function notify(){menus();const n=state.ids===null?catalog.genes.length:state.ids.length,nRNA=state.ids===null?rna.size:state.ids.filter(id=>rna.has(id)).length;$('summary').textContent=`${state.name} · ${n.toLocaleString()} genes · ${nRNA.toLocaleString()} with RNA-seq${n>nRNA?` · ${n-nRNA} without RNA-seq`:''}`;$('message').textContent=loadWarning||(!storageOK?'Browser storage unavailable. Export your collection to preserve it.':'');onChange(state);}
    function change(next){state={...next,ids:next.ids===null?null:[...next.ids]};persist();const url=new URL(location.href);url.searchParams.delete('collection');history.replaceState(null,'',url);notify();}
    function check(){try{draft=resolve(parse($('input').value,catalog.reference_id),catalog,choices);$('preview').innerHTML=`<p>${draft.ids.length} matched · ${draft.mapped.length} historical IDs mapped · ${draft.duplicates} duplicates removed · ${draft.unknown.length} unmatched · ${draft.ambiguous.length} ambiguous</p>${draft.unknown.length?`<details><summary>Unmatched IDs</summary><p>${draft.unknown.map(esc).join(', ')}</p></details>`:''}${draft.ambiguous.map(a=>`<label>${esc(a.token)}<select data-choice="${esc(a.token.toLowerCase())}"><option value="">Choose current ID</option>${a.candidates.map(id=>`<option>${esc(id)}</option>`).join('')}</select></label>`).join('')}${draft.mapped.length?`<details><summary>Historical ID mapping</summary><p>${draft.mapped.map(x=>`${esc(x.token)} → ${esc(x.id)}`).join('<br>')}</p></details>`:''}`;$('apply').disabled=$('save').disabled=!draft.ids.length;$('apply').textContent=`Use ${draft.ids.length} matched genes`;}catch(e){draft=null;$('preview').textContent=e.message;$('apply').disabled=$('save').disabled=true;}}
    $('select').onchange=()=>{const value=$('select').value;if(value.startsWith('saved:')){const x=saved[Number(value.slice(6))];change({key:value,name:x.name,ids:x.ids});}else if(presetState(value))change(presetState(value));};
    $('edit').onclick=()=>{$('editor').open=true;$('input').value=(state.ids||[]).join('\n');$('name').value=state.ids?state.name:'My collection';choices={};check();$('input').focus();};
    $('input').oninput=()=>{draft=null;choices={};$('apply').disabled=$('save').disabled=true;$('preview').textContent='Check the edited IDs before applying.';};
    $('check').onclick=check;$('preview').onchange=e=>{if(e.target.dataset.choice){choices[e.target.dataset.choice]=e.target.value;check();}};
    $('file').onchange=async()=>{try{const file=$('file').files[0];if(!file)return;if(file.size>2000000)throw Error('Please use a gene-list file smaller than 2 MB.');$('input').value=await file.text();if(file.name.endsWith('.json')){try{const loaded=JSON.parse($('input').value);if(loaded.name)$('name').value=String(loaded.name).slice(0,80);}catch{}}choices={};check();}catch(e){$('preview').textContent=e.message;}};
    const custom=()=>({key:'custom',name:$('name').value.trim()||'My collection',ids:draft.ids});
    $('apply').onclick=()=>{if(draft?.ids.length){change(custom());$('editor').open=false;}};
    $('save').onclick=()=>{if(!draft?.ids.length)return;const next=custom(),index=saved.findIndex(x=>x.name===next.name);if(index<0)saved.push({name:next.name,ids:[...next.ids]});else saved[index]={name:next.name,ids:[...next.ids]};change({...next,key:`saved:${index<0?saved.length-1:index}`});$('editor').open=false;};
    $('export').onclick=()=>download('gene_collection.json',JSON.stringify({schema_version:1,reference_id:catalog.reference_id,name:state.name,ids:state.ids||catalog.genes.map(r=>r[0])},null,2),'application/json');
    window.addEventListener('storage',e=>{if(e.key===key){read();notify();}});
    persist();notify();return {catalog,get state(){return state;},contains:id=>state.ids===null||state.ids.includes(id),select:value=>{const next=presetState(value);if(next)change(next);}};
  }
  return {parse,resolve,regionMatches,download,mount,esc};
});
