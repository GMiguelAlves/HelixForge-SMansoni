/* Track selection is independent of region filtering and never pools libraries. */
(function(root,factory){const api=factory();if(typeof module==='object'&&module.exports)module.exports=api;else root.ChipTracks=api;})(globalThis,()=>{
  'use strict';
  const key=c=>JSON.stringify([c.study,c.stage,c.factor]);
  const libraryKey=(c,library)=>JSON.stringify([c.study,c.stage,c.factor,library]);
  const matches=(c,f)=>!['study','stage','factor'].some(field=>f[field]&&c[field]!==f[field]);
  const esc=v=>String(v??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  function contexts(data){
    if(data.contexts?.length)return data.contexts;
    return data.profiles.flatMap(p=>Object.keys(p.stage_segments).map(stage=>({study:p.study_id,stage,factor:p.target,method:'replicate_support',min_libraries:2,
      libraries:p.samples.filter(s=>s.condition===stage&&s.is_control==='false').map((s,i)=>({id:s.run_accession,sample:s.sample_id,label:`Library ${i+1}`}))})));
  }
  function initial(all,params,selected,filters={}){
    const available=all.filter(c=>matches(c,filters));let requested=[];
    try{requested=JSON.parse(params.get('tracks')||'[]');}catch{}
    if(!Array.isArray(requested))requested=[];
    const chosen=available.filter(c=>requested.includes(key(c))).slice(0,4);
    if(!chosen.length){const first=available.find(c=>selected&&key(c)===key(selected))||available[0];if(first)chosen.push(first);}
    let libraries=[];try{libraries=JSON.parse(params.get('trackLibraries')||'[]');}catch{}
    if(!Array.isArray(libraries))libraries=[];
    let cores=[];try{cores=JSON.parse(params.get('trackCores')||'[]');}catch{}
    if(!Array.isArray(cores))cores=[];
    let idrs=[];try{idrs=JSON.parse(params.get('trackIDR')||'[]');}catch{}if(!Array.isArray(idrs))idrs=[];
    return {idrs:new Set(chosen.filter(c=>idrs.includes(key(c))&&c.idr_count).map(key)),cores:new Set(chosen.filter(c=>cores.includes(key(c))).map(key)),contexts:new Set(chosen.map(key)),libraries:new Set(chosen.flatMap(c=>c.libraries.filter(l=>l.coordinates_available!==false&&libraries.includes(libraryKey(c,l.id))).map(l=>libraryKey(c,l.id))).slice(0,6))};
  }
  const summaryLabel=c=>c.method==='nfcore_2_1'?'Regional candidates':`Support core · ≥${c.min_libraries||2}`;
  function mount(host,data,{params,selected,filters,onChange}){
    const all=contexts(data),state=initial(all,params,selected,filters()),byKey=new Map(all.map(c=>[key(c),c]));
    host.innerHTML='<summary>Displayed tracks <span data-tracks="count"></span></summary><div class="track-picker-body"><div class="track-picker-actions"><label>Find a context<input data-tracks="search" type="search" placeholder="Study, stage or factor"></label><button type="button" class="secondary" data-tracks="compact">Compact view</button><details class="tip"><summary aria-label="About displayed tracks">i</summary><div>Start with one study / stage / factor context. Add up to four contexts and six individual libraries for comparison. Contexts and libraries remain separate tracks. Track visibility does not change the region directory, gene associations or BED export. Library identifiers remain in tooltips and study QC. Regional candidates and atomic support cores retain their source-specific meanings.</div></details></div><div data-tracks="options"></div><div data-tracks="status" role="status"></div></div>';
    const $=name=>host.querySelector(`[data-tracks="${name}"]`);
    const visible=()=>all.filter(c=>state.contexts.has(key(c))&&matches(c,filters()));
    function sync(){
      const available=all.filter(c=>matches(c,filters()));
      for(const id of state.contexts)if(!available.some(c=>key(c)===id))state.contexts.delete(id);
      if(!state.contexts.size&&available.length)state.contexts.add(key(available[0]));
      for(const id of state.idrs)if(!state.contexts.has(id))state.idrs.delete(id);
      for(const id of state.cores)if(!state.contexts.has(id))state.cores.delete(id);
      for(const id of state.libraries)if(!available.some(c=>state.contexts.has(key(c))&&c.libraries.some(l=>libraryKey(c,l.id)===id)))state.libraries.delete(id);
      const query=$('search').value.trim().toLowerCase();
      $('count').textContent=`${visible().length} context${visible().length===1?'':'s'} · ${state.libraries.size} individual ${state.libraries.size===1?'library':'libraries'}${state.cores.size?` · ${state.cores.size} shared core`: ''}${state.idrs.size?' · IDR':''}`;
      $('options').innerHTML=available.filter(c=>[c.study,c.stage,c.factor].join(' ').toLowerCase().includes(query)).map(c=>`<div class="track-context"><label class="track-context-heading"><input type="checkbox" data-context="${esc(key(c))}" ${state.contexts.has(key(c))?'checked':''}><strong>${esc(c.stage)} · ${esc(c.factor)}</strong><small>${esc(c.study)}</small></label>${state.contexts.has(key(c))?`<div class="track-library-options">${c.method==='nfcore_2_1'?`<label><input type="checkbox" data-core="${esc(key(c))}" ${state.cores.has(key(c))?'checked':''}>Shared core · ≥2</label>`:''}${c.idr_count?`<label><input type="checkbox" data-idr="${esc(key(c))}" ${state.idrs.has(key(c))?'checked':''}>IDR ≤0.05 · ${c.idr_count.toLocaleString()} intervals</label>`:''}${c.libraries.some(l=>l.coordinates_available===false)?'<small>Individual intervals not published; library coverage is available in region details.</small>':''}${c.libraries.map(l=>`<label title="${esc(l.id)} · ${esc(l.sample)}"><input type="checkbox" data-library="${esc(libraryKey(c,l.id))}" ${state.libraries.has(libraryKey(c,l.id))?'checked':''} ${l.coordinates_available===false?'disabled':''}>${esc(l.label)}<small>${esc(l.id)}</small></label>`).join('')||'<span>No individual peak calls available.</span>'}</div>`:''}</div>`).join('')||'<div class="empty">No context matches these filters.</div>';
    }
    $('search').oninput=sync;
    $('options').onchange=e=>{
      const {context,library,core,idr}=e.target.dataset;if(!context&&!library&&!core&&!idr)return;
      const set=context?state.contexts:core?state.cores:idr?state.idrs:state.libraries,limit=context||core||idr?4:6,id=context||core||idr||library;
      if(e.target.checked&&set.size>=limit){e.target.checked=false;$('status').textContent=`Show up to ${limit} ${context?'contexts':'individual libraries'} at once. Hide one before adding another.`;return;}
      if(!e.target.checked&&context&&set.size===1){e.target.checked=true;$('status').textContent='Keep at least one context visible.';return;}
      if(e.target.checked)set.add(id);else set.delete(id);$('status').textContent='';sync();onChange();
    };
    $('compact').onclick=()=>{const first=visible()[0]||all.filter(c=>matches(c,filters()))[0];state.contexts=new Set(first?[key(first)]:[]);state.libraries.clear();state.cores.clear();state.idrs.clear();$('status').textContent='';sync();onChange();};
    function focus(region){const id=key(region),c=byKey.get(id);if(!c)return;if(!state.contexts.has(id)){state.contexts=new Set([id]);state.libraries.clear();state.cores.clear();state.idrs.clear();}if(region.library&&c.libraries.some(l=>l.id===region.library&&l.coordinates_available!==false)){const lib=libraryKey(c,region.library);if(state.libraries.size>=6&&!state.libraries.has(lib))state.libraries.delete([...state.libraries][0]);state.libraries.add(lib);}sync();}
    sync();
    return {state,sync,visible,focus,summaryLabel,coreVisible:c=>state.cores.has(key(c)),idrVisible:c=>state.idrs.has(key(c)),libraryVisible:(c,id)=>state.libraries.has(libraryKey(c,id)),writeURL:p=>{if(state.idrs.size)p.set('trackIDR',JSON.stringify([...state.idrs]));else p.delete('trackIDR');p.set('tracks',JSON.stringify([...state.contexts]));if(state.cores.size)p.set('trackCores',JSON.stringify([...state.cores]));else p.delete('trackCores');if(state.libraries.size)p.set('trackLibraries',JSON.stringify([...state.libraries]));else p.delete('trackLibraries');}};
  }
  return {key,libraryKey,matches,contexts,initial,summaryLabel,mount};
});
