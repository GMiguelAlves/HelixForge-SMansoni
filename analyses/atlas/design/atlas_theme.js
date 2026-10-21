/* Apply before paint; only a display preference is stored. Scientific data stay untouched. */
(()=>{
  'use strict';
  const key='hf-atlas-theme-v1',valid=value=>['system','light','dark'].includes(value);
  const system=window.matchMedia?window.matchMedia('(prefers-color-scheme: dark)'):null;
  let mode='system';
  try{const saved=localStorage.getItem(key);if(valid(saved))mode=saved;}catch{}
  function apply(){
    const theme=mode==='system'?(system?.matches?'dark':'light'):mode;
    document.documentElement.dataset.theme=theme;
    document.documentElement.style.colorScheme=theme;
    const select=document.getElementById('atlas-theme');
    if(select)select.value=mode;
  }
  apply();
  function bind(){
    const select=document.getElementById('atlas-theme');
    if(!select)return;
    select.value=mode;
    select.addEventListener('change',()=>{
      if(!valid(select.value))return;
      mode=select.value;apply();
      try{localStorage.setItem(key,mode);}catch{}
    });
  }
  if(document.readyState==='loading')document.addEventListener('DOMContentLoaded',bind,{once:true});else bind();
  if(system?.addEventListener)system.addEventListener('change',()=>{if(mode==='system')apply();});
  else if(system?.addListener)system.addListener(()=>{if(mode==='system')apply();});
  window.addEventListener('storage',event=>{
    if(event.key!==key&&event.key!==null)return;
    mode=valid(event.newValue)?event.newValue:'system';apply();
  });
})();
