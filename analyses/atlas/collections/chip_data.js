/* Lossless lazy tables shared by gene views and the independent ChIP browser. */
(function(root,factory){const api=factory();if(typeof module==='object'&&module.exports)module.exports=api;else root.HFChipData=api;})(typeof window==='object'?window:this,function(){
  const cache=new Map();
  const unpack=data=>data.rows.map(row=>{const r={...data.defaults,...Object.fromEntries(data.fields.map((f,i)=>[f,row[i]]))};if(data.peak_id_prefixes)r.peak_id=data.peak_id_prefixes[r.peak_id[0]]+'Interval_'+r.peak_id[1];return r;});
  const sexMatches=(sample,stage)=>!/^female_|^male_/.test(stage)||String(sample.sex||'').toLowerCase()===stage.split('_')[0];
  async function get(url){if(!cache.has(url))cache.set(url,fetch(url).then(r=>{if(!r.ok)throw Error('ChIP details unavailable; reload or consult the native downloads.');return r.json();}).catch(e=>{cache.delete(url);throw e;}));return cache.get(url);}
  async function loadGene(A,gene){
    if(!A.chip_gene_store)return [];
    const bucket=Number(gene.split('_')[1])%A.chip_gene_store.buckets;
    const data=await get(A.chip_gene_store.pattern.replace('{bucket}',bucket));
    return unpack(data).filter(r=>r.gene_id===gene).map(r=>({...r,gene_index:A.genes.indexOf(gene)}));
  }
  function evidence(A){
    const groups=new Map(),contexts=A.chip_mark_stages.map(r=>Object.fromEntries(A.chip_mark_stage_fields.map((f,i)=>[f,r[i]]))),indices=new Map(A.genes.map((g,i)=>[g,i]));
    for(const row of A.chip_summaries||[]){if(!indices.has(row[0]))continue;const [gene,ci,support,count,promoter,ambiguous,na,mask]=row,c=contexts[ci],key=gene+'|'+ci;
      if(!groups.has(key))groups.set(key,{...c,gene_index:indices.get(gene),n_peaks:0,n_promoter_peaks:0,n_single_library:0,n_multiple_libraries:0,n_ambiguous:0,n_ambiguity_not_assessed:0,mask:0,regulatory_class:'candidate_nearest_promoter_association',expected_effect:'not_inferred_from_occupancy'});
      const r=groups.get(key);r.n_peaks+=count;r.n_promoter_peaks+=promoter;r.n_single_library+=support===1?count:0;r.n_multiple_libraries+=support>=2?count:0;r.n_ambiguous+=ambiguous;r.n_ambiguity_not_assessed+=na;r.mask|=mask;
    }
    return [...groups.values()].map(r=>({...r,peak_locations:(A.chip_category_values||[]).filter((_,i)=>r.mask&(1<<i)).join(';')}));
  }
  return {unpack,get,loadGene,evidence,sexMatches};
});
