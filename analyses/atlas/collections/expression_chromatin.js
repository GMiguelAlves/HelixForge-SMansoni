/* Aligned views of frozen RNA values and native ChIP assignments, never a paired model. */
(function(root,factory){const api=factory();if(typeof module==='object'&&module.exports)module.exports=api;else root.HFExpressionChromatin=api;})(typeof window==='object'?window:this,function(){
  'use strict';
  const objects=(fields,rows)=>(rows||[]).map(r=>Object.fromEntries(fields.map((f,i)=>[f,r[i]])));
  const esc=value=>String(value??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const average=values=>values.length?values.reduce((a,b)=>a+b,0)/values.length:null;
  function create(A,api){
    const samples=objects(A.sample_fields,A.samples),contexts=objects(A.chip_mark_stage_fields,A.chip_mark_stages),loci=objects(A.chip_locus_fields,A.chip_loci),compatibility=objects(A.stage_compatibility_fields,A.stage_compatibility);
    const stages=new Map(compatibility.map(c=>[JSON.stringify([c.study,c.rna_stage,c.chip_stage]),c.match_status]));
    const byGeneContext=new Map();
    for(const row of loci){const key=JSON.stringify([row.gene_index,row.context_id]);if(!byGeneContext.has(key))byGeneContext.set(key,[]);byGeneContext.get(key).push(row);}
    const geneIndex=new Map(A.genes.map((g,i)=>[g,i]));const summary=new Map();for(const r of A.chip_summaries||[]){if(!geneIndex.has(r[0]))continue;const key=JSON.stringify([geneIndex.get(r[0]),contexts[r[1]].context_id]);if(!summary.has(key))summary.set(key,[]);summary.get(key).push(r);}
    function addLoci(records){for(const r of records){const key=JSON.stringify([r.gene_index,r.context_id]);if(!byGeneContext.has(key))byGeneContext.set(key,[]);const list=byGeneContext.get(key);if(!list.some(p=>p.peak_id===r.peak_id))list.push(r);}}
    function match(rna,chip){
      if(rna.indices.some(i=>/^female_|^male_/.test(chip.stage_or_condition)&&String(samples[i].sex).toLowerCase()!==chip.stage_or_condition.split('_')[0]))return 'NOT_COVERED';
      const statuses=rna.indices.map(i=>stages.get(JSON.stringify([rna.study,samples[i].stage,chip.stage_or_condition]))||'NOT_COVERED');
      return !statuses.length||statuses.includes('NOT_COVERED')?'NOT_COVERED':statuses.every(s=>s==='DIRECT_STAGE_CONTEXT')?'DIRECT_STAGE_CONTEXT':'CONTEXTUAL_STAGE_FAMILY';
    }
    function assigned(index,context,minSupport=1){
      // Count source regions once for this selected gene, even if an association row is duplicated.
      return [...new Map((byGeneContext.get(JSON.stringify([index,context.context_id]))||[]).filter(r=>Number(r.support)>=minSupport).map(r=>[r.peak_id,r])).values()];
    }
    function expression(index,context){
      const points=context.indices.map(i=>({sample:samples[i],tpm:A.tpm[index][i]})).filter(p=>p.tpm!==null&&p.tpm!==''&&Number.isFinite(Number(p.tpm))).map(p=>({...p,tpm:Number(p.tpm),value:Math.log2(Number(p.tpm)+1)}));
      const mean=average(points.map(p=>p.value));
      return {points,n:points.length,missing:context.indices.length-points.length,mean,meanTPM:average(points.map(p=>p.tpm))};
    }
    function build(rows,study,{context='compatible',minSupport=1,clustered=api.cluster(rows,study)}={}){
      const rna=api.contexts(study),chip=contexts.filter(c=>context==='compatible'?rna.some(r=>match(r,c)!=='NOT_COVERED'):c.context_id===context);
      const ordered=clustered.rows;
      const values=ordered.map(row=>{
        const expressionValues=rna.map(c=>expression(row.index,c)),finite=expressionValues.map(v=>v.mean).filter(v=>v!==null),center=average(finite),sd=Math.sqrt(average(finite.map(v=>(v-center)**2)))||1;
        return {row,rna:expressionValues.map(v=>({...v,z:v.mean===null?null:(v.mean-center)/sd})),chip:chip.map(c=>{const regions=assigned(row.index,c,minSupport),summed=summary.get(JSON.stringify([row.index,c.context_id]));if(summed){const selected=summed.filter(r=>r[2]>=minSupport);return {regions,count:selected.reduce((n,r)=>n+r[3],0),promoter:selected.reduce((n,r)=>n+r[4],0),assignment_states:['true','not_assessed','false'].filter(flag=>selected.some(r=>flag==='true'?r[5]>0:flag==='not_assessed'?r[6]>0:r[3]-r[5]-r[6]>0))};}return {regions,count:regions.length,promoter:regions.filter(r=>['true','1'].includes(String(r.promoter_flag).toLowerCase())).length};})};
      });
      return {rows:ordered,values,rna,chip,cluster:clustered.cluster,study,minSupport,matching:rna.map(r=>chip.map(c=>match(r,c))),linkedGenes:values.filter(v=>v.chip.some(c=>c.count>0)).length};
    }
    function chart(model){
      const {rna,chip,values}=model,left=292,cellW=36,chipW=92,gap=32,top=Math.max(170,...rna.map(c=>c.label.length*5.8*.82+72)),rowH=25,chipLeft=left+rna.length*cellW+gap,width=Math.max(740,chipLeft+chip.length*chipW+24),height=top+values.length*rowH+28;
      const title='Aligned expression and chromatin heatmap. Rows are genes clustered by RNA expression; right-hand counts are independent source-assigned ChIP regions.';
      let svg=`<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}" role="img" aria-label="${title}" class="ec-matrix"><title>${title}</title><rect width="100%" height="100%" fill="#fff"/><text x="${left}" y="25" font-size="13" font-weight="700" fill="#234451">RNA · ${esc(model.study)}</text><text x="${chipLeft}" y="25" font-size="13" font-weight="700" fill="#12665d">ChIP-seq</text>`;
      svg+=api.dendrogram(model.cluster,top,rowH,left);
      rna.forEach((c,i)=>{const x=left+i*cellW,covered=model.matching[i].some(s=>s!=='NOT_COVERED');svg+=`<text x="${x+cellW/2}" y="${top-24}" transform="rotate(-55 ${x+cellW/2} ${top-24})" font-size="10" fill="#38535f">${esc(c.label)}</text><rect x="${x}" y="${top-15}" width="${cellW}" height="12" fill="${covered?'#d7eee7':'#f1f3f5'}"/><text x="${x+cellW/2}" y="${top-5}" text-anchor="middle" font-size="9" fill="#38535f"><title>${covered?'Stage-compatible ChIP context available; independent studies':'No compatible ChIP context in this view'}</title>${covered?'ChIP':'—'}</text>`;});
      chip.forEach((c,i)=>{const x=chipLeft+i*chipW+chipW/2;svg+=`<text x="${x}" y="84" text-anchor="middle" font-size="11" font-weight="700" fill="#234451">${esc(c.stage_or_condition)}</text><text x="${x}" y="103" text-anchor="middle" font-size="11" fill="#12665d">${esc(c.mark_or_factor)}</text><text x="${x}" y="122" text-anchor="middle" font-size="9" fill="#607983">${esc(c.study_id)}</text><text x="${x}" y="151" text-anchor="middle" font-size="9" fill="#607983">≥${model.minSupport} ${model.minSupport===1?'library':'libraries'}</text>`;});
      if(!chip.length)svg+=`<text x="${chipLeft}" y="84" font-size="11" fill="#607983">No matching stage</text>`;
      values.forEach((v,i)=>{const y=top+i*rowH;svg+=`<text x="${left-94}" y="${y+17}" text-anchor="end" font-size="11" fill="#12665d" data-open-gene="${esc(v.row.gene)}" role="button" tabindex="0"><title>${esc(v.row.functional_name)}</title>${esc(v.row.gene)}</text>`;
        v.rna.forEach((r,j)=>{const label=`${v.row.gene} · ${rna[j].label} · mean log₂(TPM + 1): ${r.mean===null?'NA':r.mean.toFixed(3)} · n=${r.n}`;svg+=`<rect data-ec-rna="${i}:${j}" role="button" tabindex="0" aria-label="${esc(label)}" x="${left+j*cellW}" y="${y}" width="${cellW}" height="${rowH}" stroke="#fff" fill="${r.z===null?'#e2e8ed':api.heatColor(r.z)}"><title>${esc(label)}</title></rect>`;});
        v.chip.forEach((c,j)=>{const colors=['#f8fafb','#d7eee7','#93ccc0','#429e8e','#12665d'],label=`${v.row.gene} · ${chip[j].stage_or_condition} · ${c.count} assigned regions · ${c.promoter} native promoter regions`;svg+=`<g data-ec-chip="${i}:${j}" role="button" tabindex="0" aria-label="${esc(label)}"><title>${esc(label)}</title><rect x="${chipLeft+j*chipW}" y="${y}" width="${chipW}" height="${rowH}" stroke="#fff" fill="${colors[Math.min(c.count,4)]}"/><text x="${chipLeft+j*chipW+chipW/2}" y="${y+17}" text-anchor="middle" font-size="11" fill="${c.count>=3?'#fff':'#234451'}">${c.count}${c.promoter?' P':''}</text></g>`;});
      });
      return svg+'</svg>';
    }
    function drilldown(model,rowIndex,chipIndex,rnaIndex=null){
      const v=model.values[rowIndex];if(!v)return '';
      const chip=model.chip[chipIndex],matched=chip?model.rna.map((c,i)=>({c,i})).filter(({i})=>model.matching[i][chipIndex]!=='NOT_COVERED'):[];
      const rna=rnaIndex!==null?model.rna.map((c,i)=>({c,i})).filter(({i})=>i===rnaIndex):matched;
      const evidence=chip?v.chip[chipIndex]:null,gene=v.row.gene;
      const params=new URLSearchParams({collection:'all',gene});if(chip){params.set('study',chip.study_id);params.set('stage',chip.stage_or_condition);params.set('factor',chip.mark_or_factor);}if(model.minSupport>1)params.set('minSupport',String(model.minSupport));
      const selectedValues=rna.flatMap(({i})=>v.rna[i].points.map(p=>p.value)),max=Math.max(1,...selectedValues),left=230,plotW=320,top=32,rowH=34;
      let plot=`<svg xmlns="http://www.w3.org/2000/svg" width="580" height="${top+rna.length*rowH+28}" role="img" aria-label="Individual RNA samples and mean log2 TPM plus one"><rect width="100%" height="100%" fill="#fff"/>`;
      for(let t=0;t<=4;t++){const value=max*t/4,x=left+plotW*t/4;plot+=`<line x1="${x}" x2="${x}" y1="${top-4}" y2="${top+rna.length*rowH}" stroke="#e4ecee"/><text x="${x}" y="18" text-anchor="middle" fill="#607983" font-size="10">${value.toFixed(1)}</text>`;}
      rna.forEach(({c,i},j)=>{const r=v.rna[i],y=top+j*rowH+13;plot+=`<text x="${left-12}" y="${y+4}" text-anchor="end" fill="#38535f" font-size="10">${esc(c.label)} · n=${r.n}</text>`;r.points.forEach((p,k)=>{const x=left+p.value/max*plotW,jitter=(k-(r.n-1)/2)*Math.min(3,16/Math.max(1,r.n-1));plot+=`<circle cx="${x}" cy="${y+jitter}" r="3.5" fill="#527db5" fill-opacity=".65"><title>${esc(p.sample.sample_id)} · ${p.tpm.toFixed(3)} TPM · tissue: ${esc(p.sample.tissue)} · sex: ${esc(p.sample.sex)} · treatment: ${esc(p.sample.treatment)}</title></circle>`;});if(r.mean!==null){const x=left+r.mean/max*plotW;plot+=`<line x1="${x}" x2="${x}" y1="${y-9}" y2="${y+9}" stroke="#123a57" stroke-width="2.5"><title>Mean log₂(TPM + 1): ${r.mean.toFixed(3)}</title></line>`;}});
      plot+='</svg>';
      const matchLabel=chip&&rna.some(({i})=>model.matching[i][chipIndex]==='CONTEXTUAL_STAGE_FAMILY')?'Stage family only · ChIP age is not resolved to the RNA days':'Stage / recorded-sex context · independent studies';
      return `<div class="ec-detail-heading"><div><h3>${esc(gene)}</h3><span>${esc(v.row.functional_name)}</span></div><a href="atlas.html?gene=${encodeURIComponent(gene)}#genes">Full gene view →</a></div><div class="ec-detail-grid"><div><h4>RNA expression · ${esc(model.study)}</h4><div class="ec-sample-plot">${rna.length?plot:'<p>No RNA context matches this ChIP stage in the selected study.</p>'}</div><div class="ec-caption">log₂(TPM + 1) · dots: biological samples · line: mean</div></div><div><h4>${chip?`${esc(chip.mark_or_factor)} · ${esc(chip.stage_or_condition)}`:'ChIP-seq coverage'}</h4>${evidence?`<div class="ec-count"><strong>${evidence.count}</strong> assigned candidate region${evidence.count===1?'':'s'} <span>· ${evidence.promoter} native promoter</span></div><p class="ec-caption">${esc(chip.study_id)} · ≥${model.minSupport} calling ${model.minSupport===1?'library':'libraries'}<br>${matchLabel}<br>Gene-assignment ambiguity: ${evidence.regions.length?[...new Set(evidence.regions.map(r=>String(r.gene_assignment_ambiguous).replaceAll('_',' ')))].map(esc).join(', '):'no assigned region'}</p><a class="ec-region-link" href="chipseq.html?${esc(params)}">Inspect regions & library support →</a><details class="ec-regions"><summary>Source regions (${evidence.count})</summary>${evidence.regions.map(r=>`<div>${esc(r.peak_chrom)}:${Number(r.peak_start)+1}–${r.peak_end} · ${esc(r.peak_location)} · ${r.support} libraries · TSS ${esc(r.distance_to_tss)} bp</div>`).join('')||'<p>No source-assigned region under this support filter. This is not evidence of biological absence.</p>'}</details>`:'<p>No ChIP stage is available for this RNA context in this view.</p>'}</div></div>`;
    }
    function exportRows(model){
      const rows=[['gene_id','assay','study','context','stage','factor','n_rna_samples','mean_log2_TPM_plus_1','row_z_score','assigned_regions','native_promoter_regions','min_calling_libraries','assignment_ambiguity','matching_basis']];
      model.values.forEach(v=>{
        model.rna.forEach((c,i)=>rows.push([v.row.gene,'RNA',c.study,c.label,[...new Set(c.indices.map(j=>samples[j].stage))].join(';'),'',v.rna[i].n,v.rna[i].mean??'NA',v.rna[i].z??'NA','','','','',model.matching[i].filter(s=>s!=='NOT_COVERED').join(';')||'NOT_COVERED']));
        model.chip.forEach((c,i)=>{const r=v.chip[i];rows.push([v.row.gene,'ChIP',c.study_id,c.context_id,c.stage_or_condition,c.mark_or_factor,'','','',r.count,r.promoter,model.minSupport,(r.assignment_states||[...new Set(r.regions.map(p=>p.gene_assignment_ambiguous))]).join(';')||'NO_ASSIGNED_REGION','INDEPENDENT_STUDIES_STAGE_CONTEXT']);});
      });return rows;
    }
    return {contexts,match,assigned,expression,build,chart,drilldown,exportRows,addLoci};
  }
  return {create};
});
