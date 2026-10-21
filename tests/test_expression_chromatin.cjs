// Verify joins and plotted/exported quantities against the real frozen atlas.
const fs=require('node:fs'),path=require('node:path'),vm=require('node:vm'),assert=require('node:assert/strict');
const root=path.resolve(__dirname,'..'),bundle=process.argv[2]||path.join(root,'results/atlas/rnaseq');
const E=require('../analyses/atlas/collections/expression_chromatin.js');
const element=()=>({value:'',innerHTML:'',textContent:'',className:'',insertAdjacentHTML(){},addEventListener(){},appendChild(){},querySelectorAll(){return [];}}),elements=new Map();
const document={getElementById(id){if(!elements.has(id))elements.set(id,element());return elements.get(id);},querySelector:element,querySelectorAll(){return [];},addEventListener(){},documentElement:{dataset:{}},createElement:element};
const env=vm.createContext({window:{},document,URLSearchParams,location:{search:''},performance:require('node:perf_hooks').performance});
vm.runInContext(fs.readFileSync(path.join(bundle,'data/atlas_payload.js'),'utf8'),env);
let script=fs.readFileSync(path.join(root,'analyses/atlas/rnaseq_atlas_template.html'),'utf8').match(/<script>\s*([\s\S]*?)<\/script>/)[1];
script=script.replace('    })();','window.test={contexts:groupContexts,cluster:clusteredGroupRows,dendrogram:clusterDendrogram,heatColor};})();');
vm.runInContext(script,env);
const A=env.window.HF_ATLAS;A.chip_mark_stages=A.chip_mark_stages.filter(r=>r[A.chip_mark_stage_fields.indexOf('mark_or_factor')]==='HP1');A.chip_summaries=[];
const api=env.window.test,C=JSON.parse(fs.readFileSync(path.join(bundle,'data/gene_collections.json'))),engine=E.create(A,api);
const rows=C.presets.machinery.ids.map(gene=>({gene,index:A.genes.indexOf(gene),functional_name:''}));
const model=engine.build(rows,'PRJEB32839');
assert.equal(model.rows.length,118);assert.equal(model.rna.length,15);assert.equal(model.chip.length,2);assert.equal(model.linkedGenes,21);
assert.equal(model.values.flatMap(v=>v.chip).reduce((n,c)=>n+c.count,0),24);
assert.equal(model.matching.filter(r=>r.some(s=>s!=='NOT_COVERED')).length,8);
assert.equal(model.matching.filter(r=>r.includes('DIRECT_STAGE_CONTEXT')).length,3);
assert.equal(model.matching.filter(r=>r.includes('CONTEXTUAL_STAGE_FAMILY')).length,5);
// No signal-derived ordering or exclusion of genes without assigned regions.
assert.deepEqual(Array.from(model.rows,r=>r.gene),Array.from(api.cluster(rows,'PRJEB32839').rows,r=>r.gene));
assert.ok(model.values.some(v=>v.chip.every(c=>c.count===0)));
const samples=A.samples.map(r=>Object.fromEntries(A.sample_fields.map((f,i)=>[f,r[i]])));
const locusRows=A.chip_loci.map(r=>Object.fromEntries(A.chip_locus_fields.map((f,i)=>[f,r[i]])));
for(const v of model.values){
 for(let i=0;i<model.rna.length;i++){
  const indices=model.rna[i].indices,expected=indices.reduce((n,j)=>n+Math.log2(A.tpm[v.row.index][j]+1),0)/indices.length;
  assert.ok(Math.abs(v.rna[i].mean-expected)<1e-12);assert.equal(v.rna[i].n,indices.length);
 }
 for(let i=0;i<model.chip.length;i++)assert.equal(v.chip[i].count,new Set(locusRows.filter(r=>r.gene_index===v.row.index&&r.context_id===model.chip[i].context_id).map(r=>r.peak_id)).size);
}
const filtered=engine.build(rows,'PRJEB32839',{minSupport:2});
for(const v of filtered.values)for(const c of v.chip)assert.ok(c.regions.every(r=>Number(r.support)>=2));
assert.ok(filtered.values.flatMap(v=>v.chip).reduce((n,c)=>n+c.count,0)<24);
const promoterGene=engine.build([{gene:'Smp_244880',index:A.genes.indexOf('Smp_244880')}],'PRJEB32839');
const promoterCell=promoterGene.values[0].chip.find(c=>c.count);
assert.equal(promoterCell.regions[0].promoter_flag,'true');assert.equal(promoterCell.promoter,1);assert.match(engine.chart(promoterGene),/>1 P</);
const noMatch=engine.build(rows.slice(0,3),'PRJNA597909');assert.equal(noMatch.chip.length,0);assert.equal(noMatch.matching.every(r=>!r.length),true);assert.match(engine.chart(noMatch),/No matching stage/);
const selected=engine.build(rows,'PRJNA597909',{context:engine.contexts[0].context_id});assert.equal(selected.chip.length,1);assert.ok(selected.matching.every(r=>r[0]==='NOT_COVERED'));
// All samples in a condition must be compatible: a mixed adult/cercaria condition is uncovered.
const mixed={study:'PRJEB32839',indices:[samples.findIndex(s=>s.study==='PRJEB32839'&&s.stage==='cercaria'),samples.findIndex(s=>s.study==='PRJEB32839'&&s.stage==='egg')]};
assert.equal(engine.match(mixed,engine.contexts.find(c=>c.stage_or_condition==='cercaria')),'NOT_COVERED');
const exported=engine.exportRows(model),header=exported[0],objects=exported.slice(1).map(r=>Object.fromEntries(header.map((h,i)=>[h,r[i]])));
assert.equal(objects.length,118*17);assert.equal(objects.filter(r=>r.assay==='ChIP').reduce((n,r)=>n+r.assigned_regions,0),24);
for(const row of objects.filter(r=>r.assay==='RNA')){
 const v=model.values.find(v=>v.row.gene===row.gene_id),i=model.rna.findIndex(c=>c.label===row.context);
 assert.equal(row.mean_log2_TPM_plus_1,v.rna[i].mean);assert.equal(row.n_rna_samples,v.rna[i].n);
 assert.equal(row.assigned_regions,''); // No replicated ChIP measurement per RNA condition/sample.
}
assert.ok(objects.filter(r=>r.assay==='ChIP'&&r.assigned_regions>0).every(r=>r.assignment_ambiguity==='not_assessed'));
const svg=engine.chart(model);assert.equal((svg.match(/data-ec-rna=/g)||[]).length,118*15);assert.equal((svg.match(/data-ec-chip=/g)||[]).length,118*2);
assert.match(svg,/Rows are genes clustered by RNA expression/);
const rowIndex=model.values.findIndex(v=>v.chip.some(c=>c.count)),chipIndex=model.values[rowIndex].chip.findIndex(c=>c.count),detail=engine.drilldown(model,rowIndex,chipIndex);
assert.match(detail,/Stage family only|Stage context/);assert.match(detail,/not assessed/);assert.match(detail,/biological samples/);assert.match(detail,/collection=all/);assert.match(detail,/Source regions/);
assert.ok(!/correlation coefficient|p-value|regression/i.test(svg+detail));
console.log('PASS: native gene/stage joins; 118 genes, 21 linked, 24 regions; support filtering; no unmatched-stage pairing; exact RNA means, sample points, count exports and aligned cells.');
