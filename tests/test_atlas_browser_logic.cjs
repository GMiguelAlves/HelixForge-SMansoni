// Run with node tests/test_atlas_browser_logic.cjs [atlas bundle directory].
// Executes the real template against the bundle, with a minimal DOM for startup.
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const assert = require('node:assert/strict');
const root = path.resolve(__dirname, '..');
const bundle = process.argv[2] ? path.resolve(process.argv[2]) : path.join(root, 'results/atlas/rnaseq');
const element = () => ({value:'',innerHTML:'',textContent:'',className:'',
  insertAdjacentHTML(){},addEventListener(){},appendChild(){},querySelectorAll(){return [];}});
const elements = new Map();
const document = {getElementById(id){if(!elements.has(id))elements.set(id,element());return elements.get(id);},
  querySelector(){return element();},querySelectorAll(){return [];},addEventListener(){},documentElement:{dataset:{}},createElement:element};
const browserGlobals={URLSearchParams,location:{search:''},performance:require('node:perf_hooks').performance};
const context = vm.createContext({window:{}, document,...browserGlobals});
vm.runInContext(fs.readFileSync(path.join(bundle,'data/atlas_payload.js'),'utf8'), context);
const template = fs.readFileSync(path.join(root,'analyses/atlas/rnaseq_atlas_template.html'),'utf8');
let script = template.match(/<script>\s*([\s\S]*?)<\/script>/)[1];
script = script.replace('    })();', 'window.test={clusterDendrogram,heatColor,contextMean,groupContrasts,setCollectionScope,findGeneMatches,effectState,contrastMatchesStage,contextsForStudy,openGene,contextCatalog,navigationStage,contextMatches,expressionStats,geneConditionStats,studyConditions,studyQuality,renderStudyPage,hierarchicalCluster,clusteredGroupRows,groupExpressionHeatmap,groupEffectHeatmap,groupChipHeatmap,geneExpressionPattern,groupContexts,studyWarningsHtml,geneDistributionSvg,geneLocusSvg,geneLocusWindows,peakGeneContext};})();');
vm.runInContext(script,context);
const A=context.window.HF_ATLAS, t=context.window.test;
const modelObjects=index=>A.gene_models[index].map(row=>Object.fromEntries(A.gene_model_fields.map((field,i)=>[field,row[i]])));
const targetIndex=A.genes.indexOf('Smp_337860');
const targetAnnotation=Object.fromEntries(A.gene_annotation_fields.map((field,i)=>[field,A.gene_annotations[targetIndex][i]]));
const targetModels=modelObjects(targetIndex);
assert.equal(targetModels.length,1);
assert.equal(targetModels[0].transcript_id,'Smp_337860.1');
const targetPeaks=A.chip_loci.filter(row=>row[0]===targetIndex).map(row=>Object.fromEntries(A.chip_locus_fields.map((field,i)=>[field,row[i]])));
const targetSvg=t.geneLocusSvg(targetAnnotation,targetPeaks,targetModels);
assert.equal((targetSvg.match(/data-feature="exon"/g)||[]).length,9);
assert.equal((targetSvg.match(/data-feature="intron"/g)||[]).length,8);
assert.match(targetSvg,/data-feature="intron" data-start="6630829" data-end="6632734"/);
assert.match(targetSvg,/Smp_337860\.1 →/);
assert.equal((targetSvg.match(/data-feature="peak"/g)||[]).length,targetPeaks.length);
if(process.argv[3])fs.writeFileSync(path.resolve(process.argv[3]),targetSvg);
const syntheticAnnotation={chromosome:'chr1',start:101,end:200,strand:'+'};
const syntheticModel={transcript_id:'T1',chromosome:'chr1',start:101,end:200,strand:'+',exons:[[101,120],[181,200]]};
const nativeContext=Object.fromEntries(A.chip_mark_stage_fields.map((field,i)=>[field,A.chip_mark_stages[0][i]]));
const syntheticPeak={...nativeContext,peak_chrom:'chr1',peak_start:120,peak_end:180,peak_id:'synthetic-geometry-only',peak_location:'intron'};
const locus=t.geneLocusSvg(syntheticAnnotation,[syntheticPeak],[syntheticModel]);
assert.match(locus,/data-feature="intron" data-start="121" data-end="180"/);
assert.match(locus,/data-feature="peak" data-start="120" data-end="180"/);
const intronX=locus.match(/data-feature="intron"[^>]*x1="([^"]+)" x2="([^"]+)"/);
const peakX=locus.match(/data-feature="peak"[^>]*x="([^"]+)"[^>]*width="([^"]+)"/);
assert.equal(Number(intronX[1]),Number(peakX[1]));
assert.ok(Math.abs(Number(intronX[2])-Number(intronX[1])-Number(peakX[2]))<1e-9);
const reverse=t.geneLocusSvg({...syntheticAnnotation,strand:'-'},[],[{...syntheticModel,strand:'-'}]);
assert.match(reverse,/T1 ←/);
assert.match(reverse,/exon 2 · 101/);
const isoforms=t.geneLocusSvg(syntheticAnnotation,[],[syntheticModel,{...syntheticModel,transcript_id:'T2',exons:[[101,140],[181,200]]}]);
assert.equal((isoforms.match(/data-transcript=/g)||[]).length,2);
const unavailable=t.geneLocusSvg(syntheticAnnotation,[],[]);
assert.doesNotMatch(unavailable,/data-feature="exon"/);
assert.match(unavailable,/No exon structure annotated/);
console.log('PASS: exon boxes, intron lines, separate isoforms, both strands and GFF/BED boundary alignment.');
assert.equal(A.locus_display.reclassify_source_peaks,false);
assert.equal(A.locus_display.outside_flank_bp,2000);
assert.equal(A.locus_display.inside_gene_bp,500);
const reverseIndex=A.genes.indexOf('Smp_210650');
const reverseAnnotation=Object.fromEntries(A.gene_annotation_fields.map((field,i)=>[field,A.gene_annotations[reverseIndex][i]]));
// Geometry fixtures only: they do not assert biological occupancy in this gene.
const reversePeaks=[{...nativeContext,peak_chrom:reverseAnnotation.chromosome,peak_start:35689210,peak_end:35689729,peak_id:'synthetic-upstream',peak_location:'promoter'},
  {...nativeContext,peak_chrom:reverseAnnotation.chromosome,peak_start:35661796,peak_end:35661959,peak_id:'synthetic-downstream',peak_location:'downstream'}];
const upstreamPeak=reversePeaks.find(row=>Number(row.peak_start)===35689210&&Number(row.peak_end)===35689729);
const originalPeak=JSON.stringify(upstreamPeak),upstreamContext=t.peakGeneContext(reverseAnnotation,upstreamPeak);
assert.equal(upstreamContext.position,'upstream');
assert.equal(upstreamContext.fivePrimeDistance.join(','),'-2423,-1905');
assert.equal(upstreamContext.promoterOverlap,96);
assert.equal(upstreamPeak.peak_location,'promoter');
assert.equal(JSON.stringify(upstreamPeak),originalPeak);
const downstreamContext=t.peakGeneContext(reverseAnnotation,reversePeaks.find(row=>Number(row.peak_start)===35661796));
assert.equal(downstreamContext.position,'downstream');
assert.equal(downstreamContext.threePrimeDistance.join(','),'63,225');
const forwardGene={start:5001,end:15000,strand:'+'};
const forwardUpstream=t.peakGeneContext(forwardGene,{peak_start:3800,peak_end:4000});
assert.equal(forwardUpstream.position,'upstream');
assert.equal(forwardUpstream.fivePrimeDistance.join(','),'-1200,-1001');
const forwardDownstream=t.peakGeneContext(forwardGene,{peak_start:15500,peak_end:15600});
assert.equal(forwardDownstream.position,'downstream');
assert.equal(forwardDownstream.threePrimeDistance.join(','),'501,600');
assert.equal(forwardDownstream.downstreamOverlap,100);
assert.equal(t.peakGeneContext(forwardGene,{peak_start:4999,peak_end:5000}).position,'upstream');
assert.equal(t.peakGeneContext(forwardGene,{peak_start:5000,peak_end:5001}).position,'gene body');
assert.equal(t.peakGeneContext(forwardGene,{peak_start:4999,peak_end:5002}).position,'overlaps gene 5′ boundary');
const reverseLocus=t.geneLocusSvg(reverseAnnotation,reversePeaks,modelObjects(reverseIndex));
assert.match(reverseLocus,/data-feature="promoter-reference-window" data-start="35686806" data-end="35689306"/);
assert.ok(reverseLocus.indexOf('Downstream · 3′ flank')<reverseLocus.indexOf('Upstream · 5′ flank'));
assert.match(reverseLocus,/WBPS19 position: upstream/);
assert.match(reverseLocus,/Reference-window overlap: 5′ 96 bp/);
if(process.argv[4])fs.writeFileSync(path.resolve(process.argv[4]),reverseLocus);
console.log('PASS: strand-aware upstream/downstream positions, signed boundary distances and reference-window overlaps; source labels preserved.');
assert.equal(Object.keys(A.selected_groups).join(','),'Writers,Readers,Erasers,HP1,Chaperones,Cofactors,Remodelers');
const synthetic=t.hierarchicalCluster([[0,0],[10,10],[0.1,0],[10.1,10]]);
assert.equal(synthetic.order.join(','),'0,2,1,3');
assert.equal(t.hierarchicalCluster([]).order.length,0);
assert.equal(t.hierarchicalCluster([[3,4]]).order.join(','),'0');
const selected=A.selected_gene_status.map(row=>Object.fromEntries(A.selected_gene_fields.map((field,i)=>[field,row[i]])));
const readers=selected.filter(row=>String(row.groups).split(';').includes('Readers'));
const uniqueReaders=[...new Map(readers.map(row=>[row.resolved_gene_id,row])).values()];
const readerRows=uniqueReaders.map(row=>({gene:row.resolved_gene_id,functional_name:row.functional_name,index:A.genes.indexOf(row.resolved_gene_id)}));
assert.ok(readerRows.every(row=>row.index>=0));
const grouped=t.clusteredGroupRows(readerRows);
assert.equal(grouped.rows.length,uniqueReaders.length);
assert.equal(new Set(grouped.rows.map(row=>row.gene)).size,grouped.rows.length);
const heatmapOrder=svg=>[...svg.matchAll(/<text(?: data-open-gene="[^"]+" role="button" tabindex="0" aria-label="[^"]+")? x="260" y="[^"]+" text-anchor="end" font-size="11">(Smp_[0-9]+)/g)].map(match=>match[1]).join(',');
const expectedOrder=grouped.rows.map(row=>row.gene).join(',');
for(const svg of [t.groupExpressionHeatmap(grouped.rows,grouped.cluster),t.groupEffectHeatmap(grouped.rows,grouped.cluster)]){
  assert.equal(heatmapOrder(svg),expectedOrder);
  assert.match(svg,/Average-linkage gene dendrogram/);
}
assert.match(t.groupChipHeatmap(grouped.rows,grouped.cluster),/<svg/);
assert.match(t.groupChipHeatmap(grouped.rows,grouped.cluster),/PRJNA1090249/);
assert.equal(document.getElementById('group-study').value,'PRJEB32839');
assert.equal(t.groupContexts('PRJEB32839').length,15);
const studyGrouped=t.clusteredGroupRows(readerRows,'PRJEB32839');
const studyExpression=t.groupExpressionHeatmap(studyGrouped.rows,studyGrouped.cluster,'PRJEB32839');
const studyEffects=t.groupEffectHeatmap(studyGrouped.rows,studyGrouped.cluster,'PRJEB32839');
assert.match(studyExpression,/<svg width="1040"/);
assert.equal(heatmapOrder(studyExpression),studyGrouped.rows.map(row=>row.gene).join(','));
assert.match(studyExpression,/PRJEB32839 · egg_mixed/);
assert.doesNotMatch(studyExpression,/PRJNA602528/);
assert.match(studyEffects,/PRJEB32839/);
assert.doesNotMatch(studyEffects,/PRJNA597909/);
assert.match(t.groupEffectHeatmap(studyGrouped.rows,studyGrouped.cluster,'PRJNA602528'),/no inferential contrasts/);
const withinStudy=t.geneExpressionPattern(readerRows[0].index,'PRJEB32839');
assert.equal(withinStudy.length,15);
assert.ok(Math.abs(withinStudy.reduce((sum,value)=>sum+value,0))<1e-9);
for(const [group,genes] of Object.entries(A.selected_groups)){
  const rows=genes.map(gene=>({gene,functional_name:'',index:A.genes.indexOf(gene)}));
  assert.ok(rows.every(row=>row.index>=0),group);
  for(const study of A.study_profiles.map(profile=>profile.study)){
    const scoped=t.clusteredGroupRows(rows,study);
    const svg=t.groupExpressionHeatmap(scoped.rows,scoped.cluster,study);
    assert.equal((svg.match(/ · z=/g)||[]).length,rows.length*t.groupContexts(study).length,`${group} × ${study}`);
    assert.equal(heatmapOrder(svg),scoped.rows.map(row=>row.gene).join(','),`${group} × ${study}`);
  }
}
const activeIds=new Set(A.selected_gene_status.map(row=>row[A.selected_gene_fields.indexOf('gene_id')]));
for(const excluded of ['Smp_024380','Smp_059920','Smp_064250','Smp_128520','Smp_128530','Smp_204920','Smp_329500']){
  assert.equal(activeIds.has(excluded),false,excluded);
}
const stats=t.expressionStats([4,1,3,2]);
assert.equal(stats.median,2.5);
assert.equal(stats.q1,1.75);
assert.equal(stats.q3,3.25);
assert.ok(Math.abs(stats.sd-Math.sqrt(5/3))<1e-12);
assert.equal(t.expressionStats([2]).sd,null);
assert.equal(t.expressionStats([]).median,null);
assert.equal(t.expressionStats([0,0,0]).sd,0);
const geneValues=A.tpm[A.genes.indexOf('Smp_179650')];
assert.equal(t.geneConditionStats(geneValues,'PRJNA602528').length,10);
assert.ok(t.geneConditionStats(geneValues,'PRJNA602528').every(r=>r.n===1&&r.sd===null));
assert.equal(t.geneConditionStats(geneValues).reduce((n,r)=>n+r.n,0),268);
assert.equal(t.findGeneMatches('').length,0);
for(const [historical,current] of Object.entries({
  Smp_147870:'Smp_213920',Smp_137240:'Smp_316180',
  Smp_180990:'Smp_241580',Smp_161010:'Smp_246410',Smp_163470:'Smp_246920',
})){
  const matches=t.findGeneMatches(historical);
  assert.equal(matches.length,1,historical);
  assert.equal(matches[0],A.genes.indexOf(current),historical);
}
assert.equal(t.findGeneMatches('Smp_179650').length,1);
assert.equal(t.findGeneMatches('chromobox').length,3);
document.getElementById('gene-search').value='';t.openGene();
assert.match(document.getElementById('gene-result').textContent,/Enter a gene/);
assert.equal(t.effectState(2,null),'PADJ_UNAVAILABLE');
assert.equal(t.effectState(null,null),'NOT_TESTED');
assert.equal(t.effectState(2,.05),'NOT_SIGNIFICANT');
const samples=A.samples.map(r=>Object.fromEntries(A.sample_fields.map((f,i)=>[f,r[i]])));
assert.equal(t.navigationStage('juvenile_26d'),'juvenile');
assert.equal(t.navigationStage('adult_d21'),'adult');
assert.equal(samples.filter(r=>t.contextMatches(r,{family:'juvenile'})).length,29);
assert.equal(samples.filter(r=>t.contextMatches(r,{family:'cercaria',tissue:'tail'})).length,1);
assert.equal(samples.filter(r=>t.contextMatches(r,{family:'adult',tissue:'tail'})).length,0);
assert.equal(samples.filter(r=>t.contextMatches(r,{treatment:'control'})).length,18);
assert.equal(samples.filter(r=>r.study==='PRJNA428470'&&t.contextMatches(r,{treatment:'control'})).length,3);
assert.equal(samples.filter(r=>t.contextMatches(r,{study:'PRJEB32839',sex:'female',stage:'juvenile_26d'})).length,5);
for(const c of t.contextCatalog){
  assert.equal(new Set(c.indices.map(i=>samples[i].condition)).size,1,'Mixed conditions in '+c.label);
}
assert.equal(t.contextsForStudy('PRJEB32839').length,15);
assert.equal(t.contextsForStudy('E-MTAB-451').length,5);
assert.equal(t.contrastMatchesStage({study:'PRJEB32839',numerator:'sporocyst_1d_mixed',denominator:'miracidium_mixed'},'adult'),false);
assert.equal(t.contrastMatchesStage({study:'PRJNA597909',numerator:'adult_pzq',denominator:'adult_control'},'adult'),false); // no adult ChIP in the native package
assert.equal(t.contrastMatchesStage({study:'PRJEB32839',numerator:'sporocyst_32d_female',denominator:'sporocyst_32d_male'},'sporocyst'),true);
assert.equal(t.contrastMatchesStage({study:'PRJEB32839',numerator:'sporocyst_1d_mixed',denominator:'miracidium_mixed'},'sporocyst'),false);
assert.equal(t.contrastMatchesStage({study:'PRJNA597909',numerator:'adult_control',denominator:'juvenile_control'},'adult'),false);
// Compare all displayed significance counts with the unrounded catalog.
const lines=fs.readFileSync(path.join(bundle,'data/contrast_catalog.tsv'),'utf8').trim().split(/\r?\n/);
const fields=lines.shift().split('\t');
for(let j=0;j<lines.length;j++){
  const r=Object.fromEntries(lines[j].split('\t').map((v,i)=>[fields[i],v]));
  const n=A.genes.reduce((n,g,i)=>n+(t.effectState(A.lfc[i][j],A.padj[i][j])==='SIGNIFICANT'),0);
  assert.equal(n,Number(r.significant_genes),r.contrast_key);
}
console.log('PASS: clustered group heatmaps, unique gene rows, search, condition grouping, stage matching and all significance counts.');

assert.equal(A.study_profiles.length,12);
assert.equal(A.study_profiles.reduce((n,p)=>n+p.samples.length,0),268);
assert.equal(A.study_profiles.reduce((n,p)=>n+p.contrasts.length,0),71);
const refit=A.study_profiles.find(p=>p.study==='PRJNA312093');
assert.equal(refit.samples.length,12);
assert.equal(refit.samples.some(s=>/adult|cercaria/.test(s.condition)),false);
assert.equal(refit.contrasts.map(c=>c.significant_genes).join(','),'290,1045,437');
assert.equal(A.study_profiles.find(p=>p.study==='PRJNA1177852').contrasts.length,0);
assert.equal(A.study_profiles.find(p=>p.study==='PRJNA428470').contrasts.length,1);
for(const study of ['PRJNA312093','PRJNA1177852','PRJNA428470']){
  const profile=A.study_profiles.find(p=>p.study===study);
  assert.ok(profile.interpretation_warning);
  t.renderStudyPage(profile);
  assert.ok(document.getElementById('study-details').innerHTML.includes(profile.interpretation_warning));
  assert.ok(t.geneExpressionPattern(A.genes.indexOf('Smp_179650'),study).every(Number.isFinite));
  assert.doesNotMatch(t.geneDistributionSvg(geneValues,t.geneConditionStats(geneValues,study)),/undefined/);
  assert.ok(t.studyWarningsHtml(study).includes(profile.interpretation_warning));
}
assert.equal(A.study_profiles.some(p=>['PRJNA494090','PRJNA551487'].includes(p.study)),false);
const sampleSection=fs.readFileSync(path.join(bundle,'atlas.html'),'utf8').match(/<section id="samples">([\s\S]*?)<\/section>/)[1];
assert.match(sampleSection,/global_pca_by_study\.svg/);
assert.match(sampleSection,/global_pca_by_stage\.svg/);
assert.doesNotMatch(sampleSection,/pca_by_condition\.svg/);
const uniqueCandidates=new Set(A.selected_gene_status.map(row=>row[A.selected_gene_fields.indexOf('resolved_gene_id')]).filter(Boolean));
for(const p of A.study_profiles){
  assert.equal(t.studyConditions(p).reduce((n,c)=>n+c.n,0),p.samples.length);
  assert.equal(Object.values(t.studyQuality(p)).reduce((n,c)=>n+c,0),p.samples.length);
  assert.ok(p.sources.length>0);
  assert.ok(p.samples.filter(r=>r.qc_classification==='REVIEW').every(r=>r.qc_reason&&r.qc_reason!=='none'));
  const studyPage=fs.readFileSync(path.join(bundle,`study_${p.study}.html`),'utf8');
  assert.ok(studyPage.includes(`data-study="${p.study}"`));
  t.renderStudyPage(p);
  const studyDetails=document.getElementById('study-details').innerHTML;
  assert.ok(studyDetails.includes(`figures/${p.study}_pca_by_condition.svg`));
  assert.doesNotMatch(studyDetails,/Epigenetic-gene expression heatmap|candidate_expression_heatmap\.svg/);
  for(const heading of ['Experimental design','Quality control','Contrast catalog','Provenance & downloads']){
    assert.ok(studyDetails.includes(`<h3>${heading}</h3>`),`${p.study}: ${heading}`);
  }
  assert.deepEqual(JSON.parse(fs.readFileSync(path.join(bundle,`data/study_${p.study}.json`))),JSON.parse(JSON.stringify(p)));
  for(const src of p.sources){
    assert.equal(require('node:crypto').createHash('sha256').update(fs.readFileSync(path.join(root,src.path))).digest('hex'),src.sha256);
  }
}
const blocked=A.study_profiles.find(p=>p.study==='PRJNA602528');
assert.equal(blocked.differential_expression,false);
assert.match(blocked.blocked_reason,/no biological replication/);
assert.equal(t.studyQuality(A.study_profiles.find(p=>p.study==='PRJEB32839')).REVIEW,21);
console.log('PASS: twelve study pages, approved samples and refit contrasts, warnings, exports, QC and source checksums.');
const atlasPage=fs.readFileSync(path.join(bundle,'atlas.html'),'utf8');
const machineryPage=fs.readFileSync(path.join(bundle,'machinery.html'),'utf8');
assert.match(atlasPage,/<body class="atlas-page hf-atlas" data-portal="atlas" data-assay="rna">/);
assert.match(atlasPage,/id="gene-collection"/);
assert.match(atlasPage,/src="gene_collections.js/);
assert.doesNotMatch(atlasPage,/<a class="portal-card" href="machinery.html">/);
assert.match(machineryPage,/p.set\('collection','machinery'\)/);
assert.match(machineryPage,/location.replace/);
const machineElements=new Map();
const machineDoc={body:{dataset:{portal:'machinery'}},documentElement:{dataset:{}},addEventListener(){},querySelectorAll(){return [];},getElementById(id){if(!machineElements.has(id))machineElements.set(id,element());return machineElements.get(id);},querySelector(){return element();},createElement:element};
const machineContext=vm.createContext({window:{HF_ATLAS:A},document:machineDoc,...browserGlobals});
vm.runInContext(script,machineContext);
const machineSearch=machineContext.window.test.findGeneMatches;
const outside=A.genes.find(gene=>!uniqueCandidates.has(gene));
assert.equal(machineSearch(outside).length,0);
assert.equal(machineSearch('Smp_179650').length,1);
assert.equal(machineSearch('chromobox').length,3);
console.log('PASS: two assay portals, shared collection entry, legacy machinery redirect and scoped search.');
t.setCollectionScope(['Smp_244880']);
assert.equal(t.findGeneMatches('Smp_244880').length,1);
assert.equal(t.findGeneMatches('Smp_179650').length,0);
t.setCollectionScope([]);
assert.equal(t.findGeneMatches('Smp_244880').length,0);
t.setCollectionScope(null);
assert.equal(t.findGeneMatches('Smp_179650').length,1);
assert.equal(A.chip_scope,'all_atlas_genes');
const observedChipGenes=new Set(A.chip_loci.map(row=>A.genes[row[0]]));
assert.equal(observedChipGenes.size,2014);
assert.equal([...observedChipGenes].filter(gene=>uniqueCandidates.has(gene)).length,21);
const firstChipIndex=A.genes.indexOf('Smp_244880');
assert.ok(A.chip_loci.some(row=>row[0]===firstChipIndex));
const firstChipAnnotation=Object.fromEntries(A.gene_annotation_fields.map((field,i)=>[field,A.gene_annotations[firstChipIndex][i]]));
const firstChipLoci=A.chip_loci.filter(row=>row[0]===firstChipIndex).map(row=>Object.fromEntries(A.chip_locus_fields.map((field,i)=>[field,row[i]])));
const actualChipSvg=t.geneLocusSvg(firstChipAnnotation,firstChipLoci,modelObjects(firstChipIndex));
assert.equal((actualChipSvg.match(/data-feature="peak"/g)||[]).length,firstChipLoci.length);
assert.match(actualChipSvg,/PRJNA1090249/);
assert.doesNotMatch(atlasPage,/PRELIMINARY_REANALYSIS_REQUIRED|Preliminary ChIP evidence/);
console.log('PASS: native HF loci outside the epigenetic collection are searchable and render with study provenance.');

// Exercise the collection adapter's real render/export handlers against frozen values.
(async()=>{
  const views=new Map(),exports=[];
  const el=id=>{if(!views.has(id))views.set(id,{...element(),disabled:false,onclick:null});return views.get(id);};
  el('cv-mode').value='expression';
  const ids=['Smp_244880','Smp_179650'];
  const adapter=vm.createContext({window:{GeneCollections:{esc:String,download:(name,text)=>exports.push({name,text}),mount:async(host,options)=>{options.onChange({name:'Test collection',ids});return {};}}},document:{getElementById:el,body:{dataset:{}}},URL,URLSearchParams,location:{search:'',href:'http://localhost/atlas.html'},history:{replaceState(){}}});
  vm.runInContext(fs.readFileSync(path.join(root,'analyses/atlas/collections/expression_chromatin.js'),'utf8'),adapter);
  vm.runInContext(fs.readFileSync(path.join(root,'analyses/atlas/collections/rna_collections.js'),'utf8'),adapter);
  await adapter.window.HFRNACollections.mount({A,index:new Map(A.genes.map((id,i)=>[id.toLowerCase(),i])),annotations:A.gene_annotations.map(row=>Object.fromEntries(A.gene_annotation_fields.map((f,i)=>[f,row[i]]))),studies:A.study_profiles.map(p=>p.study),cluster:t.clusteredGroupRows,dendrogram:t.clusterDendrogram,heatColor:t.heatColor,expression:t.groupExpressionHeatmap,effects:t.groupEffectHeatmap,contexts:t.groupContexts,contrasts:t.groupContrasts,mean:t.contextMean,scope:()=>{},open:()=>{}});
  assert.match(el('cv-chart').innerHTML,/Clustered gene collection expression heatmap/);
  el('cv-export').onclick();
  const lines=exports[0].text.trim().split('\n').map(line=>line.split('\t'));
  assert.equal(lines.length,3);
  for(const row of lines.slice(1)){const gi=A.genes.indexOf(row[0]);for(const [i,c] of t.groupContexts('PRJEB32839').entries()){const expected=c.indices.reduce((n,si)=>n+Math.log2(A.tpm[gi][si]+1),0)/c.indices.length;assert.ok(Math.abs(Number(row[i+1])-expected)<1e-10);}}
  el('cv-mode').value='effects';el('cv-mode').onchange();el('cv-export').onclick();
  assert.match(exports[1].name,/contrasts/);
  const contrasts=exports[1].text.trim().split('\n').slice(1).map(line=>line.split('\t'));
  assert.equal(contrasts.length,ids.length*t.groupContrasts('PRJEB32839').length);
  for(const row of contrasts){const gi=A.genes.indexOf(row[0]),c=t.groupContrasts(row[1]).find(c=>c.numerator===row[2]&&c.denominator===row[3]);assert.equal(row[4],String(A.lfc[gi][c.index]??'NA'));assert.equal(row[5],String(A.padj[gi][c.index]??'NA'));}
  el('cv-study').value='PRJNA602528';el('cv-study').onchange();assert.match(el('cv-chart').innerHTML,/no inferential contrasts/);
  el('cv-mode').value='integrated';el('cv-study').value='PRJEB32839';el('cv-chip').value='compatible';el('cv-support').value='1';el('cv-mode').onchange();el('cv-export').onclick();
  assert.match(el('cv-chart').innerHTML,/Aligned expression and chromatin/);assert.match(exports[2].name,/expression_chromatin/);assert.match(exports[2].text,/Smp_244880\tChIP/);assert.match(exports[2].text,/not_assessed/);
  console.log('PASS: collection expression exports reproduce sample means; contrast exports retain exact original effects and padj; descriptive studies stay non-inferential.');
})().catch(error=>{console.error(error);process.exitCode=1;});
