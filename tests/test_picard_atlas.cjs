// Lossless lazy delivery and biological joins, using only the actual frozen catalogs.
const fs=require('node:fs'),path=require('node:path'),zlib=require('node:zlib'),vm=require('node:vm'),assert=require('node:assert/strict');
const D=require('../analyses/atlas/collections/chip_data.js'),E=require('../analyses/atlas/collections/expression_chromatin.js');
const root=path.resolve(__dirname,'..'),bundle=path.resolve(process.argv[2]||path.join(root,'results/atlas/rnaseq'));
const read=p=>JSON.parse(fs.readFileSync(path.join(bundle,p),'utf8')),A=read('data/chipseq/explorer.json');
assert.equal(A.total_regions,244788);assert.equal(A.consensus.length,2683);assert.equal(A.contexts.length,10);assert.equal(A.lazy_catalogs.length,8);
assert.ok(A.lazy_catalogs.every(d=>d.context.libraries.every(l=>l.coordinates_available===false)));
const nativeRoot=path.join(root,'analyses/atlas/local_inputs/hf_chipseq/PRJNA312093/regional_e583c84/results/chipseq/PRJNA312093');
const tsv=p=>{const text=zlib.gunzipSync(fs.readFileSync(p)).toString(),lines=text.trimEnd().split('\n'),fields=lines.shift().replace(/\r$/,'').split('\t');return lines.map(line=>{const values=line.replace(/\r$/,'').split('\t');return Object.fromEntries(fields.map((f,i)=>[f,values[i]]));});};
const regions=new Map(),flags={true:0,false:0,not_assessed:0};let supportCount=0,idrCount=0;
for(const d of A.lazy_catalogs){
 const data=read(d.url),records=D.unpack(data).map(r=>({...r,id:data.id_prefix+r.id})),arm=d.context.factor==='H3K4me3'?'narrow':'broad';
 const native=tsv(path.join(nativeRoot,arm,'chipseq/consensus',d.context.stage,'consolidated_peaks.tsv.gz'));
 assert.equal(records.length,d.regions);assert.equal(records.length,native.length);
 const expected=new Map(native.map(r=>[r.peak_id,r]));
 for(const r of records){const n=expected.get(r.id);assert.ok(n);assert.equal(r.chrom,n.chrom);assert.equal(r.start,Number(n.start));assert.equal(r.end,Number(n.end));assert.equal(r.support,Number(n.support));assert.ok(['true','false','not_assessed'].includes(r.ambiguity));flags[r.ambiguity]++;regions.set(r.id,r);}
 assert.equal(data.idr.length,d.context.idr_count);assert.ok(data.idr.every(r=>r[4]<=.05));idrCount+=data.idr.length;
 const files=fs.readdirSync(path.join(bundle,path.dirname(d.support_base))).filter(f=>f.startsWith(path.basename(d.support_base)));
 for(const file of files){const evidence=read(path.join(path.dirname(d.support_base),file)),byId=evidence.regions;
  for(const [id,rows] of Object.entries(byId)){assert.equal(rows.length,2);assert.ok(expected.has(evidence.id_prefix+id));supportCount+=rows.length;}
 }
}
assert.equal(regions.size,242105);assert.equal(supportCount,484210);assert.equal(idrCount,18179);assert.deepEqual(flags,{true:22,false:52644,not_assessed:189439});
let geneRows=0;
for(let b=0;b<64;b++)for(const r of D.unpack(read(`data/chipseq/genes/${b}.json`))){const n=regions.get(r.peak_id);assert.ok(n);assert.equal(r.peak_chrom,n.chrom);assert.equal(r.peak_start,n.start);assert.equal(r.peak_end,n.end);assert.equal(r.gene_assignment_ambiguous,n.ambiguity);assert.ok(n.genes.split(';').includes(r.gene_id));assert.equal(r.alternative_gene_ids,n.alternatives);geneRows++;}
assert.equal(geneRows,242105);
const context=vm.createContext({window:{}});vm.runInContext(fs.readFileSync(path.join(bundle,'data/atlas_payload.js'),'utf8'),context);const R=context.window.HF_ATLAS;
assert.equal(R.chip_loci.length,2683);assert.ok(R.chip_summaries.length>60000);const rnaGenes=new Set(R.genes);assert.equal(D.evidence(R).reduce((n,r)=>n+r.n_peaks,0),[...regions.values()].filter(r=>r.genes.split(';').every(g=>rnaGenes.has(g))).length);
const samples=R.samples.map(r=>Object.fromEntries(R.sample_fields.map((f,i)=>[f,r[i]])));
const api={contexts:study=>{const groups=new Map();samples.forEach((s,i)=>{if(s.study!==study)return;if(!groups.has(s.condition))groups.set(s.condition,{study,label:s.condition,indices:[]});groups.get(s.condition).indices.push(i);});return [...groups.values()];},cluster:rows=>({rows,cluster:{root:null,order:[]}})};
const engine=E.create(R,api),preset=read('data/gene_collections.json').presets.machinery.ids,rows=preset.map(gene=>({gene,index:R.genes.indexOf(gene)}));
const model=engine.build(rows,'PRJEB32839');assert.equal(model.chip.length,6);
const female=model.chip.find(c=>c.stage_or_condition==='female_cercaria'),rnaFemale=model.rna.find(c=>c.label==='cercaria_female'),rnaMale=model.rna.find(c=>c.label==='cercaria_male'),mixed=model.rna.find(c=>c.label==='cercaria_mixed');
assert.notEqual(engine.match(rnaFemale,female),'NOT_COVERED');assert.equal(engine.match(rnaMale,female),'NOT_COVERED');assert.equal(engine.match(mixed,female),'NOT_COVERED');
const values=engine.exportRows(model).slice(1).filter(r=>r[1]==='ChIP');assert.ok(values.some(r=>r[9]>0));assert.ok(values.filter(r=>r[9]>0).every(r=>r[12]!=='NO_ASSIGNED_REGION'));
const selected=model.values.find(v=>v.chip.some((c,i)=>model.chip[i].study_id==='PRJNA312093'&&c.count)),before=selected.chip.map(c=>c.count),gene=selected.row.gene;
const details=D.unpack(read(`data/chipseq/genes/${Number(gene.split('_')[1])%64}.json`)).filter(r=>r.gene_id===gene).map(r=>({...r,gene_index:R.genes.indexOf(gene)}));engine.addLoci(details);
assert.deepEqual(engine.build(rows,'PRJEB32839').values.find(v=>v.row.gene===gene).chip.map(c=>c.count),before);
console.log('PASS: 242,105 native regions, 484,210 support rows, 18,179 separate IDR intervals, exact lazy gene assignments/ambiguity, summary counts and sex-aware RNA joins.');
