"""Import a compact, explicitly scoped publication; never promote it to full-run QC."""
import bisect
import csv
import gzip
import hashlib
import json
import subprocess
from collections import Counter, defaultdict
from hf_chipseq_regional import rows, bed, SOURCE_FIELDS


def freeze_catalog(repository, revision, study, output):
    commit = subprocess.check_output(['git','rev-parse',revision+'^{commit}'],cwd=repository,text=True).strip()
    names = subprocess.check_output(['git','ls-tree','-r','--name-only',commit],cwd=repository,text=True).splitlines()
    prefix = f'results/chipseq/{study}/'
    selected = [p for p in names if p.startswith(prefix) and not p.endswith('.html') or p in
                {f'metadata/chipseq/{study}/design.tsv',f'config/chipseq/{study}/db_review.json',f'docs/chipseq/{study.lower()}_differential_binding.md'}]
    records=[]
    for name in selected:
        data=subprocess.check_output(['git','show',f'{commit}:{name}'],cwd=repository)
        path=output/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(data)
        records.append(dict(path=name,sha256=hashlib.sha256(data).hexdigest(),bytes=len(data)))
    manifest=dict(schema_version='1.0',study_id=study,repository='https://github.com/GMiguelAlves/HelixForge-SMansoni',commit=commit,
                  scope='published_descriptive_catalog_only',representation='nfcore_2_1_min1_homer_4_11',files=records)
    (output/'source_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    return manifest


def load_catalog(root,result,entry,source,config,annotations):
    if entry.get('admission_scope')!='published_descriptive_catalog_only':raise ValueError('Compact publication needs explicit descriptive admission scope')
    review=json.loads((result.parents[2]/f'config/chipseq/{entry["id"]}/db_review.json').read_text())
    if review['atlas_peak_inclusion']!='ACCEPTED' or review['peak_processing_status']!='COMPLETE' or review['differential_binding_executed']:raise ValueError('Catalog is not accepted or DB scope changed')
    publication=json.loads((result/'manifests/regional_publication.json').read_text())
    execution=json.loads((result/'narrow/manifests/native_execution_summary.json').read_text())
    if publication['status']!='PUBLISHED' or publication['publication_checks']!='PASS' or publication['source_workflow_validation']!='PASS':raise ValueError('Unpublished or failed catalog')
    if publication['source_direct_equivalence']!='FAIL_DIRECT_EQUIVALENCE' or publication['strict_gate_redefined']:raise ValueError('Picard comparison record changed; review source')
    failed=[r for r in publication['source_comparison'] if r['technical_validation']!='PASS']
    if {(r['condition'],r['direct_changed_peaks']) for r in failed}!={('male_adult',1),('male_cercaria',1)}:raise ValueError('Unreviewed direct comparison discrepancy')
    if execution['native_terminal_validation']!={'filesystem':'tracked_inputs_verified','schema':'valid','semantic':'valid'}:raise ValueError('Narrow native validation changed')
    if publication['differential_binding_status']!='HOLD' or execution['differential_binding']!='HOLD_NOT_EXECUTED':raise ValueError('DB status requires a new atlas review')
    reference=json.loads((root/'config/reference/chipseq_chromosomes_wbps19.json').read_text());chroms=reference['chromosomes']
    output={k:[] for k in ('metadata','evidence','positions','loci','segments','profiles','library_peaks','core_segments','regional_support','idr_segments')}
    study=entry['id'];summaries={};positions=defaultdict(set);samples={};contexts=[];all_groups=[];ambiguity_counts=Counter()
    def interval(c,s,e):
        if c not in chroms or not 0<=int(s)<int(e)<=chroms[c]:raise ValueError('Catalog interval outside reference')
    for arm,groups in [('broad',publication['groups']),('narrow',execution['groups'])]:
      for group in groups:
        stage,factor=group['condition'],group['target'];context=f'{study}|{stage}|{factor}';directory=result/f'{arm}/chipseq/consensus/{stage}'
        m=json.loads((directory/'manifest.json').read_text());all_groups.append(group)
        if m.get('reference_id',m.get('genome_id'))!=config['reference_id'] or m['reference_sha256']!=reference['reference_fasta_sha256'] or m['annotation_sha256']!='0ba2e3ff9a6328080a20f7b91bc7aa66c558eda5fdb6c382349d5e1f1a5f2d0e':raise ValueError('Catalog reference mismatch')
        if m['regional_min_libraries']!=1 or m['core_min_libraries']!=2 or m['annotation_provider']!='homer_4_11':raise ValueError('Catalog method mismatch')
        for artifact in m['artifacts'].values():
            p=result/artifact['path']
            if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=artifact['sha256']:raise ValueError('Catalog artifact checksum mismatch')
        libs={r['record_id']:r for r in m['library_metadata']}
        if len(libs)!=2 or group['total_library_count']!=2:raise ValueError('Unexpected deposited library count')
        for rid,lib in libs.items():
            samples[rid]=dict(run_accession=rid,sample_id=lib['sample_id'],condition=stage,target=factor,is_control='false',qc_classification='CATALOG_ONLY_NOT_FULL_QC')
            cid=lib['control_record_id'];samples[cid]=dict(run_accession=cid,sample_id=cid,condition=stage,target='Input',is_control='true',qc_classification='CATALOG_ONLY_NOT_FULL_QC')
        native={r['peak_id']:r for r in rows(directory/'consolidated_peaks.tsv.gz')};annotated=rows(result/f'{arm}/chipseq/peak_annotation/{stage}/annotated_peaks.tsv.gz')
        if len(native)!=group['regions'] or len(annotated)!=len(native):raise ValueError('Catalog region count mismatch')
        support=defaultdict(dict)
        for r in rows(directory/'support_by_library.tsv.gz'):
            pid,rid=r['peak_id'],r['record_id']
            if pid not in native or rid not in libs or rid in support[pid]:raise ValueError('Catalog support identities mismatch')
            length=int(native[pid]['end'])-int(native[pid]['start']);covered=int(r['covered_bp']);called=int(r['called_peak_count'])
            if int(r['region_bp'])!=length or not 0<=covered<=length or abs(float(r['covered_fraction'])-covered/length)>1e-5 or (called==0 and covered!=0):raise ValueError('Invalid native support coverage')
            support[pid][rid]=r
            output['regional_support'].append(r)
        if len(support)!=len(native) or any(set(v)!=set(libs) for v in support.values()):raise ValueError('Catalog lost zero library-support rows')
        assigned=set();seen=set();flags=Counter();categories=Counter()
        for r in annotated:
            pid=r['peak_id'];n=native[pid];interval(r['chrom'],r['start'],r['end'])
            if any(r[k]!=n[k] for k in ('chrom','start','end')):raise ValueError('Catalog coordinates changed')
            calling=sorted(rid for rid,d in support[pid].items() if int(d['called_peak_count'])>0)
            expected='single_library' if len(calling)==1 else 'all_libraries'
            if not calling or int(n['support'])!=len(calling) or r['calling_libraries']!=','.join(calling) or r['library_support_category']!=expected or int(r['calling_library_count'])!=len(calling):raise ValueError('Catalog calling libraries changed')
            for d in json.loads(r['library_support_json']):
                if any(str(d[k])!=str(support[pid][d['record_id']][k]) for k in ('record_id','called_peak_count','covered_bp','region_bp','covered_fraction','region_core_bp')):raise ValueError('Catalog annotation/support evidence differs')
            flag=r['gene_assignment_ambiguous'];flags[flag]+=1;ambiguity_counts[flag]+=1
            if flag not in {'true','false','not_assessed'} or arm=='broad' and flag!='not_assessed':raise ValueError('Catalog ambiguity state changed')
            alternatives=json.loads(r['alternative_gene_ids'])
            if flag=='true' and not alternatives or any(g not in annotations for g in alternatives):raise ValueError('Catalog ambiguity alternatives invalid')
            genes=list(filter(None,r['gene_ids'].split(';')))
            if len(genes)!=int(r['gene_count']) or any(g not in annotations or annotations[g]['chromosome']!=r['chrom'] for g in genes):raise ValueError('Catalog gene reference mismatch')
            if r['annotation_provider']!='homer_4_11' or r['association_method']!='homer_nearest_promoter' or r['occupancy_status']!='not_assessed':raise ValueError('Catalog method/status changed')
            assigned.update((pid,g) for g in genes);categories[r['category']]+=1
            extras={k:r[k] for k in SOURCE_FIELDS}
            output['segments'].append(dict(study_id=study,context_id=context,stage_or_condition=stage,mark_or_factor=factor,**{k:r[k] for k in ('peak_id','chrom','start','end','category','gene_ids')},support=len(calling),support_libraries=r['calling_libraries'],consensus_method='nfcore_2_1',**extras))
            for gene in genes:
                category=r['category'];position='promoter_TSS' if category=='promoter' else 'gene_body' if category in {'exon','intron','gene'} else 'distal_intergenic' if category=='intergenic' else 'other_annotated'
                row=dict(study_id=study,context_id=context,gene_id=gene,stage_or_condition=stage,mark_or_factor=factor,peak_id=pid,peak_chrom=r['chrom'],peak_start=int(r['start']),peak_end=int(r['end']),peak_midpoint=(int(r['start'])+int(r['end']))//2,peak_location=category,promoter_flag=str(category=='promoter').lower(),position_class=position,support=len(calling),**extras)
                output['loci'].append(row);positions[(gene,context,stage,factor,position)].add(pid)
                s=summaries.setdefault((gene,context),dict(study_id=study,context_id=context,gene_id=gene,functional_name=annotations[gene]['functional_name'],mark_or_factor=factor,stage_or_condition=stage,n_peaks=0,n_promoter_peaks=0,peak_locations=set(),regulatory_class='candidate_nearest_promoter_association',expected_effect='not_inferred_from_occupancy',n_single_library=0,n_multiple_libraries=0,n_ambiguous=0,n_ambiguity_not_assessed=0))
                s['n_peaks']+=1;s['n_promoter_peaks']+=category=='promoter';s['peak_locations'].add(category);s['n_single_library']+=len(calling)==1;s['n_multiple_libraries']+=len(calling)>1;s['n_ambiguous']+=flag=='true';s['n_ambiguity_not_assessed']+=flag=='not_assessed'
        for r in rows(result/f'{arm}/chipseq/peak_annotation/{stage}/peak_gene_associations.tsv.gz'):
            key=r['peak_id'],r['gene_id']
            if key in seen or key not in assigned:raise ValueError('Catalog association identity mismatch')
            seen.add(key)
        if assigned!=seen or len({g for _,g in assigned})!=group['associated_gene_ids'] or sum(int(n['end'])-int(n['start']) for n in native.values())!=group['region_bp']:raise ValueError('Catalog association/coverage totals changed')
        if arm=='narrow' and (flags['true']!=group['ambiguous_gene_assignments'] or flags['not_assessed']!=group['unassessed_gene_assignments']):raise ValueError('Catalog ambiguity counts changed')
        by_chrom=defaultdict(list)
        for r in native.values():by_chrom[r['chrom']].append((int(r['start']),int(r['end'])))
        for c in by_chrom:by_chrom[c].sort()
        starts={c:[r[0] for r in v] for c,v in by_chrom.items()}
        core=bed(directory/'support_core.bed.gz')
        for c,s,e,*_ in core:
            interval(c,s,e);i=bisect.bisect_right(starts[c],int(s))-1
            if i<0 or by_chrom[c][i][1]<int(e):raise ValueError('Core outside published regions')
            output['core_segments'].append(dict(study_id=study,stage_or_condition=stage,mark_or_factor=factor,chrom=c,start=int(s),end=int(e)))
        if sum(int(r[2])-int(r[1]) for r in core)!=group['core_bp']:raise ValueError('Catalog core coverage changed')
        if arm=='narrow':
            idr=rows(result/f'narrow/chipseq/idr/{stage}/consolidated_peaks.tsv.gz')
            if len(idr)!=group['idr_regions'] or any(float(r['global_idr'])>0.05 for r in idr):raise ValueError('Native IDR evidence mismatch')
            for r in idr:
                interval(r['chrom'],r['start'],r['end']);output['idr_segments'].append(dict(study_id=study,stage_or_condition=stage,mark_or_factor=factor,chrom=r['chrom'],start=int(r['start']),end=int(r['end']),peak_id=r['peak_id'],global_idr=float(r['global_idr'])))
        output['metadata'].append(dict(study_id=study,context_id=context,mark_or_factor=factor,stage_or_condition=stage,n_samples=2,n_controls=2,replicates=';'.join(r['sample_id'] for r in libs.values()),batches='review_required',qc_status='PUBLISHED_CATALOG_WITH_LIMITATIONS',consensus_method='nfcore_2_1',annotation_provider='homer_4_11',min_libraries=1))
        contexts.append(dict(study=study,stage=stage,factor=factor,method='nfcore_2_1',min_libraries=1,region_count=group['regions'],idr_count=group.get('idr_regions',0),source_comparison='NOT_PERFORMED' if arm=='narrow' else next(r['technical_validation'] for r in publication['source_comparison'] if r['condition']==stage),libraries=[dict(id=rid,sample=lib['sample_id'],label=f'Library {i+1}',coordinates_available=False) for i,(rid,lib) in enumerate(libs.items())]))
    output['evidence']=[{**r,'peak_locations':';'.join(sorted(r['peak_locations']))} for r in summaries.values()]
    output['positions']=[dict(study_id=study,context_id=c,gene_id=g,stage_or_condition=s,mark_or_factor=f,position_class=p,n_unique_peaks=len(ids)) for (g,c,s,f,p),ids in positions.items()]
    order=lambda stage,factor:(0 if stage.endswith('cercaria') else 1,stage.split('_')[0],factor)
    output['metadata'].sort(key=lambda r:order(r['stage_or_condition'],r['mark_or_factor']))
    contexts.sort(key=lambda r:order(r['stage'],r['factor']))
    output['profiles'].append(dict(study_id=study,reference_id=config['reference_id'],target='H3K27me3; H3K4me3',targets=['H3K27me3','H3K4me3'],repository=source['repository'],commit=source['commit'],admission_scope=entry['admission_scope'],qc=dict(ip_libraries=16,input_libraries=8,qc_classification='PUBLISHED_CATALOG_WITH_LIMITATIONS',differential_binding='HOLD',technical_execution='PARTIAL_ARM_REPORTS_NOT_COMBINED_CERTIFICATION'),helixforge_release='source-specific; see group manifests',helixforge_commit=execution['helixforge_commit'],regional_execution_commit='H3K27me3: cd390a0; H3K4me3: '+execution['helixforge_commit'],caller='macs3',caller_version='3.0.4',consensus_method='nfcore_2_1',annotation_provider='homer_4_11',regional_min_libraries=1,core_min_libraries=2,contexts=contexts,samples=list(samples.values()),annotation_policies=[],limitations=entry['interpretation_warning'],consensus_segments=len(output['segments']),intergenic_segments=sum(r['category']=='intergenic' for r in output['segments']),gene_associations=len(output['loci']),genes_with_associations=len({r['gene_id'] for r in output['loci']}),stage_segments=dict(Counter(r['stage_or_condition'] for r in output['segments'])),coordinate_convention='0-based_half-open',scope='all_atlas_genes',publication={**publication,'groups':all_groups},source_files=source['files'],native_execution=execution,ambiguity_counts=dict(ambiguity_counts)))
    output['profiles'][-1].update(db_review=review,peak_qc=rows(result/'narrow/qc/peak_qc_summary.tsv.gz'))
    return output


def render_catalog_page(p):
    from html import escape as e
    groups=p['publication']['groups']
    table=''.join(f'<tr><td>{e(g["target"])}</td><td>{e(g["condition"])}</td><td>{g["regions"]:,}</td><td>{g["single_library_regions"]:,}</td><td>{g["all_library_regions"]:,}</td><td>{g["associated_gene_ids"]:,}</td><td>{g.get("idr_regions","—")}</td></tr>' for g in groups)
    qc=''.join(f'<tr><td>{e(r["record_id"])}</td><td>{e(r["sample_id"])}</td><td>{100*float(r["frip"]):.2f}%</td><td>{int(r["peak_count"]):,}</td></tr>' for r in p['peak_qc'])
    sources=''.join(f'<tr><td><a href="{p["repository"]}/blob/{p["commit"]}/{r["path"]}">{e(r["path"])}</a></td><td>{r["sha256"]}</td></tr>' for r in p['source_files'])
    return f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{p['study_id']} · H3K27me3 / H3K4me3</title><style>body{{margin:0;font:15px/1.6 system-ui;color:#234451;background:#f2f7f7}}header{{padding:36px max(5%,calc((100% - 1180px)/2));background:#103f46;color:white}}h1{{margin:8px 0}}main,nav{{max-width:1180px;margin:24px auto}}nav,a{{color:#12665d}}nav{{display:flex;gap:25px}}section{{background:white;padding:24px;border:1px solid #d6e2e4;border-radius:14px;margin:20px 0}}.notice{{border-left:4px solid #bd944b}}.scroll{{overflow:auto}}table{{border-collapse:collapse;width:100%;font-size:13px}}td,th{{padding:10px;text-align:left;border-bottom:1px solid #e4ecee}}.counts{{display:flex;gap:35px;flex-wrap:wrap}}.counts strong{{font-size:27px;display:block}}details{{margin:15px 0}}summary{{cursor:pointer}}.hash td:last-child{{font:10px monospace;word-break:break-all}}</style><header><small>HELIXFORGE · PUBLISHED CHIP-SEQ CATALOGS</small><h1>{p['study_id']} · H3K27me3 / H3K4me3</h1></header><nav><a href="chipseq.html?study={p['study_id']}">Explore regions →</a><a href="atlas.html?view=integration#collection">Expression + chromatin →</a><a href="chipseq.html#studies">Browse ChIP studies</a></nav><main>
<section class="counts"><div><strong>16</strong>IP libraries</div><div><strong>8</strong>shared matched Inputs</div><div><strong>{p['consensus_segments']:,}</strong>regional candidates</div><div><strong>{p['genes_with_associations']:,}</strong>genes with source assignments</div></section>
<section class="notice"><strong>Descriptive catalogs · Differential Binding on HOLD</strong><details><summary>Scope, QC and interpretation</summary><p>{e(p['limitations'])}</p></details></section>
<section><h2>Experimental design</h2><p>Female and male adults and cercariae; two deposited IP libraries per sex / stage / mark, with matched Input identifiers. H3K27me3 uses broad calls; H3K4me3 uses narrow calls. Library labels alone do not certify independence.</p><details><summary>Processing and reference</summary><p>{e(p['reference_id'])} · MACS3 3.0.4 · nf-core/chipseq 2.1.0 minimum-one regional catalog · HOMER 4.11. H3K27me3 source: cd390a0; H3K4me3 source: {e(p['native_execution']['helixforge_commit'])}. The narrow continuation passed its declared terminal validation; this does not certify upstream alignment QC or a combined run for both marks.</p></details></section>
<section><h2>Regional candidates and separate IDR evidence</h2><div class="scroll"><table><thead><tr><th>Mark</th><th>Context</th><th>Regions</th><th>One library</th><th>Both libraries</th><th>Assigned gene IDs</th><th>IDR ≤0.05</th></tr></thead><tbody>{table}</tbody></table></div><details><summary>How to interpret support and IDR</summary><p>The annotation input is the complete minimum-one regional catalog. Shared bases and IDR-selected intervals remain separate products. IDR is ranked by signal value at 0.05; it does not certify every candidate region. Genomic coordinates and per-library covered bases are preserved. Individual-call coordinates and bigWigs are not in this compact public package; their absence is not a negative result.</p></details></section>
<section><h2>Published QC and DB review</h2><p>H3K4me3 per-library Peak QC; upstream alignment is not certified by this compact continuation report. FRiP measures the fraction of selected reads in called peaks, not biological replication.</p><div class="scroll"><table><tr><th>Library</th><th>Sample</th><th>FRiP</th><th>Called peaks</th></tr>{qc}</table></div><details><summary>Comparison-specific eligibility</summary><p>DB remains HOLD. Adult sex and sequencing instrument are completely confounded. Selected female source pools are unresolved; H3K27me3 also needs enrichment and model review. No comparison is certified for inferential DB.</p></details></section><section><h2>Gene-assignment ambiguity and source validation</h2><p>H3K27me3: not assessed. Its strict comparison retains FAIL_DIRECT_EQUIVALENCE for one nearest-gene choice in each male group (Smp_095580 / Smp_341610 share a TSS). No source assignment was changed to force agreement.</p><p>H3K4me3: native TSS snapshots preserve true, false and not_assessed states. Counts across both arms: {e(str(p['ambiguity_counts']))}.</p><details><summary>RNA integration boundaries</summary><p>The approved RNA contribution from this accession contains schistosomula only. Adult/cercaria RNA and sex contrasts remain excluded. Integration with other approved studies is by stage and recorded sex, without sample pairing; mixed-sex RNA is not matched to a sex-specific ChIP group. No regulatory direction is inferred from counts or marks.</p></details></section>
<section><h2>Provenance and downloads</h2><p>Source commit: <code>{p['commit']}</code></p><a href="data/chipseq_{p['study_id']}.json">Complete study provenance</a> · <a href="data/chipseq_regional_candidates.tsv">All regional candidates</a> · <a href="data/chipseq_gene_peak_loci.tsv">Native gene assignments</a> · <a href="data/chipseq_idr.tsv">Separate IDR intervals</a> · <a href="data/chipseq_regional_library_support.tsv">Library coverage</a> · <a href="{p['repository']}/tree/{p['commit']}/results/chipseq/{p['study_id']}">Published HF package</a><details><summary>Frozen source files and checksums</summary><div class="scroll"><table class="hash">{sources}</table></div></details></section></main></html>'''
