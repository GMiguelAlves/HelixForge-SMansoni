"""Keep large catalogs out of initial browser loads; preserve lossless native rows."""
import hashlib
import json
import os
import gzip
import shutil
from collections import defaultdict, Counter


def dump(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,separators=(',',':'),ensure_ascii=False)+'\n',encoding='utf-8')


def pack(fields,records):
    defaults={f:records[0][i] for i,f in enumerate(fields) if records and all(r[i]==records[0][i] for r in records)}
    varying=[i for i,f in enumerate(fields) if f not in defaults]
    return dict(defaults=defaults,fields=[fields[i] for i in varying],rows=[[r[i] for i in varying] for r in records])


def export_lazy(output,chip_data,fields,segments,locus_fields):
    delayed={p['study_id'] for p in chip_data['profiles'] if p.get('admission_scope')=='published_descriptive_catalog_only'}
    by_key={(c['study'],c['stage'],c['factor']):c for p in chip_data['profiles'] for c in p.get('contexts',[])}
    contexts=[by_key[(r['study_id'],r['stage_or_condition'],r['mark_or_factor'])] for r in chip_data['metadata']]
    context_index={(c['study'],c['stage'],c['factor']):i for i,c in enumerate(contexts)}
    if list(context_index) != [(r['study_id'],r['stage_or_condition'],r['mark_or_factor']) for r in chip_data['metadata']]:raise ValueError('ChIP context index differs from metadata order')
    native_segments=defaultdict(list)
    for r,native in zip(segments,chip_data['segments']):
        if native['study_id'] in delayed:native_segments[context_index[(native['study_id'],native['stage_or_condition'],native['mark_or_factor'])]].append(r)
    categories=sorted({r['peak_location'] for r in chip_data['loci']});summaries=defaultdict(lambda:[0,0,0,0,0])
    for r in chip_data['loci']:
        c=context_index[(r['study_id'],r['stage_or_condition'],r['mark_or_factor'])]
        s=summaries[(r['gene_id'],c,int(r['support']))];s[0]+=1;s[1]+=r['promoter_flag']=='true';s[2]+=r.get('gene_assignment_ambiguous')=='true';s[3]+=r.get('gene_assignment_ambiguous')=='not_assessed';s[4]|=1<<categories.index(r['peak_location'])
    summary_rows=[[g,c,s,*values] for (g,c,s),values in sorted(summaries.items())]
    descriptors=[];support_fields=['record_id','sample_id','called_peak_count','covered_bp','covered_fraction','region_bp','region_core_bp']
    support=defaultdict(list)
    for r in chip_data['regional_support']:support[r['peak_id']].append([r[k] for k in support_fields])
    cores=defaultdict(list);idr=defaultdict(list)
    for r in chip_data['core_segments']:cores[(r['study_id'],r['stage_or_condition'],r['mark_or_factor'])].append([r['chrom'],r['start'],r['end']])
    for r in chip_data.get('idr_segments',[]):idr[(r['study_id'],r['stage_or_condition'],r['mark_or_factor'])].append([r['chrom'],r['start'],r['end'],r['peak_id'],r['global_idr']])
    for i,records in native_segments.items():
        c=contexts[i];key=(c['study'],c['stage'],c['factor']);slug=f"{c['study']}_{c['factor']}_{c['stage']}";prefix=os.path.commonprefix([r[fields.index('id')] for r in records]).rsplit('Interval_',1)[0]+'Interval_'
        short=[];buckets=defaultdict(dict)
        for r in records:
            copy=list(r);pid=copy[fields.index('id')];suffix=pid[len(prefix):];copy[fields.index('id')]=suffix;short.append(copy)
            buckets[int(suffix)//1000][suffix]=support[pid]
        url=f'data/chipseq/catalogs/{slug}.json';packed=pack(fields,short)
        packed.update(id_prefix=prefix,cores=cores[key],idr=idr[key],context=c)
        dump(output/url,packed)
        base=f'data/chipseq/support/{slug}_'
        for b,content in buckets.items():dump(output/f'{base}{b}.json',dict(fields=support_fields,id_prefix=prefix,regions=content))
        descriptors.append(dict(context_index=i,context=c,url=url,regions=len(records),id_prefix=prefix,support_base=base,support_bucket_size=1000))
    gene_buckets=defaultdict(list)
    for r in chip_data['loci']:
        if r['study_id'] in delayed:gene_buckets[int(r['gene_id'].split('_')[1])%64].append([r[f] for f in locus_fields])
    for b in range(64):
        records=gene_buckets.get(b,[])
        prefixes=[];column=locus_fields.index('peak_id')
        for row in records:
            prefix,suffix=row[column].rsplit('Interval_',1)
            if prefix not in prefixes:prefixes.append(prefix)
            row[column]=[prefixes.index(prefix),suffix]
        packed=pack(locus_fields,records);packed['peak_id_prefixes']=prefixes
        dump(output/f'data/chipseq/genes/{b}.json',packed)
    return dict(delayed_studies=sorted(delayed),lazy_catalogs=descriptors,summary_fields=['gene','context','support','count','promoter','ambiguous','unassessed','category_mask'],gene_summaries=summary_rows,category_values=categories,gene_store=dict(pattern='data/chipseq/genes/{bucket}.json',buckets=64),contexts=contexts)


def compress_downloads(output):
    """Large lossless source tables are downloads, never initial browser data."""
    names=['chipseq_gene_peak_loci.tsv','chipseq_all_consensus_segments.tsv','chipseq_regional_candidates.tsv','chipseq_regional_library_support.tsv','chipseq_gene_position_mark.tsv','chipseq_idr.tsv']
    replacements={}
    for name in names:
        path=output/'data'/name
        if not path.is_file():continue
        target=path.with_suffix('.tsv.gz')
        with path.open('rb') as source,target.open('wb') as raw,gzip.GzipFile(filename='',fileobj=raw,mode='wb',mtime=0) as compressed:shutil.copyfileobj(source,compressed)
        path.unlink();replacements['data/'+name]='data/'+name+'.gz'
    for path in list(output.glob('*.html'))+list(output.glob('*.md')):
        text=path.read_text(encoding='utf-8')
        for old,new in replacements.items():text=text.replace(old+'"',new+'"').replace(old+')',new+')')
        path.write_text(text,encoding='utf-8')
