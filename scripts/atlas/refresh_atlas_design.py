"""Refresh atlas presentation from an existing validated bundle, preserving data bytes."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
from atlas_design import apply_atlas_design


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def refresh(root, output):
    manifest = json.loads((output/'manifest.json').read_text(encoding='utf-8'))
    if manifest['status'] != 'PASS':
        raise ValueError('Presentation refresh requires a previously validated atlas.')
    scientific = [r for r in manifest['outputs'] if r['path'].startswith(('data/','figures/'))]
    template = (root/'analyses/atlas/rnaseq_atlas_template.html').read_text(encoding='utf-8')
    header = re.search(r'<header>[\s\S]*?</header>',template)[0]
    for key,value in {'STUDY_COUNT':len(manifest['studies']),'SAMPLE_COUNT':manifest['biological_samples'],'GENE_COUNT':manifest['genes'],'CONTRAST_COUNT':manifest['contrasts']}.items():
        header = header.replace('{{'+key+'}}',str(value))
    if '{{' in header:
        raise ValueError('Unresolved presentation placeholder')
    style = re.search(r'<style>[\s\S]*?</style>',template)[0]
    script = re.search(r'<script>\s*[\s\S]*?</script>',template)[0]
    for page in [output/'atlas.html',*sorted(output.glob('study_*.html'))]:
        text = page.read_text(encoding='utf-8')
        study = page.stem.split('_',1)[1] if page.stem.startswith('study_') else None
        text = re.sub(r'<header>[\s\S]*?</header>',lambda m:header.replace('<h1>Gene expression</h1>',f'<h1>{study}</h1>') if study else header,text,count=1)
        text = re.sub(r'<style>[\s\S]*?</style>',lambda m:style,text,count=1)
        text = re.sub(r'<script>\s*[\s\S]*?</script>',lambda m:script,text,count=1)
        collection = re.search(r'\s*<section id="collection">[\s\S]*?</section>',text)
        if collection:
            block=collection[0];text=text[:collection.start()]+text[collection.end():]
            text=text.replace('    <section id="study-details"',block+'\n    <section id="study-details"',1)
        page.write_text(text,encoding='utf-8')
    shutil.copy2(root/'analyses/atlas/chipseq_portal/chipseq.html',output/'chipseq.html')
    asset_inputs=[]
    for directory,names in [('collections',['gene_collections.js','gene_collections.css','rna_collections.js','expression_chromatin.js','chip_data.js']),('chipseq_portal',['chipseq.css','chipseq.js','chipseq_logic.js','chipseq_tracks.js'])]:
        for name in names:
            source=root/'analyses/atlas'/directory/name
            shutil.copy2(source,output/name);asset_inputs.append(source)
    inputs=apply_atlas_design(root,output)
    for page in output.glob('*.html'):
        text=page.read_text(encoding='utf-8')
        # The canonical template names raw TSVs; this bundle uses lossless gzip downloads.
        for relative in set(re.findall(r'(?:src|href)="([^"?]+\.tsv)"',text)):
            if not (output/relative).is_file() and (output/(relative+'.gz')).is_file():
                text=text.replace('"'+relative+'"','"'+relative+'.gz"')
        def version(match):
            asset=output/match[2]
            return match[1]+match[2]+'?v='+sha(asset)[:12]+match[3] if asset.is_file() else match[0]
        text=re.sub(r'((?:src|href)=")([^"?]+\.(?:svg|js|css))(?:\?v=[a-f0-9]+)?(")',version,text)
        page.write_text(text,encoding='utf-8')
    for record in scientific:
        if sha(output/record['path'])!=record['sha256']:
            raise ValueError('Scientific data changed during visual refresh: '+record['path'])
    changed=[*inputs,*asset_inputs,root/'analyses/atlas/rnaseq_atlas_template.html',root/'analyses/atlas/chipseq_portal/chipseq.html',root/'scripts/atlas/build_rnaseq_atlas.py',Path(__file__).resolve()]
    records={r['path']:r for r in manifest['inputs']}
    for path in changed:
        relative=path.resolve().relative_to(root.resolve()).as_posix();records[relative]={'path':relative,'sha256':sha(path)}
    manifest['inputs']=sorted(records.values(),key=lambda r:r['path'])
    size=json.loads((output/'size_report.json').read_text())
    shared=['gene_collections.js','gene_collections.css','data/gene_collections.json','chip_data.js','atlas_design.css','atlas_theme.js']
    fonts=sum(p.stat().st_size for p in (output/'fonts').glob('*.woff2'))
    size['shared_collection_asset_bytes']=sum((output/name).stat().st_size for name in shared)+fonts
    size['chip_explorer_initial_asset_bytes']=sum((output/name).stat().st_size for name in ['chipseq.html','chipseq.css','chipseq_logic.js','chipseq_tracks.js','chipseq.js','data/chipseq/explorer.json',*shared])+fonts
    size['bundle_content_bytes']=sum(p.stat().st_size for p in output.rglob('*') if p.is_file() and p.name not in {'manifest.json','size_report.json'})
    (output/'size_report.json').write_text(json.dumps(size,indent=2)+'\n',encoding='utf-8')
    manifest['outputs']=[{'path':p.relative_to(output).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(output.rglob('*')) if p.is_file() and p.name!='manifest.json']
    (output/'manifest.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({'status':'PASS','unchanged_scientific_files':len(scientific),'pages':len(list(output.glob('*.html')))-1,'presentation_assets_bytes':(output/'atlas_design.css').stat().st_size+fonts},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path('.'))
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    refresh(args.root.resolve(),args.output.resolve())
