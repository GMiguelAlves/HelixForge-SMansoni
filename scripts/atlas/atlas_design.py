"""Shared presentation for both assays; no scientific data transformations."""
import hashlib
import re
import shutil


def masthead(assay):
    rna = ' aria-current="page"' if assay == 'rna' else ''
    chip = ' aria-current="page"' if assay == 'chip' else ''
    collection = 'chipseq.html#collection' if assay == 'chip' else 'atlas.html#collection'
    return f'''<nav class="atlas-masthead" aria-label="Atlas portals"><a class="atlas-brand" href="atlas.html"><svg width="28" height="26" viewBox="0 0 28 26" aria-hidden="true"><path d="M2 8h24M2 18h24" stroke="#b9ccd5" stroke-width="2"/><path d="M5 4h6v8H5zM17 4h6v8h-6z" fill="#265d8a"/><path d="M8 14h11v8H8z" fill="#14776c"/></svg><span>HelixForge atlas</span><em>Schistosoma mansoni</em></a><div class="atlas-assays"><a class="rna" href="atlas.html"{rna}>RNA-seq</a><a class="chip" href="chipseq.html"{chip}>ChIP-seq</a><a href="{collection}">Gene collections</a><label class="atlas-theme-control">Theme<select id="atlas-theme" aria-label="Color theme"><option value="system">System</option><option value="light">Light</option><option value="dark">Dark</option></select></label></div></nav>'''


def decorate(report, assay, study=None, css_version=''):
    if 'class="atlas-masthead"' in report:
        report = re.sub(r'<nav class="atlas-masthead"[\s\S]*?</nav>', lambda m:masthead(assay), report, count=1)
        return re.sub(r'href="atlas_design\.css(?:\?v=[a-f0-9]+)?"', 'href="atlas_design.css'+('?v='+css_version if css_version else '')+'"', report)
    # Older compact ChIP study exports rely on the browser's implicit head/body.
    # Give these pages explicit landmarks before adding the shared frame.
    if '<head>' not in report and '</head>' not in report:
        report = report.replace('<html lang="en">', '<html lang="en"><head>', 1)
        report = report.replace('</style>', '</style></head><body>', 1)
        report = report.replace('</html>', '</body></html>')
    def body(match):
        attrs = match[1]
        if 'class="' in attrs:
            attrs = re.sub(r'class="([^"]*)"', lambda m: 'class="'+m[1]+' hf-atlas'+(' chip-study' if assay=='chip' and study else '')+'"', attrs, count=1)
        else:
            attrs += ' class="hf-atlas'+(' chip-study' if assay=='chip' and study else '')+'"'
        return '<body'+attrs+f' data-assay="{assay}">'
    report = re.sub(r'<body([^>]*)>', body, report, count=1)
    link = '<link rel="stylesheet" href="atlas_design.css'+('?v='+css_version if css_version else '')+'">'
    report = report.replace('</head>', link+'\n</head>', 1)
    # The old assay-card landing and ChIP brand bar duplicate the common navigation.
    report = re.sub(r'\s*<nav class="mode-nav"[\s\S]*?</nav>', '', report, count=1)
    report = re.sub(r'\s*<section id="portal-choices"[\s\S]*?</section>', '', report, count=1)
    report = re.sub(r'<nav(?![^>]*class=)', '<nav class="atlas-local-nav"', report, count=1)
    if assay == 'chip' and not study:
        report = report.replace('</header>', '</header><nav class="atlas-local-nav" aria-label="ChIP-seq sections"><a href="#regions">Regions</a><a href="#collection">Gene collections</a><a href="#studies">Studies</a></nav>', 1)
    if study and assay == 'rna':
        report = re.sub(r'(<header>[\s\S]*?<h1>)[\s\S]*?(</h1>)', lambda m:m[1]+study+m[2], report, count=1)
    report = report.replace('HELIXFORGE · CHIP-SEQ STUDY', 'ChIP-seq study').replace('HELIXFORGE · PUBLISHED CHIP-SEQ CATALOGS', 'ChIP-seq study')
    report = re.sub(r'(<body[^>]*>)', lambda m:m[1]+'\n'+masthead(assay), report, count=1)
    # Link labels use native affordances; arrows are reserved for browser panning.
    report = re.sub(r'(<a\b[^>]*>)([\s\S]*?)(</a>)', lambda m:m[1]+re.sub(r'\s*[→↗]\s*$', '', m[2])+m[3], report)
    return report


def apply_atlas_design(root, output):
    source = root / 'analyses/atlas/design'
    files = [source / 'atlas_design.css', source / 'atlas_theme.js', *sorted((source/'fonts').glob('*'))]
    for path in files:
        relative = path.relative_to(source)
        target = output / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
    version = hashlib.sha256((output/'atlas_design.css').read_bytes()).hexdigest()[:12]
    for page in output.glob('*.html'):
        if page.name == 'machinery.html':
            continue  # compatibility redirect, not an independent assay interface
        assay = 'chip' if page.stem.startswith('chipseq') else 'rna'
        study = page.stem.split('_',1)[1] if page.stem.startswith(('study_','chipseq_')) else None
        text = decorate(page.read_text(encoding='utf-8'), assay, study, version)
        text = re.sub(r'<script src="atlas_theme\.js(?:\?v=[a-f0-9]+)?"></script>\s*', '', text)
        theme_version = hashlib.sha256((output/'atlas_theme.js').read_bytes()).hexdigest()[:12]
        text = text.replace('<head>', '<head>\n<script src="atlas_theme.js?v='+theme_version+'"></script>', 1)
        page.write_text(text, encoding='utf-8')
    return [root/'scripts/atlas/atlas_design.py', *files]
