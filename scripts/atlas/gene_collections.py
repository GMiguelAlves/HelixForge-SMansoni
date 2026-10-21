"""Small shared identity catalog, independent of expression matrices."""
import json
import re
import shutil


def build_gene_collections(root, output, annotations, rna_genes, groups, resolutions, reference_id):
    aliases = {}
    for gene, row in annotations.items():
        for alias in re.split(r"[;,\s]+", row.get("previous_stable_id", "")):
            if alias:
                aliases.setdefault(alias.lower(), set()).add(gene)
    for alias, gene in resolutions.items():
        if gene in annotations:
            aliases.setdefault(alias.lower(), set()).add(gene)
    payload = {
        "schema_version": 1, "reference_id": reference_id,
        "genes": [[gene, gene in rna_genes] for gene in sorted(annotations)],
        "aliases": {alias: sorted(ids) for alias, ids in sorted(aliases.items())},
        "presets": {"machinery": {"name": "Epigenetic machinery", "ids": sorted({g for ids in groups.values() for g in ids})}},
        "groups": {name: ids for name, ids in groups.items()},
    }
    (output / "data/gene_collections.json").write_text(json.dumps(payload, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")
    sources = [root / "scripts/atlas/gene_collections.py"]
    for name in ("gene_collections.js", "gene_collections.css", "rna_collections.js", "expression_chromatin.js", "chip_data.js"):
        source = root / "analyses/atlas/collections" / name
        shutil.copy2(source, output / name)
        sources.append(source)
    return sources
