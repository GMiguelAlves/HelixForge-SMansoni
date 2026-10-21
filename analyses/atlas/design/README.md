# Atlas visual system

Audience: researchers consulting S. mansoni expression, chromatin regions, gene collections and original study evidence. The first action should be a gene/locus search or selection of a biological context. Scientific values, scales, inference boundaries and downloadable results remain unchanged.

## Plan and tokens

- Surface `#F5F8FA`: cool laboratory canvas; plots remain white.
- Ink `#203B42`: body, gene identifiers and headings.
- Rule `#D7E1E5`: separates navigation, records and plotting surfaces.
- RNA `#265D8A`: RNA navigation and actions.
- Chromatin `#14776C`: ChIP navigation and actions.
- Review `#9B661F`: limitations and QC review, never decoration.
- IBM Plex Sans, locally bundled under its upstream license, for interface and headings. Tabular numerals for counts and coordinates. Existing scientific SVG encodings remain intact.
- Type scale: 12 / 14 / 16 / 20 / 26 / 36 px; restrained weights, sentence case. Text blocks generally <=76 characters. Controls at least 38 px tall; information buttons retain a compact circular shape.

The distinctive element is a paired assay navigation: RNA and ChIP share one reference, one gene collection and one interface, with assay color restricted to navigation and actions. Quantitative color scales retain their scientific meaning.

```
HelixForge atlas    Schistosoma mansoni        RNA-seq | ChIP-seq | Collections
RNA expression                              [ Gene / description ] [Find gene]
Life-cycle expression                       12 studies / samples / coverage
Contexts | Genes | Studies | Collections | Sample map | Literature | Methods

Biological contexts / chosen collection
[ study / view / support controls                     ] [export]
[               white scientific plotting surface             ]
```

ChIP uses the same frame; the region browser comes before the collection builder, with filters on the left and the locus plot plus source evidence on the right. Studies share that frame and use a compact facts strip, quiet sections for design/QC and disclosures for detailed records.

## Review before implementation

An initial idea of two large assay cards would repeat the existing landing-page treatment and postpone the search. Replace them with the paired top navigation and a compact assay switch. Avoid a decorative helix or invented track: the real heatmaps and locus browser provide the atlas's visual identity. Remove gradients, floating-card shadows and redundant uppercase labels. Use borders only to distinguish controls, selectable items, plot boundaries and evidence sections. Preserve existing small-screen behavior; this delivery focuses on desktop, as requested earlier.

## Implementation rules

`atlas_design.css` is loaded after legacy layout styles on every entry and study page. Its scope uses explicit atlas classes/attributes and avoids changing SVG scientific encodings. Shared assets are copied and hashed by `atlas_design.py` during builds. No external font request is required at runtime. CSS and navigation presentation may change; result tables, original sources and biological values may not.

## Delivery review — 2026-10-21

Applied to both entry points, twelve RNA study pages and two ChIP study pages. RNA starts with biological contexts; ChIP starts with its region browser. Search placement, local/global navigation, controls, information icons, statistics strips, table surfaces and download actions share the same tokens. The gene view retains expression, variability and original contrasts. Collection views retain the original clustered matrices and source evidence.

Visual critique removed the duplicate assay landing cards, heavy gradients/shadows, repeated study title and an overly strong secondary-button treatment. Native plots supply the visual identity. The ChIP context selector was widened to keep the study identifier readable. No mobile redesign was performed.

The complete rebuild ran out of temporary disk space before publication; its partial output was discarded. Presentation was then refreshed against the existing validated bundle with `scripts/atlas/refresh_atlas_design.py`. All 410 data and scientific figure files, including the full browser payload, were verified byte-identical. The full bundle validator, five JavaScript suites and the final 16-page landmark/link/cache-version check passed. A final browser-logic run covered the updated study heading. The shared build hook was also added for future complete builds; that full build still needs sufficient temporary space.

Manual desktop checks covered stage filtering (cercaria: 19 samples, five conditions, two studies), expandable information, Smp_244880 expression/contrasts, PRJEB32839 study design, the machinery integration with H3K4me3/female_cercaria and the HP1 locus browser. Final localhost RNA and ChIP pages showed no console errors. Scientific color scales, assignment/QC states and inference limits were preserved.

The shared CSS plus all three local font weights add approximately 217 kB, cached across pages. Font source, license and exact upstream revision are recorded in `fonts/provenance.json` and `fonts/LICENSE.txt`.

Docker packaging and remote publication remain deferred to the user's next instruction.

## Theme preference — 2026-10-21

All 16 pages now expose Theme: System / Light / Dark in the shared masthead. `atlas_theme.js` resolves the preference in the document head before the interface renders. The browser-local `hf-atlas-theme-v1` preference follows navigation and synchronizes across tabs; System follows OS changes. Unavailable storage falls back gracefully for the current page. Print styles remain light. Native scientific plot canvases, color scales and downloaded figures retain their original appearance; no inversion or recoloring filter is applied.

The dark palette uses a blue-green canvas, distinct panel/field surfaces, pale RNA and ChIP action colors and a separate amber review surface. No new framework or external request is needed. Theme-state tests cover initial system resolution, manual overrides, navigation persistence, cross-tab events, invalid stored values and blocked storage.
