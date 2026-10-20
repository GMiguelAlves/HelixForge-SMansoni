# Configuration

`reference/reference_manifest.json` freezes reference identity and download
sources. Each RNA-seq study directory under `rnaseq/` contains its HelixForge compatibility adapter,
user settings, reviewed DE specification (or explicit block), and the original
curation README. `server.env.template`, `server.p1.env.template` and
`slurm.config` are site templates; the populated environment file is ignored
by Git.

ChIP-seq-specific configuration remains under `chipseq/`. The external
`HF_RESULTS_ROOT` used by Slurm is independent of this repository layout.

Pending P1 studies freeze `HF_RNASEQ_RUN_MODE` in their tracked settings.
Inferential studies use `full`; unreplicated/descriptive studies use `import`
and carry `de_spec.blocked.json`. P1 copies its dedicated template into the
private, ignored `server.env` and verifies the HelixForge `v1.1.0` tag and
commit against the local execution plan. The release includes the native
single-end path and its reduced real-tool Slurm validation.
At the validated site, P1 uses existing Conda environments through
`HF_CERTIFIED_RUNTIME_PATH` with the `slurm` profile; it does not ask Nextflow
to solve new Conda environments or use Apptainer.
The historical `server.env.template` remains pinned to v1.0.2 for
reproducibility. This evidence does not replace the biological single-end gate
for `PRJNA312093`.
