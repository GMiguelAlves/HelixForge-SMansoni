# Configuration

`reference/reference_manifest.json` freezes reference identity and download
sources. Each study directory contains its HelixForge compatibility adapter,
user settings, reviewed DE specification (or explicit block), and the original
curation README. `server.env.template` and `slurm.config` are site templates;
the populated environment file is ignored by Git.

Pending P1 studies freeze `HF_RNASEQ_RUN_MODE` in their tracked settings.
Inferential studies use `full`; unreplicated/descriptive studies use `import`
and carry `de_spec.blocked.json`. P1 uses the private, ignored `server.env`
with the HelixForge commit pinned in the local execution plan. That commit
includes the native single-end path and its reduced real-tool Slurm validation.
The historical `server.env.template` remains pinned to v1.0.2 for
reproducibility. This evidence does not replace the biological single-end gate
for `PRJNA312093`.
