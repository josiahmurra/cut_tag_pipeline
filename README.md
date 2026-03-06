# CUT&Tag Pipeline

A Snakemake workflow for CUT&Tag data analysis, designed to run on institutional SLURM clusters. Uses environment modules (not conda) for tool management.

## Workflow overview

The pipeline runs the following steps (in dependency order):

1. **Trim** — Cutadapt on raw paired-end FASTQs (Nextera adapters).
2. **Align** — Bowtie2, then samtools sort (coordinate-sorted BAM).
3. **Picard MarkDuplicates** — Deduplicated BAM + metrics (for QC).
4. **QC** — FastQC on raw and trimmed FASTQs; MultiQC aggregates those reports. R scripts produce figures for alignment stats (alignment rate, mapped fragments, sequencing depth, duplication rate, estimated library size).
5. **Fragment BEDs and correlation** — MACS3 filterdup: aligned BAM → BEDPE (one row per fragment). BEDPE → 500 bp binned fragment count BEDs (bin500). R script merges bin500 BEDs and produces Pearson correlation on log2(counts), heatmap + dendrogram PDF (requires ≥2 samples).
6. **MACS3 peak calling** — Target vs control (IgG) pairing; always runs cutoff analysis. Declared output is the cutoff table (for elbow plots). Uses aligned BAMs (with duplicates) and `--keep-dup auto`. Broad peaks for marks listed in `macs3_broad_marks` (e.g. H3K27me3); others get narrow peaks.

## Workflow organization

**Typical usage:** Keep the pipeline in one place. Each project is a separate directory containing `config/`, `fastq/`, and (optionally) project-specific config. You submit from the pipeline directory with `PROJECT_DIR=/path/to/project`; the job runs Snakemake from the project directory and uses the pipeline Snakefile via `-s`, so you do **not** need to clone the repo per project — clone once and reuse.

```bash
# One-time setup: clone the pipeline
git clone https://github.com/josiahmurra/cut_tag_pipeline.git
```

Typical layout:

```
cut_tag_pipeline/          # Pipeline (clone once, keep central)
├── Snakefile
├── config/
├── profiles/
├── scripts/
└── ...

my_project/                # Project directory (data + config)
├── config/
│   └── config.yaml        # Project-specific settings
├── fastq/
└── ...
```

You run from the pipeline directory: set `PROJECT_DIR` and **run the script** (do not call `sbatch` yourself). The script will submit itself so that job logs go to the project directory:

```bash
cd /path/to/cut_tag_pipeline
PROJECT_DIR=/path/to/my_project ./run_snakemake_slurm.sh
```

**Important:** Use `./run_snakemake_slurm.sh`, not `sbatch run_snakemake_slurm.sh`. If you use `sbatch` directly, the script cannot set the log path and logs will go to the pipeline directory.


---

## Prerequisites

- Access to a cluster with SLURM
- Required environment modules available: `cutadapt`, `bowtie2`, `samtools`, `R` (see `config/config.yaml` for the full list and versions)
- Micromamba installed (or use conda/mamba if you prefer)

## Install Micromamba (if not already installed)

If you don't have micromamba, install it first:

```bash
# Linux/macOS
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

Then either restart your shell or source your config file (e.g. `~/.bashrc`).

## Create Snakemake environment with micromamba

Create an environment with Snakemake 8+ and the SLURM executor plugin:

```bash
micromamba create -c conda-forge -c bioconda -n snakemake snakemake snakemake-executor-plugin-slurm
```

Activate the environment before running the pipeline:

```bash
micromamba activate snakemake
```

## Pipeline structure

```
cut_tag_pipeline/          # Pipeline (central; clone once)
├── Snakefile              # Workflow definition
├── config/
│   └── config.yaml        # Default paths, tool settings
├── profiles/
│   └── slurm/
│       └── config.yaml    # SLURM settings (account, partition, resources)
├── scripts/
│   ├── alignment_qc_plots.R   # Alignment rate, mapped fragments, depth
│   ├── picard_dup_plots.R     # Duplication rate, estimated library size
│   ├── bin500_correlation.R   # Correlation heatmap + dendrogram from bin500 BEDs
│   └── make_chr1_test_fastq.sh  # Optional: create chr1-only test FASTQs from BAM
└── rmarkdowns/            # Reference Rmd workflows

project_dir/               # Project (PROJECT_DIR points here)
├── config/
│   └── config.yaml        # Project config (paths, modules, samples, MACS3 peak options)
├── fastq/                 # Raw FASTQ files ({sample}.1.fq.gz, {sample}.2.fq.gz)
├── fastq_trim/            # Trimmed FASTQs (created)
├── alignment/
│   ├── bam/               # BAMs, bowtie2 summaries, picard_summary (created)
│   └── bed/               # BEDPE, bin500 BEDs (created)
├── peakCalling/
│   └── MACS3/             # Cutoff tables, narrow/broad peaks (created)
├── qc/                    # SampleList, summary_paths, alignSummary.rds, alignDupSummary.rds (created)
├── figures/               # QC plots + correlation PDF (created)
└── job_reports/           # All rule logs: SLURM job logs + alignment_qc.log, picard_qc.log (created)
```

All rule logs live under `job_reports/`: submitted jobs (trim, align, FastQC, MultiQC, Picard, MACS3, etc.) get their own logs from SLURM, and the local rules (alignment_qc, picard_qc) write `job_reports/alignment_qc.log` and `job_reports/picard_qc.log`.

## FASTQ naming convention

The pipeline uses a dot (`.`) to separate the sample ID from the file suffix. This allows underscores in sample names for fields like genotype, target, and replicate.

**Raw input files** (place in `fastq/`):

- `{sample}.1.fq.gz` and `{sample}.2.fq.gz`
- Example: `WT_CTCF_Rep1.1.fq.gz`, `WT_CTCF_Rep1.2.fq.gz`

**Trimmed output** (created by pipeline):

- `{sample}.cut.1.fq.gz` and `{sample}.cut.2.fq.gz`

**Optional:** If you prefer not to rely on auto-discovery, list samples explicitly in `config/config.yaml`:

```yaml
samples: ["WT_CTCF_Rep1", "KO_H3K27me3_Rep2"]
```

## Configuration

### SLURM profile (`profiles/slurm/config.yaml`)

Edit this file if your cluster uses different defaults.

**YAML → SLURM (`sbatch`) mapping:**

| YAML key | SLURM equivalent | Description |
|----------|------------------|-------------|
| `slurm_account` | `#SBATCH --account` | PI/research account for billing |
| `slurm_partition` | `#SBATCH --partition` / `-p` | Partition (queue) name |
| `mem_mb_per_cpu` | `#SBATCH --mem-per-cpu` | Memory per CPU in MB. Use this (not `mem_mb`) if your cluster limits mem-per-cpu (e.g. 7.5 GB max). |
| `mem_mb` | `#SBATCH --mem` | Total memory per job in MB (alternative to mem_mb_per_cpu) |
| `runtime` | `#SBATCH --time` | Walltime in minutes (e.g. 1080 = 18 hours) |
| `slurm_logdir` | `-o`, `-e` | Directory for stdout/stderr logs |
| `jobs` | (Snakemake) | Max concurrent SLURM jobs, not an sbatch flag |

**CPUs per job:** Set by each rule's `threads:` in the Snakefile (maps to `#SBATCH --cpus-per-task`).

**Memory/CPU coordination:** If your cluster caps memory per CPU (e.g. 7.5 GB), use `mem_mb_per_cpu` and keep it ≤ 7500. Total memory = threads × mem_mb_per_cpu.

### Pipeline config (`config/config.yaml`)

- **Paths:** `fastq_dir`, `align_bam_dir`, `bed_dir`, `figures_dir`, `macs3_peak_dir`, etc. — adjust if your project layout differs.
- **`modules`:** Environment module names (e.g. `cutadapt/4.0`, `bowtie2/2.5.0`, `R/4.4.2`, `picard`, `macs3/3.0.2`) — set for your cluster.
- **`fastq_suffix_r1` / `fastq_suffix_r2`:** Default `.1.fq.gz` and `.2.fq.gz`; change if you use a different naming convention.
- **`samples`:** Optional list of sample IDs; if empty, samples are discovered from `fastq_dir` using the suffix pattern.

**MACS3 peak calling (target vs control):**

- **`control_pattern`:** Sample names containing this (e.g. `IgG`) are treated as controls; only non-control samples get peak calling.
- **Default control:** For each target, control is `{genotype}_IgG_{replicate}` (sample name split on `sample_delimiter`, default `_`, three parts).
- **`macs3_control_overrides`:** Optional map `target_id: control_id` when an IgG failed and you want to use a different control (e.g. another replicate). Only list exceptions.
- **`macs3_broad_marks`:** List of target substrings that get broad peak calling (e.g. `["H3K27me3"]`); others get narrow peaks.
- **`macs3_peak_dir`:** Output directory for cutoff tables and peak files (default `peakCalling/MACS3`).
- **`macs3_q_narrow` / `macs3_broad_cutoff`:** q-value for narrow peaks and cutoff for broad peaks.
- Cutoff analysis is always run; the workflow requests the cutoff table; narrow/broad peak files are written in the same run.

**Other:**

- **`fragment_bin_length`:** Bin size in bp for bin500 BEDs (default 500).
- **Bowtie2:** `bowtie2.index_base` and `bowtie2.threads`; index must be pre-built.

## Running the pipeline

1. Edit `run_snakemake_slurm.sh` — set `#SBATCH --mail-user` to your email (one-time).
2. From the **pipeline** directory, set `PROJECT_DIR` and **run the script** (do not use `sbatch`):
   ```bash
   cd /path/to/cut_tag_pipeline
   PROJECT_DIR=/path/to/project ./run_snakemake_slurm.sh
   ```
   Job logs (`cuttag_<jobid>.out`, `cuttag_<jobid>.err`) will be in the project directory. You get one email when the job starts, completes, or fails.

You can disconnect after submitting; you get emails when the job starts, completes, or fails. To run a project that lives in the same directory as the pipeline, use `PROJECT_DIR=.` when submitting.

**Optional dry run** (from the project directory, using the same invocation style as the batch job):
```bash
cd /path/to/project
micromamba activate snakemake
snakemake -s /path/to/cut_tag_pipeline/Snakefile --configfile config/config.yaml --profile /path/to/cut_tag_pipeline/profiles/slurm -n
```

## Useful commands

| Command | Purpose |
|---------|---------|
| From project dir: `snakemake -s /path/to/pipeline/Snakefile --configfile config/config.yaml --profile /path/to/pipeline/profiles/slurm -n` | Dry run |
| Same but add `--unlock` | Unlock if a run was interrupted |

**Optional:** To build chr1-only (or subset) test FASTQs from an existing BAM for quick pipeline tests, use `scripts/make_chr1_test_fastq.sh`; see the script for usage (requires a coordinate-sorted, indexed BAM).

## Troubleshooting

- **"Set PROJECT_DIR" / usage message** — Run the script with `PROJECT_DIR` set: `PROJECT_DIR=/path/to/project ./run_snakemake_slurm.sh`. Do not use `sbatch` yourself.
- **Job logs in the wrong place** — You used `sbatch run_snakemake_slurm.sh` instead of `./run_snakemake_slurm.sh`. Run the script directly so it can submit itself with the correct log path.

- **"module: command not found"** — The profile uses `use_envmodules: true`; the cluster's default job script may need to initialize the module system.
- **"No samples found"** — Add FASTQ files to the project's `fastq/` with the correct naming (`*.1.fq.gz`, `*.2.fq.gz`), or specify `samples` in the project's `config/config.yaml`.
- **SLURM plugin not found** — Activate the snakemake environment: `micromamba activate snakemake`
- **R / alignment QC or Picard QC errors** — The submit script runs Snakemake from the project directory with `-s` pointing at the pipeline Snakefile. If you run Snakemake manually, use the same pattern (run from the project dir and pass `-s /path/to/pipeline/Snakefile` and `--configfile config/config.yaml`) so relative paths resolve correctly.
- **Peak calling skipped** — Ensure you have both target samples (e.g. CTCF, H3K27me3) and matching control samples (e.g. IgG) in `fastq/`; control is inferred as `{genotype}_IgG_{replicate}`. Use `macs3_control_overrides` in config if a control failed and you want to use another replicate.
