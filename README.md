# CUT&Tag Pipeline

A Snakemake workflow for CUT&Tag data analysis, designed to run on institutional SLURM clusters. Uses environment modules (not conda) for tool management.

## Workflow organization

**We use a per-project layout:** clone or copy this pipeline into each project directory. The Snakefile and config live alongside your data.

```bash
git clone https://github.com/josiahmurra/cut_tag_pipeline.git
cd cut_tag_pipeline
```

```
my_project/
├── Snakefile
├── config/
├── profiles/
├── fastq/
└── ...
```

This keeps each project self-contained and makes it easy to customize config or lock to a specific pipeline version. The Snakefile stays compact (~50–250 lines), so duplicating it per project is straightforward.

**Alternative (not used here):** You could keep the pipeline in a central folder and point Snakemake at project directories with `--directory`. We prefer the per-project approach for simplicity.

---

## Prerequisites

- Access to a cluster with SLURM
- Required environment modules available: `cutadapt`, `bowtie2`, `samtools`, etc. (see `config/config.yaml` for the full list as rules are added)
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
cut_tag_pipeline/
├── Snakefile              # Workflow definition
├── config/
│   └── config.yaml        # Paths, tool settings, sample naming
├── profiles/
│   └── slurm/
│       └── config.yaml    # SLURM settings (account, partition, resources)
├── fastq/                 # Raw FASTQ files go here
├── fastq_trim/            # Trimmed output (created by pipeline)
├── job_reports/           # SLURM logs (created by pipeline)
└── rmarkdowns/            # R Markdown documentation
```

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

- `modules`: Environment module names (e.g. `cutadapt/4.0`) — adjust if your cluster uses different versions
- `fastq_suffix_r1`, `fastq_suffix_r2`: Change if using a different naming convention (default: `.1.fq.gz`, `.2.fq.gz`)

## Running the pipeline

### Dry run (no jobs submitted)

Check what would run without executing:

```bash
cd /path/to/cut_tag_pipeline
snakemake -n --profile profiles/slurm
```

### Full run (submits jobs to SLURM)

```bash
cd /path/to/cut_tag_pipeline
micromamba activate snakemake
snakemake --profile profiles/slurm
```

Snakemake will submit one SLURM job per sample per rule. Jobs run on the cluster; Snakemake monitors them from the login node until completion.

### Local run (for testing, no SLURM)

To run on your current machine without SLURM (useful for small test datasets):

```bash
snakemake --cores 4
```

Note: Use a local profile or omit the SLURM profile. The default SLURM profile submits to the cluster.

## Useful commands

| Command | Purpose |
|---------|---------|
| `snakemake -n` | Dry run — show what would be executed |
| `snakemake --dag \| dot -Tpng -o dag.png` | Generate DAG visualization |
| `snakemake -n --summary` | Summary of jobs to run |
| `snakemake --unlock` | Unlock if a previous run was interrupted |

## Troubleshooting

- **"module: command not found"** — Ensure your SLURM job script sources the module system. The profile uses `use_envmodules: true`; the cluster's default job script may need to initialize modules.
- **"No samples found"** — Add FASTQ files to `fastq/` with the correct naming (`*.1.fq.gz`, `*.2.fq.gz`), or specify `samples` in `config/config.yaml`.
- **SLURM plugin not found** — Activate the snakemake environment: `micromamba activate snakemake`
