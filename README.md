# CUT&Tag Pipeline

A Snakemake workflow for CUT&Tag data analysis, designed to run on institutional SLURM clusters. Uses environment modules (not conda) for tool management.

## Workflow organization

**Typical usage:** Keep the pipeline in one place. Each project is a separate directory containing `config/`, `fastq/`, and (optionally) project-specific config overrides. Run Snakemake by pointing `--directory` at the project. You do **not** need to clone the repo for every run — clone once (or copy) and reuse.

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

You run from the pipeline directory and pass the project path: `sbatch scripts/run_snakemake_slurm.sh /path/to/my_project`.


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
cut_tag_pipeline/          # Pipeline (central; clone once)
├── Snakefile              # Workflow definition
├── config/
│   └── config.yaml        # Default paths, tool settings
├── profiles/
│   └── slurm/
│       └── config.yaml    # SLURM settings (account, partition, resources)
├── scripts/
│   └── run_snakemake_slurm.sh
└── rmarkdowns/

project_dir/               # Project (--directory points here)
├── config/
│   └── config.yaml        # Project-specific overrides (optional)
├── fastq/                 # Raw FASTQ files
├── fastq_trim/            # Trimmed output (created)
├── alignment/             # BAMs, etc. (created)
└── job_reports/           # SLURM logs (created)
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

1. Edit `scripts/run_snakemake_slurm.sh` — set `--mail-user` to your email (one-time)
2. From the pipeline directory:
   ```bash
   cd /path/to/cut_tag_pipeline
   sbatch scripts/run_snakemake_slurm.sh /path/to/project
   ```

You can disconnect; you get emails when the job starts, completes, or fails. Per-project layout: use `.` as the project path.

**Optional dry run** (check what would run):
```bash
cd /path/to/cut_tag_pipeline
micromamba activate snakemake
snakemake -n --profile profiles/slurm --directory /path/to/project
```

## Useful commands

| Command | Purpose |
|---------|---------|
| `snakemake -n --directory /path/to/project` | Dry run |
| `snakemake --unlock --directory /path/to/project` | Unlock if a run was interrupted |

## Troubleshooting

- **"module: command not found"** — Ensure your SLURM job script sources the module system. The profile uses `use_envmodules: true`; the cluster's default job script may need to initialize modules.
- **"No samples found"** — Add FASTQ files to `fastq/` with the correct naming (`*.1.fq.gz`, `*.2.fq.gz`), or specify `samples` in `config/config.yaml`.
- **SLURM plugin not found** — Activate the snakemake environment: `micromamba activate snakemake`
