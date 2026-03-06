"""
CUT&Tag pipeline - Snakemake workflow
Uses environment modules (--use-envmodules) for tool management.
Run from the project directory (e.g. via run_snakemake_slurm.sh with -s path/to/Snakefile).
alignment_qc runs locally and passes explicit input/output file paths to R (no cwd or directory logic).
"""

configfile: "config/config.yaml"

import os

# Run alignment_qc and picard_qc on the main Snakemake process; we pass file paths so R doesn't depend on cwd.
localrules: alignment_qc, picard_qc

# Discover samples: from config if specified, otherwise from fastq directory
def get_samples():
    if config.get("samples"):
        return config["samples"]
    # Pattern: {sample_id}.{suffix} where suffix is .1.fq.gz (from config)
    pattern = f"{config['fastq_dir']}/{{sample}}{config['fastq_suffix_r1']}"
    return sorted(glob_wildcards(pattern).sample)

SAMPLES = get_samples()

# Trim rule outputs for the final target
TRIMMED = expand(
    "{dir}/{sample}.cut.{pair}.fq.gz",
    dir=config["fastq_trim_dir"],
    sample=SAMPLES,
    pair=["1", "2"],
)

# Alignment rule outputs for the final target
ALIGNED = expand(
    "{dir}/{sample}.bowtie2.bam",
    dir=config["align_bam_dir"],
    sample=SAMPLES,
)

# Alignment QC: one rule that uses all bowtie2 summary files
ALIGNMENT_QC = [
    f"{config['alignment_qc_dir']}/SampleList.txt",
    f"{config['alignment_qc_dir']}/alignSummary.rds",
    f"{config['figures_dir']}/alignment_rate.pdf",
    f"{config['figures_dir']}/mapped_fragments.pdf",
    f"{config['figures_dir']}/sequencing_depth.pdf",
]

# FastQC and MultiQC dirs (defaults so project configs without these keys still work)
# Under qc/: fastqc/raw, fastqc/trimmed, multiqc/raw, multiqc/trimmed
fastqc_raw_dir = config.get("fastqc_raw_dir", "qc/fastqc/raw")
fastqc_trimmed_dir = config.get("fastqc_trimmed_dir", "qc/fastqc/trimmed")
multiqc_raw_dir = config.get("multiqc_raw_dir", "qc/multiqc/raw")
multiqc_trimmed_dir = config.get("multiqc_trimmed_dir", "qc/multiqc/trimmed")

# FastQC outputs (raw and trimmed) for MultiQC
FASTQC_RAW_ZIPS = (
    expand(f"{fastqc_raw_dir}/{{sample}}.1_fastqc.zip", sample=SAMPLES)
    + expand(f"{fastqc_raw_dir}/{{sample}}.2_fastqc.zip", sample=SAMPLES)
)
FASTQC_TRIMMED_ZIPS = (
    expand(f"{fastqc_trimmed_dir}/{{sample}}.cut.1_fastqc.zip", sample=SAMPLES)
    + expand(f"{fastqc_trimmed_dir}/{{sample}}.cut.2_fastqc.zip", sample=SAMPLES)
)
MULTIQC_RAW_REPORT = f"{multiqc_raw_dir}/multiqc_report.html"
MULTIQC_TRIMMED_REPORT = f"{multiqc_trimmed_dir}/multiqc_report.html"

# Picard MarkDuplicates: remove duplicates, write metrics. Naming matches rmarkdowns/2-CUTandTag-RemoveDup-Picard.Rmd
picard_metrics_dir = config.get("picard_metrics_dir", "alignment/bam/picard_summary")
PICARD_RMDUP_BAMS = expand(
    f"{config['align_bam_dir']}/{{sample}}.bowtie2.rmDup.bam",
    sample=SAMPLES,
)
PICARD_METRICS = expand(
    f"{picard_metrics_dir}/{{sample}}.picard.rmDup.txt",
    sample=SAMPLES,
)

# Picard duplicate QC: duplication rate and estimated library size plots (from Picard metrics); merged with alignment summary.
PICARD_QC = [
    f"{config['alignment_qc_dir']}/alignDupSummary.rds",
    f"{config['figures_dir']}/duplication_rate.pdf",
    f"{config['figures_dir']}/estimated_library_size.pdf",
]

# MACS3 filterdup: BAMPE → BEDPE (one row per read-pair, duplicates removed). Downstream: convert BEDPE to fragment BED for correlation/peak calling.
bed_dir = config.get("bed_dir", "alignment/bed")
MACS3_BEDPE = expand(
    f"{bed_dir}/{{sample}}.bedpe",
    sample=SAMPLES,
)

# 500 bp binned fragment counts (from MACS3 .bedpe); naming matches rmarkdowns/3-CUTandTag-Make_the_bed.Rmd for correlation.
BIN500_BED = expand(
    f"{bed_dir}/{{sample}}.bowtie2.fragmentsCount.bin500.bed",
    sample=SAMPLES,
)

figures_dir = config.get("figures_dir", "figures")
BIN500_CORRELATION_PDF = f"{figures_dir}/correlation_with_dendrogram.pdf"
BIN500_CORRELATION_TARGET = [BIN500_CORRELATION_PDF] if len(SAMPLES) >= 2 else []

# MACS3 peak calling: target vs control pairing (matches rmarkdowns/4-CUTandTag-Peak_Calling.Rmd).
# Targets = samples that are not controls (name does not contain control_pattern).
# Default control = {genotype}_IgG_{replicate}; overrides in config for failed IgG.
def is_control(sample_id):
    pat = config.get("control_pattern", "IgG")
    return pat.lower() in sample_id.lower()

def default_control_for_target(sample_id):
    d = config.get("sample_delimiter", "_")
    parts = sample_id.split(d)
    if len(parts) < 3:
        return None
    geno, target, repli = parts[0], parts[1], parts[2]
    return f"{geno}_IgG_{repli}"

def get_macs3_peak_pairs():
    overrides = config.get("macs3_control_overrides", {}) or {}
    pairs = []  # list of (target, control) tuples
    for s in SAMPLES:
        if is_control(s):
            continue
        control = overrides.get(s) or default_control_for_target(s)
        if control is None or control not in SAMPLES:
            continue
        pairs.append((s, control))
    return pairs

MACS3_PEAK_PAIRS = get_macs3_peak_pairs()
TARGETS_FOR_PEAKS = [p[0] for p in MACS3_PEAK_PAIRS]
target_to_control = dict(MACS3_PEAK_PAIRS)
macs3_peak_dir = config.get("macs3_peak_dir", "peakCalling/MACS3")
# Targets whose name contains a broad mark (e.g. H3K27me3) get --broad peak calling; others get narrow only.
TARGETS_BROAD = [t for t in TARGETS_FOR_PEAKS if any(m in t for m in config.get("macs3_broad_marks", []))]
# MACS3 writes {prefix}_peaks.narrowPeak and (if --broad) {prefix}_peaks.broadPeak; -n is the prefix only.
MACS3_NARROW_PEAKS = expand(
    f"{macs3_peak_dir}/{{target}}.macs3.nolambda_peaks.narrowPeak",
    target=TARGETS_FOR_PEAKS,
)
MACS3_BROAD_PEAKS = expand(
    f"{macs3_peak_dir}/{{target}}.macs3.nolambda_peaks.broadPeak",
    target=TARGETS_BROAD,
)
MACS3_PEAK_TARGETS = list(MACS3_NARROW_PEAKS) + list(MACS3_BROAD_PEAKS) if MACS3_PEAK_PAIRS else []
# Cutoff analysis is always run; we use it as the only declared output so the rule runs when it's missing. MACS3 also writes narrow/broad peaks (side effects).
MACS3_CUTOFF_ANALYSIS = expand(
    f"{macs3_peak_dir}/{{target}}.macs3.nolambda_cutoff_analysis.txt",
    target=TARGETS_FOR_PEAKS,
) if MACS3_PEAK_PAIRS else []


rule all:
    """Build all trimmed FASTQ files, aligned BAMs, alignment QC, FastQC/MultiQC, Picard QC, MACS3 BEDPE, bin500 BED/correlation, and MACS3 cutoff analysis (peak files are produced in the same run)."""
    input:
        TRIMMED,
        ALIGNED,
        ALIGNMENT_QC,
        FASTQC_RAW_ZIPS,
        FASTQC_TRIMMED_ZIPS,
        MULTIQC_RAW_REPORT,
        MULTIQC_TRIMMED_REPORT,
        PICARD_RMDUP_BAMS,
        PICARD_METRICS,
        PICARD_QC,
        MACS3_BEDPE,
        BIN500_BED,
        BIN500_CORRELATION_TARGET,
        MACS3_CUTOFF_ANALYSIS,


rule trim:
    """Trim adapter sequences from paired-end FASTQs using cutadapt."""
    input:
        r1=f"{config['fastq_dir']}/{{sample}}{config['fastq_suffix_r1']}",
        r2=f"{config['fastq_dir']}/{{sample}}{config['fastq_suffix_r2']}",
    output:
        r1=f"{config['fastq_trim_dir']}/{{sample}}.cut.1.fq.gz",
        r2=f"{config['fastq_trim_dir']}/{{sample}}.cut.2.fq.gz",
    log:
        f"{config['fastq_trim_dir']}/logs/{{sample}}.cutadapt.log",
    threads:
        config["cutadapt"]["threads"]
    params:
        min_length=config["cutadapt"]["min_length"],
        adapter_1=config["cutadapt"]["adapter_1"],
        adapter_2=config["cutadapt"]["adapter_2"],
    envmodules:
        config["modules"]["cutadapt"],
    shell:
        """
        cutadapt -j {threads} \
            -m {params.min_length} \
            -a {params.adapter_1} \
            -A {params.adapter_2} \
            -o {output.r1} -p {output.r2} \
            {input.r1} {input.r2} \
            > {log} 2>&1
        """


rule align:
    """Align trimmed reads with Bowtie2, convert to BAM, and coordinate-sort (for Picard downstream)."""
    input:
        r1=f"{config['fastq_trim_dir']}/{{sample}}.cut.1.fq.gz",
        r2=f"{config['fastq_trim_dir']}/{{sample}}.cut.2.fq.gz",
    output:
        bam=f"{config['align_bam_dir']}/{{sample}}.bowtie2.bam",
        log=f"{config['align_summary_dir']}/{{sample}}.bowtie2.txt",
    threads:
        config["bowtie2"]["threads"]
    resources:
        mem_mb=49000  # 7 threads × 7000 MB/CPU; overrides Snakemake's input-based estimate (~1.35 GB)
    params:
        index=config["bowtie2"]["index_base"],
    envmodules:
        config["modules"]["bowtie2"],
        config["modules"]["samtools"],
    shell:
        """
        bowtie2 \
            --local --very-sensitive \
            --no-mixed --no-discordant \
            --phred33 -I 10 -X 700 \
            -p {threads} \
            -x {params.index} \
            -1 {input.r1} -2 {input.r2} \
            2> {output.log} \
        | samtools view -@ {threads} -bS - \
        | samtools sort -@ {threads} -o {output.bam} -
        """


rule picard_markdup:
    """Run Picard MarkDuplicates on sorted BAMs; remove duplicates and write deduplicated BAM plus metrics.
    Output naming matches rmarkdowns/2-CUTandTag-RemoveDup-Picard.Rmd: .bowtie2.rmDup.bam, picard_summary/{sample}.picard.rmDup.txt
    ASSUME_SORT_ORDER=coordinate: input is coordinate-sorted (align rule); tells Picard to assume this order (see Picard MarkDuplicates docs)."""
    input:
        bam=f"{config['align_bam_dir']}/{{sample}}.bowtie2.bam",
    output:
        bam=f"{config['align_bam_dir']}/{{sample}}.bowtie2.rmDup.bam",
        metrics=f"{picard_metrics_dir}/{{sample}}.picard.rmDup.txt",
    log:
        "job_reports/picard_markdup/{{sample}}.markdup.log",
    params:
        logdir="job_reports/picard_markdup",
    threads: 4
    resources:
        mem_mb=16000,
        runtime=120,
    envmodules:
        config.get("modules", {}).get("picard", "picard"),
    shell:
        """
        mkdir -p {picard_metrics_dir} {params.logdir} && \
        java -jar $PICARD MarkDuplicates \
            I={input.bam} \
            O={output.bam} \
            M={output.metrics} \
            ASSUME_SORT_ORDER=coordinate \
            TAG_DUPLICATE_SET_MEMBERS=true \
            REMOVE_DUPLICATES=true \
            >> {log} 2>&1
        """


rule macs3_filterdup:
    """Run MACS3 filterdup on aligned BAM to produce BEDPE (one row per read-pair). MACS3 removes duplicates; use coordinate-sorted BAM from align rule."""
    input:
        bam=f"{config['align_bam_dir']}/{{sample}}.bowtie2.bam",
    output:
        bedpe=f"{bed_dir}/{{sample}}.bedpe",
    log:
        "job_reports/macs3_filterdup/{{sample}}.filterdup.log",
    params:
        logdir="job_reports/macs3_filterdup",
        gsize=config.get("macs3_gsize", "hs"),
        keep_dup=config.get("macs3_keep_dup", "auto"),
    resources:
        mem_mb=16000,
        runtime=120,
    envmodules:
        config.get("modules", {}).get("macs3", "macs3/3.0.2"),
    shell:
        """
        mkdir -p {bed_dir} {params.logdir} && \
        macs3 filterdup \
            -i {input.bam} \
            -f BAMPE \
            -o {output.bedpe} \
            -g {params.gsize} \
            --keep-dup {params.keep_dup} \
            >> {log} 2>&1
        """


rule bin500_bed:
    """Bin MACS3 fragment (bedpe) into 500 bp windows and count fragments per bin. Output naming matches rmarkdowns/3-CUTandTag-Make_the_bed.Rmd for correlation (chrom, bin, count)."""
    input:
        bedpe=f"{bed_dir}/{{sample}}.bedpe",
    output:
        bin500=f"{bed_dir}/{{sample}}.bowtie2.fragmentsCount.bin500.bed",
    params:
        bin_len=config.get("fragment_bin_length", 500),
    shell:
        """
        awk -v w={params.bin_len} '{{print $1, int(($2 + $3)/(2*w))*w + w/2}}' {input.bedpe} \
        | sort -k1,1V -k2,2n \
        | uniq -c \
        | awk -v OFS="\\t" '{{print $2, $3, $1}}' \
        | sort -k1,1V -k2,2n \
        > {output.bin500}
        """


rule bin500_correlation:
    """Compute Pearson correlation across samples using 500 bp binned fragment counts and save heatmap + dendrogram (matches rmarkdowns/3-CUTandTag-Make_the_bed.Rmd)."""
    input:
        script=workflow.basedir + "/scripts/bin500_correlation.R",
        beds=BIN500_BED,
    output:
        pdf=BIN500_CORRELATION_PDF,
    log:
        "job_reports/bin500_correlation.log",
    params:
        figures_dir=figures_dir,
    envmodules:
        config["modules"]["r"],
    shell:
        """
        mkdir -p job_reports {params.figures_dir} && \
        Rscript {input.script} {output.pdf} {input.beds} >> {log} 2>&1
        """


rule macs3_callpeak:
    """MACS3 peak calling: one rule triggered when cutoff file is missing. Always runs --cutoff-analysis; narrow/broad peaks are written by MACS3 as side effects (same dir)."""
    input:
        target_bam=lambda w: f"{config['align_bam_dir']}/{w.target}.bowtie2.bam",
        control_bam=lambda w: f"{config['align_bam_dir']}/{target_to_control[w.target]}.bowtie2.bam",
    output:
        cutoff=f"{macs3_peak_dir}/{{target}}.macs3.nolambda_cutoff_analysis.txt",
    log:
        "job_reports/macs3_callpeak/{target}.callpeak.log",
    params:
        logdir="job_reports/macs3_callpeak",
        peak_dir=macs3_peak_dir,
        gsize=config.get("macs3_gsize", "hs"),
        q_narrow=config.get("macs3_q_narrow", 0.00001),
        broad_cutoff=config.get("macs3_broad_cutoff", 0.0001),
        is_broad=lambda w: w.target in set(TARGETS_BROAD),
    resources:
        mem_mb=16000,
        runtime=90,
    envmodules:
        config.get("modules", {}).get("macs3", "macs3/3.0.2"),
    shell:
        """
        mkdir -p {params.peak_dir} {params.logdir} && \
        BROAD="" && \
        if [ "{params.is_broad}" = "True" ]; then \
          BROAD="--broad --broad-cutoff {params.broad_cutoff}"; \
        else \
          BROAD="-q {params.q_narrow}"; \
        fi && \
        macs3 callpeak \
          -t {input.target_bam} \
          -c {input.control_bam} \
          -f BAMPE $BROAD --cutoff-analysis \
          -g {params.gsize} -B --nomodel \
          --extsize 200 --shift 0 \
          --nolambda --scale-to large \
          --keep-dup auto \
          -n {params.peak_dir}/{wildcards.target}.macs3.nolambda \
          >> {log} 2>&1
        """


rule fastqc_raw:
    """Run FastQC on raw (pre-trim) paired-end FASTQs; one job per sample."""
    input:
        r1=f"{config['fastq_dir']}/{{sample}}{config['fastq_suffix_r1']}",
        r2=f"{config['fastq_dir']}/{{sample}}{config['fastq_suffix_r2']}",
    output:
        r1_zip=f"{fastqc_raw_dir}/{{sample}}.1_fastqc.zip",
        r2_zip=f"{fastqc_raw_dir}/{{sample}}.2_fastqc.zip",
    log:
        f"{fastqc_raw_dir}/logs/{{sample}}.fastqc.log",
    params:
        outdir=fastqc_raw_dir,
    threads: 1
    resources:
        mem_mb=2048,  # FastQC needs very little; 2 GB is plenty per sample
        runtime=30,   # minutes; FastQC is quick per sample
    envmodules:
        config.get("modules", {}).get("fastqc", "fastqc"),
    shell:
        """
        mkdir -p {params.outdir}/logs && \
        fastqc -o {params.outdir} -f fastq {input.r1} {input.r2} >> {log} 2>&1
        """


rule fastqc_trimmed:
    """Run FastQC on trimmed paired-end FASTQs; one job per sample."""
    input:
        r1=f"{config['fastq_trim_dir']}/{{sample}}.cut.1.fq.gz",
        r2=f"{config['fastq_trim_dir']}/{{sample}}.cut.2.fq.gz",
    output:
        r1_zip=f"{fastqc_trimmed_dir}/{{sample}}.cut.1_fastqc.zip",
        r2_zip=f"{fastqc_trimmed_dir}/{{sample}}.cut.2_fastqc.zip",
    log:
        f"{fastqc_trimmed_dir}/logs/{{sample}}.fastqc.log",
    params:
        outdir=fastqc_trimmed_dir,
    threads: 1
    resources:
        mem_mb=2048,  # FastQC needs very little; 2 GB is plenty per sample
        runtime=30,   # minutes; FastQC is quick per sample
    envmodules:
        config.get("modules", {}).get("fastqc", "fastqc"),
    shell:
        """
        mkdir -p {params.outdir}/logs && \
        fastqc -o {params.outdir} -f fastq {input.r1} {input.r2} >> {log} 2>&1
        """


rule multiqc_raw:
    """Aggregate FastQC results for raw FASTQs into one MultiQC report."""
    input:
        FASTQC_RAW_ZIPS,
    output:
        report=MULTIQC_RAW_REPORT,
    log:
        f"{multiqc_raw_dir}/multiqc.log",
    params:
        outdir=multiqc_raw_dir,
    threads: 1
    resources:
        mem_mb=2048,  # MultiQC is lightweight; 2 GB is plenty
        runtime=30,   # minutes; MultiQC aggregates quickly
    envmodules:
        config.get("modules", {}).get("multiqc", "multiqc"),
    shell:
        """
        mkdir -p {params.outdir} && \
        multiqc --force {fastqc_raw_dir} -o {params.outdir} >> {log} 2>&1
        """


rule multiqc_trimmed:
    """Aggregate FastQC results for trimmed FASTQs into one MultiQC report."""
    input:
        FASTQC_TRIMMED_ZIPS,
    output:
        report=MULTIQC_TRIMMED_REPORT,
    log:
        f"{multiqc_trimmed_dir}/multiqc.log",
    params:
        outdir=multiqc_trimmed_dir,
    threads: 1
    resources:
        mem_mb=2048,  # MultiQC is lightweight; 2 GB is plenty
        runtime=30,   # minutes; MultiQC aggregates quickly
    envmodules:
        config.get("modules", {}).get("multiqc", "multiqc"),
    shell:
        """
        mkdir -p {params.outdir} && \
        multiqc --force {fastqc_trimmed_dir} -o {params.outdir} >> {log} 2>&1
        """


rule alignment_qc:
    """Plot alignment QC (rate, mapped fragments, depth) from all bowtie2 summary files.
    Summary manifest is written with Python so paths are absolute; R runs via shell so envmodules (R) is loaded."""
    input:
        script=workflow.basedir + "/scripts/alignment_qc_plots.R",
        summary_files=expand(
            "{dir}/{sample}.bowtie2.txt",
            dir=config["align_summary_dir"],
            sample=SAMPLES,
        ),
    output:
        sample_list=f"{config['alignment_qc_dir']}/SampleList.txt",
        summary_manifest=f"{config['alignment_qc_dir']}/summary_paths.txt",
        rds=f"{config['alignment_qc_dir']}/alignSummary.rds",
        fig_align_rate=f"{config['figures_dir']}/alignment_rate.pdf",
        fig_mapped=f"{config['figures_dir']}/mapped_fragments.pdf",
        fig_depth=f"{config['figures_dir']}/sequencing_depth.pdf",
    log:
        "job_reports/alignment_qc.log",
    params:
        alignment_qc_dir=config["alignment_qc_dir"],
        figures_dir=config["figures_dir"],
    envmodules:
        config["modules"]["r"],
    shell:
        """
        mkdir -p job_reports {params.alignment_qc_dir} {params.figures_dir} && \
        for f in {input.summary_files}; do basename "$$f" .bowtie2.txt; done | sort -u > {output.sample_list} && \
        python3 -c "import os, sys; [print(os.path.abspath(p)) for p in sys.argv[1:]]" {input.summary_files} > {output.summary_manifest} && \
        Rscript {input.script} {output.sample_list} {output.rds} {output.fig_align_rate} {output.fig_mapped} {output.fig_depth} {output.summary_manifest} >> {log} 2>&1
        """


rule picard_qc:
    """Plot duplication rate and estimated library size from Picard MarkDuplicates metrics (same pattern as alignment_qc)."""
    input:
        script=workflow.basedir + "/scripts/picard_dup_plots.R",
        metrics_files=PICARD_METRICS,
        align_summary_rds=f"{config['alignment_qc_dir']}/alignSummary.rds",
    output:
        sample_list=f"{config['alignment_qc_dir']}/picard_sample_list.txt",
        metrics_manifest=f"{config['alignment_qc_dir']}/picard_metrics_paths.txt",
        rds=f"{config['alignment_qc_dir']}/alignDupSummary.rds",
        fig_dup_rate=f"{config['figures_dir']}/duplication_rate.pdf",
        fig_lib_size=f"{config['figures_dir']}/estimated_library_size.pdf",
    log:
        "job_reports/picard_qc.log",
    params:
        alignment_qc_dir=config["alignment_qc_dir"],
        figures_dir=config["figures_dir"],
    envmodules:
        config["modules"]["r"],
    shell:
        """
        mkdir -p job_reports {params.alignment_qc_dir} {params.figures_dir} && \
        for f in {input.metrics_files}; do basename "$$f" .picard.rmDup.txt; done > {output.sample_list} && \
        python3 -c "import os, sys; [print(os.path.abspath(p)) for p in sys.argv[1:]]" {input.metrics_files} > {output.metrics_manifest} && \
        Rscript {input.script} {output.sample_list} {input.align_summary_rds} {output.rds} {output.fig_dup_rate} {output.fig_lib_size} {output.metrics_manifest} >> {log} 2>&1
        """
