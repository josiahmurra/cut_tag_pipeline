"""
CUT&Tag pipeline - Snakemake workflow
Uses environment modules (--use-envmodules) for tool management.
"""

configfile: "config/config.yaml"

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


rule all:
    """Build all trimmed FASTQ files and aligned BAMs."""
    input:
        TRIMMED,
        ALIGNED,


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
    """Align trimmed reads with Bowtie2 and convert to BAM."""
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
        | samtools view -@ {threads} -bS -o {output.bam} -
        """
