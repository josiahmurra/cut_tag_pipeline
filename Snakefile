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


rule all:
    """Build all trimmed FASTQ files."""
    input:
        TRIMMED


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
