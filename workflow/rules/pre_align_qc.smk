rule fastqc_se:
    input:
        lambda wildcards: samplesheet[wildcards.se_sample]["r1"],
    output:
        html="results/qc/fastqc/{se_sample}_fastqc.html",
        zip="results/qc/fastqc/{se_sample}_fastqc.zip",
    log:
        "logs/qc/fastqc/{se_sample}.log",
    threads: config_threads
    wrapper:
        "v1.24.0/bio/fastqc"


rule fastqc_pe:
    input:
        lambda wildcards: samplesheet[wildcards.pe_sample][wildcards.group],
    output:
        html="results/qc/fastqc/{pe_sample}_{group}_fastqc.html",
        zip="results/qc/fastqc/{pe_sample}_{group}_fastqc.zip",
    log:
        "logs/qc/fastqc/{pe_sample}_{group}.log",
    threads: config_threads
    wrapper:
        "v1.24.0/bio/fastqc"


rule trimgalore_se:
    input:
        get_fastq_se,
    output:
        "results/qc/trimgalore/{sample}_trimmed.fq.gz",
    conda:
        "../envs/trimgalore.yaml"
    log:
        "logs/qc/trimgalore/{sample}.log",
    params:
        quality=config["trimgalore"]["quality"],
        stringency=config["trimgalore"]["stringency"],
        e=config["trimgalore"]["e"],
        #trimgalore does not recommend using more than 4 threads
        threads_actual=config_threads if config_threads < 4 else 4,
    # setting it to mid ensures not no many trimgalore jobs launched
    threads: config_threads
    shell:
        "trim_galore "
        "--quality {params.quality} "
        "--stringency {params.stringency} "
        "-e {params.e} "
        "--gzip "
        "--output_dir results/qc/trimgalore/ "
        "--cores {params.threads_actual} "
        "--basename {wildcards.sample} "
        "--no_report_file "
        "{input} 2>&1 | tee {log}"


rule trimgalore_pe:
    input:
        get_fastq_pe,
    output:
        r1="results/qc/trimgalore/{sample}_val_1.fq.gz",
        r2="results/qc/trimgalore/{sample}_val_2.fq.gz",
    conda:
        "../envs/trimgalore.yaml"
    log:
        "logs/qc/trimgalore/{sample}.log",
    params:
        quality=config["trimgalore"]["quality"],
        stringency=config["trimgalore"]["stringency"],
        e=config["trimgalore"]["e"],
        #trimgalore does not recommend using more than 4 threads
        threads_actual=config_threads if config_threads < 4 else 4,
    # setting it to mid ensures not no many trimgalore jobs launched
    threads: config_threads
    shell:
        "trim_galore "
        "--quality {params.quality} "
        "--stringency {params.stringency} "
        "-e {params.e} "
        "--gzip "
        "--output_dir results/qc/trimgalore/ "
        "--cores {params.threads_actual} "
        "--basename {wildcards.sample} "
        "--paired --no_report_file "
        "{input[0]} {input[1]} 2>&1 | tee {log}"
