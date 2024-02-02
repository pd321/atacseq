rule bowtie:
    input:
        rules.trimgalore_se.output,
    output:
        temp("results/bam/{sample}.bowtie.bam"),
    conda:
        "../envs/bowtie.yaml"
    threads: config_threads
    log:
        "logs/bowtie/{sample}.log",
    params:
        idx=config["bowtie"]["idx"],
        max_multiple_aln=config["bowtie"]["max_multiple_aln"],
        max_aln_per_read=config["bowtie"]["max_aln_per_read"],
        chunkmbs=config["bowtie"]["chunkmbs"],
        platform=config["platform"],
        labname="--sam-RG CN:{}".format(config["labname"]) if config["labname"] else "",
        date=time.strftime("%Y-%m-%d"),
    shell:
        "zcat {input} | "
        "bowtie "
        "-m {params.max_multiple_aln} "
        "-k {params.max_aln_per_read} "
        "--best --strata "
        "--sam "
        "--chunkmbs {params.chunkmbs} "
        "--time "
        "--threads {threads} "
        "--sam-RG SM:{wildcards.sample} "
        "--sam-RG LB:{wildcards.sample} "
        "--sam-RG ID:{wildcards.sample} "
        "--sam-RG PL:{params.platform} "
        "{params.labname} "
        "--sam-RG DT:{params.date} "
        "{params.idx} - 2> {log} | "
        "samtools view -@ {threads} -bS - | "
        "samtools sort -@ {threads} -o {output} -"


rule bowtie2:
    input:
        sample=[rules.trimgalore_pe.output.r1, rules.trimgalore_pe.output.r2],
        idx=multiext(
            config["bowtie2"]["idx"],
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    output:
        temp("results/bam/{sample}.bowtie2.bam"),
    log:
        "logs/bowtie2/{sample}.log",
    params:
        extra="--very-sensitive --no-mixed --no-discordant --time --maxins {maxins}".format(
            maxins=config["bowtie2"]["maxins"]
        ),
    threads: config_threads
    wrapper:
        "v1.24.0/bio/bowtie2/align"


rule sortbam:
    input:
        get_sortbam_input,
    output:
        bam=temp("results/bam/{sample}.sorted.bam"),
        idx=temp("results/bam/{sample}.sorted.bai"),
    log:
        "logs/samtools/sortbam/{sample}.log",
    threads: config_threads
    wrapper:
        "v1.23.5-48-gf27313f0/bio/samtools/sort"
