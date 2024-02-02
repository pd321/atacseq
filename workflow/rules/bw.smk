rule bamcoverage_bw:
    input:
        remdup_bam=rules.re_sort_bam.output.bam,
    output:
        "results/bw/{sample}.bw",
    conda:
        "../envs/deeptools.yaml"
    log:
        "logs/bamcoverage_bw/{sample}.log",
    threads: config_threads
    params:
        bin_size=config["bamcoverage_bw"]["bin_size"],
        smooth_length=config["bamcoverage_bw"]["smooth_length"],
        normalize_using=config["bamcoverage_bw"]["normalize_using"],
        read_extension_length=get_read_extension_length,
    shell:
        "bamCoverage --bam {input.remdup_bam} "
        "--outFileName {output} "
        "--outFileFormat bigwig "
        "--numberOfProcessors {threads} "
        "--extendReads {params.read_extension_length} "
        "--normalizeUsing {params.normalize_using} "
        "--binSize {params.bin_size} "
        "--smoothLength {params.smooth_length} 2>&1 | tee {log}"
