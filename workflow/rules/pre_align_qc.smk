rule fastqc:
	input:
		lambda wildcards: samplesheet[wildcards.sample][wildcards.group]
	output:
		html="results/qc/fastqc/{sample}_{group}_fastqc.html",
		zip="results/qc/fastqc/{sample}_{group}_fastqc.zip"
	log:
		"logs/qc/fastqc/{sample}_{group}_fastqc.log"
	threads: config_threads
	wrapper:
		"v1.24.0/bio/fastqc"

rule trimgalore:
	input:
		get_fastq
	output:
		r1 = temp("results/bam/{sample}_R1_val_1.fq.gz"),
		r2 = temp("results/bam/{sample}_R2_val_2.fq.gz")
	conda:
		"envs/trimgalore.yaml"
	log:
		"logs/trimgalore/{sample}.log"
	params:
		quality = config['trimgalore']['quality'],
		stringency = config['trimgalore']['stringency'],
		e = config['trimgalore']['e']
	threads: threads_mid
	shell:
		'trim_galore '
		'--quality {params.quality} '
		'--stringency {params.stringency} '
		'-e {params.e} '
		'--gzip '
		'--output_dir results/bam/ '
		'--cores 4 '
		'--basename {wildcards.sample} '
		'--paired --no_report_file '
		'{input[0]} {input[1]} 2>&1 | tee {log}'
