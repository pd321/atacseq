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
		r1 = temp("results/qc/trimgalore/{sample}_val_1.fq.gz"),
		r2 = temp("results/qc/trimgalore/{sample}_val_2.fq.gz")
	conda:
		"../envs/trimgalore.yaml"
	log:
		"logs/qc/trimgalore/{sample}.log"
	params:
		threads_actual = config_threads if config_threads < 4 else 4
	threads: config_threads
	shell:
		'trim_galore '
		'--gzip '
		'--output_dir results/qc/trimgalore '
		'--cores {params.threads_actual} '
		'--basename {wildcards.sample} '
		'--paired --no_report_file '
		'{input[0]} {input[1]} 2>&1 | tee {log}'
