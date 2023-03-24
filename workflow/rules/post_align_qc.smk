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

rule remdup:
	input:
		rules.bowtie2.output
	output:
		bam = temp("results/bam/{sample}_remdup.bam"),
		bai = temp("results/bam/{sample}_remdup.bai")
	conda:
		"envs/picard.yaml"
	log: 
		metrics = "results/qc/picard/{sample}.txt",
		logfile = "logs/picard/{sample}.log"
	threads: threads_mid
	shell:
		'picard MarkDuplicates '
		'INPUT={input} '
		'OUTPUT={output.bam} '
		'REMOVE_DUPLICATES=true '
		'ASSUME_SORTED=true '
		'CREATE_INDEX=true '
		'METRICS_FILE={log.metrics} '
		'2> {log.logfile}'


rule bamfilter:
	input:
		bam = rules.remdup.output.bam,
		bai = rules.remdup.output.bai
	output:
		"results/bam/{sample}_filt.bam"
	conda:
		"envs/samtools.yaml"
	threads: threads_mid
	params:	
		remove = config['bamfilter']['remove'],
		keep = config['bamfilter']['keep'],
		mapqual = config['bamfilter']['mapqual']
	shell:
		'samtools idxstats '
		'{input.bam} | cut -f 1 | grep -v chrM | xargs samtools view '
		'--threads {threads} '
		'-F {params.remove} '
		'-f {params.keep} '
		'-q {params.mapqual} '
		'-o {output} {input.bam}'

rule flagstat:
	input:
		rules.bamfilter.output
	output:
		"results/qc/flagstat/{sample}.txt"
	conda:
		"envs/samtools.yaml"
	threads: threads_high
	shell:
		'samtools flagstat '
		'--threads {threads} '
		'{input} > {output}'

rule index_bam:
	input:
		rules.bamfilter.output
	output:
		"results/bam/{sample}_filt.bam.bai"
	conda:
		"envs/samtools.yaml"
	threads: threads_mid
	shell:
		'samtools index -@ {threads} {input}'

rule phantompeakqual:
	input:
		rules.bamfilter.output
	output:
		stats = "results/qc/phantompeakqual/{sample}.spp.out",
		crosscorrplot = report("results/qc/phantompeakqual/{sample}_filt.pdf", caption="report/phantompeakquals.rst", category="Quality control")
	conda:
		"../envs/phantompeakqualtools.yaml"
	log:
		"logs/phantompeakqual/{sample}.log"
	threads: config_threads	
	shell:
		'run_spp.R '
		'-c={input} '
		'-p={threads} '
		'-savp '
		'-odir=results/qc/phantompeakqual '
		'-out={output.stats} 2>&1 | tee {log}'

rule ataqv:
	input:
		bam = rules.bamfilter.output,
		bam_index = rules.index_bam.output,
		peaks = rules.blklist_filt.output
	output:
		temp("results/qc/ataqv/{sample}.ataqv.json")
	conda:
		"../envs/ataqv.yaml"
	log:
		"logs/ataqv/{sample}.log"
	params:
		organism = config['ataqv']['organism'],
		tssfile = config['ataqv']['tssfile'],
		excludedregionfile = config['ataqv']['excludedregionfile']
	threads: config_threads
	shell:
		'ataqv '
		'--peak-file {input.peaks} '
		'--threads {threads} '
		'--metrics-file {output} '
		'--name {wildcards.sample} '
		'--tss-file {params.tssfile} '
		'--excluded-region-file {params.excludedregionfile} '
		'{params.organism} {input.bam} 2>&1 | tee {log}'

rule ataqv_merge:
	input:
		expand("results/qc/ataqv/{sample}.ataqv.json", sample = samples)
	output:
		"results/qc/ataqv/report/index.html"
	conda:
		"../envs/ataqv.yaml"
	threads: config_threads
	shell:
		'mkarv '
		'--concurrency {threads} '
		'--force '
		'results/qc/ataqv/report {input}'

rule multiqc:
	input:
		expand(["results/qc/fastqc/{sample}_{group}_fastqc.html", "results/qc/flagstat/{sample}.txt", 
			"results/qc/phantompeakqual/{sample}.spp.out"], group=["R1", "R2"], sample = samples)
	output:
		report("results/qc/multiqc/multiqc_report.html", caption="report/multiqc.rst", category="Quality control")
	conda:
		"../envs/multiqc.yaml"
	threads: 1
	log:
		"logs/multiqc/multiqc.log"
	shell:
		'multiqc '
		'--force '
		'--outdir results/qc/multiqc '
		'--zip-data-dir '
		'. 2>&1 | tee {log}'
