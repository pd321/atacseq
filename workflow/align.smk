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

rule bowtie2:
	input:
		r1 = rules.trimgalore.output.r1,
		r2 = rules.trimgalore.output.r2
	output:
		temp("results/bam/{sample}_bowtie2.bam")
	conda:
		"envs/bowtie2.yaml"
	log:
		"logs/bowtie2/{sample}.log"
	params:	
		idx = config['bowtie2']['idx'],
		maxins = config['bowtie2']['maxins']
	threads: threads_high
	shell:
		'bowtie2 '
		'--very-sensitive '
		'--maxins {params.maxins} '
		'--no-mixed '
		'--no-discordant '
		'--time '
		'--threads {threads} '
		'-x {params.idx} '
		'-1 {input[0]} '
		'-2 {input[1]} '
		'2> {log} | '
		'samtools view '
		'--threads {threads} '
		'-b -| '
		'samtools sort '
		'--threads {threads} '
		'-o {output}'

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