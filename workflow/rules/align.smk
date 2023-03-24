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

