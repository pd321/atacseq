rule bowtie2:
	input:
		sample = [rules.trimgalore.output.r1, rules.trimgalore.output.r2],
		idx=multiext(
			config['bowtie2']['idx'],
			".1.bt2",
			".2.bt2",
			".3.bt2",
			".4.bt2",
			".rev.1.bt2",
			".rev.2.bt2",
		),
	output:
		temp("results/bam/{sample}.bowtie2.bam")
	log:
		"logs/bowtie2/{sample}.log"
	params:	
		extra = "--very-sensitive --no-mixed --no-discordant --time --maxins {maxins}".format(maxins = config['bowtie2']['maxins'])
	threads: config_threads
	wrapper:
		"v1.24.0/bio/bowtie2/align"

rule sortbam:
	input: rules.bowtie2.output
	output: 
		bam = temp("results/bam/{sample}.sorted.bam"),
		idx = temp("results/bam/{sample}.sorted.bai")
	log:
		"logs/sortbam/{sample}.log"
	threads: config_threads
	wrapper:
		"v1.23.5-48-gf27313f0/bio/samtools/sort"
