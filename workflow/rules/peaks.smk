rule macs2:
	input:
		treatment = rules.re_sort_bam.output.bam
	output:
		multiext("results/peaks/{sample}/{sample}",
				 "_peaks.xls",
				 "_peaks.narrowPeak",
				 "_summits.bed"
				 )
	log:
		"logs/macs2/{sample}.log"
	threads: config_threads
	params:	
		extra = "--format BAMPE --gsize {macs2_genome} --keep-dup all --nomodel".format(macs2_genome = config['macs2']['gsize'])
	wrapper:
		"v1.24.0/bio/macs2/callpeak"
