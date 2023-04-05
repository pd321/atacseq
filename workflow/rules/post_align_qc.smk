rule remdup:
	input:
		bams = rules.sortbam.output.bam
	output:
		bam = temp("results/bam/{sample}.sorted.remdup.bam"),
		bai = temp("results/bam/{sample}.sorted.remdup.bai"),
		metrics = "results/qc/remdup/{sample}.metrics.txt"
	threads: 
		config_threads
	log:
		"logs/qc/remdup/{sample}.log"
	params:
		extra="--REMOVE_DUPLICATES true --CREATE_INDEX true",
	wrapper:
		"v1.23.5-48-gf27313f0/bio/picard/markduplicates"

rule remove_blacklist_reads:
	input:
		left = rules.remdup.output.bam,
		right = config['blklist_regions']
	output:
		temp("results/bam/{sample}.sorted.remdup.nonblklst.bam")
	log:
		"logs/qc/remove_blacklist_reads/{sample}.log"
	params:
		extra="-v"
	wrapper:
		"v1.24.0/bio/bedtools/intersect"

rule index_bam:
	input:
		rules.remove_blacklist_reads.output
	output:
		temp("results/bam/{sample}.sorted.remdup.nonblklst.bai")
	log:
		"logs/qc/samtools_index/{sample}.log"
	threads: config_threads
	wrapper:
		"v1.25.0-52-g79342b73/bio/samtools/index"

rule bamfilter:
	input:
		bam = rules.remove_blacklist_reads.output,
		bai = rules.index_bam.output
	output:
		bam = temp("results/bam/{sample}.sorted.remdup.nonblklst.filt.bam"),
		idx = temp("results/bam/{sample}.sorted.remdup.nonblklst.filt.bai")
	conda:
		"../envs/samtools.yaml"
	threads: config_threads
	params:	
		remove_reads = config['bamfilter']['remove'],
		keep_reads = config['bamfilter']['keep'],
		mapqual = config['bamfilter']['mapqual']
	shell:
		'samtools idxstats '
		'{input.bam} | cut -f 1 | grep -v chrM | xargs samtools view '
		'--threads {threads} '
		'--write-index '
		'-F {params.remove_reads} '
		'-f {params.keep_reads} '
		'-q {params.mapqual} '
		'-o {output.bam}##idx##{output.idx} {input.bam}'

rule re_sort_bam:
	input: 
		bam = rules.bamfilter.output.bam,
		bai = rules.bamfilter.output.idx
	output: 
		bam = "results/bam/{sample}.sorted.remdup.nonblklst.filt.resort.bam",
		idx = "results/bam/{sample}.sorted.remdup.nonblklst.filt.resort.bai"
	log:
		"logs/samtools/resortbam/{sample}.log"
	threads: config_threads
	wrapper:
		"v1.23.5-48-gf27313f0/bio/samtools/sort"

rule flagstat:
	input:
		rules.re_sort_bam.output.bam
	output:
		"results/qc/flagstat/{sample}.txt"
	log:
		"logs/qc/flagstat/{sample}.log",
	threads: config_threads
	wrapper:
		"v1.24.0/bio/samtools/flagstat"

rule phantompeakqual:
	input:
		rules.re_sort_bam.output.bam
	output:
		stats = "results/qc/phantompeakqual/{sample}.sorted.remdup.nonblklst.filt.resort.spp.out",
		crosscorrplot = report("results/qc/phantompeakqual/{sample}.sorted.remdup.nonblklst.filt.resort.pdf", caption="report/phantompeakquals.rst", category="Quality control")
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
