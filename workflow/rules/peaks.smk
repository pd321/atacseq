rule macs2:
	input:
		rules.bamfilter.output
	output:
		peaks = temp("results/peaks/{sample}/{sample}_peaks.narrowPeak"),
		trtbdg = temp("results/peaks/{sample}/{sample}_treat_pileup.bdg"),
		cntbdg = temp("results/peaks/{sample}/{sample}_control_lambda.bdg")
	conda:
		"envs/macs2.yaml"
	log:
		"logs/macs2/{sample}.log"
	threads: threads_mid
	params:	
		qvalue = config['macs2']['qvalue'],
		gsize = config['macs2']['gsize']
	shell:
		'macs2 '
		'callpeak '
		'--treatment {input} '
		'--format BAMPE '
		'--qvalue {params.qvalue} '
		'--nomodel '
		'--gsize {params.gsize} '
		'--keep-dup all '
		'--bdg --SPMR '
		'--outdir results/peaks/{wildcards.sample} '
		'--name {wildcards.sample} 2> {log}'

rule bdg2bw:
	"""
	Convert the bdg files to bw files
	"""
	input:
		rules.macs2.output.trtbdg
	output:
		bw = "results/bw/{sample}.bw",
		clip = temp("results/peaks/{sample}/{sample}_treat_pileup.bdg.clip"),
		sort_clip = temp("results/peaks/{sample}/{sample}_treat_pileup.bdg.sort.clip")
	conda:
		"envs/bdg2bw.yaml"
	threads: threads_mid
	params:
		genome_sizes = config['bdg2bw']['genome_sizes']
	shell:
		'../scripts/bdg2bw '
		'{input} '
		'{params.genome_sizes} '
		'{output.bw}'

rule blklist_filt:
	input:
		rules.macs2.output.peaks
	output:
		"results/peaks/{sample}/{sample}_peaks_filt.narrowPeak"
	conda:
		"envs/bedtools.yaml"
	threads: threads_low
	params:
		blklist_regions = config['blklist_filt']['blklist_regions']
	shell:
		'bedtools intersect '
		'-v -a {input} -b {params.blklist_regions} > {output}'

rule homer_annotate:
	input:
		rules.blklist_filt.output
	output:
		temp("results/peaks_annot/{sample}.bed")
	log:
		"logs/homer/{sample}.log"
	threads: threads_low
	params:
		genome = config['homer']['genome']
	shell:
		'annotatePeaks.pl '
		'{input} '
		'{params.genome} 1> {output} 2> {log}'

rule clean_merge_homer:
	input:
		macs2_out = rules.macs2.output.peaks,
		homer_out = rules.homer_annotate.output
	output:
		"results/peaks_annot/{sample}.xls"
	threads: threads_low
	log:
		"logs/clean_merge_homer/{sample}.log"
	shell:
		'Rscript --vanilla ../scripts/homerCleanAndMerge.R '
		'--macs2Peaks {input.macs2_out} '
		'--homerOutFile {input.homer_out} '
		'--outFile {output} 2>&1 | tee {log}'

