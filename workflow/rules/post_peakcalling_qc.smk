rule ataqv:
	input:
		bam = rules.re_sort_bam.output.bam,
		peaks = "results/peaks/{sample}/{sample}_peaks.narrowPeak"
	output:
		temp("results/qc/ataqv/{sample}.ataqv.json")
	conda:
		"../envs/ataqv.yaml"
	log:
		"logs/ataqv/{sample}.log"
	params:
		organism = config['ataqv']['organism'],
		tssfile = config['ataqv']['tssfile'],
		excludedregionfile = config['blklist_regions']
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
			"results/qc/phantompeakqual/{sample}.spp.out"], group=["r1", "r2"], sample = samples)
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
