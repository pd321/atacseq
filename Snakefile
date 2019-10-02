include: "src/common.smk"
include: "src/align.smk"
include: "src/peaks.smk"
include: "src/qc.smk"

rule all:
	input:
		"results/qc/multiqc/multiqc_report.html",
		"results/qc/ataqv/report/index.html",
		expand("results/peaks_annot/{sample}.xls", sample = samples),
		expand("results/bw/{sample}.bw", sample = samples)