import os
import pandas as pd
from snakemake.utils import validate
import time


report: "../report/workflow.rst"


configfile: "config/config.yaml"


validate(config, schema="../schemas/config.schema.yaml")

# Setup vars
config_threads = int(config["threads"])
run_name = os.path.basename(os.getcwd())

# Load in metadata
metadata_file = "config/metadata.tsv"
metadata_df = pd.read_csv(metadata_file, sep="\t").set_index("sample_name", drop=False)

# Setup samplesheet
samplesheet = metadata_df.to_dict(orient="index")
samples = sorted(samplesheet.keys())
se_samples = [
    sample for sample in samplesheet if samplesheet[sample]["end_type"] == "se"
]
pe_samples = [
    sample for sample in samplesheet if samplesheet[sample]["end_type"] == "pe"
]


# Some of the required functions
def get_fastq_se(wildcards):
    return metadata_df.loc[wildcards.sample, "r1"]


def get_fastq_pe(wildcards):
    return {
        "r1": metadata_df.loc[wildcards.sample, "r1"],
        "r2": metadata_df.loc[wildcards.sample, "r2"],
    }

def get_sortbam_input(wildcards):
    if metadata_df.loc[wildcards.sample, "end_type"] == "se":
        return "results/bam/{}.bowtie.bam".format(wildcards.sample)
    elif metadata_df.loc[wildcards.sample, "end_type"] == "pe":
        return "results/bam/{}.bowtie2.bam".format(wildcards.sample)


def get_bamfilter_params(wildcards):
    if metadata_df.loc[wildcards.sample, "end_type"] == "se":
        return " -F {remove} -q {mapqual} ".format(
            remove=config["filter_bam"]["se_remove_reads_flags"],
            mapqual=config["filter_bam"]["mapqual_filter"],
        )
    elif metadata_df.loc[wildcards.sample, "end_type"] == "pe":
        return " -F {remove} -f {keep} -q {mapqual} ".format(
            remove=config["filter_bam"]["pe_remove_reads_flags"],
            keep=config["filter_bam"]["pe_keep_reads_flags"],
            mapqual=config["filter_bam"]["mapqual_filter"],
        )


def get_macs2_params(wildcards):
    if metadata_df.loc[wildcards.sample, "end_type"] == "se":
        return "--format BAM --gsize {macs2_genome} --keep-dup all --nomodel --extsize 200".format(
            macs2_genome=config["macs2"]["gsize"]
        )
    elif metadata_df.loc[wildcards.sample, "end_type"] == "pe":
        return "--format BAMPE --gsize {macs2_genome} --keep-dup all --nomodel".format(
            macs2_genome=config["macs2"]["gsize"]
        )


def get_read_extension_length(wildcards):
    if metadata_df.loc[wildcards.sample, "end_type"] == "se":
        return config["bamcoverage_bw"]["read_extension_length"]
    elif metadata_df.loc[wildcards.sample, "end_type"] == "pe":
        return ""
