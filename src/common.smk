import pandas as pd

report: "report/workflow.rst"

configfile: 'config.yaml'

# Setup vars
threads_high = int(config['general']['threads'])
threads_mid = int(threads_high/2)
threads_low = int(threads_high/4)

# Load in metadata
metadata_file = config['general']['metadata_file']
metadata_df = pd.read_csv(metadata_file, sep = "\t", index_col = "SampleName")

# Setup samplesheet
samplesheet =  metadata_df.to_dict(orient='index')
samples = sorted(samplesheet.keys())

# Some of the required functions
def get_fastq(wildcards):
	return metadata_df.loc[(wildcards.sample), ["R1", "R2"]]

