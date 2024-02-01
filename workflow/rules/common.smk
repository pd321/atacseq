import pandas as pd
from snakemake.utils import validate

report: "../report/workflow.rst"

configfile: 'config/config.yaml'

validate(config, schema="../schemas/config.schema.yaml")

# Setup vars
config_threads = int(config["threads"])
run_name = os.path.basename(os.getcwd())

# Load in metadata
metadata_file = "config/metadata.tsv"
metadata_df = pd.read_csv(metadata_file, sep = "\t").set_index('sample_name', drop=False)

# Setup samplesheet
samplesheet =  metadata_df.to_dict(orient='index')
samples = sorted(samplesheet.keys())

# Some of the required functions
def get_fastq(wildcards):
    return metadata_df.loc[(wildcards.sample), ["r1", "r2"]]

