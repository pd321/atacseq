#!/usr/bin/env bash
set -euxo pipefail 

# Set base outdir
base_out_dir="data/external"
mkdir -p $base_out_dir/{blacklists,chrom_sizes,ataqv_tss}

# Fetch blacklists
wget -O $base_out_dir/blacklists/hg19_blacklist.bed.gz https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.bed.gz
gunzip $base_out_dir/blacklists/hg19_blacklist.bed.gz

# Fetch chrom sizes
wget -O $base_out_dir/chrom_sizes/hg19.chrom.sizes http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes

# Get tss
wget -O $base_out_dir/ataqv_tss/hg19.tss.refseq.bed.gz https://github.com/ParkerLab/ataqv/raw/master/data/tss/hg19.tss.refseq.bed.gz
gunzip $base_out_dir/ataqv_tss/hg19.tss.refseq.bed.gz