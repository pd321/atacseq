# Libload
suppressPackageStartupMessages(library(optparse))
library(dplyr)
library(readr)
library(futile.logger)

# Getting in the command line options
option_list <- list(
  make_option(c("-i", "--macs2Peaks"), 
              help="A narrowPeak/broad file, REQUIRED",
              metavar="FILEPATH"),
  make_option(c("-m", "--homerOutFile"), 
              help="A homer annotated file of the narrowPeak file, REQUIRED",
              metavar="FILEPATH"),
  make_option(c("-o","--outFile"),
              help="File to write sanitized and merge output to, REQUIRED")
)
opt <- parse_args(OptionParser(option_list=option_list))
flog.info("Recieved the following options")
opt

# Read in the summary file
summary.df <- readr::read_tsv(opt$macs2Peaks, col_names = FALSE)
# Failure is observed when there are no peaks called. Handling it with this check. If no peaks called no need to adjust.
# Just write a blank output
# Handling it at macs end because homer adds a header
if(nrow(summary.df) > 0){
  # Taking only first 9 columns to keep handling same for narrow and broad peaks
  summary.df <- summary.df[,1:9]
  colnames(summary.df) <- c("chrom", "chromStart", "chromEnd", "peakID", "score", "strand", "foldChange", "pValue", "qValue")
  
  
  # Read in homer file
  homer.df <- readr::read_tsv(opt$homerOutFile)
  # Rename the Peak ID column where homer writes the command too
  colnames(homer.df)[1] <- "peakID"
  # Remove the feature details
  homer.df <- homer.df %>% mutate(Annotation = sapply(Annotation, function(x){strsplit(as.character(x), split = " (", fixed = TRUE)[[1]][1]}))
  # Keep only required columns
  homer.df <- homer.df %>% select(c("peakID", "Annotation", "Gene Name", "Gene Type", "Distance to TSS", "Gene Alias"	,"Gene Description", "Nearest PromoterID", "Nearest Refseq", "Nearest Ensembl"))
  
  # Now lets merge these two dfs
  summary_homer_merged <- dplyr::left_join(summary.df, homer.df, by.y = "peakID")
  
  # Remove the peak id
  summary_homer_merged <- summary_homer_merged %>% select(-peakID, -score, -strand, -pValue) %>% arrange(desc(qValue))
  
  # Write output
  write.table(summary_homer_merged, file = opt$outFile, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
} else {
  out_write_cmd <- paste0("touch ", opt$outFile)
  system(out_write_cmd)
}
sessionInfo()
