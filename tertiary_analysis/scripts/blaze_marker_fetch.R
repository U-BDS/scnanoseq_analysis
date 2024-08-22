# A brief script to fetch gencode v31 transcript IDs, gene names and transcript names.
# The input IDs are derived from the BLAZE Marker gene analysis

### Libraries ###
library(rtracklayer)
library(dplyr)

### Import GTF from GENCODE v31 (BLAZE original analysis) ###

# import GTFs
hg38_gtf_31 <- rtracklayer::import("./input/references/gencode.v31.annotation.gtf")

# specify genome
genome(hg38_gtf_31) <- "GRCh38"

hg38_gtf_31

### Import BLAZE gene markers ###

markers <- read.csv("./input/references/blaze_gene_markers_gencode_v31.csv")

### Select relevant metadata from annotations ###

hg38_gtf_31_meta <- as.data.frame(hg38_gtf_31) %>%
  filter(type == "transcript") %>%
  select(gene_id, gene_name, transcript_id, transcript_name)


### Merge selected metadata to acquire variables of interest to BLAZE markers ###

## v31 ##

blaze_marker_31 <- na.omit(merge(x = markers, 
                                 y = hg38_gtf_31_meta,
                                 by = "gene_id"))

write.csv(blaze_marker_31, file = "./input/references/blaze_gene_markers_gencode_v31_full.csv", row.names = FALSE)
