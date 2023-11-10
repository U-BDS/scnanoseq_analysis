# A brief script to fetch gencode v31 transcript IDs, gene names and transcript names.
# It will also convert the data to gencode v40
# The input IDs are derived from the BLAZE Marker gene analysis

### Libraries ###
library(rtracklayer)
library(dplyr)

### Import GTF from GENCODE v31 (BLAZE original) and v40 (scnanoseq) ###

# import GTFs
hg38_gtf_31 <- rtracklayer::import("./input/references/gencode.v31.annotation.gtf")
hg38_gtf_40 <- rtracklayer::import("./input/references/gencode.v40.annotation.gtf")

# specify genome
genome(hg38_gtf_31) <- "GRCh38"
genome(hg38_gtf_40) <- "GRCh38"

hg38_gtf_31
hg38_gtf_40

### Import BLAZE gene markers ###

markers <- read.csv("./input/references/blaze_gene_markers_gencode_v31.csv")

### Select relevant metadata from annotations ###

hg38_gtf_31_meta <- as.data.frame(hg38_gtf_31) %>%
  filter(type == "transcript") %>%
  select(gene_id, gene_name, transcript_id, transcript_name)


hg38_gtf_40_meta <- as.data.frame(hg38_gtf_40) %>%
  filter(type == "transcript") %>%
  select(gene_id, gene_name, transcript_id, transcript_name)


### Merge selected metadata to acquire variables of interest to BLAZE markers ###

## v31 ##

blaze_marker_31 <- na.omit(merge(x = markers, 
                                 y = hg38_gtf_31_meta,
                                 by = "gene_id"))

write.csv(blaze_marker_31, file = "./input/references/blaze_gene_markers_gencode_v31_full.csv", row.names = FALSE)

## v40 ##

# strip the version from gene_id from reference and query and look-up based on that

markers$gene_id_no_version <- gsub("\\..+$", "", markers$gene_id)
hg38_gtf_40_meta$gene_id_no_version <- gsub("\\..+$", "", hg38_gtf_40_meta$gene_id)

# merge

blaze_marker_40 <- na.omit(merge(x = markers,
                                 y = hg38_gtf_40_meta,
                                 by = "gene_id_no_version"))

# clean-up: remove un-needed colum and rename gene_id containing v40

blaze_marker_40 <- blaze_marker_40 %>%
  dplyr::rename("gene_id" = "gene_id.y") %>%
  dplyr::select(-gene_id_no_version, -gene_id.x)

write.csv(blaze_marker_40, file = "./input/references/blaze_gene_markers_gencode_v40_full.csv", row.names = FALSE)
