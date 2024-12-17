library(Seurat)

oarfish_input_dir <- "./input/oarfish/base_minimap/"
output_mtx <- "./input/oarfish/SC3pv3_GEX_Human_PBMC_ONT.transcript_counts.tsv.base_minimap"

oarfish_output <- Read10X(oarfish_input_dir,
                          gene.column = 1,
                          cell.column = 2)

rows <- rownames(oarfish_output)
transcript_names <- do.call(rbind, strsplit(rows, split = '|', fixed = TRUE))
row.names(oarfish_output) <- transcript_names[,1]

write.table(oarfish_output, output_mtx, sep='\t',quote = FALSE )
