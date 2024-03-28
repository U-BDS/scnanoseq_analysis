# make vector of all gene names mapped to chrM
hg38_mito_gene_ids <- as.data.frame(hg38_annot) %>%
dplyr::filter(seqnames == "MT", type == "gene") %>%
dplyr::select(seqnames, gene_id, gene_name, gene_biotype)

# make vector of all transcript names mapped to chrM
hg38_mito_transcript_ids <- as.data.frame(hg38_annot) %>%
dplyr::filter(seqnames == "MT", type == "transcript") %>%
dplyr::select(seqnames, transcript_id, transcript_name, gene_biotype)
hg38_mito_gene_ids
hg38_mito_transcript_ids



feature_type <- c(feature_type_scnanoseq, feature_type_CR)
# pipeline input due to small differences in feature names for now as scnanoseq outputs Ensembl IDs
pipeline <- c("scnanoseq", "scnanoseq","scnanoseq", "scnanoseq","CR", "CR")

### TODO: just need to re-add the section where this object already contains the gene/transcript names gsub for . for scnanoseq sections
ids_per_obj <- mapply(FUN = function(x,y, p) {
  if (y == "gene") {
    obj_id <- hg38_mito_gene_ids
    # which MT genes are in the object
    if (p == "scnanoseq") {
      #gsub_version
      gene_names_scnano <- gsub("-",".",hg38_mito_gene_ids$gene_name)
      obj_id$in_data <- gene_names_scnano %in% rownames(x@assays$RNA@counts)
      
      #TODO: see note above (same below)
      
    } else {
      obj_id$in_data <- hg38_mito_gene_ids$gene_name %in% rownames(x@assays$RNA@counts)
    }
    # keep only those which are in the obj
    obj_id <- dplyr::filter(obj_id, in_data==TRUE)
  } else {
    obj_id <- hg38_mito_transcript_ids
    if (p == "scnanoseq") {
      trans_name_scnano <- gsub("-",".",hg38_mito_transcript_ids$transcript_name)
      obj_id$in_data <- trans_name_scnano %in% rownames(x@assays$RNA@counts)
    } else {
      # which MT transcripts are in the object
      obj_id$in_data <- hg38_mito_transcript_ids$transcript_name %in% rownames(x@assays$RNA@counts)
    }
    # keep only those which are in the obj
    obj_id <- dplyr::filter(obj_id, in_data==TRUE)
  }
  return(obj_id)
}, x=seurat_objs, y=feature_type, p=pipeline, SIMPLIFY = FALSE)

ids_per_obj # notes Illumina had less genes overall than oxford data
