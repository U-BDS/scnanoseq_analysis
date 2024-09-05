PREFIX=/data/project/U_BDS/pipeline_dev_space/pipeline_testing/nanopore_dev/UBDS_nf_core_scnanoseq/lianov/PRE_RELEASE_PREPRINT_DATA/5_prime/scnanoseq/results/SC5pv2_GEX_Human_Lung_Carcinoma_DTC_ONT

##############i
### RAW BC ###
##############

#IN_RAW_BC=$PREFIX/blaze/SC5pv2_GEX_Human_Lung_Carcinoma_DTC_ONT.putative_bc.no_header.csv

#cut -f2 -d',' $IN_RAW_BC | awk '{if ($0 != "") {print $0}}' | sort -u > scnanoseq_oxford.raw.tsv

####################
### CORRECTED BC ###
####################

IN_CORRECT_BC=$PREFIX/bam/barcode_tagged/SC5pv2_GEX_Human_Lung_Carcinoma_DTC_ONT.tagged.bam
samtools view -d 'CB' $IN_CORRECT_BC | awk '{print $NF}' | cut -f3 -d':' | sort -u > scnanoseq_oxford.correct.tsv

################
### DEDUP BC ###
################

IN_DEDUP_BC=$PREFIX/bam/dedup/SC5pv2_GEX_Human_Lung_Carcinoma_DTC_ONT.dedup.sorted.bam
samtools view -d 'CB' $IN_DEDUP_BC | awk '{print $NF}' | cut -f3 -d':' | sort -u > scnanoseq_oxford.dedup.tsv
