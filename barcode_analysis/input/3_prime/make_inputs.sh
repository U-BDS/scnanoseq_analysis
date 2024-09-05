##############
### RAW BC ###
##############

#IN_RAW_BC=/data/project/U_BDS/pipeline_dev_space/pipeline_testing/nanopore_dev/UBDS_nf_core_scnanoseq/lianov/PRE_RELEASE_PREPRINT_DATA/3_prime/scnanoseq/results/SC3pv3_GEX_Human_PBMC_ONT/blaze/SC3pv3_GEX_Human_PBMC_ONT.putative_bc.no_header.csv
#
#cut -f2 -d',' $IN_RAW_BC | awk '{if ($0 != "") {print $0}}' | sort -u > scnanoseq_oxford.raw.tsv

####################
### CORRECTED BC ###
####################

IN_CORRECT_BC=/data/project/U_BDS/pipeline_dev_space/pipeline_testing/nanopore_dev/UBDS_nf_core_scnanoseq/lianov/PRE_RELEASE_PREPRINT_DATA/3_prime/scnanoseq/results/SC3pv3_GEX_Human_PBMC_ONT/bam/barcode_tagged/SC3pv3_GEX_Human_PBMC_ONT.tagged.bam
samtools view -d 'CB' $IN_CORRECT_BC | awk '{print $NF}' | cut -f3 -d':' | sort -u > scnanoseq_oxford.correct.tsv

################
### DEDUP BC ###
################

IN_DEDUP_BC=/data/project/U_BDS/pipeline_dev_space/pipeline_testing/nanopore_dev/UBDS_nf_core_scnanoseq/lianov/PRE_RELEASE_PREPRINT_DATA/3_prime/scnanoseq/results/SC3pv3_GEX_Human_PBMC_ONT/bam/dedup/SC3pv3_GEX_Human_PBMC_ONT.dedup.sorted.bam
samtools view -d 'CB' $IN_DEDUP_BC | awk '{print $NF}' | cut -f3 -d':' | sort -u > scnanoseq_oxford.dedup.tsv
