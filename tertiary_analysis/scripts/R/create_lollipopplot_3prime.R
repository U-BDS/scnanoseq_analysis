#################
### LIBRARIES ###
#################

library(lubridate)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(viridis)

#################
### CONSTANTS ###
#################

SPLIT_PROCESSES <- c(
  "FASTQC",
  "NANOPLOT",
  "TOULLIGQC",
  "PREEXTRACT_FASTQ",
  "NANOFILT",
  "CORRECT_BARCODES",
  "UMITOOLS_DEDUP (GENOME)",
  "UMITOOLS_DEDUP (TRANSCRIPTOME)",
  "ISOQUANT (GENOME)",
  "SEURAT (GENOME)"
)

GENOME_WORKFLOW <- "PROCESS_LONGREAD_SCRNA_GENOME"
TRANSCRIPT_WORKFLOW <- "PROCESS_LONGREAD_SCRNA_TRANSCRIPT"

PROCESSES_ORDERED <- c(
  "FASTQC",
  "NANOPLOT",
  "TOULLIGQC",
  "NANOFILT",
  "BLAZE",
  "PREEXTRACT_FASTQ",
  "CORRECT_BARCODES",
  "MINIMAP2_ALIGN (GENOME)",
  "RSEQC_READDISTRIBUTION (GENOME)",
  "TAG_BARCODES (GENOME)",
  "UMITOOLS_DEDUP (GENOME)",
  "ISOQUANT (GENOME)",
  "SEURAT (GENOME)",
  "MINIMAP2_ALIGN (TRANSCRIPTOME)",
  "TAG_BARCODES (TRANSCRIPTOME)",
  "UMITOOLS_DEDUP (TRANSCRIPTOME)",
  "OARFISH (TRANSCRIPTOME)",
  "SEURAT (TRANSCRIPTOME)",
  "MULTIQC_FINALQC"
)

#################
### FUNCTIONS ###
#################

clean_stats_dataframe <- function(stats_path){
  stats_df <- read.csv(stats_path, header = FALSE)
  colnames(stats_df) <- c("name","cpus","rss","vmem","duration","realtime")
  
  ######################
  ### CLEAN REALTIME ###
  ######################
  
  stats_df["realtime_mod"] <- stats_df["realtime"]
  
  ### Remove time suffix
  stats_df <- stats_df %>%
    mutate(
      realtime_mod = ifelse(grepl("m", realtime_mod), gsub("m", "", realtime_mod), paste0("0 ",realtime_mod))
    ) %>%
    mutate(
      realtime_mod = ifelse(grepl("h", realtime_mod), gsub("h", "", realtime_mod), paste0("0 ", realtime_mod))
    ) %>%
    mutate(
      realtime_mod = ifelse(grepl("s", realtime_mod), gsub("s", "", realtime_mod), paste0(realtime_mod, " 0"))
    )
  
  data_realtime <- data.frame(str_split_fixed(stats_df$realtime_mod, ' ', 3))
  colnames(data_realtime) <- c("time_h","time_m","time_s")
  
  ### Calculate total
  data_realtime <- data_realtime %>%
    mutate(
      realtime_s = (as.integer(time_h) * 3600) + (as.integer(time_m) * 60) + as.integer(time_s)
    )
  
  stats_df["realtime_s"] <- data_realtime["realtime_s"]
  stats_df["realtime_m"] <- data_realtime["realtime_s"] / 60
  stats_df["realtime_h"] <- data_realtime["realtime_s"] / 3600
  
  ##################
  ### CLEAN VMEM ###
  ##################
  
  ### Clean up the vmem
  stats_df <- stats_df %>% mutate( 
    vmem_gb = ifelse(grepl("GB", vmem),as.integer(gsub(" GB", "", vmem)),as.integer(gsub(" MB", "", vmem)) / 1000)
  )
  
  ##################
  ### CLEAN NAME ###
  ##################
  
  ### Clean up the name and add the source (i.e. GENOME or TRANSCRIPTOME only process)
  stats_df <- 
    stats_df %>%
    mutate (name_base = gsub(".*:","", name)) %>%
    mutate(name_base = gsub(" \\(.*\\)", "", name_base)) %>%
    mutate (process_source = ifelse(
      grepl(GENOME_WORKFLOW, name), 
      "GENOME", 
      ifelse(
        grepl(TRANSCRIPT_WORKFLOW, name),
        "TRANSCRIPTOME", "")
    )
    ) %>%
    mutate (process_with_src = ifelse(
      process_source == "",
      name_base,
      paste0(name_base, " (", process_source,")")
    )
    )
  return(stats_df)
}

get_longest_process <- function(process_name, process_df) {
  single_process_df <- process_df[process_df$process_with_src == process_name,]
  return(single_process_df[which.max(single_process_df$realtime_s),])
}

###############
### DO WORK ###
###############
in_file <- "input/3_prime.all_stats.csv"
out_dir <- "output/lollipop_plots/"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

all_df <- clean_stats_dataframe(in_file)

### Get the Unique Processes

filtered_df <- all_df[all_df$process_with_src %in% PROCESSES_ORDERED, ]

# For all the split processes apart of this path, grab the max time
split_process_df <- do.call(
  rbind,
  lapply(
    SPLIT_PROCESSES,
    get_longest_process,
    process_df = filtered_df
  )
)

# Combine the max_time split process df with the nonsplit df
main_process_df <- rbind(
  split_process_df,
  filtered_df[! filtered_df$name_base %in% unique(split_process_df$name_base),]
)

# Order the processes based on order of the nf-core/scnanoseq pipeline
 main_process_df$process_with_src <- factor(
   main_process_df$process_with_src, 
   levels = rev(PROCESSES_ORDERED)
 )
 
 ############
 ### PLOT ###
 #############

lollipop_plot <- main_process_df %>%
  ggplot( aes(x=realtime_h, y=process_with_src)) +
  geom_segment( aes(x=0, xend=realtime_h, y=process_with_src, yend=process_with_src)) + 
  geom_point(aes(size = as.factor(cpus), color = vmem_gb)) +
  scale_y_discrete(
    labels=c(
      "PREEXTRACT_FASTQ" = "PREEXTRACT_FASTQ*",
      "NANOFILT" = "NANOFILT*",
      "CORRECT_BARCODES" = "CORRECT_BARCODES*",
      "UMITOOLS_DEDUP (GENOME)" = "UMITOOLS_DEDUP (GENOME)*",
      "UMITOOLS_DEDUP (TRANSCRIPTOME)" = "UMITOOLS_DEDUP (TRANSCRIPTOME)*",
      "ISOQUANT (GENOME)" = "ISOQUANT (GENOME)*"
    ) 
  ) + 
  scale_colour_viridis(option="D", name = "Memory (GB)") +
  labs(
    title = "nf-core/scnanoseq Process Requirements",
    x="Duration (hours)",
    y="Process Name",
    size="CPUs"
  ) +
  theme(
    axis.title = element_text(size = 20),
    axis.text.x = element_text(size=10),
    axis.text.y = element_text(size=10),
    plot.title = element_text(size = 20, color = "#009E73", hjust = 0.5),
    plot.subtitle = element_text(size = 15, face = "bold", hjust = 0.5),
  ) + theme_bw()

pdf(paste0(out_dir, "LollipopPlot_3prime.pdf"), height = 5, width = 8)
plot(lollipop_plot)
dev.off()