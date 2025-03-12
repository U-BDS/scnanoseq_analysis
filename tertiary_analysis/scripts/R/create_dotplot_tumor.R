# This script will generate the DotPlot present in Supplementary Figure 1
# Providing a way to visualize the runtimes for analayzing all samples from
#   Shiau et al. datasets

#################
### LIBRARIES ###
#################

options(scipen = 100, digits = 2)

library(lubridate)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(viridis)
library(scales)

#################
### CONSTANTS ###
#################

BOX_OUTDIR <- "output/dot_plots/"

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
  "SEURAT (TRANSCRIPTOME)"
)

#################
### FUNCTIONS ###
#################

clean_stats_dataframe <- function(stats_path){
  stats_df <- read.csv(stats_path, header = FALSE)
  colnames(stats_df) <- c("name","cpus","rss","vmem","duration","realtime")
  
  ######################
  ### CLEAN realtime ###
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

make_dot_plot <- function(all_process_df){
  
  # Set jitter to help with dot visibility
  jitter = position_jitter(
    width = 0.1,
    height = 0.1,
    seed = 12345
  )

  # Dot Plot
  dot_plot <- 
    ggplot(
      all_process_df,
      aes(
        x=realtime_h,
        y=process_with_src,
        color=sample_name,
        fill=sample_name
      )
    ) +
    geom_dotplot(
      position = jitter,
      dotsize=0.5,
      stackdir="center",
      stackratio=0,
      binaxis="x",
      alpha=0.4
    ) +
    scale_x_continuous(
      trans = trans_new(
        name = "custom",
        transform = function(x) ifelse(x < 1, x*10, x+9),
        inverse = function(x) ifelse(x < 1, x / 10, x-9)
      ),
      breaks=c(0:10)
    ) +
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
    labs(
      title = "nf-core/scnanoseq Runtimes",
      x="Duration (hours)",
      y="Process Name"
    ) +
    theme(
      axis.title = element_text(size = 20),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(size = 20, color = "#009E73", hjust = 0.5),
      plot.subtitle = element_text(size = 15, face = "bold", hjust = 0.5),
    ) + 
    theme_bw()
  
  pdf(paste0(BOX_OUTDIR, "DotPlot_tumor_datasets.pdf"), height = 5, width = 8)
  plot(dot_plot)
  
  dev.off()
}

###############
### DO WORK ###
###############

# Make output directories
dir.create(BOX_OUTDIR, showWarnings = FALSE, recursive = TRUE)

all_stats_df <- clean_stats_dataframe("input/cancer_datasets.all_stats.csv")

# Only parse out the 'important' processes
main_process_df <- all_stats_df[
  all_stats_df$process_with_src %in% PROCESSES_ORDERED,
]

# Parse out the sample name into its own column
sample_name <- data.frame(
  do.call(
    'rbind',
    strsplit(as.character(main_process_df$name), ' ')
  )
)

colnames(sample_name) <- c('process_name', 'sample')
sample_name <- sample_name %>%
  mutate(
    cleaned_sample = gsub("\\(", "", 
                     gsub("\\)", "", 
                     gsub("\\..*", "", sample)))
)

main_process_df$sample_name <- sample_name$cleaned_sample

# Order the processes
main_process_df$process_with_src <- factor(
  main_process_df$process_with_src, 
  levels = rev(PROCESSES_ORDERED)
)

############
### PLOT ###
############

result <- make_dot_plot(
  main_process_df
)
