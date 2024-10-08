library(tidyverse)
library(devtools)
library(data.table)
library(cowplot)
library(yaml)
library(glue)
library(here)
config <- read_yaml(here("config.yaml"))

library(signals)
#devtools::load_all(config$signals)

font_size <- config$font_size
gg_font_size <- font_size / 1 
gg_theme <- theme_cowplot(font_size = gg_font_size, font_family = "sans", line_size = 0.25)
theme_set(gg_theme)

removexaxis <- theme(axis.line.x=element_blank(),
                     axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank())

removeyaxis <- theme(axis.line.y=element_blank(),
                     axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())

get_data <- function(read_all = FALSE, filtbins = TRUE){
  library(glue)
  cnaneuploid <- fread(file.path(config$basedir, "zenododirectory/normal_breast/cn_aneuploid.csv.gz"))
  cnarm <- fread(file.path(config$basedir,"zenododirectory/normal_breast/cn_arm.csv.gz"))
  metrics <- fread(file.path(config$basedir, "zenododirectory/normal_breast/metrics.csv.gz"))
  
  cells_to_keep <- metrics %>% filter(keep_cell == TRUE) %>% pull(cell_id)
  
  arm_frequencies <- fread(file.path(config$basedir, "zenododirectory/normal_breast/chromosome_frequencies.csv.gz"))
  multiarm_frequencies <- fread(file.path(config$basedir, "zenododirectory/normal_breast/chromosomepair_frequencies.csv.gz"))
  
  #clinical data
  clinical <- xlsx::read.xlsx(config$supplementary_tables, sheetName = "TableS1") %>% distinct()
  clinical <- rename(clinical, current_past_cancer = `current_past_cancer.in.contralateral.breast.`)
  clinical <- mutate(clinical, genotype = ifelse(patient_id == "B1-6548", "BRCA1", genotype))
  
  cnaneuploid <- filter(cnaneuploid, cell_id %in% cells_to_keep)
  cnarm <- filter(cnarm, cell_id %in% cells_to_keep)
  
  metrics <- metrics %>% filter(sample %in% clinical$patient_id)
  arm_frequencies <- arm_frequencies %>% filter(sample %in% clinical$patient_id)
  multiarm_frequencies <- multiarm_frequencies %>% filter(sample %in% clinical$patient_id)
  cnarm <- cnarm %>% filter(sample %in% clinical$patient_id)
  cnaneuploid <- cnaneuploid %>% filter(sample %in% clinical$patient_id)
  
  
  if (read_all == TRUE){
    cn <- fread(file.path(config$basedir, "zenododirectory/normal_breast/cn.csv.gz"))
    if (filtbins == TRUE){
      cn <- cn[filt == FALSE]
    }
    cn <- cn %>% filter(sample %in% clinical$patient_id)
    return(list(cells_to_keep = cells_to_keep, cnarm = cnarm, clinical = clinical,  cn = cn, cnaneuploid = cnaneuploid, metrics = metrics, arm_frequencies = arm_frequencies, multiarm_frequencies = multiarm_frequencies))
  } else{
    return(list(cells_to_keep = cells_to_keep, cnarm = cnarm, clinical = clinical, cnaneuploid = cnaneuploid, metrics = metrics, arm_frequencies = arm_frequencies, multiarm_frequencies = multiarm_frequencies))
  }
}
