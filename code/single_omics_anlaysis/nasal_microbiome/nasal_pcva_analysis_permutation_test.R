#################
#####don't use this one, use the nasal_pcva_anlaysis2.R

##Principal Variance Component Analysis
no_source()

# set work directory
setwd(masstools::get_project_wd())
library(tidyverse)
library(phyloseq)
rm(list = ls())

source("code/tools.R")

{
  load("data_analysis/nasal_microbiome/data_preparation/sample_info")
  load("data_analysis/nasal_microbiome/data_preparation/expression_data")
  load("data_analysis/nasal_microbiome/data_preparation/variable_info")
  
  dim(variable_info)
  
  ns_microbiome_sample_info = sample_info
  ns_microbiome_expression_data = expression_data
  ns_microbiome_variable_info = variable_info
  
  ns_microbiome_sample_info =
    ns_microbiome_sample_info %>%
    dplyr::mutate(SSPG = as.numeric(SSPG)) %>%
    dplyr::mutate(iris = case_when(SSPG > 125 ~ "IR",
                                   SSPG < 125 ~ "IS"))
  
  library(lubridate)
  
  ns_microbiome_sample_info$days =
    ns_microbiome_sample_info$Date %>%
    lubridate::as_date() %>%
    yday()
  
  ns_microbiome_sample_info =
    ns_microbiome_sample_info %>%
    dplyr::mutate(
      days =  yday(lubridate::as_date(Date)),
      months = lubridate::month(lubridate::as_date(ns_microbiome_sample_info$Date)),
      weeks = lubridate::week(lubridate::as_date(ns_microbiome_sample_info$Date))
    )
}

#######work directory
setwd(masstools::get_project_wd())
dir.create("data_analysis/nasal_microbiome/season_analysis/permutation_test",
           recursive = TRUE)
setwd("data_analysis/nasal_microbiome/season_analysis/permutation_test")

zero_percent =
  ns_microbiome_expression_data %>%
  apply(1, function(x) {
    sum(x == 0) / ncol(ns_microbiome_expression_data)
  })

sum(zero_percent > 0.99) / nrow(ns_microbiome_variable_info)

##here we remove the genus with 0 > 99%
remain_idx = which(zero_percent < 0.99)

length(remain_idx)

ns_microbiome_expression_data =
  ns_microbiome_expression_data[remain_idx, ]

ns_microbiome_variable_info =
  ns_microbiome_variable_info[remain_idx, ]

##remove duplicated sample
which(duplicated(ns_microbiome_sample_info$sample_id))
# ns_microbiome_sample_info$sample_id[910]
# ns_microbiome_sample_info = ns_microbiome_sample_info[-910,]
# ns_microbiome_expression_data = ns_microbiome_expression_data[,-910]

library(plyr)
library(lubridate)
library(mgcv)
library(gratia)
library(pvca)

temp_data <-
  apply(ns_microbiome_expression_data, 1, function(x) {
    (x - mean(x)) / sd(x)
  })

rownames(temp_data) == ns_microbiome_sample_info$sample_id

# ####pvca analysis
library(pvca)
pvca_object_list <-
  purrr::map(1:100, function(i) {
    idx <-
      sample(1:nrow(temp_data), nrow(temp_data), replace = TRUE) %>%
      unique()
    
    counts <- as.matrix(t(temp_data))
    
    meta <- ns_microbiome_sample_info %>%
      dplyr::select(sample_id,
                    subject_id,
                    days) %>%
      as.data.frame()
    
    rownames(meta) = NULL
    
    meta =
      meta %>%
      tibble::column_to_rownames(var = "sample_id")
    
    pvca_object <-
      PVCA(
        counts = counts[, idx],
        meta = meta[idx, ],
        threshold = 0.6,
        inter = FALSE
      )
    
    pvca_object
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

save(pvca_object_list, file = "pvca_object_list")
