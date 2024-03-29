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
setwd("data_analysis/nasal_microbiome/season_analysis")

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

####use the wilcox test to find the different taxa
###set it to summer (121 - 300) and winter (0-121, 301-365)

summer_idx <-
  which(ns_microbiome_sample_info$days >= 121 &
          ns_microbiome_sample_info$days <= 300)

winter_idx <-
  which(ns_microbiome_sample_info$days < 121 |
          ns_microbiome_sample_info$days > 300)

all_p_values <-
  seq_len(ncol(temp_data)) %>%
  purrr::map(function(i) {
    test <-
      wilcox.test(temp_data[summer_idx, i],
                  temp_data[winter_idx, i])
    test$p.value
  })

season_p_value <-
  data.frame(
    variable_id = colnames(temp_data),
    p_value = unlist(all_p_values),
    fdr = p.adjust(unlist(all_p_values), method = "fdr")
  )

save(season_p_value, file = "season_p_value")
load("season_p_value")
season_p_value$variable_id == variable_info$variable_id

season_p_value <-
  season_p_value %>%
  dplyr::left_join(variable_info, by = "variable_id")

openxlsx::write.xlsx(season_p_value, file = "season_p_value.xlsx", asTable = TRUE)
