##Principal Variance Component Analysis

no_source()

# set work directory
setwd(masstools::get_project_wd())
library(tidyverse)
library(phyloseq)
rm(lioral = ls())

source("code/tools.R")

{
  load("data_analysis/oral_microbiome/data_preparation/sample_info")
  load("data_analysis/oral_microbiome/data_preparation/expression_data")
  load("data_analysis/oral_microbiome/data_preparation/variable_info")
  
  dim(variable_info)
  
  oral_microbiome_sample_info = sample_info
  oral_microbiome_expression_data = expression_data
  oral_microbiome_variable_info = variable_info
  
  oral_microbiome_sample_info =
    oral_microbiome_sample_info %>%
    dplyr::mutate(SSPG = as.numeric(SSPG)) %>%
    dplyr::mutate(iris = case_when(SSPG > 125 ~ "IR",
                                   SSPG < 125 ~ "IS"))
  
  library(lubridate)
  
  oral_microbiome_sample_info$days =
    oral_microbiome_sample_info$Date %>%
    lubridate::as_date() %>%
    yday()
  
  oral_microbiome_sample_info =
    oral_microbiome_sample_info %>%
    dplyr::mutate(
      days =  yday(lubridate::as_date(Date)),
      months = lubridate::month(lubridate::as_date(oral_microbiome_sample_info$Date)),
      weeks = lubridate::week(lubridate::as_date(oral_microbiome_sample_info$Date))
    )
}

#######work directory
setwd(masstools::get_project_wd())
setwd("data_analysis/oral_microbiome/season_analysis")

zero_percent =
  oral_microbiome_expression_data %>%
  apply(1, function(x) {
    sum(x == 0) / ncol(oral_microbiome_expression_data)
  })

sum(zero_percent > 0.99) / nrow(oral_microbiome_variable_info)

##here we remove the genus with 0 > 99%
remain_idx = which(zero_percent < 0.99)

length(remain_idx)

oral_microbiome_expression_data =
  oral_microbiome_expression_data[remain_idx, ]

oral_microbiome_variable_info =
  oral_microbiome_variable_info[remain_idx, ]

##remove duplicated sample
which(duplicated(oral_microbiome_sample_info$sample_id))
# oral_microbiome_sample_info$sample_id[910]
# oral_microbiome_sample_info = oral_microbiome_sample_info[-910,]
# oral_microbiome_expression_data = oral_microbiome_expression_data[,-910]

library(plyr)
library(lubridate)
library(mgcv)
library(gratia)
library(pvca)

# temp_data <-
#   apply(oral_microbiome_expression_data, 1, function(x) {
#     (x - mean(x)) / sd(x)
#   })

temp_data <-
  oral_microbiome_expression_data %>% 
  t()

rownames(temp_data) == oral_microbiome_sample_info$sample_id

####use the wilcox test to find the different taxa
###set it to summer (121 - 300) and winter (0-121, 301-365)

summer_idx <-
  which(oral_microbiome_sample_info$days >= 121 &
          oral_microbiome_sample_info$days <= 300)

winter_idx <-
  which(oral_microbiome_sample_info$days < 121 |
          oral_microbiome_sample_info$days > 300)

all_p_values <-
  seq_len(ncol(temp_data)) %>%
  purrr::map(function(i) {
    test <-
      wilcox.test(temp_data[summer_idx, i],
                  temp_data[winter_idx, i])
    data.frame(p = test$p.value, 
               winter_mean = mean(temp_data[winter_idx, i], na.rm = TRUE),
               summer_mean = mean(temp_data[summer_idx, i], na.rm = TRUE),
               winter_summer = mean(temp_data[winter_idx, i], na.rm = TRUE) - mean(temp_data[summer_idx, i], na.rm = TRUE),
               fc = mean(temp_data[winter_idx, i], na.rm = TRUE)/mean(temp_data[summer_idx, i], na.rm = TRUE))
    
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

season_p_value <-
  data.frame(
    variable_id = colnames(temp_data),
    p_value = all_p_values$p,
    fdr = p.adjust(all_p_values$p, method = "fdr"),
    fc = all_p_values$fc,
    winter_mean = all_p_values$winter_mean,
    summer_mean = all_p_values$summer_mean,
    winter_summer = all_p_values$winter_summer
  )

sum(season_p_value$p_value < 0.05)
sum(season_p_value$fdr < 0.05)
save(season_p_value, file = "season_p_value")
load("season_p_value")
season_p_value <-
  season_p_value %>%
  dplyr::left_join(variable_info, by = "variable_id")

openxlsx::write.xlsx(season_p_value, file = "season_p_value.xlsx", asTable = TRUE)
