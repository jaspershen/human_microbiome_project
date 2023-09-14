#################
#################
##Principal Variance Component Analysis
no_source()

# set work directory
setwd(masstools::get_project_wd())
library(tidyverse)
library(phyloseq)
rm(list = ls())
library(microbiomedataset)

source("code/tools.R")

{
  load("data_analysis/skin_microbiome/data_preparation/sample_info")
  load("data_analysis/skin_microbiome/data_preparation/expression_data")
  load("data_analysis/skin_microbiome/data_preparation/variable_info")
  
  object <-
    create_microbiome_dataset(
      expression_data = expression_data,
      sample_info = data.frame(sample_info, class = "Subject"),
      variable_info = variable_info
    )
  
  ###genus level
  object <-
    object %>%
    microbiomedataset::summarise_variables(what = "sum_intensity",
                                           group_by = "Genus")
  
  ###relative
  object <-
    microbiomedataset::transform2relative_intensity(object)
  
  skin_microbiome_sample_info = object@sample_info
  skin_microbiome_expression_data = object@expression_data
  skin_microbiome_variable_info = object@variable_info
  
  skin_microbiome_sample_info =
    skin_microbiome_sample_info %>%
    dplyr::mutate(SSPG = as.numeric(SSPG)) %>%
    dplyr::mutate(iris = case_when(SSPG > 125 ~ "IR",
                                   SSPG < 125 ~ "IS"))
  
  library(lubridate)
  
  skin_microbiome_sample_info$days =
    skin_microbiome_sample_info$Date %>%
    lubridate::as_date() %>%
    yday()
  
  skin_microbiome_sample_info =
    skin_microbiome_sample_info %>%
    dplyr::mutate(
      days =  yday(lubridate::as_date(Date)),
      months = lubridate::month(lubridate::as_date(skin_microbiome_sample_info$Date)),
      weeks = lubridate::week(lubridate::as_date(skin_microbiome_sample_info$Date))
    )
}

#######work directory
setwd(masstools::get_project_wd())
dir.create("data_analysis/skin_microbiome/season_analysis/genus_relative",
           recursive = TRUE)
setwd("data_analysis/skin_microbiome/season_analysis/genus_relative")

zero_percent =
  skin_microbiome_expression_data %>%
  apply(1, function(x) {
    sum(x == 0) / ncol(skin_microbiome_expression_data)
  })

sum(zero_percent > 0.99) / nrow(skin_microbiome_variable_info)

##here we remove the genus with 0 > 99%
remain_idx = which(zero_percent < 0.99)

length(remain_idx)

skin_microbiome_expression_data =
  skin_microbiome_expression_data[remain_idx,]

skin_microbiome_variable_info =
  skin_microbiome_variable_info[remain_idx,]

##remove duplicated sample
which(duplicated(skin_microbiome_sample_info$sample_id))
# skin_microbiome_sample_info$sample_id[910]
# skin_microbiome_sample_info = skin_microbiome_sample_info[-910,]
# skin_microbiome_expression_data = skin_microbiome_expression_data[,-910]

library(plyr)
library(lubridate)
library(mgcv)
library(gratia)
library(pvca)

# temp_data <-
#   apply(skin_microbiome_expression_data, 1, function(x) {
#     (x - mean(x)) / sd(x)
#   })

temp_data <-
  skin_microbiome_expression_data %>%
  t()

rownames(temp_data) == skin_microbiome_sample_info$sample_id

####use the wilcox test to find the different taxa
###set it to summer (121 - 300) and winter (0-121, 301-365)

summer_idx <-
  which(skin_microbiome_sample_info$days >= 121 &
          skin_microbiome_sample_info$days <= 300)

winter_idx <-
  which(skin_microbiome_sample_info$days < 121 |
          skin_microbiome_sample_info$days > 300)

all_p_values <-
  seq_len(ncol(temp_data)) %>%
  purrr::map(function(i) {
    test <-
      wilcox.test(temp_data[summer_idx, i],
                  temp_data[winter_idx, i])
    data.frame(
      p = test$p.value,
      winter_mean = mean(temp_data[winter_idx, i], na.rm = TRUE),
      summer_mean = mean(temp_data[summer_idx, i], na.rm = TRUE),
      winter_summer = mean(temp_data[winter_idx, i], na.rm = TRUE) - mean(temp_data[summer_idx, i], na.rm = TRUE),
      fc = mean(temp_data[winter_idx, i], na.rm = TRUE) / mean(temp_data[summer_idx, i], na.rm = TRUE)
    )
    
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

save(season_p_value, file = "season_p_value")
load("season_p_value")
season_p_value$variable_id == variable_info$variable_id

season_p_value <-
  season_p_value %>%
  dplyr::left_join(variable_info, by = "variable_id")

openxlsx::write.xlsx(season_p_value, file = "season_p_value.xlsx", asTable = TRUE)
