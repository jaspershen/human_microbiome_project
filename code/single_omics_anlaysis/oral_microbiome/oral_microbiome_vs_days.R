###https://biofam.github.io/MOFA2/index.html
##https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-020-02015-1.pdf
no_source()

# set work directory
setwd(masstools::get_project_wd())
library(tidyverse)
library(phyloseq)
rm(list = ls())

source("code/tools.R")

###load oral_microbiome
{
  load("data_analysis/oral_microbiome/data_preparation/sample_info")
  load("data_analysis/oral_microbiome/data_preparation/expression_data")
  load("data_analysis/oral_microbiome/data_preparation/variable_info")
  
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
  
  ###only remain mike's samples
  oral_microbiome_sample_info =
    oral_microbiome_sample_info %>%
    dplyr::filter(subject_id == "69-001")
  
  oral_microbiome_expression_data =
    oral_microbiome_expression_data[, oral_microbiome_sample_info$sample_id]
  
  
  zero_percent =
    oral_microbiome_expression_data %>%
    apply(1, function(x) {
      sum(x == 0) / ncol(oral_microbiome_expression_data)
    })
  
  sum(zero_percent > 0.99) / nrow(oral_microbiome_variable_info)
  
  ##here we remove the genus with 0 > 99%
  remain_idx = which(zero_percent < 0.99)
  
  oral_microbiome_expression_data = oral_microbiome_expression_data[remain_idx, ]
  oral_microbiome_variable_info = oral_microbiome_variable_info[remain_idx, ]
  
}

#######work directory
setwd(masstools::get_project_wd())
dir.create("data_analysis/oral_microbiome/season_analysis/days_vs_microbiome")
setwd("data_analysis/oral_microbiome/season_analysis/days_vs_microbiome")

temp_expression_data <-
  log(oral_microbiome_expression_data + 1, 2)

temp_expression_data <-
  temp_expression_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

rownames(temp_expression_data)

oral_microbiome_expression_data2 = temp_expression_data

#####just use the pca to do dimension reduction for exposome_chemical
days <- matrix(oral_microbiome_sample_info$days, nrow = 1) %>%
  as.data.frame()

colnames(days) <- colnames(oral_microbiome_expression_data2)

days <-
  days %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

#-------------------------------------------------------------------------------
###prepare data
dim(oral_microbiome_expression_data)
dim(oral_microbiome_sample_info)

####multiple linear regression
total_r2 <-
  oral_microbiome_expression_data2 %>%
  t() %>%
  as.data.frame() %>%
  purrr::map(
    .f = function(x) {
      temp_data <-
        rbind(y = x,
              days) %>%
        t() %>%
        as.data.frame()
      
      colnames(temp_data) <-
        c("y", "days")
      
      lm_object <-
        lm(formula = y ~ .,
           data = temp_data)
      data.frame(r2 = summary(lm_object)$r.squared,
                 p = summary(lm_object)$coefficients[2, 4])
    }
  ) %>%
  do.call(rbind, .) %>%
  as.data.frame()

save(total_r2, file = "total_r2")

oral_microbiome_variable_info$variable_id == rownames(total_r2)

total_r2 <-
  cbind(oral_microbiome_variable_info, total_r2)

write.csv(total_r2, file = "total_r2.csv", row.names = FALSE)
