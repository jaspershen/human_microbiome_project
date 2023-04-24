###https://biofam.github.io/MOFA2/index.html
##https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-020-02015-1.pdf
no_source()

# set work directory
setwd(masstools::get_project_wd())
library(tidyverse)
library(phyloseq)
rm(list = ls())

source("code/tools.R")

###load nasal_microbiome
{
  load("data_analysis/nasal_microbiome/data_preparation/sample_info")
  load("data_analysis/nasal_microbiome/data_preparation/expression_data")
  load("data_analysis/nasal_microbiome/data_preparation/variable_info")
  
  nasal_microbiome_sample_info = sample_info
  nasal_microbiome_expression_data = expression_data
  nasal_microbiome_variable_info = variable_info
  
  nasal_microbiome_sample_info =
    nasal_microbiome_sample_info %>%
    dplyr::mutate(SSPG = as.numeric(SSPG)) %>%
    dplyr::mutate(iris = case_when(SSPG > 125 ~ "IR",
                                   SSPG < 125 ~ "IS"))
  
  library(lubridate)
  
  nasal_microbiome_sample_info$days =
    nasal_microbiome_sample_info$Date %>%
    lubridate::as_date() %>%
    yday()
  
  nasal_microbiome_sample_info =
    nasal_microbiome_sample_info %>%
    dplyr::mutate(
      days =  yday(lubridate::as_date(Date)),
      months = lubridate::month(lubridate::as_date(nasal_microbiome_sample_info$Date)),
      weeks = lubridate::week(lubridate::as_date(nasal_microbiome_sample_info$Date))
    )
  
  ###only remain mike's samples
  nasal_microbiome_sample_info =
    nasal_microbiome_sample_info %>%
    dplyr::filter(subject_id == "69-001")
  
  nasal_microbiome_expression_data =
    nasal_microbiome_expression_data[, nasal_microbiome_sample_info$sample_id]
  
  
  zero_percent =
    nasal_microbiome_expression_data %>%
    apply(1, function(x) {
      sum(x == 0) / ncol(nasal_microbiome_expression_data)
    })
  
  sum(zero_percent > 0.99) / nrow(nasal_microbiome_variable_info)
  
  ##here we remove the genus with 0 > 99%
  remain_idx = which(zero_percent < 0.99)
  
  nasal_microbiome_expression_data = nasal_microbiome_expression_data[remain_idx,]
  nasal_microbiome_variable_info = nasal_microbiome_variable_info[remain_idx,]
  
}

#######work directory
setwd(masstools::get_project_wd())
dir.create("data_analysis/nasal_microbiome/season_analysis/days_vs_microbiome")
setwd("data_analysis/nasal_microbiome/season_analysis/days_vs_microbiome")

temp_expression_data <-
  log(nasal_microbiome_expression_data + 1, 2)

temp_expression_data <-
  temp_expression_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

rownames(temp_expression_data)

nasal_microbiome_expression_data2 = temp_expression_data

#####just use the pca to do dimension reduction for exposome_chemical
days <- matrix(nasal_microbiome_sample_info$days, nrow = 1) %>% 
  as.data.frame()

colnames(days) <- colnames(nasal_microbiome_expression_data2)

days <-
  days %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()
#-------------------------------------------------------------------------------
###prepare data
dim(nasal_microbiome_expression_data)
dim(nasal_microbiome_sample_info)

####multiple linear regression
total_r2 <-
  nasal_microbiome_expression_data2 %>%
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
                 p = summary(lm_object)$coefficients[2,4])
    }
  ) %>%
  do.call(rbind, .) %>% 
  as.data.frame()

save(total_r2, file = "total_r2")

nasal_microbiome_variable_info$variable_id == rownames(total_r2)

total_r2 <-
  cbind(nasal_microbiome_variable_info, total_r2)

write.csv(total_r2, file = "total_r2.csv", row.names = FALSE)



