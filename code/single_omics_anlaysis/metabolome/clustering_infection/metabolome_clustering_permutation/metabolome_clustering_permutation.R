###https://biofam.github.io/MOFA2/index.html
##https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-020-02015-1.pdf
no_source()

# set work directory
setwd(masstools::get_project_wd())
library(tidyverse)
library(phyloseq)
rm(list = ls())

source("code/tools.R")

###plasma metabolome
{
  load(here::here(
    "data_analysis/metabolome/data_preparation/expression_data"
  ))
  load(here::here("data_analysis/metabolome/data_preparation/sample_info"))
  load(here::here(
    "data_analysis/metabolome/data_preparation/variable_info"
  ))
}

metabolome_expression_data = expression_data
metabolome_sample_info = sample_info
metabolome_variable_info = variable_info

metabolome_expression_data <-
  log(metabolome_expression_data + 1, 2)

metabolome_sample_info <-
  metabolome_sample_info %>% 
  dplyr::mutate(Date = CollectionDate)

metabolome_sample_info$Date =
  as.Date(metabolome_sample_info$Date, "%m/%d/%y")

dim(metabolome_expression_data)
length(unique(metabolome_sample_info$subject_id))

#######work directory
setwd(masstools::get_project_wd())
dir.create("data_analysis/metabolome/clustering_permutation",
           recursive = TRUE)
setwd("data_analysis/metabolome/clustering_permutation")

head(metabolome_sample_info)
dim(metabolome_expression_data)

###only remain the infectional and pro and post two samples
infection_idx <-
  which(metabolome_sample_info$CL4 == "Infection")

#pre-healthy (–H) state (healthy baselines within 186 days before the event’s onset),
#event early (EE) state (visits on days 1–6 of the event),
#event late (EL) state (visits on days 7–14 since the event started),
#recovery (RE) state (visits within days 15–40 since the event started),
#and finally post-healthy (+H) state (visits within 186 days after the event;

all_index <-
  infection_idx %>%
  purrr::map(function(i) {
    temp_info <-
      metabolome_sample_info[i,]
    idx <-
      which(metabolome_sample_info$subject_id == temp_info$subject_id)
    if (length(idx) < 5) {
      return(
        list(
          negtaive_h_id = character(0),
          ee_id = character(0),
          el_id = character(0),
          re_id = character(0),
          positive_h_id = character(0)
        )
      )
    }
    temp_sample_info <- metabolome_sample_info[idx,]
    day_diff <-
      as.numeric(temp_sample_info$Date - metabolome_sample_info[i,]$Date)
    
    negtaive_h_idx <-
      which(day_diff < -6 & day_diff >= -180)
    
    ee_idx <-
      which(day_diff <= 0 & day_diff >= -6)
    
    el_idx <-
      which(day_diff >= 1 & day_diff <= 7)
    
    re_idx <-
      which(day_diff >= 8 & day_diff <= 32)
    
    positive_h_idx <-
      which(day_diff >= 33 & day_diff <= 180)
    
    list(
      negtaive_h_id = temp_sample_info$sample_id[negtaive_h_idx],
      ee_id = temp_sample_info$sample_id[ee_idx],
      el_id = temp_sample_info$sample_id[el_idx],
      re_id = temp_sample_info$sample_id[re_idx],
      positive_h_id = temp_sample_info$sample_id[positive_h_idx]
    )
  })

####remove some index
temp <-
  all_index %>%
  purrr::map(function(x) {
    unlist(lapply(x, length))
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  apply(1, function(x) {
    all(x > 0)
  }) %>%
  which()

all_index <-
  all_index[temp]

###remove some duplicated index
all_ee_id <-
  all_index %>%
  lapply(function(x) {
    x$ee_id
  }) %>%
  unlist()

duplicated_ee_id <-
  unique(all_ee_id[which(duplicated(all_ee_id))])

if (length(duplicated_ee_id) > 0) {
  remove_idx <-
    duplicated_ee_id %>%
    purrr::map(function(id) {
      all_index %>%
        lapply(function(x) {
          id %in% x$ee_id
        }) %>%
        unlist() %>%
        which() %>%
        `[`(1)
    }) %>%
    unlist() %>%
    unique()
  
  all_index <-
    all_index[-remove_idx]
  
}


#####permutation
for (i in 1:100) {
  cat(i, " ")
  all_index_new <-
    all_index %>%
    lapply(function(x) {
      x %>%
        lapply(function(y) {
          sample(y, size = length(y), replace = TRUE) %>%
            unique()
        })
    })
  
  new_expression_data <-
    all_index_new %>%
    purrr::map(function(idx) {
      negtaive_h <-
       metabolome_expression_data[, idx$negtaive_h_id, drop = FALSE] %>%
        apply(1, mean)
      
      ee <-
       metabolome_expression_data[, idx$ee_id, drop = FALSE] %>%
        apply(1, mean)
      # `-`(negtaive_h)
      
      el <-
       metabolome_expression_data[, idx$el_id, drop = FALSE] %>%
        apply(1, mean)
      # `-`(negtaive_h)
      
      re <-
       metabolome_expression_data[, idx$re_id, drop = FALSE] %>%
        apply(1, mean)
      # `-`(negtaive_h)
      
      positive_h <-
       metabolome_expression_data[, idx$positive_h_id, drop = FALSE] %>%
        apply(1, mean)
      # `-`(negtaive_h)
      
      temp <-
        cbind(negtaive_h, ee, el, re, positive_h) %>%
        as.data.frame()
      
      subject_id <-
        idx$ee_id[1]
      
      colnames(temp) <-
        paste(colnames(temp), subject_id, sep = ".")
      temp
    }) %>%
    do.call(cbind, .) %>%
    as.data.frame()
  
  colnames(new_expression_data)
  
  new_sample_info <-
    data.frame(
      sample_id = colnames(new_expression_data),
      class = stringr::str_split(colnames(new_expression_data), "\\.") %>%
        lapply(function(x)
          x[1]) %>%
        unlist,
      sample_id = stringr::str_split(colnames(new_expression_data), "\\.") %>%
        lapply(function(x)
          x[2]) %>%
        unlist
    )
  
  new_variable_info <-
   metabolome_variable_info
  
  save(new_expression_data,
       file = paste("new_expression_data", i, sep = "_"))
  save(new_variable_info, file = paste("new_variable_info", i, sep = "_"))
  save(new_sample_info, file = paste("new_sample_info", i, sep = "_"))
}
