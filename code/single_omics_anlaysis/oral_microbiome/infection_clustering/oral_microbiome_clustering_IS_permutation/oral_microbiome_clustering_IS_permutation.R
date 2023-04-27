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
  
  expression_data =
    data.table::fread(here::here("data/from_xin/Genus Table/OR/Genus_OR.csv")) %>%
    as.data.frame() %>%
    tibble::column_to_rownames(var = "SampleID") %>%
    dplyr::select(-c(V1:SubjectID)) %>%
    t() %>%
    as.data.frame()
  
  oral_microbiome_variable_info =
    oral_microbiome_variable_info[match(rownames(expression_data),
                                        oral_microbiome_variable_info$Genus),]
  
  oral_microbiome_variable_info$Genus == rownames(expression_data)
  
  ###remove the variables which Genus are NA
  remove_idx = which(is.na(oral_microbiome_variable_info$Genus))
  remove_idx
  if (length(remove_idx) > 0) {
    oral_microbiome_variable_info = oral_microbiome_variable_info[-remove_idx,]
    expression_data = expression_data[-remove_idx,]
  }
  
  rownames(expression_data) = oral_microbiome_variable_info$variable_id
  
  oral_microbiome_expression_data <-
    expression_data
  
  oral_microbiome_variable_info =
    oral_microbiome_variable_info %>%
    dplyr::filter(!stringr::str_detect(Genus, "Unclassified_Bacteria"))
  
  oral_microbiome_expression_data =
    oral_microbiome_expression_data[oral_microbiome_variable_info$variable_id,]
  
  dim(oral_microbiome_sample_info)
  dim(oral_microbiome_variable_info)
  
  rownames(oral_microbiome_expression_data) == oral_microbiome_variable_info$variable_id
  colnames(oral_microbiome_expression_data) == oral_microbiome_sample_info$sample_id
  
  load("data_analysis/phenotype_data_sample_wise/data_preparation/sample_info")
  
  sample_info <-
    sample_info %>%
    dplyr::distinct(sample_id, .keep_all = TRUE)
  
  oral_microbiome_sample_info <-
    oral_microbiome_sample_info %>%
    dplyr::left_join(sample_info[, c("sample_id", "CL4")], by = "sample_id")
  
}

#######work directory
setwd(masstools::get_project_wd())
dir.create("data_analysis/oral_microbiome/clustering/phylum_IS_permutation",
           recursive = TRUE)
setwd("data_analysis/oral_microbiome/clustering/phylum_IS_permutation")

head(oral_microbiome_sample_info)
dim(oral_microbiome_expression_data)

###only remain IS
oral_microbiome_sample_info <-
  oral_microbiome_sample_info %>%
  # dplyr::left_join(metadata[, c("SubjectID", "IRIS")],
  #                  by = c("subject_id" = "SubjectID")) %>%
  dplyr::filter(!is.na(IRIS)) %>%
  dplyr::filter(IRIS == "IS")

oral_microbiome_expression_data <-
  oral_microbiome_expression_data[, oral_microbiome_sample_info$sample_id]

###only remain the infectional and pro and post two samples
infection_idx <-
  which(oral_microbiome_sample_info$CL4 == "Infection")

#pre-healthy (–H) state (healthy baselines within 186 days before the event’s onset),
#event early (EE) state (visits on days 1–6 of the event),
#event late (EL) state (visits on days 7–14 since the event started),
#recovery (RE) state (visits within days 15–40 since the event started),
#and finally post-healthy (+H) state (visits within 186 days after the event;

all_index <-
  infection_idx %>%
  purrr::map(function(i) {
    temp_info <-
      oral_microbiome_sample_info[i, ]
    idx <-
      which(oral_microbiome_sample_info$subject_id == temp_info$subject_id)
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
    temp_sample_info <- oral_microbiome_sample_info[idx, ]
    day_diff <-
      as.numeric(temp_sample_info$Date - oral_microbiome_sample_info[i, ]$Date)
    
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
        oral_microbiome_expression_data[, idx$negtaive_h_id, drop = FALSE] %>%
        apply(1, mean)
      
      ee <-
        oral_microbiome_expression_data[, idx$ee_id, drop = FALSE] %>%
        apply(1, mean)
      # `-`(negtaive_h)
      
      el <-
        oral_microbiome_expression_data[, idx$el_id, drop = FALSE] %>%
        apply(1, mean)
      # `-`(negtaive_h)
      
      re <-
        oral_microbiome_expression_data[, idx$re_id, drop = FALSE] %>%
        apply(1, mean)
      # `-`(negtaive_h)
      
      positive_h <-
        oral_microbiome_expression_data[, idx$positive_h_id, drop = FALSE] %>%
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
    oral_microbiome_variable_info
  
  save(new_expression_data,
       file = paste("new_expression_data", i, sep = "_"))
  save(new_variable_info, file = paste("new_variable_info", i, sep = "_"))
  save(new_sample_info, file = paste("new_sample_info", i, sep = "_"))
}
