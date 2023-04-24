#' ---
#' title: "oral microbiome oral_microbiome correlation"
#' author:
#'   - name: "Xiaotao Shen"
#'     url: https://www.shenxt.info/
#'     affiliation: Stanford School of Medicine
#' date: "`r Sys.Date()`"
#' site: distill::distill_website
#' output:
#'   distill::distill_article:
#'     code_folding: false
#' ---

#+ r setup, echo=TRUE, eval = TRUE, include = TRUE

no_function()
# set work directory

setwd(masstools::get_project_wd())
library(tidyverse)
rm(list = ls())

source("code/tools.R")

######work directory
setwd(masstools::get_project_wd())

##oral
edge_percentage_ratio <- vector(mode = "list", length = 20)
phylum_percentage_ratio <- vector(mode = "list", length = 20)

for (i in 1:20) {
  cat(i, " ")
  load(
    paste0(
      "data_analysis/correlation_network/whole_data_set/intra_oral_microbiome/permutation/intra_oral_microbiome_lm_adjusted_cor_",
      i
    )
  )
  
  oral_sample_cor <-
    intra_oral_microbiome_lm_adjusted_cor
  
  load(
    "data_analysis/correlation_network/whole_data_set/intra_oral_microbiome/permutation/intra_oral_sample_wise_dim"
  )
  
  #####sample wise
  oral_sample_cor <-
    oral_sample_cor %>%
    dplyr::rename(from = microbiome, to = metabolite) %>%
    dplyr::mutate(from_class = "oral", to_class = "oral") %>%
    dplyr::mutate(from = paste("oral", from, sep = ""),
                  to = paste("oral", to, sep = ""),
    ) %>%
    dplyr::select(-name)
  
  sample_cor <-
    rbind(oral_sample_cor) %>%
    dplyr::filter(p_adjust < 0.05 * 0.005)
  
  name = apply(sample_cor, 1, function(x) {
    paste(sort(as.character(x[1:2])), collapse = "_")
  })
  
  sample_cor =
    sample_cor %>%
    dplyr::mutate(name = name) %>%
    dplyr::distinct(name, .keep_all = TRUE)
  
  #####intra-stool microbiome
  sample_cor$class_name <-
    purrr::map(as.data.frame(t(sample_cor)), function(x) {
      paste(sort(as.character(x)[7:8]), collapse = "_")
    })  %>%
    unlist() %>%
    unname()
  
  sample_wise_number <-
    purrr::map(as.data.frame(t(sample_cor)), function(x) {
      c(sort(as.character(x)[7:8]), name = paste(sort(as.character(x)[7:8]), collapse = "_"))
    }) %>%
    do.call(rbind, .) %>%
    as.data.frame() %>%
    dplyr::group_by(name) %>%
    dplyr::summarise(significant_n = n()) %>%
    dplyr::arrange(name) %>%
    dplyr::mutate(theoretical_n = NA)
  
  sample_wise_number$theoretical_n[sample_wise_number$name == "oral_oral"] =
    (intra_oral_sample_wise_dim[1] * (intra_oral_sample_wise_dim[1] - 1)) /
    2
  
  sample_wise_number <-
    sample_wise_number %>%
    dplyr::mutate(percentage = significant_n * 100 / theoretical_n)
  
  new_info =
    sample_wise_number$name %>%
    purrr::map(function(x) {
      significant_n_pos =
        sample_cor %>%
        dplyr::filter(class_name == x & cor > 0) %>%
        nrow()
      
      significant_n_neg =
        sample_cor %>%
        dplyr::filter(class_name == x & cor < 0) %>%
        nrow()
      
      data.frame(significant_n_pos = significant_n_pos,
                 significant_n_neg = significant_n_neg)
      
    }) %>%
    dplyr::bind_rows()
  
  sample_wise_number =
    cbind(sample_wise_number, new_info)
  
  sample_wise_number =
    sample_wise_number %>%
    dplyr::mutate(
      percentage_pos = significant_n_pos * 100 / theoretical_n,
      percentage_neg = significant_n_neg * 100 / theoretical_n
    )
  
  load(here::here(
    "data_analysis/oral_microbiome/data_preparation/variable_info"
  ))
  oral_variable_info = variable_info
  
  from <-
    sample_cor %>%
    dplyr::filter(class_name == "oral_oral") %>%
    pull(from) %>%
    stringr::str_replace("oral", "") %>%
    unique()
  
  temp <-
    oral_variable_info %>%
    dplyr::filter(variable_id %in% from)
  
  temp <-
    temp %>%
    dplyr::count(Phylum) %>%
    dplyr::mutate(n = n * 100 / sum(n))
  
  phylum_percentage_ratio[[i]] <-
    temp
  
  edge_data =
    sample_wise_number %>%
    dplyr::select(name, percentage) %>%
    tidyr::separate(col = name,
                    sep = "_",
                    into = c("from", "to")) %>%
    dplyr::mutate(from = stringr::str_to_title(from)) %>%
    dplyr::mutate(to = stringr::str_to_title(to))
  
  edge_percentage_ratio[[i]] <-
    edge_data
}

temp <-
  edge_percentage_ratio %>%
  do.call(rbind, .) %>%
  as.data.frame()

temp_oral <- temp


##skin
edge_percentage_ratio <- vector(mode = "list", length = 20)
phylum_percentage_ratio <- vector(mode = "list", length = 20)

for (i in 1:20) {
  cat(i, " ")
  load(
    paste0(
      "data_analysis/correlation_network/whole_data_set/intra_skin_microbiome/permutation/intra_skin_microbiome_lm_adjusted_cor_",
      i
    )
  )
  
  skin_sample_cor <-
    intra_skin_microbiome_lm_adjusted_cor
  
  load(
    "data_analysis/correlation_network/whole_data_set/intra_skin_microbiome/permutation/intra_skin_sample_wise_dim"
  )
  
  #####sample wise
  skin_sample_cor <-
    skin_sample_cor %>%
    dplyr::rename(from = microbiome, to = metabolite) %>%
    dplyr::mutate(from_class = "skin", to_class = "skin") %>%
    dplyr::mutate(from = paste("skin", from, sep = ""),
                  to = paste("skin", to, sep = ""),
    ) %>%
    dplyr::select(-name)
  
  sample_cor <-
    rbind(skin_sample_cor) %>%
    dplyr::filter(p_adjust < 0.05 * 0.005)
  
  name = apply(sample_cor, 1, function(x) {
    paste(sort(as.character(x[1:2])), collapse = "_")
  })
  
  sample_cor =
    sample_cor %>%
    dplyr::mutate(name = name) %>%
    dplyr::distinct(name, .keep_all = TRUE)
  
  #####intra-stool microbiome
  sample_cor$class_name <-
    purrr::map(as.data.frame(t(sample_cor)), function(x) {
      paste(sort(as.character(x)[7:8]), collapse = "_")
    })  %>%
    unlist() %>%
    unname()
  
  sample_wise_number <-
    purrr::map(as.data.frame(t(sample_cor)), function(x) {
      c(sort(as.character(x)[7:8]), name = paste(sort(as.character(x)[7:8]), collapse = "_"))
    }) %>%
    do.call(rbind, .) %>%
    as.data.frame() %>%
    dplyr::group_by(name) %>%
    dplyr::summarise(significant_n = n()) %>%
    dplyr::arrange(name) %>%
    dplyr::mutate(theoretical_n = NA)
  
  sample_wise_number$theoretical_n[sample_wise_number$name == "skin_skin"] =
    (intra_skin_sample_wise_dim[1] * (intra_skin_sample_wise_dim[1] - 1)) /
    2
  
  sample_wise_number <-
    sample_wise_number %>%
    dplyr::mutate(percentage = significant_n * 100 / theoretical_n)
  
  new_info =
    sample_wise_number$name %>%
    purrr::map(function(x) {
      significant_n_pos =
        sample_cor %>%
        dplyr::filter(class_name == x & cor > 0) %>%
        nrow()
      
      significant_n_neg =
        sample_cor %>%
        dplyr::filter(class_name == x & cor < 0) %>%
        nrow()
      
      data.frame(significant_n_pos = significant_n_pos,
                 significant_n_neg = significant_n_neg)
      
    }) %>%
    dplyr::bind_rows()
  
  sample_wise_number =
    cbind(sample_wise_number, new_info)
  
  sample_wise_number =
    sample_wise_number %>%
    dplyr::mutate(
      percentage_pos = significant_n_pos * 100 / theoretical_n,
      percentage_neg = significant_n_neg * 100 / theoretical_n
    )
  
  load(here::here(
    "data_analysis/skin_microbiome/data_preparation/variable_info"
  ))
  skin_variable_info = variable_info
  
  from <-
    sample_cor %>%
    dplyr::filter(class_name == "skin_skin") %>%
    pull(from) %>%
    stringr::str_replace("skin", "") %>%
    unique()
  
  temp <-
    skin_variable_info %>%
    dplyr::filter(variable_id %in% from)
  
  temp <-
    temp %>%
    dplyr::count(Phylum) %>%
    dplyr::mutate(n = n * 100 / sum(n))
  
  phylum_percentage_ratio[[i]] <-
    temp
  
  edge_data =
    sample_wise_number %>%
    dplyr::select(name, percentage) %>%
    tidyr::separate(col = name,
                    sep = "_",
                    into = c("from", "to")) %>%
    dplyr::mutate(from = stringr::str_to_title(from)) %>%
    dplyr::mutate(to = stringr::str_to_title(to))
  
  edge_percentage_ratio[[i]] <-
    edge_data
}

temp <-
  edge_percentage_ratio %>%
  do.call(rbind, .) %>%
  as.data.frame()

temp_skin <- temp

##nasal
edge_percentage_ratio <- vector(mode = "list", length = 20)
phylum_percentage_ratio <- vector(mode = "list", length = 20)

for (i in 1:20) {
  cat(i, " ")
  load(
    paste0(
      "data_analysis/correlation_network/whole_data_set/intra_nasal_microbiome/permutation/intra_nasal_microbiome_lm_adjusted_cor_",
      i
    )
  )
  
  nasal_sample_cor <-
    intra_nasal_microbiome_lm_adjusted_cor
  
  load(
    "data_analysis/correlation_network/whole_data_set/intra_nasal_microbiome/permutation/intra_nasal_sample_wise_dim"
  )
  
  #####sample wise
  nasal_sample_cor <-
    nasal_sample_cor %>%
    dplyr::rename(from = microbiome, to = metabolite) %>%
    dplyr::mutate(from_class = "nasal", to_class = "nasal") %>%
    dplyr::mutate(from = paste("nasal", from, sep = ""),
                  to = paste("nasal", to, sep = ""),
    ) %>%
    dplyr::select(-name)
  
  sample_cor <-
    rbind(nasal_sample_cor) %>%
    dplyr::filter(p_adjust < 0.05 * 0.005)
  
  name = apply(sample_cor, 1, function(x) {
    paste(sort(as.character(x[1:2])), collapse = "_")
  })
  
  sample_cor =
    sample_cor %>%
    dplyr::mutate(name = name) %>%
    dplyr::distinct(name, .keep_all = TRUE)
  
  #####intra-stool microbiome
  sample_cor$class_name <-
    purrr::map(as.data.frame(t(sample_cor)), function(x) {
      paste(sort(as.character(x)[7:8]), collapse = "_")
    })  %>%
    unlist() %>%
    unname()
  
  sample_wise_number <-
    purrr::map(as.data.frame(t(sample_cor)), function(x) {
      c(sort(as.character(x)[7:8]), name = paste(sort(as.character(x)[7:8]), collapse = "_"))
    }) %>%
    do.call(rbind, .) %>%
    as.data.frame() %>%
    dplyr::group_by(name) %>%
    dplyr::summarise(significant_n = n()) %>%
    dplyr::arrange(name) %>%
    dplyr::mutate(theoretical_n = NA)
  
  sample_wise_number$theoretical_n[sample_wise_number$name == "nasal_nasal"] =
    (intra_nasal_sample_wise_dim[1] * (intra_nasal_sample_wise_dim[1] - 1)) /
    2
  
  sample_wise_number <-
    sample_wise_number %>%
    dplyr::mutate(percentage = significant_n * 100 / theoretical_n)
  
  new_info =
    sample_wise_number$name %>%
    purrr::map(function(x) {
      significant_n_pos =
        sample_cor %>%
        dplyr::filter(class_name == x & cor > 0) %>%
        nrow()
      
      significant_n_neg =
        sample_cor %>%
        dplyr::filter(class_name == x & cor < 0) %>%
        nrow()
      
      data.frame(significant_n_pos = significant_n_pos,
                 significant_n_neg = significant_n_neg)
      
    }) %>%
    dplyr::bind_rows()
  
  sample_wise_number =
    cbind(sample_wise_number, new_info)
  
  sample_wise_number =
    sample_wise_number %>%
    dplyr::mutate(
      percentage_pos = significant_n_pos * 100 / theoretical_n,
      percentage_neg = significant_n_neg * 100 / theoretical_n
    )
  
  load(here::here(
    "data_analysis/nasal_microbiome/data_preparation/variable_info"
  ))
  nasal_variable_info = variable_info
  
  from <-
    sample_cor %>%
    dplyr::filter(class_name == "nasal_nasal") %>%
    pull(from) %>%
    stringr::str_replace("nasal", "") %>%
    unique()
  
  temp <-
    nasal_variable_info %>%
    dplyr::filter(variable_id %in% from)
  
  temp <-
    temp %>%
    dplyr::count(Phylum) %>%
    dplyr::mutate(n = n * 100 / sum(n))
  
  phylum_percentage_ratio[[i]] <-
    temp
  
  edge_data =
    sample_wise_number %>%
    dplyr::select(name, percentage) %>%
    tidyr::separate(col = name,
                    sep = "_",
                    into = c("from", "to")) %>%
    dplyr::mutate(from = stringr::str_to_title(from)) %>%
    dplyr::mutate(to = stringr::str_to_title(to))
  
  edge_percentage_ratio[[i]] <-
    edge_data
}

temp <-
  edge_percentage_ratio %>%
  do.call(rbind, .) %>%
  as.data.frame()

temp_nasal <- temp


##stool
edge_percentage_ratio <- vector(mode = "list", length = 20)
phylum_percentage_ratio <- vector(mode = "list", length = 20)

for (i in 1:20) {
  cat(i, " ")
  load(
    paste0(
      "data_analysis/correlation_network/whole_data_set/intra_stool_microbiome/permutation/intra_stool_microbiome_lm_adjusted_cor_",
      i
    )
  )
  
  stool_sample_cor <-
    intra_stool_microbiome_lm_adjusted_cor
  
  load(
    "data_analysis/correlation_network/whole_data_set/intra_stool_microbiome/permutation/intra_stool_sample_wise_dim"
  )
  
  #####sample wise
  stool_sample_cor <-
    stool_sample_cor %>%
    dplyr::rename(from = microbiome, to = metabolite) %>%
    dplyr::mutate(from_class = "stool", to_class = "stool") %>%
    dplyr::mutate(from = paste("stool", from, sep = ""),
                  to = paste("stool", to, sep = ""),
    ) %>%
    dplyr::select(-name)
  
  sample_cor <-
    rbind(stool_sample_cor) %>%
    dplyr::filter(p_adjust < 0.05 * 0.005)
  
  name = apply(sample_cor, 1, function(x) {
    paste(sort(as.character(x[1:2])), collapse = "_")
  })
  
  sample_cor =
    sample_cor %>%
    dplyr::mutate(name = name) %>%
    dplyr::distinct(name, .keep_all = TRUE)
  
  #####intra-stool microbiome
  sample_cor$class_name <-
    purrr::map(as.data.frame(t(sample_cor)), function(x) {
      paste(sort(as.character(x)[7:8]), collapse = "_")
    })  %>%
    unlist() %>%
    unname()
  
  sample_wise_number <-
    purrr::map(as.data.frame(t(sample_cor)), function(x) {
      c(sort(as.character(x)[7:8]), name = paste(sort(as.character(x)[7:8]), collapse = "_"))
    }) %>%
    do.call(rbind, .) %>%
    as.data.frame() %>%
    dplyr::group_by(name) %>%
    dplyr::summarise(significant_n = n()) %>%
    dplyr::arrange(name) %>%
    dplyr::mutate(theoretical_n = NA)
  
  sample_wise_number$theoretical_n[sample_wise_number$name == "stool_stool"] =
    (intra_stool_sample_wise_dim[1] * (intra_stool_sample_wise_dim[1] - 1)) /
    2
  
  sample_wise_number <-
    sample_wise_number %>%
    dplyr::mutate(percentage = significant_n * 100 / theoretical_n)
  
  new_info =
    sample_wise_number$name %>%
    purrr::map(function(x) {
      significant_n_pos =
        sample_cor %>%
        dplyr::filter(class_name == x & cor > 0) %>%
        nrow()
      
      significant_n_neg =
        sample_cor %>%
        dplyr::filter(class_name == x & cor < 0) %>%
        nrow()
      
      data.frame(significant_n_pos = significant_n_pos,
                 significant_n_neg = significant_n_neg)
      
    }) %>%
    dplyr::bind_rows()
  
  sample_wise_number =
    cbind(sample_wise_number, new_info)
  
  sample_wise_number =
    sample_wise_number %>%
    dplyr::mutate(
      percentage_pos = significant_n_pos * 100 / theoretical_n,
      percentage_neg = significant_n_neg * 100 / theoretical_n
    )
  
  load(here::here(
    "data_analysis/stool_microbiome/data_preparation/variable_info"
  ))
  stool_variable_info = variable_info
  
  from <-
    sample_cor %>%
    dplyr::filter(class_name == "stool_stool") %>%
    pull(from) %>%
    stringr::str_replace("stool", "") %>%
    unique()
  
  temp <-
    stool_variable_info %>%
    dplyr::filter(variable_id %in% from)
  
  temp <-
    temp %>%
    dplyr::count(Phylum) %>%
    dplyr::mutate(n = n * 100 / sum(n))
  
  phylum_percentage_ratio[[i]] <-
    temp
  
  edge_data =
    sample_wise_number %>%
    dplyr::select(name, percentage) %>%
    tidyr::separate(col = name,
                    sep = "_",
                    into = c("from", "to")) %>%
    dplyr::mutate(from = stringr::str_to_title(from)) %>%
    dplyr::mutate(to = stringr::str_to_title(to))
  
  edge_percentage_ratio[[i]] <-
    edge_data
}

temp <-
  edge_percentage_ratio %>%
  do.call(rbind, .) %>%
  as.data.frame()

temp_stool <- temp

###wilcox test
wilcox.test(temp_oral$percentage, temp_skin$percentage)
wilcox.test(temp_oral$percentage, temp_nasal$percentage)
wilcox.test(temp_oral$percentage, temp_stool$percentage)

wilcox.test(temp_skin$percentage, temp_nasal$percentage)
wilcox.test(temp_skin$percentage, temp_stool$percentage)

wilcox.test(temp_nasal$percentage, temp_stool$percentage)
