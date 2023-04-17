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

masstools::setwd_project()
library(tidyverse)
rm(list = ls())

source("code/tools.R")

######work directory
masstools::setwd_project()

edge_percentage_ratio <- vector(mode = "list", length = 50)
phylum_percentage_ratio <- vector(mode = "list", length = 50)

for (i in 1:50) {
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

mean(temp$percentage)
sd(temp$percentage)
