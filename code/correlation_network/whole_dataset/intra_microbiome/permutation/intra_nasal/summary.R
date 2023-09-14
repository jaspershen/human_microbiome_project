#' ---
#' title: "oral microbiome nasal_microbiome correlation"
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
    intra_nasal_microbiome_lm_adjusted_cor %>%
    dplyr::filter(microbiome != metabolite)
  
  load(
    "data_analysis/correlation_network/whole_data_set/intra_nasal_microbiome/permutation/intra_nasal_sample_wise_dim"
  )
  
  #####sample wise
  nasal_sample_cor <-
    nasal_sample_cor %>%
    dplyr::rename(from = microbiome, to = metabolite) %>%
    dplyr::mutate(from_class = "nasal", to_class = "nasal") %>%
    dplyr::mutate(from = paste("nasal", from, sep = ""),
                  to = paste("nasal", to, sep = ""),) %>%
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

mean(temp$percentage)
sd(temp$percentage)
