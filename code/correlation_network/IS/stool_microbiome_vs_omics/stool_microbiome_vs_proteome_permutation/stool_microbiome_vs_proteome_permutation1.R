#' ---
#' title: "Stool microbiome proteome correlation"
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
library(ggraph)
library(igraph)
library(tidygraph)

###example network
library(plyr)
rm(list = ls())

source("code/tools.R")

######work directory
masstools::setwd_project()
dir.create(
  "data_analysis/correlation_network/whole_data_set_IS/stool_microbiome_vs_proteome_IS_permutation"
)
setwd(
  "data_analysis/correlation_network/whole_data_set_IS/stool_microbiome_vs_proteome_IS_permutation"
)

####load data
{
  load("../stool_microbiome_vs_proteome_IS/stool_microbiome_expression_data")
  load("../stool_microbiome_vs_proteome_IS/stool_microbiome_variable_info")
  load("../stool_microbiome_vs_proteome_IS/stool_microbiome_sample_info")
  
  load("../stool_microbiome_vs_proteome_IS/proteome_expression_data")
  load("../stool_microbiome_vs_proteome_IS/proteome_variable_info")
  load("../stool_microbiome_vs_proteome_IS/proteome_sample_info")
}

dim(stool_microbiome_expression_data)
dim(proteome_expression_data)

######--------------------------------------------------------------------------
##proteome data have been preprocessed
library(plyr)

###because our microbiome are percentage data, so here we use the CTL method
library(compositions)
stool_microbiome_expression_data =
  stool_microbiome_expression_data %>%
  purrr::map(function(x) {
    x = compositions::clr(x) %>%
      as.numeric()
    x
  }) %>%
  do.call(cbind, .) %>%
  as.data.frame()

rownames(stool_microbiome_expression_data) = stool_microbiome_variable_info$variable_id

#####
dim(proteome_expression_data)
dim(stool_microbiome_expression_data)

stool_microbiome_sample_info$subject_id == proteome_sample_info$subject_id

##https://rpubs.com/DKCH2020/578881
##https://ourcodingclub.github.io/tutorials/mixed-models/
##https://zhuanlan.zhihu.com/p/63092231
##https://www.linglab.cn/knowledge/10

###step 1
###linear mixed model to adjust the subject ID random effect, and then use the partial correlation
###to get the correlation between microbiome and proteome
library(lme4)
library(rmcorr)

library(future)
library(furrr)

###permutation

for (i in 1:7) {
  cat(i, " ")
  idx <-
    sample(
      1:ncol(stool_microbiome_expression_data),
      ncol(stool_microbiome_expression_data),
      replace = TRUE
    ) %>%
    unique() %>%
    sort()
  
  stool_microbiome_proteome_lm_adjusted_cor =
    lm_adjusted_cor(
      data_set1 = stool_microbiome_expression_data[, idx],
      data_set2 = proteome_expression_data[, idx],
      sample_info = stool_microbiome_sample_info[idx, ],
      method = "all",
      threads = 8
    )
  
  stool_microbiome_proteome_lm_adjusted_cor_spearman <-
    stool_microbiome_proteome_lm_adjusted_cor[[1]]
  
  cor_data =
    stool_microbiome_proteome_lm_adjusted_cor_spearman %>%
    dplyr::filter(p_adjust < 0.2)
  
  edge_data =
    cor_data %>%
    dplyr::mutate(p = -log(p_adjust, 10)) %>%
    dplyr::select(microbiome, metabolite, p, p_adjust, cor) %>%
    dplyr::rename(from = microbiome,
                  to = metabolite,
                  cor = cor)
  
  proteome_variable_info$Lipid_Name[is.na(proteome_variable_info$Lipid_Name)] =
    proteome_variable_info$variable_id[is.na(proteome_variable_info$Lipid_Name)]
  
  node_data =
    data.frame(node = unique(c(edge_data$from, edge_data$to))) %>%
    dplyr::left_join(proteome_variable_info[, c("variable_id", "annotation")],
                     by = c("node" = "variable_id")) %>%
    dplyr::rename(true_name = annotation) %>%
    dplyr::mutate(class = case_when(
      is.na(true_name) ~ "Stool microbiome",
      !is.na(true_name) ~ "Lipidome"
    )) %>%
    dplyr::left_join(stool_microbiome_variable_info[, c("variable_id", "Genus")],
                     by = c("node" = "variable_id")) %>%
    dplyr::mutate(true_name =
                    case_when(is.na(true_name) ~ Genus,
                              !is.na(true_name) ~ true_name)) %>%
    dplyr::select(-Genus)
  
  node_data =
    node_data %>%
    dplyr::filter(!stringr::str_detect(true_name, "C[0-9]{1,2}H[0-9]{1,2}"))
  
  node_data$true_name = stringr::str_replace_all(node_data$true_name, '\"', "")
  
  edge_data =
    edge_data %>%
    dplyr::filter(from %in% node_data$node & to %in% node_data$node)
  
  node_data =
    node_data %>%
    dplyr::filter(node %in% edge_data$from | node %in% edge_data$to)
  
  ####output data
  edge_data$from_true_name =
    node_data$true_name[match(edge_data$from, node_data$node)]
  
  edge_data$to_true_name =
    node_data$true_name[match(edge_data$to, node_data$node)]
  
  edge_data =
    edge_data %>%
    dplyr::mutate(significance = case_when(
      p_adjust < 0.05 ~ "p.adj<0.05",
      p_adjust >= 0.05 ~ "0.05<p.adj<0.2",
    ))
  
  save(node_data, file = paste("node_data", i, sep = "_"))
  save(edge_data, file = paste("edge_data", i, sep = "_"))
}


# for (i in 1:12) {
#   load(paste("edge_data", i, sep = "_"))
#   cat(nrow(edge_data), " ")
# }
# 
# for (i in 1:12) {
#   load(paste("node_data", i, sep = "_"))
#   cat(nrow(node_data), " ")
# }
