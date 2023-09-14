#' ---
#' title: "stool microbiome oral_microbiome correlation"
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
dir.create(
  "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_oral_microbiome/permutation",
  recursive = TRUE
)
setwd(
  "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_oral_microbiome/permutation"
)

{
  load("../stool_microbiome_expression_data")
  load("../stool_microbiome_variable_info")
  load("../stool_microbiome_sample_info")
  
  load("../oral_microbiome_expression_data")
  load("../oral_microbiome_variable_info")
  load("../oral_microbiome_sample_info")
}

dim(stool_microbiome_expression_data)
dim(oral_microbiome_expression_data)

stool_oral_sample_wise_stool_dim =
  dim(stool_microbiome_expression_data)

stool_oral_sample_wise_oral_dim =
  dim(oral_microbiome_expression_data)

# save(stool_oral_sample_wise_stool_dim, file = "stool_oral_sample_wise_stool_dim")
# save(stool_oral_sample_wise_oral_dim, file = "stool_oral_sample_wise_oral_dim")


###finally, for stool microbiome, 106 genus, for oral_microbiome, 76 genus
######--------------------------------------------------------------------------
library(plyr)

#####
dim(oral_microbiome_expression_data)
dim(stool_microbiome_expression_data)

stool_microbiome_sample_info$subject_id == oral_microbiome_sample_info$subject_id

##https://rpubs.com/DKCH2020/578881
##https://ourcodingclub.github.io/tutorials/mixed-models/
##https://zhuanlan.zhihu.com/p/63092231
##https://www.linglab.cn/knowledge/10

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

library(compositions)
oral_microbiome_expression_data =
  oral_microbiome_expression_data %>%
  purrr::map(function(x) {
    x = compositions::clr(x) %>%
      as.numeric()
    x
  }) %>%
  do.call(cbind, .) %>%
  as.data.frame()

rownames(oral_microbiome_expression_data) = oral_microbiome_variable_info$variable_id


##step 1
###linear mixed model to adjust the subject ID random effect, and then use the partial correlation
###to get the correlation between microbiome and oral_microbiome
library(lme4)
library(rmcorr)

library(future)
library(furrr)

###permutation
for (i in 1:20) {
  cat(i, " ")
  sample_idx <-
    sample(
      1:ncol(stool_microbiome_expression_data),
      ncol(stool_microbiome_expression_data),
      replace = TRUE
    ) %>%
    unique() %>%
    sort()
  
  stool_microbiome_oral_microbiome_lm_adjusted_cor =
    lm_adjusted_cor(
      data_set1 = stool_microbiome_expression_data[, sample_idx],
      data_set2 = oral_microbiome_expression_data[, sample_idx],
      sample_info = stool_microbiome_sample_info[sample_idx,],
      method = "all",
      threads = 8
    )
  
  stool_microbiome_oral_microbiome_lm_adjusted_cor_spearman =
    stool_microbiome_oral_microbiome_lm_adjusted_cor[[1]]
  
  save(
    stool_microbiome_oral_microbiome_lm_adjusted_cor_spearman,
    file = "stool_microbiome_oral_microbiome_lm_adjusted_cor_spearman",
    compress = "xz"
  )
  
  save(
    stool_microbiome_oral_microbiome_lm_adjusted_cor_spearman,
    file = paste0(
      "stool_microbiome_oral_microbiome_lm_adjusted_cor_spearman_",
      i
    ),
    compress = "xz"
  )
  
}
