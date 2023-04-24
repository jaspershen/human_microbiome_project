#' ---
#' title: "oral microbiome skin_microbiome correlation"
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
  "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_skin_microbiome/permutation",
  recursive = TRUE
)
setwd(
  "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_skin_microbiome/permutation"
)

{
  load("../oral_microbiome_expression_data")
  load("../oral_microbiome_variable_info")
  load("../oral_microbiome_sample_info")
  
  load("../skin_microbiome_expression_data")
  load("../skin_microbiome_variable_info")
  load("../skin_microbiome_sample_info")
}

dim(oral_microbiome_expression_data)
dim(skin_microbiome_expression_data)

oral_skin_sample_wise_oral_dim =
  dim(oral_microbiome_expression_data)

oral_skin_sample_wise_skin_dim =
  dim(skin_microbiome_expression_data)

save(oral_skin_sample_wise_oral_dim, file = "oral_skin_sample_wise_oral_dim")
save(oral_skin_sample_wise_skin_dim, file = "oral_skin_sample_wise_skin_dim")


###finally, for oral microbiome, 106 genus, for skin_microbiome, 76 genus
######--------------------------------------------------------------------------
library(plyr)

#####
dim(skin_microbiome_expression_data)
dim(oral_microbiome_expression_data)

oral_microbiome_sample_info$subject_id == skin_microbiome_sample_info$subject_id

##https://rpubs.com/DKCH2020/578881
##https://ourcodingclub.github.io/tutorials/mixed-models/
##https://zhuanlan.zhihu.com/p/63092231
##https://www.linglab.cn/knowledge/10

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

library(compositions)
skin_microbiome_expression_data =
  skin_microbiome_expression_data %>%
  purrr::map(function(x) {
    x = compositions::clr(x) %>%
      as.numeric()
    x
  }) %>%
  do.call(cbind, .) %>%
  as.data.frame()

rownames(skin_microbiome_expression_data) = skin_microbiome_variable_info$variable_id


##step 1
###linear mixed model to adjust the subject ID random effect, and then use the partial correlation
###to get the correlation between microbiome and skin_microbiome
library(lme4)
library(rmcorr)

library(future)
library(furrr)

###permutation
for (i in 1:20) {
  cat(i, " ")
  sample_idx <-
    sample(
      1:ncol(oral_microbiome_expression_data),
      ncol(oral_microbiome_expression_data),
      replace = TRUE
    ) %>%
    unique() %>%
    sort()
  
  oral_microbiome_skin_microbiome_lm_adjusted_cor =
    lm_adjusted_cor(
      data_set1 = oral_microbiome_expression_data[, sample_idx],
      data_set2 = skin_microbiome_expression_data[, sample_idx],
      sample_info = oral_microbiome_sample_info[sample_idx,],
      method = "all",
      threads = 8
    )
  
  oral_microbiome_skin_microbiome_lm_adjusted_cor_spearman =
    oral_microbiome_skin_microbiome_lm_adjusted_cor[[1]]
  
  save(
    oral_microbiome_skin_microbiome_lm_adjusted_cor_spearman,
    file = "oral_microbiome_skin_microbiome_lm_adjusted_cor_spearman",
    compress = "xz"
  )
  
  save(
    oral_microbiome_skin_microbiome_lm_adjusted_cor_spearman,
    file = paste0(
      "oral_microbiome_skin_microbiome_lm_adjusted_cor_spearman_",
      i
    ),
    compress = "xz"
  )
  
}
