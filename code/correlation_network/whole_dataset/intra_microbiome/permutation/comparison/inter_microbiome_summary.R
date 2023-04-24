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

##nasal vs skin
load(
  "data_analysis/correlation_network/whole_data_set/skin_microbiome_vs_nasal_microbiome/permutation/temp"
)
nasal_skin_temp <- temp

##nasal vs stool
load(
  "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_nasal_microbiome/permutation/temp"
)
nasal_stool_temp <- temp

##nasal vs stool
load(
  "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_nasal_microbiome/permutation/temp"
)
nasal_stool_temp <- temp

##skin vs stool
load(
  "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_skin_microbiome/permutation/temp"
)
skin_stool_temp <- temp

##nasal vs oral
load(
  "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_nasal_microbiome/permutation/temp"
)
nasal_oral_temp <- temp

##oral vs skin
load(
  "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_skin_microbiome/permutation/temp"
)
oral_skin_temp <- temp

##oral vs stool
load(
  "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_oral_microbiome/permutation/temp"
)
oral_stool_temp <- temp

wilcox.test(nasal_skin_temp$percentage,
            nasal_stool_temp$percentage)

wilcox.test(nasal_skin_temp$percentage,
            skin_stool_temp$percentage)

wilcox.test(nasal_skin_temp$percentage,
            nasal_oral_temp$percentage)

wilcox.test(nasal_skin_temp$percentage,
            oral_skin_temp$percentage)

wilcox.test(nasal_skin_temp$percentage,
            oral_stool_temp$percentage)



wilcox.test(nasal_stool_temp$percentage,
            skin_stool_temp$percentage)

wilcox.test(nasal_stool_temp$percentage,
            nasal_oral_temp$percentage)

wilcox.test(nasal_stool_temp$percentage,
            oral_skin_temp$percentage)

wilcox.test(nasal_stool_temp$percentage,
            oral_stool_temp$percentage)





wilcox.test(skin_stool_temp$percentage,
            nasal_oral_temp$percentage)

wilcox.test(skin_stool_temp$percentage,
            oral_skin_temp$percentage)

wilcox.test(skin_stool_temp$percentage,
            oral_stool_temp$percentage)





wilcox.test(nasal_oral_temp$percentage,
            oral_skin_temp$percentage)

wilcox.test(nasal_oral_temp$percentage,
            oral_stool_temp$percentage)


wilcox.test(oral_skin_temp$percentage,
            oral_stool_temp$percentage)
