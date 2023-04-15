###this is used for cluster running

### data_analysis to the cluster labs/mpsynder/shenxt/human_microbiome_project/
####copy skin_microbiome to data_analysis folder
 
##2. this r script should be put in labs/mpsynder/shenxt/human_microbiome_project/
 
###library packages
 
list.of.packages <- c("tidyverse", "skimr", "brms")
library(tidyverse)
library(skimr)
library(brms)
library(masstools, lib.loc="./R_library")
source("code/tools.R")
 
###load data
###skin microbiome genus level
{
  load("data_analysis/skin_microbiome/data_preparation/sample_info")
 
  skin_microbiome_sample_info = sample_info
 
  skin_microbiome_expression_data =
    data.table::fread("data/from_xin/Genus Table/SK/Genus_SK.csv") %>%
    tibble::column_to_rownames(var = "SampleID") %>%
    dplyr::select(-c(V1:SubjectID)) %>%
    t() %>%
    as.data.frame()
 
  skin_microbiome_expression_data =
    skin_microbiome_expression_data[,skin_microbiome_sample_info$sample_id]
 
  skin_microbiome_variable_info =
    data.frame(
      variable_id = rownames(skin_microbiome_expression_data)
    )
 
  #####from relative to count.
  skin_microbiome_expression_data =
    round(skin_microbiome_expression_data * 20000)
}
 
skin_microbiome_variable_info$variable_id =
  stringr::str_replace(skin_microbiome_variable_info$variable_id, "\\/", "_")
 
rownames(skin_microbiome_expression_data) = skin_microbiome_variable_info$variable_id
 
cat(skin_microbiome_variable_info$variable_id[87])
cat("\n")
skin_microbiome_variable_info =
  skin_microbiome_variable_info %>%
  dplyr::filter(!stringr::str_detect(variable_id, "Unclassified_Bacteria"))
 
skin_microbiome_expression_data =
  skin_microbiome_expression_data[skin_microbiome_variable_info$variable_id,]
 
####load cytokine data
load("data_analysis/cytokine/data_preparation/expression_data")
load("data_analysis/cytokine/data_preparation/sample_info")
load("data_analysis/cytokine/data_preparation/variable_info")
 
cytokine_expression_data = t(expression_data)
cytokine_sample_info = sample_info
cytokine_variable_info = variable_info
 
cytokine_sample_info =
cytokine_sample_info %>%
  dplyr::rename(sample_id = SampleID,
                subject_id = SubjectID)
 
cytokine_sample_info =
  cytokine_sample_info %>%
  dplyr::filter(!stringr::str_detect(subject_id, "\\?"))
 
cytokine_expression_data =
  cytokine_expression_data[,cytokine_sample_info$sample_id]
 
###add days
library(plyr)
cytokine_sample_info =
cytokine_sample_info %>%
  dplyr::mutate(CollectionDate = as.Date(CollectionDate, format = "%m/%d/%y")) %>%
  plyr::dlply(.variables = .(subject_id)) %>%
  purrr::map(function(x){
    x =
    x %>%
      dplyr::arrange(CollectionDate)
   
    x %>%  
    dplyr::mutate(days = as.numeric(CollectionDate - x$CollectionDate[1]))
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()
 
cytokine_expression_data =
  cytokine_expression_data[,cytokine_sample_info$sample_id] 
 
###add days
library(plyr)
skin_microbiome_sample_info =
  skin_microbiome_sample_info %>%
  dplyr::rename(CollectionDate = Date) %>%
  dplyr::mutate(CollectionDate = as.Date(CollectionDate, format = "%m/%d/%y")) %>%
  plyr::dlply(.variables = .(subject_id)) %>%
  purrr::map(function(x){
    x =
      x %>%
      dplyr::arrange(CollectionDate)
   
    x %>%  
      dplyr::mutate(days = as.numeric(CollectionDate - x$CollectionDate[1]))
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()
 
skin_microbiome_expression_data =
  skin_microbiome_expression_data[,skin_microbiome_sample_info$sample_id] 
 
####match sample
dim(cytokine_sample_info)
dim(skin_microbiome_sample_info)
 
intersect_sample_id =
  intersect(cytokine_sample_info$sample_id,
            skin_microbiome_sample_info$sample_id)
 
length(intersect_sample_id)
 
cytokine_sample_info =
  cytokine_sample_info[match(intersect_sample_id, cytokine_sample_info$sample_id),]
 
cytokine_expression_data =
  cytokine_expression_data[,cytokine_sample_info$sample_id]
 
skin_microbiome_sample_info =
  skin_microbiome_sample_info[match(intersect_sample_id, skin_microbiome_sample_info$sample_id),]
 
skin_microbiome_expression_data =
  skin_microbiome_expression_data[,skin_microbiome_sample_info$sample_id]
 
# colnames(skin_microbiome_expression_data) == colnames(cytokine_expression_data)
 
##only remain the genus at least in 10% subjects
remain_idx =
  which(rowSums(skin_microbiome_expression_data) > 0)
 
skin_microbiome_expression_data = skin_microbiome_expression_data[remain_idx,]
skin_microbiome_variable_info = skin_microbiome_variable_info[remain_idx,,drop = FALSE]
 
remain_idx =
  skin_microbiome_expression_data %>%
  apply(1, function(x){
    sum(as.numeric(x) == 0) / ncol(skin_microbiome_expression_data)
  }) %>%
  `<`(0.9) %>%
  which()
 
length(remain_idx)
 
skin_microbiome_expression_data = skin_microbiome_expression_data[remain_idx,]
skin_microbiome_variable_info = skin_microbiome_variable_info[remain_idx,,drop = FALSE]
 
dim(skin_microbiome_sample_info)
length(unique(skin_microbiome_sample_info$subject_id))
dim(skin_microbiome_variable_info)
 
skin_microbiome_sample_info =
  skin_microbiome_sample_info %>%
  dplyr::filter(!stringr::str_detect(subject_id, "\\?"))
 
skin_microbiome_expression_data =
  skin_microbiome_expression_data[,skin_microbiome_sample_info$sample_id]
 
###transformation
cytokine_expression_data_arcsine =
  cytokine_expression_data %>%
  t() %>%
  as.data.frame() %>%
  purrr::map(function(x) {
    x = normalize(x)
    x = ifelse(x == 0, x + (min(x[x > 0]) / 2), x)
    x = ifelse(x == 1, x - (min(x[x > 0]) / 2), x)
    arcsine.trans = arcsine(x)
    arcsine.trans
    # logit.trans = logit.trans(x)
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()
 
cytokine_expression_data_logit =
  cytokine_expression_data %>%
  t() %>%
  as.data.frame() %>%
  purrr::map(function(x) {
    x = normalize(x)
    x = ifelse(x == 0, x + (min(x[x > 0]) / 2), x)
    x = ifelse(x == 1, x - (min(x[x > 0]) / 2), x)
    # arcsine.trans = arcsine(x)
    # arcsine.trans
    logit.trans = logit(x)
    logit.trans
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()
 
colnames(cytokine_expression_data_logit) =
  colnames(cytokine_expression_data_arcsine) =
  colnames(cytokine_expression_data)
 
skin_microbiome_expression_data_arcsine =
  skin_microbiome_expression_data %>%
  t() %>%
  as.data.frame() %>%
  purrr::map(function(x) {
    x = normalize(x)
    x = ifelse(x == 0, x + (min(x[x > 0]) / 2), x)
    x = ifelse(x == 1, x - (min(x[x > 0]) / 2), x)
    arcsine.trans = arcsine(x)
    arcsine.trans
    # logit.trans = logit.trans(x)
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()
 
skin_microbiome_expression_data_logit =
  skin_microbiome_expression_data %>%
  t() %>%
  as.data.frame() %>%
  purrr::map(function(x) {
    x = normalize(x)
    x = ifelse(x == 0, x + (min(x[x > 0]) / 2), x)
    x = ifelse(x == 1, x - (min(x[x > 0]) / 2), x)
    # arcsine.trans = arcsine(x)
    # arcsine.trans
    logit.trans = logit(x)
    logit.trans
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()
 
colnames(skin_microbiome_expression_data_logit) =
  colnames(skin_microbiome_expression_data_arcsine) =
  colnames(skin_microbiome_expression_data)
 
rownames(skin_microbiome_expression_data) = skin_microbiome_variable_info$variable_id
 
####----------------------------------------------------------------------------
####2. model
# Model
# This should be modified to include all of the cytokines and genera of interest
# The format of the data frame should have each cytokine and genus as a column
cytokines <- rownames(cytokine_expression_data)
##remove chex
remain_idx = which(!stringr::str_detect(cytokines, "CHEX"))
cytokines = cytokines[remain_idx]
rownames(skin_microbiome_expression_data) =
  rownames(skin_microbiome_expression_data) %>%
  stringr::str_replace("\\/", "_")
genera <- rownames(skin_microbiome_expression_data)
 
################################
# skin.cytokine.model.output = vector(mode = "list", length = length(cytokines))
# names(skin.cytokine.model.output) = cytokines
for(i in 37:42) {
  cat(i, " ")
  c = cytokines[i]
  result =
    purrr::map(genera, function(g) {
      f <- as.formula(paste(g, "~ IRIS *", c, "+ days + (1|subject_id)",
                            sep = " "))
     
      x =
        data.frame(
          x = as.numeric(cytokine_expression_data[c, ]),
          y = as.numeric(skin_microbiome_expression_data[g, ]),
          skin_microbiome_sample_info
        ) %>% 
        dplyr::filter(IRIS != "Unknown")
 
      colnames(x)[1:2] = c(c, g)
 
      brm <-
        brms::brm(f,
                  data = x,
                  sparse = TRUE,
                  family = negbinomial())
      save(brm, file = file.path(
        "data_analysis/brms/skin_IRIS",
        paste("brm_skin",c,g, sep = "_")
      ))
    })
}
