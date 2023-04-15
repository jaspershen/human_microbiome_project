######--------------------------------------------------------------------------
######this code is used for cluster running

####1. copy the code, data, data_analysis to the cluster labs/mpsynder/shenxt/human_microbiome_project/
####copy oral_microbiome to data_analysis folder

##2. this r script should be put in labs/mpsynder/shenxt/human_microbiome_project/

###library packages
###

list.of.packages <- c("tidyverse", "skimr", "brms")
library(tidyverse)
library(skimr)
library(brms)
masstools::setwd_project()
source("code/tools.R")

###load data
###oral microbiome genus level
{
  load("data_analysis/oral_microbiome/data_preparation/sample_info")
  
  oral_microbiome_sample_info = sample_info
  
  oral_microbiome_expression_data =
    data.table::fread("data/from_xin/Genus Table/OR/Genus_OR.csv") %>%
    tibble::column_to_rownames(var = "SampleID") %>%
    dplyr::select(-c(V1:SubjectID)) %>%
    t() %>%
    as.data.frame()
  
  oral_microbiome_expression_data =
    oral_microbiome_expression_data[, oral_microbiome_sample_info$sample_id]
  
  oral_microbiome_variable_info =
    data.frame(variable_id = rownames(oral_microbiome_expression_data))
  
  #####from relative to count.
  oral_microbiome_expression_data =
    round(oral_microbiome_expression_data * 20000)
}

oral_microbiome_variable_info$variable_id =
  stringr::str_replace(oral_microbiome_variable_info$variable_id, "\\/", "_")

rownames(oral_microbiome_expression_data) = oral_microbiome_variable_info$variable_id

cat(oral_microbiome_variable_info$variable_id[87])
cat("\n")
oral_microbiome_variable_info =
  oral_microbiome_variable_info %>%
  dplyr::filter(!stringr::str_detect(variable_id, "Unclassified_Bacteria"))

oral_microbiome_expression_data =
  oral_microbiome_expression_data[oral_microbiome_variable_info$variable_id, ]

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
  cytokine_expression_data[, cytokine_sample_info$sample_id]

###add days
library(plyr)
cytokine_sample_info =
  cytokine_sample_info %>%
  dplyr::mutate(CollectionDate = as.Date(CollectionDate, format = "%m/%d/%y")) %>%
  plyr::dlply(.variables = .(subject_id)) %>%
  purrr::map(function(x) {
    x =
      x %>%
      dplyr::arrange(CollectionDate)
    
    x %>%
      dplyr::mutate(days = as.numeric(CollectionDate - x$CollectionDate[1]))
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

cytokine_expression_data =
  cytokine_expression_data[, cytokine_sample_info$sample_id]

###add days
library(plyr)
oral_microbiome_sample_info =
  oral_microbiome_sample_info %>%
  dplyr::rename(CollectionDate = Date) %>%
  dplyr::mutate(CollectionDate = as.Date(CollectionDate, format = "%m/%d/%y")) %>%
  plyr::dlply(.variables = .(subject_id)) %>%
  purrr::map(function(x) {
    x =
      x %>%
      dplyr::arrange(CollectionDate)
    
    x %>%
      dplyr::mutate(days = as.numeric(CollectionDate - x$CollectionDate[1]))
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

oral_microbiome_expression_data =
  oral_microbiome_expression_data[, oral_microbiome_sample_info$sample_id]

####match sample
dim(cytokine_sample_info)
dim(oral_microbiome_sample_info)

intersect_sample_id =
  intersect(cytokine_sample_info$sample_id,
            oral_microbiome_sample_info$sample_id)

length(intersect_sample_id)

cytokine_sample_info =
  cytokine_sample_info[match(intersect_sample_id, cytokine_sample_info$sample_id), ]

cytokine_expression_data =
  cytokine_expression_data[, cytokine_sample_info$sample_id]

oral_microbiome_sample_info =
  oral_microbiome_sample_info[match(intersect_sample_id, oral_microbiome_sample_info$sample_id), ]

oral_microbiome_expression_data =
  oral_microbiome_expression_data[, oral_microbiome_sample_info$sample_id]

# colnames(oral_microbiome_expression_data) == colnames(cytokine_expression_data)

##only remain the genus at least in 10% subjects
remain_idx =
  which(rowSums(oral_microbiome_expression_data) > 0)

oral_microbiome_expression_data = oral_microbiome_expression_data[remain_idx, ]
oral_microbiome_variable_info = oral_microbiome_variable_info[remain_idx, , drop = FALSE]

remain_idx =
  oral_microbiome_expression_data %>%
  apply(1, function(x) {
    sum(as.numeric(x) == 0) / ncol(oral_microbiome_expression_data)
  }) %>%
  `<`(0.9) %>%
  which()

length(remain_idx)

oral_microbiome_expression_data = oral_microbiome_expression_data[remain_idx, ]
oral_microbiome_variable_info = oral_microbiome_variable_info[remain_idx, , drop = FALSE]

dim(oral_microbiome_sample_info)
length(unique(oral_microbiome_sample_info$subject_id))
dim(oral_microbiome_variable_info)

oral_microbiome_sample_info =
  oral_microbiome_sample_info %>%
  dplyr::filter(!stringr::str_detect(subject_id, "\\?"))

oral_microbiome_expression_data =
  oral_microbiome_expression_data[, oral_microbiome_sample_info$sample_id]

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

oral_microbiome_expression_data_arcsine =
  oral_microbiome_expression_data %>%
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

oral_microbiome_expression_data_logit =
  oral_microbiome_expression_data %>%
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

colnames(oral_microbiome_expression_data_logit) =
  colnames(oral_microbiome_expression_data_arcsine) =
  colnames(oral_microbiome_expression_data)

rownames(oral_microbiome_expression_data) = oral_microbiome_variable_info$variable_id

####----------------------------------------------------------------------------
####2. model
# Model
# This should be modified to include all of the cytokines and genera of interest
# The format of the data frame should have each cytokine and genus as a column
cytokines <- rownames(cytokine_expression_data)
##remove chex
remain_idx = which(!stringr::str_detect(cytokines, "CHEX"))
cytokines = cytokines[remain_idx]
rownames(oral_microbiome_expression_data) =
  rownames(oral_microbiome_expression_data) %>%
  stringr::str_replace("\\/", "_")
genera <- rownames(oral_microbiome_expression_data)

length(cytokines)
length(genera)

62 * 64

all_name =
  lapply(cytokines, function(x) {
    paste(x, genera, sep = "_")
  }) %>%
  unlist()

all_name1 = all_name[1:400]
all_name2 = all_name[401:800]
all_name3 = all_name[801:1200]
all_name4 = all_name[1201:1600]
all_name5 = all_name[1601:2000]
all_name6 = all_name[2001:2400]
all_name7 = all_name[2401:2800]
all_name8 = all_name[2801:3200]
all_name9 = all_name[3201:3600]
all_name10 = all_name[3601:length(all_name)]

for (i in 1:length(all_name7)) {
  cat(i, " ")
  save(i, file = "index_for_7")
  c = stringr::str_split(all_name7[i], "_", n = 2)[[1]][1]
  g = stringr::str_split(all_name7[i], "_", n = 2)[[1]][2]
  
  f <- as.formula(paste(g, "~", c, "+ days + (1|subject_id)",
                        sep = " "))
  
  x =
    data.frame(
      x = as.numeric(cytokine_expression_data[c,]),
      y = as.numeric(oral_microbiome_expression_data[g,]),
      cytokine_sample_info
    )
  
  colnames(x)[1:2] = c(c, g)
  
  brm <-
    brms::brm(f,
              data = x,
              sparse = TRUE,
              family = negbinomial())
  
  save(brm, file = file.path("data_analysis/brms/oral",
                             paste("brm_oral", c, g, sep = "_")))
}
