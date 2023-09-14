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

# ###load data
##stool vs skin
{
  load(
    "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_skin_microbiome/subject_wise/cor_data"
  )
  
  stool_skin_subject_cor = cor_data
  
  ##stool vs oral
  load(
    "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_oral_microbiome/subject_wise/cor_data"
  )
  stool_oral_subject_cor = cor_data
  
  ##stool vs nasal
  load(
    "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_nasal_microbiome/subject_wise/cor_data"
  )
  stool_nasal_subject_cor = cor_data
  
  ##oral vs nasal
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_nasal_microbiome/subject_wise/cor_data"
  )
  oral_nasal_subject_cor = cor_data
  
  ##oral vs skin
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_skin_microbiome/subject_wise/cor_data"
  )
  oral_skin_subject_cor = cor_data
  
  ##skin vs nasal
  load(
    "data_analysis/correlation_network/whole_data_set/skin_microbiome_vs_nasal_microbiome/subject_wise/cor_data"
  )
  skin_nasal_subject_cor = cor_data
}


#
#
####load the dims
####stool with skin
{
  load(
    "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_skin_microbiome/subject_wise/stool_skin_subject_wise_skin_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_skin_microbiome/subject_wise/stool_skin_subject_wise_stool_dim"
  )
  
  ####stool with oral
  load(
    "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_oral_microbiome/subject_wise/stool_oral_subject_wise_oral_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_oral_microbiome/subject_wise/stool_oral_subject_wise_stool_dim"
  )
  
  ####stool with nasal
  load(
    "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_nasal_microbiome/subject_wise/stool_nasal_subject_wise_nasal_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_nasal_microbiome/subject_wise/stool_nasal_subject_wise_stool_dim"
  )
  
  
  ####oral with nasal
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_nasal_microbiome/subject_wise/oral_nasal_subject_wise_nasal_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_nasal_microbiome/subject_wise/oral_nasal_subject_wise_oral_dim"
  )
  
  ####oral with skin
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_skin_microbiome/subject_wise/oral_skin_subject_wise_skin_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_skin_microbiome/subject_wise/oral_skin_subject_wise_oral_dim"
  )
  
  ####oral with nasal
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_nasal_microbiome/subject_wise/oral_nasal_subject_wise_nasal_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_nasal_microbiome/subject_wise/oral_nasal_subject_wise_oral_dim"
  )
  
  ####skin with nasal
  load(
    "data_analysis/correlation_network/whole_data_set/skin_microbiome_vs_nasal_microbiome/subject_wise/skin_nasal_subject_wise_nasal_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/skin_microbiome_vs_nasal_microbiome/subject_wise/skin_nasal_subject_wise_skin_dim"
  )
  setwd(
    "data_analysis/correlation_network/whole_data_set/intra_microbiome_correlation"
  )
}

# ######subject wise
# subject_cor =
#   rbind(stool_skin_subject_cor,
#         stool_nasal_subject_cor,
#         stool_oral_subject_cor,
#         oral_nasal_subject_cor,
#         oral_skin_subject_cor,
#         skin_nasal_subject_cor) %>% 
#   dplyr::filter(p_value_adjust < 0.05)
# 
# name = apply(subject_cor, 1, function(x){
#   paste(sort(as.character(x[1:2])), collapse = "_")
# })
# 
# subject_cor =
# subject_cor %>%
#   dplyr::mutate(name = name) %>%
#   dplyr::distinct(name, .keep_all = TRUE)
# 
# 
# #
# save(subject_cor, file = "subject_cor")

load("subject_cor")

sum(subject_cor$p_value_adjust < 0.05)

subject_cor <-
  subject_cor %>% 
  dplyr::filter(p_value_adjust < 0.05)

sum(subject_cor$from_class == subject_cor$to_class)
sum(subject_cor$from_class != subject_cor$to_class)

subject_cor$class_name =
  purrr::map(as.data.frame(t(subject_cor)), function(x) {
    paste(sort(as.character(x)[4:5]), collapse = "_")
  })  %>%
  unlist() %>%
  unname()

subject_wise_number =
  purrr::map(as.data.frame(t(subject_cor)), function(x) {
    c(sort(as.character(x)[4:5]), name = paste(sort(as.character(x)[4:5]), collapse = "_"))
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::group_by(name) %>%
  dplyr::summarise(significant_n = n()) %>%
  dplyr::arrange(name) %>%
  dplyr::mutate(theoretical_n = NA)

subject_wise_number$theoretical_n[subject_wise_number$name == "nasal_nasal"] =
  (stool_nasal_subject_wise_nasal_dim[1] * (stool_nasal_subject_wise_nasal_dim[1] -
                                              1)) / 2

subject_wise_number$theoretical_n[subject_wise_number$name == "nasal_oral"] =
  oral_nasal_subject_wise_oral_dim[1] * oral_nasal_subject_wise_nasal_dim[1]

subject_wise_number$theoretical_n[subject_wise_number$name == "nasal_skin"] =
  skin_nasal_subject_wise_skin_dim[1] * skin_nasal_subject_wise_nasal_dim[1]

subject_wise_number$theoretical_n[subject_wise_number$name == "nasal_stool"] =
  stool_nasal_subject_wise_stool_dim[1] * stool_nasal_subject_wise_nasal_dim[1]

subject_wise_number$theoretical_n[subject_wise_number$name == "oral_oral"] =
  (stool_oral_subject_wise_oral_dim[1] * (stool_oral_subject_wise_oral_dim[1] -
                                            1)) / 2

subject_wise_number$theoretical_n[subject_wise_number$name == "oral_skin"] =
  oral_skin_subject_wise_oral_dim[1] * oral_skin_subject_wise_skin_dim[1]

subject_wise_number$theoretical_n[subject_wise_number$name == "oral_stool"] =
  stool_oral_subject_wise_oral_dim[1] * stool_oral_subject_wise_stool_dim[1]

subject_wise_number$theoretical_n[subject_wise_number$name == "skin_skin"] =
  (stool_skin_subject_wise_skin_dim[1] * (stool_skin_subject_wise_skin_dim[1] - 1)) /
  2

subject_wise_number$theoretical_n[subject_wise_number$name == "skin_stool"] =
  stool_skin_subject_wise_stool_dim[1] * stool_skin_subject_wise_skin_dim[1]

subject_wise_number$theoretical_n[subject_wise_number$name == "stool_stool"] =
  (stool_skin_subject_wise_stool_dim[1] * (stool_skin_subject_wise_stool_dim[1] - 1)) /
  2

subject_wise_number <-
  subject_wise_number %>%
  dplyr::mutate(percentage = significant_n * 100 / theoretical_n)

new_info =
  subject_wise_number$name %>%
  purrr::map(function(x) {
    significant_n_pos =
      subject_cor %>%
      dplyr::filter(class_name == x & cor > 0) %>%
      nrow()
    
    significant_n_neg =
      subject_cor %>%
      dplyr::filter(class_name == x & cor < 0) %>%
      nrow()
    
    data.frame(significant_n_pos = significant_n_pos,
               significant_n_neg = significant_n_neg)
    
  }) %>%
  dplyr::bind_rows()

subject_wise_number =
  cbind(subject_wise_number, new_info)

subject_wise_number <-
  subject_wise_number %>%
  dplyr::mutate(
    percentage_pos = significant_n_pos * 100 / theoretical_n,
    percentage_neg = significant_n_neg * 100 / theoretical_n
  )

load(here::here(
  "data_analysis/stool_microbiome/data_preparation/variable_info"
))
stool_variable_info = variable_info

load(here::here(
  "data_analysis/skin_microbiome/data_preparation/variable_info"
))
skin_variable_info = variable_info

load(here::here(
  "data_analysis/nasal_microbiome/data_preparation/variable_info"
))
nasal_variable_info = variable_info

load(here::here(
  "data_analysis/oral_microbiome/data_preparation/variable_info"
))
oral_variable_info = variable_info

######phylum class distributation
###stool
dim(stool_variable_info)
###nasal
dim(nasal_variable_info)

new_stool_variable_info <-
  stool_variable_info %>%
  dplyr::mutate(variable_id = paste0("stool", variable_id),
                class = "stool")

new_skin_variable_info <-
  skin_variable_info %>%
  dplyr::mutate(variable_id = paste0("skin", variable_id),
                class = "skin")

new_oral_variable_info <-
  oral_variable_info %>%
  dplyr::mutate(variable_id = paste0("oral", variable_id),
                class = "oral")

new_nasal_variable_info <-
  nasal_variable_info %>%
  dplyr::mutate(variable_id = paste0("nasal", variable_id),
                class = "nasal")

new_variable_info <-
  rbind(
    new_stool_variable_info,
    new_skin_variable_info,
    new_oral_variable_info,
    new_nasal_variable_info
  )

subject_cor_output <-
  subject_cor %>%
  dplyr::filter(from != to)

temp <-
  unique(c(subject_cor_output$from, subject_cor_output$to))

temp_stool <- grep("stool", temp, value = TRUE)
temp_skin <- grep("skin", temp, value = TRUE)
temp_oral <- grep("oral", temp, value = TRUE)
temp_nasal <- grep("nasal", temp, value = TRUE)

new_variable_info$Genus[match(temp_stool, new_variable_info$variable_id)]
unique(new_variable_info$Genus[match(temp_stool, new_variable_info$variable_id)])
new_variable_info$Genus[match(temp_skin, new_variable_info$variable_id)]
unique(new_variable_info$Genus[match(temp_skin, new_variable_info$variable_id)])
new_variable_info$Genus[match(temp_oral, new_variable_info$variable_id)]
unique(new_variable_info$Genus[match(temp_oral, new_variable_info$variable_id)])
new_variable_info$Genus[match(temp_nasal, new_variable_info$variable_id)]
unique(new_variable_info$Genus[match(temp_nasal, new_variable_info$variable_id)])

subject_cor_output <-
  subject_cor_output %>%
  dplyr::left_join(new_variable_info %>% dplyr::select(-class),
                   by = c("from" = "variable_id")) %>%
  dplyr::rename(
    from_Kingdom = Kingdom,
    from_Phylum = Phylum,
    from_Class = Class,
    from_Order = Order,
    from_Family = Family,
    from_Genus = Genus,
    from_Species = Species
  ) %>%
  dplyr::left_join(new_variable_info %>% dplyr::select(-class),
                   by = c("to" = "variable_id")) %>%
  dplyr::rename(
    to_Kingdom = Kingdom,
    to_Phylum = Phylum,
    to_Class = Class,
    to_Order = Order,
    to_Family = Family,
    to_Genus = Genus,
    to_Species = Species
  )

subject_inter <-
  subject_cor_output %>%
  dplyr::filter(from_class != to_class)

subject_intra <-
  subject_cor_output %>%
  dplyr::filter(from_class == to_class)

write.csv(subject_inter, file = "subject_inter.csv", row.names = FALSE)
write.csv(subject_intra, file = "subject_intra.csv", row.names = FALSE)

sum(subject_inter$from_Genus == subject_inter$to_Genus)
sum(subject_intra$from_Genus == subject_intra$to_Genus)

