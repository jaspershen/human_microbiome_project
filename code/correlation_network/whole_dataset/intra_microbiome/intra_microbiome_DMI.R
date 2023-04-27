


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

####load the dims
####stool with skin
{
  load(
    "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_skin_microbiome/stool_skin_sample_wise_skin_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_skin_microbiome/stool_skin_sample_wise_stool_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_skin_microbiome/subject_wise/stool_skin_subject_wise_skin_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_skin_microbiome/subject_wise/stool_skin_subject_wise_stool_dim"
  )
  
  ####stool with oral
  load(
    "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_oral_microbiome/stool_oral_sample_wise_stool_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_oral_microbiome/stool_oral_sample_wise_oral_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_oral_microbiome/subject_wise/stool_oral_subject_wise_oral_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_oral_microbiome/subject_wise/stool_oral_subject_wise_stool_dim"
  )
  
  ####stool with nasal
  load(
    "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_nasal_microbiome/stool_nasal_sample_wise_stool_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_nasal_microbiome/stool_nasal_sample_wise_nasal_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_nasal_microbiome/subject_wise/stool_nasal_subject_wise_nasal_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_nasal_microbiome/subject_wise/stool_nasal_subject_wise_stool_dim"
  )
  
  
  ####oral with nasal
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_nasal_microbiome/oral_nasal_sample_wise_oral_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_nasal_microbiome/oral_nasal_sample_wise_nasal_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_nasal_microbiome/subject_wise/oral_nasal_subject_wise_nasal_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_nasal_microbiome/subject_wise/oral_nasal_subject_wise_oral_dim"
  )
  
  ####oral with skin
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_skin_microbiome/oral_skin_sample_wise_oral_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_skin_microbiome/oral_skin_sample_wise_skin_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_skin_microbiome/subject_wise/oral_skin_subject_wise_skin_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_skin_microbiome/subject_wise/oral_skin_subject_wise_oral_dim"
  )
  
  ####oral with nasal
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_nasal_microbiome/oral_nasal_sample_wise_oral_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_nasal_microbiome/oral_nasal_sample_wise_nasal_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_nasal_microbiome/subject_wise/oral_nasal_subject_wise_nasal_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_nasal_microbiome/subject_wise/oral_nasal_subject_wise_oral_dim"
  )
  
  ####skin with nasal
  load(
    "data_analysis/correlation_network/whole_data_set/skin_microbiome_vs_nasal_microbiome/skin_nasal_sample_wise_skin_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/skin_microbiome_vs_nasal_microbiome/skin_nasal_sample_wise_nasal_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/skin_microbiome_vs_nasal_microbiome/subject_wise/skin_nasal_subject_wise_nasal_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/skin_microbiome_vs_nasal_microbiome/subject_wise/skin_nasal_subject_wise_skin_dim"
  )
  
  
  #####intra-stool microbiome
  load(
    "data_analysis/correlation_network/whole_data_set/intra_stool_microbiome/intra_stool_sample_wise_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/intra_skin_microbiome/intra_skin_sample_wise_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/intra_nasal_microbiome/intra_nasal_sample_wise_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/intra_oral_microbiome/intra_oral_sample_wise_dim"
  )
}

load(
  "data_analysis/combine_microbiome/distance/stool/personalized_score_permutation"
)
stool_personalized_score_permutation <-
  personalized_score_permutation

load(
  "data_analysis/combine_microbiome/distance/skin/personalized_score_permutation"
)
skin_personalized_score_permutation <-
  personalized_score_permutation

load(
  "data_analysis/combine_microbiome/distance/nasal/personalized_score_permutation"
)
nasal_personalized_score_permutation <-
  personalized_score_permutation

load(
  "data_analysis/combine_microbiome/distance/oral/personalized_score_permutation"
)
oral_personalized_score_permutation <-
  personalized_score_permutation

setwd(
  "data_analysis/correlation_network/whole_data_set/intra_microbiome_correlation"
)

load("sample_cor")

sum(sample_cor$from_class == sample_cor$to_class)
sum(sample_cor$from_class != sample_cor$to_class)

sample_cor$class_name =
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

sample_wise_number$theoretical_n[sample_wise_number$name == "nasal_oral"] =
  oral_nasal_sample_wise_oral_dim[1] * oral_nasal_sample_wise_nasal_dim[1]

sample_wise_number$theoretical_n[sample_wise_number$name == "nasal_skin"] =
  skin_nasal_sample_wise_skin_dim[1] * skin_nasal_sample_wise_nasal_dim[1]

sample_wise_number$theoretical_n[sample_wise_number$name == "nasal_stool"] =
  stool_nasal_sample_wise_stool_dim[1] * stool_nasal_sample_wise_nasal_dim[1]

sample_wise_number$theoretical_n[sample_wise_number$name == "oral_oral"] =
  (intra_oral_sample_wise_dim[1] * (intra_oral_sample_wise_dim[1] - 1)) /
  2

sample_wise_number$theoretical_n[sample_wise_number$name == "oral_skin"] =
  oral_skin_sample_wise_oral_dim[1] * oral_skin_sample_wise_skin_dim[1]

sample_wise_number$theoretical_n[sample_wise_number$name == "oral_stool"] =
  stool_oral_sample_wise_oral_dim[1] * stool_oral_sample_wise_stool_dim[1]

sample_wise_number$theoretical_n[sample_wise_number$name == "skin_skin"] =
  (intra_skin_sample_wise_dim[1] * (intra_skin_sample_wise_dim[1] - 1)) /
  2

sample_wise_number$theoretical_n[sample_wise_number$name == "skin_stool"] =
  stool_skin_sample_wise_stool_dim[1] * stool_skin_sample_wise_skin_dim[1]

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

{
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
  }


######phylum class distributation
###stool
dim(stool_variable_info)

from <-
  sample_cor %>%
  dplyr::filter(class_name == "stool_stool") %>%
  pull(from) %>%
  stringr::str_replace("stool", "") %>%
  unique()

to <-
  sample_cor %>%
  dplyr::filter(class_name == "stool_stool") %>%
  pull(to) %>%
  stringr::str_replace("stool", "") %>%
  unique()

temp <-
  stool_variable_info %>%
  dplyr::filter(variable_id %in% c(from, to))

temp %>%
  dplyr::count(Phylum) %>%
  ggplot(aes(x = "", y = n, fill = Phylum)) +
  geom_bar(stat = "identity",
           width = 1,
           color = "black") +
  scale_fill_manual(values = phylum_color) +
  coord_polar("y", start = 0) +
  theme_void()

dmi_in <-
  stool_personalized_score_permutation %>%
  dplyr::filter(genus %in% temp$Genus) %>%
  dplyr::filter(fc1 > 15 * fc1_sd)

dmi_out <-
  stool_personalized_score_permutation %>%
  dplyr::filter(!genus %in% temp$Genus) %>%
  dplyr::filter(fc1 > 15 * fc1_sd)

dmi_in$fc1
dmi_out$fc1

library(ggsignif)
library(ggbeeswarm)
plot <-
  rbind(data.frame(dmi_in, class = "In"),
        data.frame(dmi_out, class = "Out")) %>%
  ggplot(aes(class, fc1)) +
  geom_violin() +
  geom_beeswarm() +
  theme_bw() +
  labs(x = "", y = "DMI") +
  geom_signif(
    comparisons = list(c("In", "Out")),
    map_signif_level = TRUE,
    test = "t.test"
  )
plot
ggsave(plot,
       filename = "stool_stool_dmi.pdf",
       width = 7,
       height = 7)

###skin
dim(skin_variable_info)

from <-
  sample_cor %>%
  dplyr::filter(class_name == "skin_skin") %>%
  pull(from) %>%
  stringr::str_replace("skin", "") %>%
  unique()

to <-
  sample_cor %>%
  dplyr::filter(class_name == "skin_skin") %>%
  pull(to) %>%
  stringr::str_replace("skin", "") %>%
  unique()

temp <-
  skin_variable_info %>%
  dplyr::filter(variable_id %in% c(from, to))

plot <-
  temp %>%
  dplyr::count(Phylum) %>%
  ggplot(aes(x = "", y = n, fill = Phylum)) +
  geom_bar(stat = "identity",
           width = 1,
           color = "black") +
  scale_fill_manual(values = phylum_color) +
  coord_polar("y", start = 0) +
  theme_void()
plot

dmi_in <-
  skin_personalized_score_permutation %>%
  dplyr::filter(genus %in% temp$Genus) %>%
  dplyr::filter(fc1 > 15 * fc1_sd)

dmi_out <-
  skin_personalized_score_permutation %>%
  dplyr::filter(!genus %in% temp$Genus) %>%
  dplyr::filter(fc1 > 15 * fc1_sd)

dmi_in$fc1
dmi_out$fc1

library(ggsignif)

plot <-
  rbind(data.frame(dmi_in, class = "In"),
        data.frame(dmi_out, class = "Out")) %>%
  ggplot(aes(class, fc1)) +
  geom_violin() +
  geom_beeswarm() +
  theme_bw() +
  labs(x = "", y = "DMI") +
  geom_signif(
    comparisons = list(c("In", "Out")),
    map_signif_level = TRUE,
    test = "t.test"
  )
plot
ggsave(plot,
       filename = "skin_skin_dmi.pdf",
       width = 7,
       height = 7)


###oral
dim(oral_variable_info)

from <-
  sample_cor %>%
  dplyr::filter(class_name == "oral_oral") %>%
  pull(from) %>%
  stringr::str_replace("oral", "") %>%
  unique()

to <-
  sample_cor %>%
  dplyr::filter(class_name == "oral_oral") %>%
  pull(to) %>%
  stringr::str_replace("oral", "") %>%
  unique()

temp <-
  oral_variable_info %>%
  dplyr::filter(variable_id %in% c(from, to))

plot <-
  temp %>%
  dplyr::count(Phylum) %>%
  ggplot(aes(x = "", y = n, fill = Phylum)) +
  geom_bar(stat = "identity",
           width = 1,
           color = "black") +
  scale_fill_manual(values = phylum_color) +
  coord_polar("y", start = 0) +
  theme_void()
plot

dmi_in <-
  oral_personalized_score_permutation %>%
  dplyr::filter(genus %in% temp$Genus) %>%
  dplyr::filter(fc1 > 15 * fc1_sd)

dmi_out <-
  oral_personalized_score_permutation %>%
  dplyr::filter(!genus %in% temp$Genus) %>%
  dplyr::filter(fc1 > 15 * fc1_sd)

dmi_in$fc1
dmi_out$fc1

library(ggsignif)

plot <-
  rbind(data.frame(dmi_in, class = "In"),
        data.frame(dmi_out, class = "Out")) %>%
  ggplot(aes(class, fc1)) +
  geom_violin() +
  geom_beeswarm() +
  theme_bw() +
  labs(x = "", y = "DMI") +
  geom_signif(
    comparisons = list(c("In", "Out")),
    map_signif_level = TRUE,
    test = "t.test"
  )
plot
ggsave(plot,
       filename = "oral_oral_dmi.pdf",
       width = 7,
       height = 7)

###nasal
dim(nasal_variable_info)

from <-
  sample_cor %>%
  dplyr::filter(class_name == "nasal_nasal") %>%
  pull(from) %>%
  stringr::str_replace("nasal", "") %>%
  unique()

to <-
  sample_cor %>%
  dplyr::filter(class_name == "nasal_nasal") %>%
  pull(to) %>%
  stringr::str_replace("nasal", "") %>%
  unique()

temp <-
  nasal_variable_info %>%
  dplyr::filter(variable_id %in% c(from, to))

plot <-
  temp %>%
  dplyr::count(Phylum) %>%
  ggplot(aes(x = "", y = n, fill = Phylum)) +
  geom_bar(stat = "identity",
           width = 1,
           color = "black") +
  scale_fill_manual(values = phylum_color) +
  coord_polar("y", start = 0) +
  theme_void()
plot

dmi_in <-
  nasal_personalized_score_permutation %>%
  dplyr::filter(genus %in% temp$Genus) %>%
  dplyr::filter(fc1 > 15 * fc1_sd)

dmi_out <-
  nasal_personalized_score_permutation %>%
  dplyr::filter(!genus %in% temp$Genus) %>%
  dplyr::filter(fc1 > 15 * fc1_sd)

dmi_in$fc1
dmi_out$fc1

library(ggsignif)

plot <-
  rbind(data.frame(dmi_in, class = "In"),
        data.frame(dmi_out, class = "Out")) %>%
  ggplot(aes(class, fc1)) +
  geom_violin() +
  geom_beeswarm() +
  theme_bw() +
  labs(x = "", y = "DMI") +
  geom_signif(
    comparisons = list(c("In", "Out")),
    map_signif_level = TRUE,
    test = "t.test"
  )
plot
ggsave(plot,
       filename = "nasal_nasal_dmi.pdf",
       width = 7,
       height = 7)
