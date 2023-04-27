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

{
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
}

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

#######between network
library(ggraph)
library(tidygraph)
library(igraph)

#####ggraph
#####sample wise
library(ggraph)
sample_wise_number

edge_data =
  sample_wise_number %>%
  dplyr::select(name, percentage) %>%
  tidyr::separate(col = name,
                  sep = "_",
                  into = c("from", "to")) %>%
  dplyr::mutate(from = stringr::str_to_title(from)) %>%
  dplyr::mutate(to = stringr::str_to_title(to))

node_data =
  data.frame(
    node = c("stool", "skin", "nasal", "oral"),
    true_name = c("stool", "skin", "nasal", "oral")
  ) %>%
  dplyr::mutate(node = stringr::str_to_title(node)) %>%
  dplyr::mutate(true_name = stringr::str_to_title(true_name))

#######within network
library(ggraph)
library(tidygraph)
library(igraph)

library(patchwork)

#######between network
library(ggraph)
library(tidygraph)
library(igraph)

edge_data1 =
  edge_data %>%
  dplyr::filter(from != to) %>%
  dplyr::mutate(from = paste(from, "from", sep = '_')) %>%
  dplyr::mutate(to = paste(to, "to", sep = '_')) %>%
  dplyr::arrange(desc(percentage))

edge_data1$from = masstools::name_duplicated(edge_data1$from)
edge_data1$to = masstools::name_duplicated(edge_data1$to)

node_data1 =
  data.frame(node = unique(c(edge_data1$from, edge_data1$to))) %>%
  dplyr::mutate(true_name = stringr::str_extract(node, "Stool|Skin|Nasal|Oral"))

temp_graph <-
  tidygraph::tbl_graph(nodes = node_data1,
                       edges = edge_data1,
                       directed = FALSE)
#####up-down
g <- temp_graph
library(igraph)
library(ggraph)

V(g)$type <- bipartite_mapping(g)$type

coords <-
  create_layout(g, layout = "bipartite") %>%
  dplyr::select(node, true_name, x, y)

coords$x[grep("from", coords$node)] = 0
coords$x[grep("to", coords$node)] = 1

coords$y[coords$node == "Nasal_from_1"] = 6
coords$y[coords$node == "Skin_to_1"] = 6

coords$y[coords$node == "Nasal_from_2"] = 5
coords$y[coords$node == "Stool_to_1"] = 5

coords$y[coords$node == "Skin_from"] = 4
coords$y[coords$node == "Stool_to_2"] = 4

coords$y[coords$node == "Nasal_from_3"] = 3
coords$y[coords$node == "Oral_to"] = 3

coords$y[coords$node == "Oral_from_1"] = 2
coords$y[coords$node == "Skin_to_2"] = 2

coords$y[coords$node == "Oral_from_2"] = 1
coords$y[coords$node == "Stool_to_3"] = 1

my_graph <-
  create_layout(
    graph = g,
    layout = "manual",
    x = coords$x,
    y = coords$y
    # node.position = coords
  )

plot =
  ggraph(my_graph,
         layout = 'bipartite') +
  geom_edge_link(
    aes(edge_width = percentage,
        label = paste(round(percentage, 2), "%", sep = "")),
    angle_calc = 'along',
    label_dodge = unit(5, 'mm'),
    show.legend = TRUE
  ) +
  ggraph::geom_node_point(
    aes(fill = true_name),
    shape = 21,
    size = 8,
    show.legend = TRUE
  ) +
  guides(fill = guide_legend(title = "Body site"),
         size = FALSE) +
  scale_fill_manual(values = body_site_color) +
  scale_edge_width(limits = c(1, 100), range = c(1, 20)) +
  ggraph::theme_graph() +
  theme(legend.position = "bottom")

plot

sample_cor2 =
  sample_cor %>%
  dplyr::mutate(
    from1 = to,
    to1 = from,
    from_class1 = to_class,
    to_class1 = from_class
  ) %>%
  dplyr::select(-c(from, to, from_class, to_class)) %>%
  dplyr::rename(
    from = from1,
    to = to1,
    from_class = from_class1,
    to_class = to_class1
  ) %>%
  dplyr::select(colnames(sample_cor))

sample_cor2$class_name =
  sample_cor2$class_name %>%
  stringr::str_split("_") %>%
  purrr::map(function(x) {
    paste(x[2:1], collapse = "_")
  }) %>%
  unlist()

temp_sample_cor = rbind(sample_cor,
                        sample_cor2)

#####nasal skin
nasal_skin <-
  temp_sample_cor %>%
  dplyr::filter(from_class == "nasal" & to_class == "skin") %>%
  dplyr::mutate(from = stringr::str_replace(from, "nasal", "")) %>%
  dplyr::mutate(to = stringr::str_replace(to, "skin", "")) %>%
  dplyr::left_join(nasal_variable_info[, c("variable_id", "Genus")],
                   by = c("from" = "variable_id")) %>%
  dplyr::rename(from_Genus = Genus) %>%
  dplyr::left_join(skin_variable_info[, c("variable_id", "Genus")],
                   by = c("to" = "variable_id")) %>%
  dplyr::rename(to_Genus = Genus)

dmi_in <-
  nasal_personalized_score_permutation %>%
  dplyr::filter(genus %in% nasal_skin$from_Genus) %>%
  dplyr::filter(fc1 > 15 * fc1_sd)

dmi_out <-
  nasal_personalized_score_permutation %>%
  dplyr::filter(!genus %in% nasal_skin$from_Genus) %>%
  dplyr::filter(fc1 > 15 * fc1_sd)

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
       filename = "nasal_vs_skin_nasal_dmi.pdf",
       width = 7,
       height = 7)





dmi_in <-
  nasal_personalized_score_permutation %>%
  dplyr::filter(genus %in% nasal_skin$to_Genus) %>%
  dplyr::filter(fc1 > 15 * fc1_sd)

dmi_out <-
  nasal_personalized_score_permutation %>%
  dplyr::filter(!genus %in% nasal_skin$to_Genus) %>%
  dplyr::filter(fc1 > 15 * fc1_sd)

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
       filename = "nasal_vs_skin_skin_dmi.pdf",
       width = 7,
       height = 7)










#####nasal stool
nasal_stool <-
  temp_sample_cor %>%
  dplyr::filter(from_class == "nasal" & to_class == "stool") %>%
  dplyr::mutate(from = stringr::str_replace(from, "nasal", "")) %>%
  dplyr::mutate(to = stringr::str_replace(to, "stool", "")) %>%
  dplyr::left_join(nasal_variable_info[, c("variable_id", "Genus")],
                   by = c("from" = "variable_id")) %>%
  dplyr::rename(from_Genus = Genus) %>%
  dplyr::left_join(stool_variable_info[, c("variable_id", "Genus")],
                   by = c("to" = "variable_id")) %>%
  dplyr::rename(to_Genus = Genus)

dmi_in <-
  nasal_personalized_score_permutation %>%
  dplyr::filter(genus %in% nasal_stool$from_Genus) %>%
  dplyr::filter(fc1 > 15 * fc1_sd)

dmi_out <-
  nasal_personalized_score_permutation %>%
  dplyr::filter(!genus %in% nasal_stool$from_Genus) %>%
  dplyr::filter(fc1 > 15 * fc1_sd)

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
       filename = "nasal_vs_stool_nasal_dmi.pdf",
       width = 7,
       height = 7)



dmi_in <-
  nasal_personalized_score_permutation %>%
  dplyr::filter(genus %in% nasal_stool$to_Genus) %>%
  dplyr::filter(fc1 > 15 * fc1_sd)

dmi_out <-
  nasal_personalized_score_permutation %>%
  dplyr::filter(!genus %in% nasal_stool$to_Genus) %>%
  dplyr::filter(fc1 > 15 * fc1_sd)

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
       filename = "nasal_vs_stool_stool_dmi.pdf",
       width = 7,
       height = 7)










#####skin stool
skin_stool <-
  temp_sample_cor %>%
  dplyr::filter(from_class == "skin" & to_class == "stool") %>%
  dplyr::mutate(from = stringr::str_replace(from, "skin", "")) %>%
  dplyr::mutate(to = stringr::str_replace(to, "stool", "")) %>%
  dplyr::left_join(skin_variable_info[, c("variable_id", "Genus")],
                   by = c("from" = "variable_id")) %>%
  dplyr::rename(from_Genus = Genus) %>%
  dplyr::left_join(stool_variable_info[, c("variable_id", "Genus")],
                   by = c("to" = "variable_id")) %>%
  dplyr::rename(to_Genus = Genus)

dmi_in <-
  skin_personalized_score_permutation %>%
  dplyr::filter(genus %in% skin_stool$from_Genus) %>%
  dplyr::filter(fc1 > 15 * fc1_sd)

dmi_out <-
  skin_personalized_score_permutation %>%
  dplyr::filter(!genus %in% skin_stool$from_Genus) %>%
  dplyr::filter(fc1 > 15 * fc1_sd)

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
       filename = "skin_vs_stool_skin_dmi.pdf",
       width = 7,
       height = 7)



dmi_in <-
  skin_personalized_score_permutation %>%
  dplyr::filter(genus %in% skin_stool$to_Genus) %>%
  dplyr::filter(fc1 > 15 * fc1_sd)

dmi_out <-
  skin_personalized_score_permutation %>%
  dplyr::filter(!genus %in% skin_stool$to_Genus) %>%
  dplyr::filter(fc1 > 15 * fc1_sd)

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
       filename = "skin_vs_stool_stool_dmi.pdf",
       width = 7,
       height = 7)









#####nasal oral
nasal_oral <-
  temp_sample_cor %>%
  dplyr::filter(from_class == "nasal" & to_class == "oral") %>%
  dplyr::mutate(from = stringr::str_replace(from, "nasal", "")) %>%
  dplyr::mutate(to = stringr::str_replace(to, "oral", "")) %>%
  dplyr::left_join(nasal_variable_info[, c("variable_id", "Genus")],
                   by = c("from" = "variable_id")) %>%
  dplyr::rename(from_Genus = Genus) %>%
  dplyr::left_join(oral_variable_info[, c("variable_id", "Genus")],
                   by = c("to" = "variable_id")) %>%
  dplyr::rename(to_Genus = Genus)

dmi_in <-
  nasal_personalized_score_permutation %>%
  dplyr::filter(genus %in% nasal_oral$from_Genus) %>%
  dplyr::filter(fc1 > 15 * fc1_sd)

dmi_out <-
  nasal_personalized_score_permutation %>%
  dplyr::filter(!genus %in% nasal_oral$from_Genus) %>%
  dplyr::filter(fc1 > 15 * fc1_sd)

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
       filename = "nasal_vs_oral_nasal_dmi.pdf",
       width = 7,
       height = 7)



dmi_in <-
  nasal_personalized_score_permutation %>%
  dplyr::filter(genus %in% nasal_oral$to_Genus) %>%
  dplyr::filter(fc1 > 15 * fc1_sd)

dmi_out <-
  nasal_personalized_score_permutation %>%
  dplyr::filter(!genus %in% nasal_oral$to_Genus) %>%
  dplyr::filter(fc1 > 15 * fc1_sd)

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
       filename = "nasal_vs_oral_oral_dmi.pdf",
       width = 7,
       height = 7)













#####oral skin
oral_skin <-
  temp_sample_cor %>%
  dplyr::filter(from_class == "oral" & to_class == "skin") %>%
  dplyr::mutate(from = stringr::str_replace(from, "oral", "")) %>%
  dplyr::mutate(to = stringr::str_replace(to, "skin", "")) %>%
  dplyr::left_join(oral_variable_info[, c("variable_id", "Genus")],
                   by = c("from" = "variable_id")) %>%
  dplyr::rename(from_Genus = Genus) %>%
  dplyr::left_join(skin_variable_info[, c("variable_id", "Genus")],
                   by = c("to" = "variable_id")) %>%
  dplyr::rename(to_Genus = Genus)

dmi_in <-
  oral_personalized_score_permutation %>%
  dplyr::filter(genus %in% oral_skin$from_Genus) %>%
  dplyr::filter(fc1 > 15 * fc1_sd)

dmi_out <-
  oral_personalized_score_permutation %>%
  dplyr::filter(!genus %in% oral_skin$from_Genus) %>%
  dplyr::filter(fc1 > 15 * fc1_sd)

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
       filename = "oral_vs_skin_oral_dmi.pdf",
       width = 7,
       height = 7)





dmi_in <-
  oral_personalized_score_permutation %>%
  dplyr::filter(genus %in% oral_skin$to_Genus) %>%
  dplyr::filter(fc1 > 15 * fc1_sd)

dmi_out <-
  oral_personalized_score_permutation %>%
  dplyr::filter(!genus %in% oral_skin$to_Genus) %>%
  dplyr::filter(fc1 > 15 * fc1_sd)

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
       filename = "oral_vs_skin_skin_dmi.pdf",
       width = 7,
       height = 7)













#####oral stool
oral_stool <-
  temp_sample_cor %>%
  dplyr::filter(from_class == "oral" & to_class == "stool") %>%
  dplyr::mutate(from = stringr::str_replace(from, "oral", "")) %>%
  dplyr::mutate(to = stringr::str_replace(to, "stool", "")) %>%
  dplyr::left_join(oral_variable_info[, c("variable_id", "Genus")],
                   by = c("from" = "variable_id")) %>%
  dplyr::rename(from_Genus = Genus) %>%
  dplyr::left_join(stool_variable_info[, c("variable_id", "Genus")],
                   by = c("to" = "variable_id")) %>%
  dplyr::rename(to_Genus = Genus)

dmi_in <-
  oral_personalized_score_permutation %>%
  dplyr::filter(genus %in% oral_stool$from_Genus) %>%
  dplyr::filter(fc1 > 15 * fc1_sd)

dmi_out <-
  oral_personalized_score_permutation %>%
  dplyr::filter(!genus %in% oral_stool$from_Genus) %>%
  dplyr::filter(fc1 > 15 * fc1_sd)

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
       filename = "oral_vs_stool_oral_dmi.pdf",
       width = 7,
       height = 7)





dmi_in <-
  oral_personalized_score_permutation %>%
  dplyr::filter(genus %in% oral_stool$to_Genus) %>%
  dplyr::filter(fc1 > 15 * fc1_sd)

dmi_out <-
  oral_personalized_score_permutation %>%
  dplyr::filter(!genus %in% oral_stool$to_Genus) %>%
  dplyr::filter(fc1 > 15 * fc1_sd)

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
       filename = "oral_vs_stool_stool_dmi.pdf",
       width = 7,
       height = 7)
