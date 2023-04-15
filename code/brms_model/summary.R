no_function()

masstools::setwd_project()
rm(list = ls())
source(here::here("code/tools.R"))

load("data_analysis/stool_microbiome/data_preparation/variable_info")
stool_microbiome_variable_info <-
  variable_info

load("data_analysis/skin_microbiome/data_preparation/variable_info")
skin_microbiome_variable_info <-
  variable_info

load("data_analysis/oral_microbiome/data_preparation/variable_info")
oral_microbiome_variable_info <-
  variable_info

load("data_analysis/nasal_microbiome/data_preparation/variable_info")
nasal_microbiome_variable_info <-
  variable_info

setwd("data_analysis/brms/")

stool_cytokine_beta <-
  readxl::read_xlsx("stool/result/stool_cytokine_beta.xlsx")

skin_cytokine_beta <-
  readxl::read_xlsx("skin/result/skin_cytokine_beta.xlsx")

oral_cytokine_beta <-
  readxl::read_xlsx("oral/result/oral_cytokine_beta.xlsx")

nasal_cytokine_beta <-
  readxl::read_xlsx("nasal/result/nasal_cytokine_beta.xlsx")

library(tidygraph)
library(ggraph)

edge_data <-
  rbind(stool_cytokine_beta,
        skin_cytokine_beta,
        oral_cytokine_beta,
        nasal_cytokine_beta)

edge_data <-
  edge_data %>%
  dplyr::rename(from = cytokine, to = genus)

edge_data$to <-
  paste(edge_data$class, edge_data$to, sep = "_")

edge_data <-
  edge_data %>%
  dplyr::filter(abs(b) > 0.005)

node_data <-
  data.frame(node = c(unique(edge_data$from),
                      unique(edge_data$to)),
             class = c(
               rep("Cytokine", length(unique(edge_data$from))),
               stringr::str_extract(unique(edge_data$to), "Stool|Skin|Oral|Nasal")
             )) %>%
  dplyr::mutate(true_name = stringr::str_replace(node, "Stool\\_|Skin\\_|Oral\\_|Nasal\\_", ""))


library(plyr)

node_data <-
  node_data %>%
  plyr::dlply(.variables = .(class)) %>%
  purrr::map(function(x) {
    if (x$class[1] == "Cytokine") {
      x$Phylum <- NA
    }
    
    if (x$class[1] == "Stool") {
      x <-
        x %>%
        dplyr::left_join(
          stool_microbiome_variable_info[, c("Genus", "Phylum")] %>%
            dplyr::distinct(.keep_all = TRUE),
          by = c("true_name" = "Genus")
        )
    }
    
    if (x$class[1] == "Skin") {
      x <-
        x %>%
        dplyr::left_join(
          skin_microbiome_variable_info[, c("Genus", "Phylum")] %>%
            dplyr::distinct(.keep_all = TRUE),
          by = c("true_name" = "Genus")
        )
    }
    
    if (x$class[1] == "Nasal") {
      x <-
        x %>%
        dplyr::left_join(
          nasal_microbiome_variable_info[, c("Genus", "Phylum")] %>%
            dplyr::distinct(.keep_all = TRUE),
          by = c("true_name" = "Genus")
        )
    }
    
    if (x$class[1] == "Oral") {
      x <-
        x %>%
        dplyr::left_join(
          oral_microbiome_variable_info[, c("Genus", "Phylum")] %>%
            dplyr::distinct(.keep_all = TRUE),
          by = c("true_name" = "Genus")
        )
    }
    x
    
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

node_data$Phylum[is.na(node_data$Phylum)] <-
  node_data$class[is.na(node_data$Phylum)]

##remove Archaea
remove_name =
  c(
    "Zea",
    stool_microbiome_variable_info %>%
      dplyr::filter(Kingdom == "Archaea") %>%
      dplyr::pull(Genus),
    skin_microbiome_variable_info %>%
      dplyr::filter(Kingdom == "Archaea") %>%
      dplyr::pull(Genus),
    oral_microbiome_variable_info %>%
      dplyr::filter(Kingdom == "Archaea") %>%
      dplyr::pull(Genus),
    nasal_microbiome_variable_info %>%
      dplyr::filter(Kingdom == "Archaea") %>%
      dplyr::pull(Genus)
  ) %>%
  unique()

if (length(remove_name) > 0) {
  node_data =
    node_data %>%
    dplyr::filter(!true_name %in% remove_name)
}

edge_data =
  edge_data %>%
  dplyr::filter(from %in% node_data$node & to %in% node_data$node)

node_data$Phylum[!node_data$Phylum %in% c("Actinobacteria",
                                          "Bacteroidetes",
                                          "Firmicutes",
                                          "Proteobacteria",
                                          "Cytokine")] <-
  "xOther"

range(edge_data$b)

library(tidyverse)

edge_data <-
  edge_data %>%
  dplyr::mutate(cor = case_when(b > 0 ~ "1",
                                b <= 0 ~ "-1"))

dim(edge_data)
dim(node_data)

# openxlsx::write.xlsx(node_data,
#                      file = "node_data.xlsx",
#                      asTable = TRUE,
#                      overwrite = TRUE)
#
# openxlsx::write.xlsx(edge_data,
#                      file = "node_data.xlsx",
#                      asTable = TRUE,
#                      overwrite = TRUE)


total_graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

library(igraph)
g <- total_graph
V(g)$type <- bipartite_mapping(g)$type

coords <-
  create_layout(g, layout = "bipartite") %>%
  dplyr::select(node, class, true_name, Phylum, x, y)

coords$y[coords$class == "Cytokine"] = 0
coords$y[coords$class == "Stool"] = -1
coords$y[coords$class == "Skin"] = -1
coords$y[coords$class == "Oral"] = 1
coords$y[coords$class == "Nasal"] = 1
table(coords$class)

coords <-
  coords %>%
  dplyr::rename(x = y, y = x)

range(coords$y)

table(coords$class)

###Stool
lm_trans <-
  lm(seq(
    from = 1,
    to = 105,
    length.out = 100
  ) ~ seq(
    from = 0,
    to = 275,
    length.out = 100
  ))
coefficients(lm_trans)[1]
coords$y[coords$class == "Stool"] <-
  coords$y[coords$class == "Stool"] * coefficients(lm_trans)[2] + coefficients(lm_trans)[1]

##Skin
lm_trans <-
  lm(seq(
    from = -60,
    to = -1,
    length.out = 100
  ) ~ seq(
    from = 0,
    to = 275,
    length.out = 100
  ))
coefficients(lm_trans)[1]
coords$y[coords$class == "Skin"] <-
  coords$y[coords$class == "Skin"] * coefficients(lm_trans)[2] + coefficients(lm_trans)[1]

###Oral
lm_trans <-
  lm(seq(
    from = 11,
    to = 105,
    length.out = 100
  ) ~ seq(
    from = 0,
    to = 275,
    length.out = 100
  ))
coefficients(lm_trans)[1]
coords$y[coords$class == "Oral"] <-
  coords$y[coords$class == "Oral"] * coefficients(lm_trans)[2] + coefficients(lm_trans)[1]

##Nasal
lm_trans <-
  lm(seq(
    from = -60,
    to = 10,
    length.out = 100
  ) ~ seq(
    from = 0,
    to = 275,
    length.out = 100
  ))
coefficients(lm_trans)[1]
coords$y[coords$class == "Nasal"] <-
  coords$y[coords$class == "Nasal"] * coefficients(lm_trans)[2] + coefficients(lm_trans)[1]

###cytokine
lm_trans <-
  lm(seq(
    from = -50 ,
    to = 90,
    length.out = 100
  ) ~ seq(
    from = 0,
    to = 275,
    length.out = 100
  ))
coefficients(lm_trans)[1]
coords$y[coords$class == "Cytokine"] <-
  coords$y[coords$class == "Cytokine"] * coefficients(lm_trans)[2] + coefficients(lm_trans)[1]

library(plyr)
coords2 = coords
coords =
  coords %>%
  plyr::dlply(.variables = .(class)) %>%
  purrr::map(function(x) {
    x =
      x %>%
      dplyr::arrange(desc(Phylum), y)
    x$y = seq(from = range(x$y)[1],
              to = range(x$y)[2],
              length.out = nrow(x))
    x
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

coords2$y <-
  coords$y[match(coords2$node, coords$node)]

# coords3 = coords2
# coords3 <-
#   coords3 %>%
#   plyr::dlply(.variables = .(class)) %>%
#   purrr::map(function(x) {
#     x <-
#       x %>%
#       dplyr::arrange(Phylum, y)
#     x$y = seq(from = range(x$y)[1],
#               to = range(x$y)[2],
#               length.out = nrow(x))
#     x
#   }) %>%
#   dplyr::bind_rows() %>%
#   as.data.frame()
#
# coords3$y <-
#   coords2$y[match(coords3$node, coords2$node)]

my_graph <-
  create_layout(
    graph = g,
    layout = "manual",
    x = coords2$x,
    y = coords2$y
    # node.position = coords
  )

plot1 <-
  ggraph(my_graph,
         layout = 'bipartite') +
  geom_edge_diagonal(
    aes(color = cor),
    show.legend = TRUE,
    width = 0.3,
    alpha = 0.3
  ) +
  geom_node_text(aes(
    x = x * 1.05,
    y = y * 1,
    hjust = ifelse(class == "Stool" |
                     class == "Skin", 1, 0),
    label = ifelse(class == "Cytokine", NA, true_name)
  ),
  size = 2) +
  shadowtext::geom_shadowtext(
    aes(
      x = x * 1.05,
      y = y * 1,
      hjust = ifelse(class == "Stool" | class == "Skin", 1, 0),
      label = ifelse(class == "Cytokine", true_name, NA)
    ) ,
    size = 3,
    bg.colour = 'white',
    colour = "black"
  ) +
  geom_node_point(aes(
    color = Phylum,
    size = Degree,
    show.legend = TRUE,
  ),
  alpha = 0.8) +
  scale_color_manual(
    values = c(
      phylum_color["Actinobacteria"],
      phylum_color["Bacteroidetes"],
      phylum_color["Firmicutes"],
      phylum_color["Proteobacteria"],
      omics_color["Cytokine"],
      "xOther" = "grey"
    )
  ) +
  guides(
    edge_color = ggraph::guide_edge_colorbar(title = "Beta r"),
    fill = guide_legend(
      title = "Class",
      override.aes = list(size = 4, linetype = "blank")
    ),
    size = guide_legend(title = "Degree", override.aes = list(linetype = 0))
  ) +
  ggraph::scale_edge_color_manual(values = c("1" = "red", "-1" = "darkblue")) +
  ggraph::scale_edge_width(range = c(0.2, 2)) +
  scale_size_continuous(range = c(1, 5)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

plot1

# ggsave(plot1, filename = "brms_network.pdf", width = 5, height = 13)
