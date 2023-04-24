no_function()

list.of.packages <- c("tidyverse", "skimr", "brms")

##############################next code should be run in scg
###dir is /labs/mpsnyder/shenxt/human_microbiome_project/data_analysis/brms/
library(tidyverse)
library(skimr)
library(brms)

# x = dir() %>%
#   grep("brm_skin_", ., value = TRUE)
# 
# result = vector(mode = "list", length = length(x))
# 
# for(i in i:length(x)){
#   cat(i, " ")
#   name = x[i]
#   cytokine =
#     name %>%
#     stringr::str_replace("brm_skin_", "") %>%
#     stringr::str_split(pattern = "_", n = 2) %>%
#     unlist() %>%
#     `[`(1)
#   genus =
#     name %>%
#     stringr::str_replace("brm_skin_", "") %>%
#     stringr::str_split(pattern = "_", n = 2) %>%
#     unlist() %>%
#     `[`(2)
#   if(stringr::str_detect(name, "box admin")){
#     result[[i]] =
#       data.frame(cytokine = cytokine, genus = genus, estimate = NA, low = NA, high = NA)
#   }else{
#     load(name)
#     temp =
#       tryCatch(expr = broom.mixed::tidy(brm) %>%
#                  as.data.frame(), error = function(e){
#                    return(NULL)
#                  })
#     if(is.null(temp)){
#       estimate = NA
#       low = NA
#       high = NA
#       result[[i]] =
#         data.frame(cytokine = cytokine, genus = genus, estimate = unname(estimate), low = unname(low), high = unname(high))
#     }else{
#       estimate = temp$estimate[2]
#       low = temp$conf.low[2]
#       high = temp$conf.high[2]
#       result[[i]] =
#         data.frame(cytokine = cytokine, genus = genus, estimate = unname(estimate), low = unname(low), high = unname(high))
#     }
#   }
# 
# }
# 
# save(result, file = "result")

#####just copy result to local
##/Users/xiaotaoshen/Box/Xiaotao Shen's Files/human_microbiome_project/data_analysis/brms/skin/result

####next code are in local PC
setwd(masstools::get_project_wd())
setwd("data_analysis/brms/skin/result/")
rm(list = ls())

load("result")

names(result)
length(result)

###load variable information for variables
load(here::here("data_analysis/skin_microbiome/data_preparation/variable_info"))

final_table = 
  result %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

final_table =
  result %>%
  do.call(rbind, .) %>%
  as.data.frame()

####output result
library(openxlsx)
# openxlsx::write.xlsx(
#   final_table,
#   "skin_cytokine_result.xlsx",
#   asTable = TRUE,
#   overwrite = TRUE
# )
head(final_table)

final_table$estimate

sum(is.na(final_table$estimate))

final_table = 
final_table %>% 
  dplyr::filter(!is.na(estimate)) 

final_table1 =
  final_table %>%
  dplyr::select(-c(low, high)) %>%
  tidyr::pivot_wider(names_from = genus,
                     values_from =  estimate) %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var = "cytokine")

final_table2 =
  final_table %>%
  dplyr::mutate(sig = low * high) %>% 
  dplyr::mutate(significant = case_when(
    sig > 0 ~ 1,
    sig < 0 ~ -1
  )) %>% 
  dplyr::select(c(cytokine, genus, significant)) %>%
  tidyr::pivot_wider(names_from = genus,
                     values_from =  significant) %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var = "cytokine")

remove_idx = 
final_table1 %>% 
  apply(1, function(x){
    any(is.na(x))
  }) %>% 
  which()

if(length(remove_idx) > 0){
  final_table1 =
    final_table1[-remove_idx,]
  final_table2 =
    final_table2[-remove_idx,]
}

library(ComplexHeatmap)

# robust_dist = function(x, y) {
#   qx = quantile(x, c(0.1, 0.9))
#   qy = quantile(y, c(0.1, 0.9))
#   l = x > qx[1] & x < qx[2] & y > qy[1] & y < qy[2]
#   x = x[l]
#   y = y[l]
#   sqrt(sum((x - y)^2))
# }

phylum = 
variable_info$Phylum[match(colnames(final_table1), variable_info$Genus)]

phylum_color = ggsci::pal_nejm()(n=8)

names(phylum_color) = unique(phylum)

plot = 
Heatmap(
  matrix = as.matrix(final_table1),
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_column_dend = TRUE,
  show_row_dend = TRUE, 
  name = "Beta coefficient",
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D",
  clustering_distance_columns = "euclidean",
  clustering_method_columns = "ward.D",
  column_names_rot = 45,
  show_column_names = FALSE,
  cell_fun = function(j, i, x, y, w, h, fill){
    gb = textGrob("*")
    gb_w = convertWidth(grobWidth(gb), "mm")
    gb_h = convertHeight(grobHeight(gb), "mm")
    if(final_table2[i, j] == 1) {
      grid.text("*", x, y - gb_h*0.5 + gb_w*0.4)
    }
  },
  column_names_gp = gpar(cex = 0.7),
  row_names_gp = gpar(cex = 0.5),
  bottom_annotation = columnAnnotation(
    bacteria = anno_text(
      x = colnames(final_table1),
      show_name = FALSE,
      gp = gpar(col = phylum_color[phylum], cex = 0.5),
      rot = 45
    )
  )
) 

plot = ggplotify::as.ggplot(plot)
plot

# ggsave(plot,
#        filename = "brms_skin_microbiome_vs_cytokine.pdf",
#        width = 14,
#        height = 7)


result1 <- 
  final_table1 %>% 
  tibble::rownames_to_column(var = "cytokine") %>% 
  tidyr::pivot_longer(cols = -cytokine, names_to = "genus", values_to = "b")

result2 <- 
  final_table2 %>% 
  tibble::rownames_to_column(var = "cytokine") %>% 
  tidyr::pivot_longer(cols = -cytokine, names_to = "genus", values_to = "significant")

result = 
  result1 %>% 
  dplyr::left_join(result2, by = c("cytokine", "genus")) %>% 
  dplyr::filter(significant == 1)

result <- 
  data.frame(result, class = "Skin")

# openxlsx::write.xlsx(
#   result,
#   "skin_cytokine_beta.xlsx",
#   asTable = TRUE,
#   overwrite = TRUE
# )


