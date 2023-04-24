###https://biofam.github.io/MOFA2/index.html
##https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-020-02015-1.pdf
# no_source()

# set work directory
setwd(masstools::get_project_wd())
library(tidyverse)
library(phyloseq)
rm(list = ls())

source("code/tools.R")

metadata <-
  readr::read_csv("Figures/Figure6/metadata/metadata.subject.csv")

###plasma proteome
{
  load(here::here(
    "data_analysis/proteome/data_preparation/expression_data"
  ))
  load(here::here("data_analysis/proteome/data_preparation/sample_info"))
  load(here::here(
    "data_analysis/proteome/data_preparation/variable_info"
  ))
}

proteome_expression_data = expression_data
proteome_sample_info = sample_info
proteome_variable_info = variable_info

proteome_sample_info <-
  proteome_sample_info %>% 
  dplyr::mutate(Date = CollectionDate)

proteome_sample_info$Date =
  as.Date(proteome_sample_info$Date, "%m/%d/%y")

dim(proteome_expression_data)
length(unique(proteome_sample_info$subject_id))

#######work directory
setwd(masstools::get_project_wd())
dir.create("data_analysis/proteome/Ant_clustering_IR",
           recursive = TRUE)
setwd("data_analysis/proteome/Ant_clustering_IR")

head(proteome_sample_info)
dim(proteome_expression_data)

###only remain IR
proteome_sample_info <-
  proteome_sample_info %>%
  # dplyr::left_join(metadata[, c("SubjectID", "IRIS")],
  #                  by = c("subject_id" = "SubjectID")) %>%
  dplyr::filter(!is.na(IRIS)) %>%
  dplyr::filter(IRIS == "IR")

proteome_expression_data <-
  proteome_expression_data[, proteome_sample_info$sample_id]

###only remain the infectional and pro and post two samples
infection_idx <-
  which(proteome_sample_info$CL4 == "Ant")

#pre-healthy (–H) state (healthy baselines within 186 days before the event’s onset),
#event early (EE) state (visits on days 1–6 of the event),
#event late (EL) state (visits on days 7–14 since the event started),
#recovery (RE) state (visits within days 15–40 since the event started),
#and finally post-healthy (+H) state (visits within 186 days after the event;

all_index <-
  infection_idx %>%
  purrr::map(function(i) {
    temp_info <-
      proteome_sample_info[i,]
    idx <-
      which(proteome_sample_info$subject_id == temp_info$subject_id)
    if (length(idx) < 5) {
      return(
        list(
          negtaive_h_id = character(0),
          ee_id = character(0),
          el_id = character(0),
          re_id = character(0),
          positive_h_id = character(0)
        )
      )
    }
    temp_sample_info <- proteome_sample_info[idx,]
    day_diff <-
      as.numeric(temp_sample_info$Date - proteome_sample_info[i,]$Date)
    
    negtaive_h_idx <-
      which(day_diff < -6 & day_diff >= -180)
    
    ee_idx <-
      which(day_diff <= 0 & day_diff >= -6)
    
    el_idx <-
      which(day_diff >= 1 & day_diff <= 7)
    
    re_idx <-
      which(day_diff >= 8 & day_diff <= 32)
    
    positive_h_idx <-
      which(day_diff >= 33 & day_diff <= 180)
    
    list(
      negtaive_h_id = temp_sample_info$sample_id[negtaive_h_idx],
      ee_id = temp_sample_info$sample_id[ee_idx],
      el_id = temp_sample_info$sample_id[el_idx],
      re_id = temp_sample_info$sample_id[re_idx],
      positive_h_id = temp_sample_info$sample_id[positive_h_idx]
    )
  })

####remove some index
temp <-
  all_index %>%
  purrr::map(function(x) {
    unlist(lapply(x, length))
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  apply(1, function(x) {
    all(x > 0)
  }) %>%
  which()

all_index <-
  all_index[temp]

###remove some duplicated index
all_ee_id <-
  all_index %>%
  lapply(function(x) {
    x$ee_id
  }) %>%
  unlist()

duplicated_ee_id <-
  unique(all_ee_id[which(duplicated(all_ee_id))])

if (length(duplicated_ee_id) > 0) {
  remove_idx <-
    duplicated_ee_id %>%
    purrr::map(function(id) {
      all_index %>%
        lapply(function(x) {
          id %in% x$ee_id
        }) %>%
        unlist() %>%
        which() %>%
        `[`(1)
    }) %>%
    unlist() %>%
    unique()
  
  all_index <-
    all_index[-remove_idx]
  
}

new_expression_data <-
  all_index %>%
  purrr::map(function(idx) {
    negtaive_h <-
      proteome_expression_data[, idx$negtaive_h_id, drop = FALSE] %>%
      apply(1, mean)
    
    ee <-
      proteome_expression_data[, idx$ee_id, drop = FALSE] %>%
      apply(1, mean)
      # `-`(negtaive_h)
    
    el <-
      proteome_expression_data[, idx$el_id, drop = FALSE] %>%
      apply(1, mean)
      # `-`(negtaive_h)
    
    re <-
      proteome_expression_data[, idx$re_id, drop = FALSE] %>%
      apply(1, mean)
      # `-`(negtaive_h)
    
    positive_h <-
      proteome_expression_data[, idx$positive_h_id, drop = FALSE] %>%
      apply(1, mean)
      # `-`(negtaive_h)
    
    temp <-
      cbind(negtaive_h, ee, el, re, positive_h) %>%
      as.data.frame()
    
    subject_id <-
      idx$ee_id[1]
    
    colnames(temp) <-
      paste(colnames(temp), subject_id, sep = ".")
    temp
  }) %>%
  do.call(cbind, .) %>%
  as.data.frame()

colnames(new_expression_data)

new_sample_info <-
  data.frame(
    sample_id = colnames(new_expression_data),
    class = stringr::str_split(colnames(new_expression_data), "\\.") %>%
      lapply(function(x)
        x[1]) %>%
      unlist,
    sample_id = stringr::str_split(colnames(new_expression_data), "\\.") %>%
      lapply(function(x)
        x[2]) %>%
      unlist
  )

new_variable_info <-
  proteome_variable_info

sum(is.na(new_expression_data))

save(new_expression_data, file = "new_expression_data")
save(new_variable_info, file = "new_variable_info")
save(new_sample_info, file = "new_sample_info")
# 
# temp_data <-
#   unique(new_sample_info$class) %>%
#   purrr::map(function(x) {
#     new_expression_data[, which(new_sample_info$class == x)] %>%
#       apply(1, mean)
#   }) %>%
#   do.call(cbind, .) %>%
#   as.data.frame()
# 
# colnames(temp_data) <- unique(new_sample_info$class)
# 
# ###remove the tax that all 0 in all the samples
# idx <-
#   apply(temp_data, 1, function(x) {
#     sum(x == 0) / ncol(temp_data)
#   }) %>%
#   `<`(0.5) %>%
#   which()
# 
# temp_data <-
#   temp_data[idx, ]
# 
# new_variable_info <-
#   new_variable_info[idx, ]
# 
# rownames(temp_data) == new_variable_info$variable_id
# 
# ###clustering
# library(Mfuzz)
# #first get the time point data together:
# time <- c(1:5)
# 
# temp_data2 <-
#   apply(temp_data, 1, function(x) {
#     (x - mean(x)) / sd(x)
#   }) %>%
#   t() %>%
#   as.data.frame()
# 
# temp_data <- rbind(time, temp_data)
# 
# row.names(temp_data)[1] <- "time"
# 
# #save it to a temp file so ti doesnt clutter up my blog directory
# write.table(
#   temp_data,
#   file = "temp_data.txt",
#   sep = '\t',
#   quote = FALSE,
#   col.names = NA
# )
# 
# data <- table2eset(filename = "temp_data.txt")
# 
# data.s <- standardise(data)
# 
# m1 <- mestimate(data.s)
# m1
# 
# plot <-
#   Dmin(
#     data.s,
#     m = m1,
#     crange = seq(2, 20, 1),
#     repeats = 3,
#     visu = TRUE
#   )
# 
# plot <-
#   plot %>%
#   data.frame(distance = plot,
#              k = seq(2, 20, 1)) %>%
#   ggplot(aes(k, distance)) +
#   geom_point(shape = 21,
#              size = 4,
#              fill = "black") +
#   geom_smooth() +
#   geom_segment(aes(
#     x = k,
#     y = 0,
#     xend = k,
#     yend = distance
#   )) +
#   theme_bw() +
#   theme(
#     # legend.position = c(0, 1),
#     # legend.justification = c(0, 1),
#     panel.grid = element_blank(),
#     axis.title = element_text(size = 13),
#     axis.text = element_text(size = 12),
#     panel.background = element_rect(fill = "transparent", color = NA),
#     plot.background = element_rect(fill = "transparent", color = NA),
#     legend.background = element_rect(fill = "transparent", color = NA)
#   ) +
#   labs(x = "Cluster number",
#        y = "Min. centroid distance") +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
# 
# plot
# 
# ggsave(plot,
#        filename = "distance_k_number.pdf",
#        width = 7,
#        height = 7)
# 
# clust = 2
# 
# c <- mfuzz(data.s, c = clust, m = m1)
# 
# mfuzz.plot(
#   eset = data.s,
#   min.mem = 0.5,
#   cl = c,
#   mfrow = c(4, 4),
#   time.labels = time,
#   new.window = FALSE
# )
# 
# # names(c$cluster) <- rownames(temp_data2)[-1]
# # rownames(c$membership) <- rownames(temp_data2)[-1]
# # save(c, file = "c")
# load("c")
# 
# center <- c$centers
# rownames(center) <- paste("Cluster", rownames(center), sep = ' ')
# corrplot::corrplot(
#   corr = cor(t(center)),
#   type = "full",
#   diag = TRUE,
#   order = "hclust",
#   hclust.method = "ward.D",
#   # addrect = 5,
#   col = colorRampPalette(colors = rev(
#     RColorBrewer::brewer.pal(n = 11, name = "Spectral")
#   ))(n = 100),
#   number.cex = .7,
#   addCoef.col = "black"
# )
# 
# acore <- acore(data.s, c, min.acore = 0)
# acore
# 
# centers <- c$centers
# names(c$cluster) == rownames(c$membership)
# 
# cluster_info <-
#   data.frame(
#     variable_id = names(c$cluster),
#     c$membership,
#     cluster = c$cluster,
#     stringsAsFactors = FALSE
#   ) %>%
#   arrange(cluster)
# 
# openxlsx::write.xlsx(x = cluster_info,
#                      file = "cluster_info.xlsx", asTable = TRUE)
# 
# #####output the expression data of different clusters
# 
# ##plot for each cluster
# for (cluster_idx in 1:clust) {
#   cat(cluster_idx, " ")
#   dir.create(paste("cluster", cluster_idx, sep = "_"))
#   cluster_data <-
#     cluster_info %>%
#     dplyr::filter(cluster_idx == cluster_idx) %>%
#     dplyr::select(1, 1 + cluster_idx)
#   
#   colnames(cluster_data) <- c("variable_id", "membership")
#   
#   cluster_data <-
#     cluster_data %>%
#     dplyr::filter(membership > 0.5) %>%
#     dplyr::left_join(proteome_variable_info, by = "variable_id")
#   
#   openxlsx::write.xlsx(
#     x = cluster_data,
#     file = file.path(
#       paste("cluster", cluster_idx, sep = "_"),
#       paste("cluster", cluster_idx, ".xlsx", sep = "")
#     ),
#     asTable = TRUE,
#     overwrite = TRUE
#   )
#   
#   ###cluster plot
#   
#   temp =
#     temp_data2[cluster_data$variable_id,] %>%
#     data.frame(
#       membership = cluster_data$membership,
#       .,
#       stringsAsFactors = FALSE,
#       check.names = FALSE
#     ) %>%
#     tibble::rownames_to_column(var = "variable_id") %>%
#     tidyr::pivot_longer(
#       cols = -c(variable_id, membership),
#       names_to = "sample_id",
#       values_to = "value"
#     ) %>%
#     dplyr::mutate(sample_id = factor(
#       sample_id,
#       levels = c("negtaive_h", "ee", "el", "re", "positive_h")
#     ))
#   # dplyr::left_join(sample_info[, c("sample_id", "accurate_time")], by = "sample_id")
#   
#   plot <-
#     ggplot() +
#     geom_hline(yintercept = 0) +
#     geom_line(aes(sample_id, value, group = variable_id, color = membership),
#               data = temp) +
#     scale_color_gradientn(colours = c(RColorBrewer::brewer.pal(n = 9, name = "OrRd")[c(1:7)])) +
#     theme_bw() +
#     theme(
#       panel.grid = element_blank(),
#       legend.position = c(0, 1),
#       legend.justification = c(0, 1),
#       panel.grid.minor = element_blank(),
#       axis.title = element_text(size = 13),
#       axis.text = element_text(size = 12),
#       # axis.text.x = element_blank(),
#       axis.ticks.x = element_blank(),
#       panel.background = element_rect(fill = "transparent", color = NA),
#       plot.background = element_rect(fill = "transparent", color = NA),
#       legend.background = element_rect(fill = "transparent", color = NA)
#     ) +
#     labs(
#       x = "",
#       y = "Z-score",
#       title = paste(
#         "Cluster ",
#         cluster_idx,
#         " (",
#         nrow(cluster_data),
#         " taxa)",
#         sep = ""
#       )
#     )
#   
#   plot
#   
#   ggsave(
#     plot,
#     filename = file.path(
#       paste("cluster", cluster_idx, sep = "_"),
#       paste("cluster", cluster_idx, ".pdf", sep = "")
#     ),
#     width = 8,
#     height = 6
#   )
# }
