no_source()

# set work directory
setwd(masstools::get_project_wd())
library(tidyverse)
library(phyloseq)
rm(list = ls())

source("code/tools.R")

# ####load data
# ###stool microsample geneus level data
# load("data/from_xin/Genus Table/Genus.RData")
# physeqGenus_ST
# 
# expression_data = phyloseq::otu_table(object = physeqGenus_ST)
# expression_data = 
#   expression_data@.Data %>% 
#   t() %>% 
#   as.data.frame()
# 
# sample_info = phyloseq::sample_data(object = physeqGenus_ST)
# sample_info = sample_info@.Data %>% 
#   do.call(cbind, .) %>% 
#   as.data.frame()
# 
# colnames(sample_info) = colnames(phyloseq::sample_data(object = physeqGenus_ST))
# 
# variable_info = phyloseq::tax_table(object = physeqGenus_ST)
# variable_info = variable_info@.Data %>% 
#   as.data.frame()
# 
# st_microbiome_expression_data = expression_data
# st_microbiome_sample_info = sample_info
# st_microbiome_variable_info = variable_info
# 
# colnames(st_microbiome_expression_data) == st_microbiome_sample_info$RandomID
# colnames(st_microbiome_expression_data) = st_microbiome_sample_info$SampleID
# rownames(st_microbiome_expression_data) == rownames(st_microbiome_variable_info)
# 
# st_microbiome_variable_info = 
# st_microbiome_variable_info %>% 
#   tibble::rownames_to_column(var = "variable_id")
# 
# ###st microbiome
# load("data_analysis/st_microbiome/data_preparation/sample_info")
# 
# dim(sample_info)
# dim(st_microbiome_sample_info)
# sample_info$RandomID == st_microbiome_sample_info$RandomID
# st_microbiome_sample_info = sample_info

load("data_analysis/st_microbiome/data_preparation/sample_info")
load("data_analysis/st_microbiome/data_preparation/expression_data")
load("data_analysis/st_microbiome/data_preparation/variable_info")

dim(variable_info)

st_microbiome_sample_info = sample_info
st_microbiome_expression_data = expression_data
st_microbiome_variable_info = variable_info

st_microbiome_sample_info = 
st_microbiome_sample_info %>% 
  dplyr::mutate(SSPG = as.numeric(SSPG)) %>% 
  dplyr::mutate(iris = case_when(
    SSPG > 125 ~ "IR",
    SSPG < 125 ~ "IS"
  )) 

library(lubridate)

st_microbiome_sample_info$days = 
  st_microbiome_sample_info$Date %>% 
  lubridate::as_date() %>% 
  yday() 

st_microbiome_sample_info = 
st_microbiome_sample_info %>% 
  dplyr::mutate(
    days =  yday(lubridate::as_date(Date)),
    months = lubridate::month(lubridate::as_date(st_microbiome_sample_info$Date)),
    weeks = lubridate::week(lubridate::as_date(st_microbiome_sample_info$Date))
  ) 

#######work directory
setwd(masstools::get_project_wd())
setwd("data_analysis/st_microbiome/season_analysis")

zero_percent = 
st_microbiome_expression_data %>% 
  apply(1, function(x){
    sum(x == 0)/ncol(st_microbiome_expression_data)
  })

sum(zero_percent > 0.99)/nrow(st_microbiome_variable_info)

##here we remove the genus with 0 > 99%
remain_idx = which(zero_percent < 0.95)

length(remain_idx)

save(remain_idx, file = "remain_idx")

library(plyr)
library(lubridate)
library(mgcv)
library(gratia)

#https://fromthebottomoftheheap.net/2014/05/09/modelling-seasonal-data-with-gam/
###GAMM to find the seasonality of microbiomes
ctrl <- lmeControl(opt = 'optim')

p_value_all = NULL
coeffs_all = NULL

dir.create("smooth_plot")

# for(idx in remain_idx){
#   cat(idx, " ")
#   x = as.numeric(st_microbiome_expression_data[idx, ])
#   x = (x - mean(x))/sd(x)
#   temp_data =
#     data.frame(days = st_microbiome_sample_info$days,
#                weeks = st_microbiome_sample_info$weeks,
#                months = st_microbiome_sample_info$months,
#                x = x,
#                subject_id = st_microbiome_sample_info$subject_id,
#                bmi = st_microbiome_sample_info$BMI,
#                SSPG = st_microbiome_sample_info$SSPG,
#                iris = st_microbiome_sample_info$iris) %>%
#     dplyr::filter(!is.na(SSPG))
# 
#   mod <- try(gamm(
#     formula = x ~ bmi + iris + s(days, bs = "cc"),
#     data = temp_data,
#     method = "REML",
#     control = ctrl,
#     random = list(subject_id = ~ 1),
#     knots = list(TimeOfYear = c(0, 356))
#   ), silent = TRUE)
# 
#   if(class(mod) == "try-error"){
#     p_value = 1
#     coeffs = rep(NA, 11)
#   }else{
#     p_value = summary(mod$gam)$s.table[4]
#     coeffs = mod$gam$coefficients
#   }
# 
#   p_value_all = c(p_value_all, p_value)
#   coeffs_all = rbind(coeffs_all, coeffs)
# 
#   if(p_value < 0.05){
#     temp_plot =
#       gratia::draw(
#         object = mod$gam,
#         residuals = FALSE,
#         ci_level = 0.95,
#         select = 1,
#         parametric = FALSE,
#         contour = TRUE,
#         contour_col = "red"
#       ) +
#       annotate(geom = "rect", xmin = 0, xmax = 90,
#                ymin = -Inf, ymax = Inf, fill = season_color[1], alpha = 0.1) +
#       annotate(geom = "rect", xmin = 90, xmax = 180,
#                ymin = -Inf, ymax = Inf, fill = season_color[2], alpha = 0.1) +
#       annotate(geom = "rect", xmin = 180, xmax = 270,
#                ymin = -Inf, ymax = Inf, fill = season_color[3], alpha = 0.1) +
#       annotate(geom = "rect", xmin = 270, xmax = 360,
#                ymin = -Inf, ymax = Inf, fill = season_color[4], alpha = 0.1) +
#       theme_bw() +
#       labs(x = "Days of year", ylab = "Effect",
#            title = paste(st_microbiome_variable_info$Genus[idx],"/p value:", round(p_value, 4))) +
#       theme(axis.title = element_text(size = 13),
#             axis.text = element_text(size = 12)) +
#       scale_x_continuous(expand = expansion(mult = c(0.05,0.05)),
#                          breaks = c(0, 90, 180, 270, 356),
#                          labels = c(0, 90, 180, 270, 356))
#     ggsave(
#       temp_plot,
#       filename = file.path(
#         "smooth_plot",
#         paste(st_microbiome_variable_info$Genus[idx], "pdf", sep = ".")
#       ),
#       width = 9,
#       height = 7
#     )
#   }
# }
# 
# names(p_value_all) = st_microbiome_variable_info$variable_id[remain_idx]
# coeffs_all =
# coeffs_all %>%
# as.data.frame()
# 
# rownames(coeffs_all) = st_microbiome_variable_info$variable_id[remain_idx]
# 
# save(p_value_all, file = "p_value_all")
# save(coeffs_all, file = "coeffs_all")
# 
# p_value_adj_all = p.adjust(p_value_all, method = "BH")
# save(p_value_adj_all, file = "p_value_adj_all")

load("p_value_all")
load("coeffs_all")
load("p_value_adj_all")

sum(p_value_all < 0.05)
sum(p_value_adj_all < 0.05)














  


