#mediation cases
library(ggplot2)
library(dplyr)

setwd("~/Library/CloudStorage/Box-Box/human_microbiome_project/data_analysis/mediation_analysis/")

IS_mediation <- read.csv("./Previous_Version/sample_wise_IS/summary/mediation_result_all_IS.csv")
IR_mediation <- read.csv("./Previous_Version/sample_wise_IR/summary/mediation_result_all_IR.csv")
all_mediation <- read.csv("./Previous_Version/sample_wise/summary/mediation_result_all.csv")

table(IS_mediation$treat_mediator_p_adjust > 0.05)
table(IR_mediation$treat_mediator_p_adjust > 0.05)

IS_mediation %>% 
  select(treat_mediator_cor, treat_mediator_p,mediator_phenotype_cor, treat_phenotype_cor,mediator_phenotype_p, treat_phenotype_p)


IS_mediation %>% sort()

hist(IS_mediation$prop_mediate)

IS_mediation %>% filter(Species == "uniformis")


IR_mediation %>% filter(phenotype == "LDL") #%>%
select(treat_mediator_cor, treat_mediator_p,mediator_phenotype_cor, treat_phenotype_cor,mediator_phenotype_p, treat_phenotype_p,mediator_true_name, treat_true_name, phenotype_true_name)
table(IS_mediation$phenotype)


IR_mediation %>% filter(treat_mediator_p < 0.05 & mediator_phenotype_p < 0.05) %>%
select(treat_mediator_cor, treat_mediator_p,mediator_phenotype_cor, treat_phenotype_cor,mediator_phenotype_p, treat_phenotype_p,mediator_true_name, treat_true_name, phenotype_true_name)

IR_mediation %>% filter(treat_class == "Oral microbiome" & mediator_class =="Cytokine")

intersect(paste(IS_mediation$treat_true_name, IS_mediation$phenotype_true_name, sep="_"), paste(IR_mediation$treat_true_name, IR_mediation$phenotype_true_name, sep="_"))
intersect(paste(all_mediation$treat_true_name, all_mediation$phenotype_true_name, sep="_"), paste(IS_mediation$treat_true_name, IS_mediation$phenotype_true_name, sep="_"))
intersect(paste(all_mediation$treat_true_name, all_mediation$phenotype_true_name, sep="_"), paste(IR_mediation$treat_true_name, IR_mediation$phenotype_true_name, sep="_"))


all_mediation %>% filter(treat_true_name == "Haemophilus" & phenotype_true_name == "TGL") %>% 
  select(treat_mediator_cor, treat_mediator_p,mediator_phenotype_cor, treat_phenotype_cor,mediator_phenotype_p, treat_phenotype_p,mediator_true_name, treat_true_name, phenotype_true_name, phenotype_true_name,acme_p, acme_p_inverse,prop_mediate)

IS_mediation %>% filter(treat_true_name == "Haemophilus" & phenotype_true_name == "TGL") %>% 
  select(treat_mediator_cor, treat_mediator_p,mediator_phenotype_cor, treat_phenotype_cor,mediator_phenotype_p, treat_phenotype_p,mediator_true_name, treat_true_name, phenotype_true_name, phenotype_true_name,acme_p, acme_p_inverse,prop_mediate)



all_mediation %>% filter(treat_true_name == "Akkermansia" & phenotype_true_name == "A1C" & mediator_true_name == "IL15") %>% 
  select(treat_mediator_cor, treat_mediator_p,mediator_phenotype_cor, treat_phenotype_cor,mediator_phenotype_p, treat_phenotype_p,mediator_true_name, treat_true_name, phenotype_true_name,acme_p, acme_p_inverse,prop_mediate)

IS_mediation %>% filter(treat_true_name == "Akkermansia" & phenotype_true_name == "A1C" & mediator_true_name == "IL15") %>% 
  select(treat_mediator_cor, treat_mediator_p,mediator_phenotype_cor, treat_phenotype_cor,mediator_phenotype_p, treat_phenotype_p,mediator_true_name, treat_true_name, phenotype_true_name,acme_p, acme_p_inverse,prop_mediate)

IR_mediation %>% filter(treat_true_name == "Haemophilus") %>% 
  select(treat_mediator_cor, treat_mediator_p,mediator_phenotype_cor, treat_phenotype_cor,mediator_phenotype_p, treat_phenotype_p,mediator_true_name, treat_true_name, phenotype_true_name,acme_p, acme_p_inverse,prop_mediate)

IS_mediation %>% filter(treat_true_name == "Akkermansia")



filter(all_mediation, phenotype_true_name == "GLU")
filter(all_mediation, phenotype_true_name == "fpg_mg_ml")
fpg_mg_ml

table(all_mediation$phenotype_true_name)


