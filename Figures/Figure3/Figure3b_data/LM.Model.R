library(ggplot2)
library(ggstatsplot)
library(lme4)
library(lmerTest)
library(patchwork)
library(dplyr)
library(jtools)
library(ggstance)
library(broom.mixed)

setwd("~/Library/CloudStorage/Box-Box/human_microbiome_project/Figures/Figure3/Figure3b_data/")

load("../../../../XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/Robject/Revision_MultiOmes_0509.RData")
source("../../../code/tools.R")
load("./rawdata/raw.dist.RData")

#load("./significant.corre.only/temp_data_stool")
#load("./significant.corre.only/temp_data_skin")
#load("./significant.corre.only/temp_data_oral")
#load("./significant.corre.only/temp_data_nasal")
 
ggplot(stool_braydist_by_sample, aes(x=diffdays, y=dist)) + geom_point()

scatter_stool <- ggscatterstats(stool_braydist_by_sample, 
               x = diffdays,
               y = dist, 
               title = "Stool")
scatter_stool

scatter_skin <-ggscatterstats(skin_braydist_by_sample, 
               x = diffdays,
               y = dist, 
               title = "Skin")

scatter_oral <-ggscatterstats(oral_braydist_by_sample, 
               x = diffdays,
               y = dist, 
               title = "Oral")

scatter_nasal <-ggscatterstats(nasal_braydist_by_sample, 
               x = diffdays,
               y = dist, 
               title = "Nasal")

scatter.all <- scatter_stool + scatter_skin + scatter_oral + scatter_nasal
scatter.all
#ggsave("./Overall.Scatter.Plot.pdf", width = 14, height = 12, dpi = 300)

##############################################################################################################

all <- bind_rows(nasal_braydist_by_sample,
                 oral_braydist_by_sample,
                 skin_braydist_by_sample,
                 stool_braydist_by_sample)
head(all)
#if load processed data, then you do not need to load other data above
#load("./rawdata/processed.data.RData")
##############################################################################################################

all.plot <- ggscatterstats(all, 
               x = diffdays,
               y = dist, 
               color=dataset)
all.plot
 
all
all$IRIS <- "Unknown"
all$IRIS[all$subject_id1 %in% filter(sc, IRIS == "IS")$SubjectID] <- "IS"
all$IRIS[all$subject_id1 %in% filter(sc, IRIS == "IR")$SubjectID] <- "IR"

all$dataset <- factor(all$dataset, levels = c("stool", "skin", "oral","nasal"))
all$IRIS <- factor(all$IRIS, levels = c("IS", "IR", "Unknown"))

all.mixed <- lmer(dist ~ 1 + diffdays + (1 | subject_id1), data = all)
summary(all.mixed)

ggcoefstats(all.mixed)

stool.mixed <- lmer(dist ~ diffdays + (1 | subject_id1), data = filter(all, dataset=="stool"))
skin.mixed <- lmer(dist ~ diffdays + (1 | subject_id1), data = filter(all, dataset=="skin"))
oral.mixed <- lmer(dist ~  diffdays + (1 | subject_id1), data = filter(all, dataset=="oral"))
nasal.mixed <- lmer(dist ~  diffdays + (1 | subject_id1), data = filter(all, dataset=="nasal"))

ALL.mixed <-lmer(dist ~  dataset:diffdays + (1 | subject_id1), data = all)
summary(ALL.mixed)

summary(stool.mixed)
summary(skin.mixed)
summary(oral.mixed)
summary(nasal.mixed)

stool.table <-ggcoefstats(stool.mixed,
                          exclude.intercept = F,
                          output = "tidy")
stool.table

skin.table <-ggcoefstats(skin.mixed,
                          exclude.intercept = F,
                          output = "tidy")
skin.table
oral.table <-ggcoefstats(oral.mixed,
                          exclude.intercept = F,
                          output = "tidy")
oral.table
nasal.table <-ggcoefstats(nasal.mixed,
                          exclude.intercept = F,
                          output = "tidy")
nasal.table

stool.table$bodysite <- "Stool"
skin.table$bodysite <- "Skin"
oral.table$bodysite <- "Oral"
nasal.table$bodysite <- "Nasal"

modelresult <- bind_rows(stool.table,skin.table,oral.table,nasal.table) %>% select(-expression)
#write.csv(file = "./table/model.result.csv", modelresult)

stool.mixed.iris <- lmer(dist ~  IRIS *diffdays + (1 | subject_id1), data = filter(all, dataset=="stool") %>% filter(IRIS != "Unknown"))
skin.mixed.iris <- lmer(dist ~  IRIS *diffdays + (1 | subject_id1), data = filter(all, dataset=="skin") %>% filter(IRIS != "Unknown"))
oral.mixed.iris<- lmer(dist ~   IRIS *diffdays + (1 | subject_id1), data = filter(all, dataset=="oral") %>% filter(IRIS != "Unknown"))
nasal.mixed.iris <- lmer(dist ~   IRIS *diffdays + (1 | subject_id1), data = filter(all, dataset=="nasal") %>% filter(IRIS != "Unknown"))

summary(stool.mixed.iris)
summary(skin.mixed.iris)
summary(oral.mixed.iris)
summary(nasal.mixed.iris)

stool.table.iris <- ggcoefstats(stool.mixed.iris,exclude.intercept = F,output = "tidy")
skin.table.iris <- ggcoefstats(skin.mixed.iris,exclude.intercept = F,output = "tidy")
oral.table.iris <- ggcoefstats(oral.mixed.iris,exclude.intercept = F,output = "tidy")
nasal.table.iris <- ggcoefstats(nasal.mixed.iris,exclude.intercept = F,output = "tidy")

stool.table.iris$bodysite <- "Stool"
skin.table.iris$bodysite <- "Skin"
oral.table.iris$bodysite <- "Oral"
nasal.table.iris$bodysite <- "Nasal"

modelresult.iris <- bind_rows(stool.table.iris,skin.table.iris,oral.table.iris,nasal.table.iris) %>% select(-expression)
#write.csv(file = "./table/model.result.iris.csv", modelresult.iris)

#################################################
stool.mixed.is <- lmer(dist ~  diffdays + (1 | subject_id1), data = filter(all, dataset=="stool") %>% filter(IRIS == "IS"))
stool.mixed.ir <- lmer(dist ~  diffdays + (1 | subject_id1), data = filter(all, dataset=="stool") %>% filter(IRIS == "IR"))
summary(stool.mixed.is)
summary(stool.mixed.ir)

filter(all, dataset=="stool") %>% filter(IRIS != "Unknown") %>% ggplot(aes(x=diffdays, y= dist, color=IRIS)) + geom_point(size=0) + geom_smooth(method = "lm")
filter(all, dataset=="skin") %>% filter(IRIS != "Unknown") %>% ggplot(aes(x=diffdays, y= dist, color=IRIS)) + geom_point(size=0) + geom_smooth(method = "lm")
filter(all, dataset=="oral") %>% filter(IRIS != "Unknown") %>% ggplot(aes(x=diffdays, y= dist, color=IRIS)) + geom_point(size=0) + geom_smooth(method = "lm")
filter(all, dataset=="nasal") %>% filter(IRIS != "Unknown") %>% ggplot(aes(x=diffdays, y= dist, color=IRIS)) + geom_point(size=0) + geom_smooth(method = "lm")

#################################################

p.effi <- plot_summs(stool.mixed, skin.mixed, oral.mixed, nasal.mixed,
                     model.names = c("Stool", "Skin", "Oral", "Nasal"), plot.distributions = TRUE)
p.effi <- p.effi + scale_color_manual(values = body_site_color) + ylab("Sample Collection Date Interval") + scale_fill_manual(values = body_site_color)
p.effi
#ggsave(filename = "./Model.Results.pdf", p.effi, width = 5, height = 4, dpi = 300)

################################################################
#adding data set as an interaction term, to estimate the difference in their beta coef between body sites

all$dataset <- factor(all$dataset, levels = c("stool", "skin", "oral", "nasal"))
stool.mixed.comare <- lmer(dist ~ -1 + diffdays*dataset + (1 | subject_id1), data = all)

all$dataset <- factor(all$dataset, levels = c("skin", "oral","nasal","stool"))
skin.mixed.comare  <- lmer(dist ~ -1 + diffdays*dataset + (1 | subject_id1), data = all)

all$dataset <- factor(all$dataset, levels = c("oral","nasal","stool", "skin"))
oral.mixed.comare  <- lmer(dist ~  -1 + diffdays*dataset + (1 | subject_id1), data = all)

all$dataset <- factor(all$dataset, levels = c("nasal","stool", "skin", "oral"))
nasal.mixed.comare  <- lmer(dist ~  -1 + diffdays*dataset + (1 | subject_id1), data = all)

all$dataset <- factor(all$dataset, levels = c("stool", "skin", "oral","nasal"))

summary(stool.mixed.comare)
summary(skin.mixed.comare)
summary(oral.mixed.comare)
summary(nasal.mixed.comare)

#################################################################
summary_info <- bind_rows(stool_summary_info,
                          skin_summary_info,
                          oral_summary_info,
                          nasal_summary_info)

summary_info$IRIS <- "Unknown"
summary_info$IRIS[summary_info$subject_id %in% filter(sc, IRIS == "IS")$SubjectID] <- "IS"
summary_info$IRIS[summary_info$subject_id %in% filter(sc, IRIS == "IR")$SubjectID] <- "IR"
summary_info$IRIS <- factor(summary_info$IRIS, levels = c("IS", "IR", "Unknown"))

between.stool <- summary_info %>% filter(IRIS != "Unknown") %>% filter(class=="stool") %>% 
  ggbetweenstats(x=IRIS, y=cor)

between.skin <-summary_info %>% filter(IRIS != "Unknown") %>% filter(class=="skin") %>% 
  ggbetweenstats(x=IRIS, y=cor)

between.oral <-summary_info %>% filter(IRIS != "Unknown") %>% filter(class=="oral") %>% 
  ggbetweenstats(x=IRIS, y=cor)

between.nasal <-summary_info %>% filter(IRIS != "Unknown") %>% filter(class=="nasal") %>% 
  ggbetweenstats(x=IRIS, y=cor)

between.stats <- between.stool + between.skin + between.oral + between.nasal
between.stats
#ggsave(filename = "./between_corr.pdf", between.stats, width = 10, height = 7, dpi = 300)


#################################################################################################################################
#-log(1-distance)
all.trans <- all 
all.trans$dist[all.trans$dist == 1] <- 0.9999999999999999
all.trans <- all.trans %>% mutate(dist_trans = log10(-log10(1-dist)))
pafter <- gghistostats(all.trans,dist_trans, title="After log10(-log10(1-dist)) transformation")
pbefore <- gghistostats(all.trans,dist, title="BC distance before transformation")
ptransformation <- pbefore + pafter
ptransformation
#ggsave("./Meng_Suggestion/Tranformation.pdf",ptransformation, width = 8, height = 6, dpi=300)

all.trans %>% ggscatterstats(x= dist, y = dist_trans, xlab="", ylab="", title = "")

stool.trans.mixed <- lmer(dist_trans ~ diffdays + (1 | subject_id1), data = filter(all.trans, dataset=="stool"))
skin.trans.mixed <- lmer(dist_trans ~ diffdays + (1 | subject_id1), data = filter(all.trans, dataset=="skin"))
oral.trans.mixed <- lmer(dist_trans ~  diffdays + (1 | subject_id1), data = filter(all.trans, dataset=="oral"))
nasal.trans.mixed <- lmer(dist_trans ~  diffdays + (1 | subject_id1), data = filter(all.trans, dataset=="nasal"))

ALL.trans.mixed <-lmer(dist ~ -1 + dataset:diffdays + (1 | subject_id1) + ( 1 + dataset | subject_id1), data = all.trans)
summary(ALL.trans.mixed)

summary(stool.trans.mixed)
summary(skin.trans.mixed)
summary(oral.trans.mixed)
summary(nasal.trans.mixed)


##################################################################################################################################
stool.trans.mixed.iris <- lmer(dist_trans ~ -1 + IRIS * diffdays + (1 | subject_id1) + (1 + IRIS | subject_id1), data = filter(all.trans, dataset=="stool") %>% filter(IRIS != "Unknown"))
summary(stool.trans.mixed.iris)

skin.trans.mixed.iris <- lmer(dist_trans ~ -1 + IRIS * diffdays + (1 | subject_id1) + (1 + IRIS | subject_id1), data = filter(all.trans, dataset=="skin") %>% filter(IRIS != "Unknown"))
summary(skin.trans.mixed.iris)

oral.trans.mixed.iris <- lmer(dist_trans ~ -1 + IRIS * diffdays + (1 | subject_id1) + (1 + IRIS | subject_id1), data = filter(all.trans, dataset=="oral") %>% filter(IRIS != "Unknown"))
summary(oral.trans.mixed.iris)

nasal.trans.mixed.iris <- lmer(dist_trans ~ -1 + IRIS * diffdays + (1 | subject_id1) + ( 1+ IRIS | subject_id1), data = filter(all.trans, dataset=="nasal") %>% filter(IRIS != "Unknown"))
summary(nasal.trans.mixed.iris)

stool.iris.smooth <- filter(all.trans, dataset=="stool") %>% filter(IRIS != "Unknown") %>% 
  ggplot(aes(x=diffdays, y=dist_trans, color=IRIS))  + geom_smooth(method = "lm") + ggtitle("stool")+ base_theme + geom_point(size=0.01)

skin.iris.smooth <- filter(all.trans, dataset=="skin") %>% filter(IRIS != "Unknown") %>% 
  ggplot(aes(x=diffdays, y=dist_trans, color=IRIS)) + geom_point(size=0.01) + geom_smooth(method = "lm") + ggtitle("skin")+ base_theme+ geom_point(size=0.01)

oral.iris.smooth <- filter(all.trans, dataset=="oral") %>% filter(IRIS != "Unknown") %>% 
  ggplot(aes(x=diffdays, y=dist_trans, color=IRIS)) + geom_point(size=0.01) + geom_smooth(method = "lm") + ggtitle("oral")+ base_theme+ geom_point(size=0.01)

nasal.iris.smooth <- filter(all.trans, dataset=="nasal") %>% filter(IRIS != "Unknown") %>% 
  ggplot(aes(x=diffdays, y=dist_trans, color=IRIS)) + geom_point(size=0.01) + geom_smooth(method = "lm") + ggtitle("nasal")+ base_theme+ geom_point(size=0.01)

iris.smooth <- stool.iris.smooth + skin.iris.smooth + oral.iris.smooth + nasal.iris.smooth
iris.smooth 
#ggsave(filename = "./IRIS.diffdaysbydist.smooth.pdf", iris.smooth, width = 8, height = 6, dpi = 300)

#################################################################################################################################
#no significant difference on number of observation
summary_info %>% dplyr::filter(IRIS != "Unknown") %>% dplyr::filter(class=="stool") %>% ggstatsplot::ggbetweenstats(x= IRIS, y= number) + ggplot2::stat_summary(fun.y = median)

################
#check if the correlation was similar within two body sites
st.sk <- left_join((summary_info %>%  filter(class=="stool")), 
                      (summary_info  %>% filter(class=="skin")), by="subject_id", all = T)
st.sk.or <- left_join(st.sk, (summary_info  %>% filter(class=="oral")),by="subject_id", all = T)
st.sk.or.ns <-left_join(st.sk.or,(summary_info %>% filter(class=="nasal")),by="subject_id", all = T)
st.sk.or.ns <- st.sk.or.ns %>% rename(stool.cor = cor.x) %>% rename(skin.cor = cor.y)%>% rename(oral.cor = cor.x.x)%>% rename(nasal.cor = cor.y.y)

shapiro.test(st.sk.or.ns$stool.cor)
shapiro.test(st.sk.or.ns$skin.cor)
shapiro.test(st.sk.or.ns$oral.cor)
shapiro.test(st.sk.or.ns$nasal.cor)

cor.test(st.sk.or.ns$stool.cor, st.sk.or.ns$skin.cor, method = "pearson")

p.st.sk.all <- ggscatterstats(st.sk.or.ns, x=stool.cor, y=skin.cor, title = "All Stool Skin", type = "noonparametric")
p.st.or.all <- ggscatterstats(st.sk.or.ns, x=stool.cor, y=oral.cor, title = "All Stool Oral", type = "noonparametric")
p.st.ns.all <- ggscatterstats(st.sk.or.ns, x=stool.cor, y=nasal.cor, title = "All Stool Nasal", type = "noonparametric")
p.sk.or.all <- ggscatterstats(st.sk.or.ns, x=skin.cor, y=oral.cor, title = "All Skin Oral", type = "noonparametric")
p.sk.ns.all <- ggscatterstats(st.sk.or.ns, x=skin.cor, y=nasal.cor, title = "All Skin Nasal", type = "noonparametric")
p.or.ns.all <- ggscatterstats(st.sk.or.ns, x=oral.cor, y=nasal.cor, title = "All Oral Nasal", type = "noonparametric")

p.all <- p.st.sk.all + p.st.or.all +p.st.ns.all + p.sk.or.all + p.sk.ns.all + p.or.ns.all
p.all
#ggsave(filename = "./All_between_bodysite.pdf", p.all, height = 10, width = 17, dpi = 300)

#################
#seperate IR and IS
st.sk.is <- left_join((summary_info %>% filter(IRIS == "IS") %>% filter(class=="stool")), 
      (summary_info %>% filter(IRIS == "IS") %>% filter(class=="skin")), by="subject_id", all = T)
st.sk.or.is <- left_join(st.sk.is, (summary_info %>% filter(IRIS == "IS") %>% filter(class=="oral")),by="subject_id", all = T)
st.sk.or.ns.is <-left_join(st.sk.or.is,(summary_info %>% filter(IRIS == "IS") %>% filter(class=="nasal")),by="subject_id", all = T)
st.sk.or.ns.is <- st.sk.or.ns.is %>% rename(stool.cor = cor.x) %>% rename(skin.cor = cor.y)%>% rename(oral.cor = cor.x.x)%>% rename(nasal.cor = cor.y.y)

shapiro.test(st.sk.or.ns.is$stool.cor)

cor.test(st.sk.or.ns.is$stool.cor, st.sk.or.ns.is$skin.cor, method = "pearson")

p.st.sk.is <- ggscatterstats(st.sk.or.ns.is, x=stool.cor, y=skin.cor, title = "IS Stool Skin")
p.st.or.is <- ggscatterstats(st.sk.or.ns.is, x=stool.cor, y=oral.cor, title = "IS Stool Oral")
p.st.ns.is <- ggscatterstats(st.sk.or.ns.is, x=stool.cor, y=nasal.cor, title = "IS Stool Nasal")
p.sk.or.is <- ggscatterstats(st.sk.or.ns.is, x=skin.cor, y=oral.cor, title = "IS Skin Oral")
p.sk.ns.is <- ggscatterstats(st.sk.or.ns.is, x=skin.cor, y=nasal.cor, title = "IS Skin Nasal")
p.or.ns.is <- ggscatterstats(st.sk.or.ns.is, x=oral.cor, y=nasal.cor, title = "IS Oral Nasal")

p.is <- p.st.sk.is + p.st.or.is +p.st.ns.is + p.sk.or.is + p.sk.ns.is + p.or.ns.is
p.is
#ggsave(filename = "./IS_between_bodysite.pdf", p.is, height = 10, width = 17, dpi = 300)

st.sk.ir <- left_join((summary_info %>% filter(IRIS == "IR") %>% filter(class=="stool")), 
                      (summary_info %>% filter(IRIS == "IR") %>% filter(class=="skin")), by="subject_id", all = T)
st.sk.or.ir <- left_join(st.sk.ir, (summary_info %>% filter(IRIS == "IR") %>% filter(class=="oral")),by="subject_id", all = T)
st.sk.or.ns.ir <-left_join(st.sk.or.ir,(summary_info %>% filter(IRIS == "IR") %>% filter(class=="nasal")),by="subject_id", all = T)
st.sk.or.ns.ir <- st.sk.or.ns.ir %>% rename(stool.cor = cor.x) %>% rename(skin.cor = cor.y)%>% rename(oral.cor = cor.x.x)%>% rename(nasal.cor = cor.y.y)

shapiro.test(st.sk.or.ns.ir$stool.cor)

cor.test(st.sk.or.ns.ir$stool.cor, st.sk.or.ns.ir$skin.cor, method = "pearson")

p.st.sk.ir <- ggscatterstats(st.sk.or.ns.ir, x=stool.cor, y=skin.cor, title = "IR Stool Skin")
p.st.or.ir <- ggscatterstats(st.sk.or.ns.ir, x=stool.cor, y=oral.cor, title = "IR Stool Oral")
p.st.ns.ir <- ggscatterstats(st.sk.or.ns.ir, x=stool.cor, y=nasal.cor, title = "IR Stool Nasal")
p.sk.or.ir <- ggscatterstats(st.sk.or.ns.ir, x=skin.cor, y=oral.cor, title = "IR Skin Oral")
p.sk.ns.ir <- ggscatterstats(st.sk.or.ns.ir, x=skin.cor, y=nasal.cor, title = "IR Skin Nasal")
p.or.ns.ir <- ggscatterstats(st.sk.or.ns.ir, x=oral.cor, y=nasal.cor, title = "IR Oral Nasal")

p.ir <- p.st.sk.ir + p.st.or.ir + p.st.ns.ir + p.sk.or.ir + p.sk.ns.ir + p.or.ns.ir
p.ir
#ggsave(filename = "./IR_between_bodysite.pdf", p.ir, height = 10, width = 17, dpi = 300)

save(file = "./rawdata/processed.data.RData", all,summary_info)


















