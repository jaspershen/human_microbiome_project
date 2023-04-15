#Figure2 Statistics

library(ggstatsplot)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(patchwork)

setwd("~/Library/CloudStorage/Box-Box/human_microbiome_project/Figures/Figure2/")
source("../../code/tools.R")
phylum_color2 <- phylum_color[c(1:4, 6)]
names(phylum_color2)[3] <- "Other"

familyscore <- read.csv(file = "../../Supplementary_data/family_score_based_on_phylum.csv", header = T)
DMIscore <- read.csv(file = "../../Supplementary_data/dmi_based_on_phylum.csv", header = T)

ST.pre <- read.csv("../../../XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/tables/ST.pr.csv", header = T, row.names = 1)
SK.pre <- read.csv("../../../XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/tables/SK.pr.csv", header = T, row.names = 1)
OR.pre <- read.csv("../../../XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/tables/OR.pr.csv", header = T, row.names = 1)
NS.pre <- read.csv("../../../XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/tables/NS.pr.csv", header = T, row.names = 1)

ST.pre.mean <- colMeans(ST.pre) %>% data.frame() 
SK.pre.mean <- colMeans(SK.pre) %>% data.frame()
OR.pre.mean <- colMeans(OR.pre) %>% data.frame()
NS.pre.mean <- colMeans(NS.pre) %>% data.frame()

colnames(ST.pre.mean) <- "Prev"
colnames(SK.pre.mean) <- "Prev"
colnames(OR.pre.mean) <- "Prev"
colnames(NS.pre.mean) <- "Prev"

table(familyscore$class)

familyscore$genus[familyscore$genus == "Escherichia/Shigella"] <- "Escherichia.Shigella"

ST.pre.mean.fam <- merge(ST.pre.mean, filter(familyscore, class== "Stool"),by.x="row.names",by.y = "genus")
SK.pre.mean.fam <- merge(SK.pre.mean, filter(familyscore, class== "Skin"),by.x="row.names",by.y = "genus")
OR.pre.mean.fam <- merge(OR.pre.mean, filter(familyscore, class== "Oral"),by.x="row.names",by.y = "genus")
NS.pre.mean.fam <- merge(NS.pre.mean, filter(familyscore, class== "Nasal"),by.x="row.names",by.y = "genus")


stool.pre.fm <- ggscatterstats(data=ST.pre.mean.fam,
               x=Prev,
               y=family_score,
               ggtheme=base_theme,
               ylab="Family Score",
               xlab="Prevalence",
               title="Stool Microbiome")
stool.pre.fm

skin.pre.fm <- ggscatterstats(data=SK.pre.mean.fam,
                               x=Prev,
                               y=family_score,
                               ggtheme=base_theme,
                               ylab="Family Score",
                               xlab="Prevalence",
                              title="Skin Microbiome")
skin.pre.fm

oral.pre.fm <- ggscatterstats(data=OR.pre.mean.fam,
                               x=Prev,
                               y=family_score,
                               ggtheme=base_theme,
                               ylab="Family Score",
                               xlab="Prevalence",
                              title="Oral Microbiome")
oral.pre.fm

nasal.pre.fm <- ggscatterstats(data=NS.pre.mean.fam,
                               x=Prev,
                               y=family_score,
                               ggtheme=base_theme,
                               ylab="Family Score",
                               xlab="Prevalence",
                               title="Nasal Microbiome")
nasal.pre.fm

p1 <- stool.pre.fm + skin.pre.fm + oral.pre.fm + nasal.pre.fm
#ggsave(filename = "./Figure2Stat_Fam_Pre.pdf",p1, width = 14, height = 10, dpi = 300)

ST.pre.mean.fam[ST.pre.mean.fam$family_score ==1,]

###########perform same analysis on DMI
table(DMIscore$class)
DMIscore$genus[DMIscore$genus == "Escherichia/Shigella"] <- "Escherichia.Shigella"
ST.pre.mean.dmi <- merge(ST.pre.mean, filter(DMIscore, class== "Stool"),by.x="row.names",by.y = "genus")
SK.pre.mean.dmi <- merge(SK.pre.mean, filter(DMIscore, class== "Skin"),by.x="row.names",by.y = "genus")
OR.pre.mean.dmi <- merge(OR.pre.mean, filter(DMIscore, class== "Oral"),by.x="row.names",by.y = "genus")
NS.pre.mean.dmi <- merge(NS.pre.mean, filter(DMIscore, class== "Nasal"),by.x="row.names",by.y = "genus")

stool.pre.DMI <- ggscatterstats(data=ST.pre.mean.dmi,
                               x=Prev,
                               y=dmi,
                               ggtheme=base_theme,
                               ylab="DMI Score",
                               xlab="Prevalence",
                               title="Stool Microbiome")
stool.pre.DMI

skin.pre.DMI <- ggscatterstats(data=SK.pre.mean.dmi,
                              x=Prev,
                              y=dmi,
                              ggtheme=base_theme,
                              ylab="DMI Score",
                              xlab="Prevalence",
                              title="Skin Microbiome")
skin.pre.DMI

oral.pre.DMI <- ggscatterstats(data=OR.pre.mean.dmi,
                              x=Prev,
                              y=dmi,
                              ggtheme=base_theme,
                              ylab="DMI Score",
                              xlab="Prevalence",
                              title="Oral Microbiome")
oral.pre.DMI

nasal.pre.DMI <- ggscatterstats(data=NS.pre.mean.dmi,
                               x=Prev,
                               y=dmi,
                               ggtheme=base_theme,
                               ylab="DMI Score",
                               xlab="Prevalence",
                               title="Nasal Microbiome")
nasal.pre.DMI

p2 <- stool.pre.DMI + skin.pre.DMI + oral.pre.DMI + nasal.pre.DMI
#ggsave(filename = "./Figure2Stat_DMI_Pre.pdf",p2, width = 14, height = 10, dpi = 300)

####################################################################################################
familyscore$class <- factor(familyscore$class, levels = c("Stool", "Skin", "Oral", "Nasal"))
#Shapiro wilk test indicate this data is significantly different from normal distribution (p-value < 2.2e-16)
shapiro.test(familyscore$family_score)

p3 <- ggbetweenstats(
  data  = familyscore,
  x     = class,
  y     = family_score,
  type = "noneparametric",
  p.adjust.method = "BH",
  title = "Faimly Score Cross Body Sites")
p3 <- p3 + scale_color_manual(values = body_site_color) + base_theme
p3
#ggsave(filename = "./FamilyScore.Bodysite.pdf", p3, width = 6, height = 5, dpi=300)

DMIscore$class <- factor(DMIscore$class, levels = c("Stool", "Skin", "Oral", "Nasal"))
#Shapiro wilk test indicate this data is significantly different from normal distribution (p-value < 2.2e-16)
shapiro.test(DMIscore$dmi)

p4 <- ggbetweenstats(
  data  = DMIscore,
  x     = class,
  y     = dmi,
  type = "noneparametric",
  p.adjust.method = "BH",
  title = "DMI Score Cross Body Sites")
p4 <- p4  + scale_color_manual(values = body_site_color) + base_theme
p4
#ggsave(filename = "./DMIScore.Bodysite.pdf", p4, width = 6, height = 5, dpi=300)


##################################################################################################################################
familyscore$Phylum <- factor(familyscore$Phylum, levels=c("Actinobacteria","Bacteroidetes","Firmicutes","Proteobacteria","Other"))

p5 <- ggbetweenstats(
  data  = filter(familyscore, class == "Stool"),
  x = Phylum,
  y = family_score,
  type = "noneparametric",
  p.adjust.method = "BH",
  title = "Faimly Score | Stool Microbiome")
p5 <- p5 + scale_color_manual(values = phylum_color2) + base_theme
p5

#ggsave(filename = "./FamilyScore_Stool.pdf",p5, width = 6, height = 5, dpi=300)

p6 <- ggbetweenstats(
  data  = filter(familyscore, class == "Skin"),
  x = Phylum,
  y = family_score,
  type = "noneparametric",
  p.adjust.method = "BH",
  title = "Faimly Score | Skin Microbiome")
p6 <- p6 + scale_color_manual(values = phylum_color2) + base_theme
p6

ggsave(filename = "./FamilyScore_Skin.pdf",p6, width = 6, height = 5, dpi=300)

p7 <- ggbetweenstats(
  data  = filter(familyscore, class == "Oral"),
  x = Phylum,
  y = family_score,
  type = "noneparametric",
  p.adjust.method = "BH",
  title = "Faimly Score | Oral Microbiome")
p7 <- p7 + scale_color_manual(values = phylum_color2) + base_theme
p7

ggsave(filename = "./FamilyScore_Oral.pdf",p7, width = 6, height = 5, dpi=300)


p8 <- ggbetweenstats(
  data  = filter(familyscore, class == "Nasal"),
  x = Phylum,
  y = family_score,
  type = "noneparametric",
  p.adjust.method = "BH",
  title = "Faimly Score | Nasal Microbiome")
p8 <- p8 + scale_color_manual(values = phylum_color2) + base_theme
p8

ggsave(filename = "./FamilyScore_Nasal.pdf",p8, width = 6, height = 5, dpi=300)

filter(familyscore, class == "Oral") %>% arrange(family_score)
filter(familyscore, class == "Nasal") %>% arrange(family_score)

##################################################################################################################################
DMIscore$Phylum <- factor(DMIscore$Phylum, levels=c("Actinobacteria","Bacteroidetes","Firmicutes","Proteobacteria","Other"))

p9 <- ggbetweenstats(
  data  = filter(DMIscore, class == "Stool"),
  x = Phylum,
  y = dmi,
  type = "noneparametric",
  p.adjust.method = "BH",
  title = "DMI Score | Stool Microbiome")
p9 <- p9 + scale_color_manual(values = phylum_color2) + base_theme
p9
#ggsave(filename = "./DMIScore_Stool.pdf",p9, width = 6, height = 5, dpi=300)
wilcox.test(filter(DMIscore,class == "Stool" & Phylum == "Bacteroidetes")$dmi, filter(DMIscore, class == "Stool" & Phylum == "Firmicutes")$dmi)
p.adjust(1.871e-07,method = "BH", 10)

t.test(filter(DMIscore,class == "Skin" & Phylum == "Actinobacteria")$dmi, filter(DMIscore, class == "Nasal" & Phylum == "Actinobacteria")$dmi)


p10 <- ggbetweenstats(
  data  = filter(DMIscore, class == "Skin"),
  x = Phylum,
  y = dmi,
  type = "noneparametric",
  p.adjust.method = "BH",
  title = "DMI Score | Skin Microbiome")
p10 <- p10 + scale_color_manual(values = phylum_color2) + base_theme
p10
ggsave(filename = "./DMIScore_Skin.pdf",p10, width = 6, height = 5, dpi=300)

p11 <- ggbetweenstats(
  data  = filter(DMIscore, class == "Oral"),
  x = Phylum,
  y = dmi,
  type = "noneparametric",
  p.adjust.method = "BH",
  title = "DMI Score | Oral Microbiome")
p11 <- p11 + scale_color_manual(values = phylum_color2) + base_theme
p11
ggsave(filename = "./DMIScore_Oral.pdf",p11, width = 6, height = 5, dpi=300)

p12 <- ggbetweenstats(
  data  = filter(DMIscore, class == "Nasal"),
  x = Phylum,
  y = dmi,
  type = "noneparametric",
  p.adjust.method = "BH",
  title = "DMI Score | Nasal Microbiome")
p12 <- p12 + scale_color_manual(values = phylum_color2) + base_theme
p12
ggsave(filename = "./DMIScore_Nasal.pdf",p12, width = 6, height = 5, dpi=300)

############################

p13 <- ggbetweenstats(
  data  = filter(familyscore, Phylum== "Actinobacteria"),
  x     = class,
  y     = family_score,
  type = "noneparametric",
  p.adjust.method = "BH",
  title = "Faimly Score of Actinobacteria Cross Body Sites")
p13 <- p13 + scale_color_manual(values = body_site_color) + base_theme
p13
ggsave(filename = "./FamilyScore_Actinobacteria.Bodysite.pdf", p13, width = 6, height = 5, dpi=300)

p14 <- ggbetweenstats(
  data  = filter(familyscore, Phylum== "Bacteroidetes"),
  x     = class,
  y     = family_score,
  type = "noneparametric",
  p.adjust.method = "BH",
  title = "Faimly Score of Bacteroidetes Cross Body Sites")
p14 <- p14 + scale_color_manual(values = body_site_color) + base_theme
p14
ggsave(filename = "./FamilyScore_Bacteroidetes.Bodysite.pdf", p14, width = 6, height = 5, dpi=300)

p15 <- ggbetweenstats(
  data  = filter(familyscore, Phylum== "Firmicutes"),
  x     = class,
  y     = family_score,
  type = "noneparametric",
  p.adjust.method = "BH",
  title = "Faimly Score of Firmicutes Cross Body Sites")
p15 <- p15 + scale_color_manual(values = body_site_color) + base_theme
p15
ggsave(filename = "./FamilyScore_Firmicutes.Bodysite.pdf", p15, width = 6, height = 5, dpi=300)

p16 <- ggbetweenstats(
  data  = filter(familyscore, Phylum== "Proteobacteria"),
  x     = class,
  y     = family_score,
  type = "noneparametric",
  p.adjust.method = "BH",
  title = "Faimly Score of Proteobacteria Cross Body Sites")
p16 <- p16 + scale_color_manual(values = body_site_color) + base_theme
p16
ggsave(filename = "./FamilyScore_Proteobacteria.Bodysite.pdf", p16, width = 6, height = 5, dpi=300)

p17 <- ggbetweenstats(
  data  = filter(familyscore, Phylum== "Other"),
  x     = class,
  y     = family_score,
  type = "noneparametric",
  p.adjust.method = "BH",
  title = "Faimly Score of Other Phylum Cross Body Sites")
p17 <- p17 + scale_color_manual(values = body_site_color) + base_theme
p17
ggsave(filename = "./FamilyScore_Other.Bodysite.pdf", p17, width = 6, height = 5, dpi=300)

p18 <- ggbetweenstats(
  data  = filter(DMIscore, Phylum== "Actinobacteria"),
  x     = class,
  y     = dmi,
  type = "noneparametric",
  p.adjust.method = "BH",
  title = "DMI Score of Actinobacteria Cross Body Sites")
p18 <- p18 + scale_color_manual(values = body_site_color) + base_theme
p18
ggsave(filename = "./DMIScore_Actinobacteria.Bodysite.pdf", p18, width = 6, height = 5, dpi=300)

p19 <- ggbetweenstats(
  data  = filter(DMIscore, Phylum== "Bacteroidetes"),
  x     = class,
  y     = dmi,
  type = "noneparametric",
  p.adjust.method = "BH",
  title = "DMI Score of Bacteroidetes Cross Body Sites")
p19 <- p19 + scale_color_manual(values = body_site_color) + base_theme
p19
ggsave(filename = "./DMIScore_Bacteroidetes.Bodysite.pdf", p19, width = 6, height = 5, dpi=300)

p20 <- ggbetweenstats(
  data  = filter(DMIscore, Phylum== "Firmicutes"),
  x     = class,
  y     = dmi,
  type = "noneparametric",
  p.adjust.method = "BH",
  title = "DMI Score of Firmicutes Cross Body Sites")
p20 <- p20 + scale_color_manual(values = body_site_color) + base_theme
p20
ggsave(filename = "./DMIScore_Firmicutes.Bodysite.pdf", p20, width = 6, height = 5, dpi=300)

p21 <- ggbetweenstats(
  data  = filter(DMIscore, Phylum== "Proteobacteria"),
  x     = class,
  y     = dmi,
  type = "noneparametric",
  p.adjust.method = "BH",
  title = "DMI Score of Proteobacteria Cross Body Sites")
p21 <- p21 + scale_color_manual(values = body_site_color) + base_theme
p21
ggsave(filename = "./DMIScore_Proteobacteria.Bodysite.pdf", p21, width = 6, height = 5, dpi=300)

p22 <- ggbetweenstats(
  data  = filter(DMIscore, Phylum== "Other"),
  x     = class,
  y     = dmi,
  type = "noneparametric",
  p.adjust.method = "BH",
  title = "DMI Score of Other Phylum Cross Body Sites")
p22 <- p22 + scale_color_manual(values = body_site_color) + base_theme
p22
ggsave(filename = "./DMIScore_other.Bodysite.pdf", p22, width = 6, height = 5, dpi=300)

#familyscore %>% 
#DMIscore %>% 

dim(familyscore)
dim(DMIscore)

table(DMIscore$class)
table(familyscore$class)
