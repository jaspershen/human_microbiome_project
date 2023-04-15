
if (!require("install.load")) {
  install.packages("install.load")
  library(install.load)
}

# CRAN packages
list.of.packages <- c("tidyverse", "brms", "broom")

install.load::install_load(list.of.packages)

masstools::setwd_project()
setwd("data_analysis/brms/")

load("model-df.Rdata")
load("model.trans.Rdata")

model.trans$cluster <- relevel(model.trans$cluster, ref = 3)

###what is cluster means and why remove cluster 2
x <- model.trans %>%
  filter(cluster != 2) %>%
  mutate(cluster = fct_drop(cluster))

table(x$cluster)

cluster.model.output <- list()

genus = c("Akkermansia")

for (g in genus) {
  brm <- brms::brm(Akkermansia ~ IL17F + days + cluster + (1|SubjectID),
                   data = x, sparse = TRUE, family = negbinomial())
  
  brm.df <- broom.mixed::tidy(brm)
  
  cluster.model.output <- 
    brm.df %>%
    filter()
}

# brm00 <- brms::brm(Akkermansia ~ IL17F + days + cluster + (1|SubjectID),
#                    data = x, sparse = TRUE, family = negbinomial())
# summary(brm05)
# 
# mdf <- tidy(brm05)
# 
# ci.eval <- !mdf$lower <= 0 & mdf$upper >=0 
# mdf$term[ci.eval]
# 
# plot(brm05)

# brm06 <- brms::brm(Akkermansia ~ IL17F + days + cluster + (1|SubjectID),
#                    data = x, sparse = TRUE, family = negbinomial())
# 
# summary(brm06)

brm03 <- brms::brm(Akkermansia ~ IL17F + (1 | cluster:SubjectID),
                   data = x, sparse = TRUE, family = negbinomial())
plot(brm03)

summary(brm03)

plot(marginal_effects(brm03))
zz <- broom.mixed::tidy(brm03)
zz

# brm08 <- brms::brm(Akkermansia ~ IL17F + cluster + (1 | SubjectID),
#                    data = x, sparse = TRUE, family = negbinomial())
# 
# summary(brm08)
# 
# zz <- broom.mixed::tidy(brm08)
# 
# lme00 <- lmer(Akkermansia ~ IL17F + cluster + days + (1 | SubjectID),
#               data = x)
# 
# summary(lme00)



