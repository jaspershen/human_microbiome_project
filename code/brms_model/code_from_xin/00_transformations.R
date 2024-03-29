####this code is used to transform all the microbiome and cytokine data

if (!require("install.load")) {
  install.packages("install.load")
  library(install.load)
}

# CRAN packages
list.of.packages <- c("tidyverse", "skimr", "broom")

install.load::install_load(list.of.packages)

setwd(masstools::get_project_wd())
setwd("data_analysis/brms/")
load("model-df.Rdata")

# Transformation functions
normalize <- function(x){
  (x - min(x)) / (max(x) - min(x))
} 

arcsine <- function(p) asin(sqrt(p))

logit <- function(y) log(y / (1 - y))

# Variable Summary
model.df %>%
  skimr::skim()


# Cytokine IL17F distributation
model.df %>%
  mutate(no.trans = normalize(IL17F)) %>%
  dplyr::select(no.trans) %>%
  mutate(no.trans = ifelse(no.trans == 0, no.trans + (min(no.trans[no.trans > 0])/2), no.trans)) %>%
  mutate(no.trans = ifelse(no.trans == 1, no.trans - (min(no.trans[no.trans > 0])/2), no.trans)) %>%
  mutate(arcsine.trans = arcsine(no.trans)) %>%
  mutate(logit.trans = logit(no.trans)) %>%
  gather %>%
  ggplot(aes(value)) +
  geom_density() +
  facet_wrap(~key, scales = "free") +
  theme_bw()


# Akkermansia
model.df %>%
  mutate(no.trans = normalize(Akkermansia)) %>%
  dplyr::select(no.trans) %>%
  mutate(no.trans = ifelse(no.trans == 0, no.trans + (min(no.trans[no.trans > 0])/2), no.trans)) %>%
  mutate(no.trans = ifelse(no.trans == 1, no.trans - (min(no.trans[no.trans > 0])/2), no.trans)) %>%
  mutate(arcsine = arcsine(no.trans)) %>%
  mutate(logit = logit(no.trans)) %>%
  mutate(log = log(no.trans)) %>%
  rename("Untransformed" = "no.trans") %>%
  gather(key, value) %>%
  mutate(key = fct_relevel(key, "Untransformed")) %>%
  ggplot(aes(value)) +
  geom_density() +
  facet_wrap(~key, scales = "free") +
  theme_bw()


## Output final transformations
model.trans <- 
  model.df %>%
  mutate(Akk.norm = normalize(Akkermansia)) %>%
  mutate(Akk.norm = ifelse(Akk.norm == 0, Akk.norm + (min(Akk.norm[Akk.norm > 0])/2), Akk.norm)) %>%
  mutate(Akk.norm = ifelse(Akk.norm == 1, Akk.norm - (min(Akk.norm[Akk.norm > 0])/2), Akk.norm)) %>%
  mutate(Akk.logit = logit(Akk.norm)) %>%
  dplyr::select(-Akk.norm) %>%
  mutate(il17.norm = normalize(IL17F)) %>%
  mutate(il17.norm = ifelse(il17.norm == 0, il17.norm + (min(il17.norm[il17.norm > 0])/2), il17.norm)) %>%
  mutate(il17.norm = ifelse(il17.norm == 1, il17.norm - (min(il17.norm[il17.norm > 0])/2), il17.norm)) %>%
  mutate(il17.logit = logit(il17.norm)) %>%
  dplyr::select(-il17.norm)


p1 <- 
  model.df %>%
  ggplot(aes(IL17F)) +
  geom_density(size = 1) +
  theme_bw(base_size = 12) +
  labs(title = "All individuals")

p2 <- model.df %>%
  filter(SubjectID == "69-001") %>%
  ggplot(aes(IL17F)) +
  geom_density(size = 1) +
  theme_bw(base_size = 12) +
  labs(title = "69-001 only")

gridExtra::grid.arrange(p1, p2, ncol = 2)


save(model.trans, file = "model.trans.Rdata")
