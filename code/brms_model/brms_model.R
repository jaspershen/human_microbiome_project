###library packages
if (!require("install.load")) {
  install.packages("install.load")
  library(install.load)
}

list.of.packages <- c("tidyverse", "skimr", "brms")
install.load::install_load(list.of.packages)

####----------------------------------------------------------------------------
####transformation data
##Continuous data (cytokine, microbe) were analyzed using means and standard deviations and
##categorical data (infection state, season) using counts and percents.
##The normality distribution assumption was evaluated visually with histograms.
##Both continuous data variables showed significant deviation from normality. 
##A and logit transformations were used in models that assume a 
##Gaussian distribution of errors (See Supplemental Figures). 
##To avoid infinite values after transformation, 
##zeros and ones were adjusted by adding or 
##subtracting half of the smallest observed value, as described previously.


masstools::setwd_project()
setwd("data_analysis/brms/")
rm(list = ls())

load("model-df.Rdata")

head(model.df)

###three Transformation functions
normalize <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

arcsine <- function(p) {
  asin(sqrt(p))
}

logit <- function(y) {
  log(y / (1 - y))
}


# Variable Summary
model.df %>%
  skimr::skim()


###Cytokine IL17F and Akkermansia as example
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
  dplyr::filter(SubjectID == "69-001") %>%
  ggplot(aes(IL17F)) +
  geom_density(size = 1) +
  theme_bw(base_size = 12) +
  labs(title = "69-001 only")

gridExtra::grid.arrange(p1, p2, ncol = 2)

save(model.trans, file = "model.trans.Rdata")

####----------------------------------------------------------------------------
####2. model
load("model.trans.Rdata")

###what is the cluster?
x <- model.trans %>%
  dplyr::filter(cluster != 2) %>%
  dplyr::mutate(cluster = fct_drop(cluster))

table(x$cluster)

# Model
# This should be modified to include all of the cytokines and genera of interest
# The format of the data frame should have each cytokine and genus as a column
cytokines <- c("IL17F")
genera <- c("Akkermansia")
################################


cytokine.model.output <- list()

for (c in cytokines) {
  genera.model.output <- list()
  for (g in genera) {
    f <- as.formula(paste(g, "~ days + (1|SubjectID) + cluster *", c,
                          sep = " "))
    brm <-
      brms::brm(f,
                data = x,
                sparse = TRUE,
                family = negbinomial())

    # Reformat output
    # brm.df <- broom::tidy(brm)
    brm.df <- broom.mixed::tidy(brm)

    # Create column to track the genus being modelled
    brm.df$genus <- g

    # # Create list of outputs for each genus
    # genera.model.output[[g]] <-
    #   brm.df %>%
    #   mutate(significant = !lower <= 0 & upper >= 0) %>%
    #   filter(grepl(c, brm.df$term))

    genera.model.output[[g]] <-
      brm.df %>%
      mutate(significant = !conf.low <= 0 & conf.high >= 0) %>%
      filter(grepl(c, brm.df$term))

  }

  cytokine.model.output[[c]] <- bind_rows(genera.model.output)
}


##The hypothesis is that interaction term will be significant more frequently than the cytokine fixed term alone.
save(brm, file = "model-output.RData")
save(cytokine.model.output, file = "cytokine.model.output")
