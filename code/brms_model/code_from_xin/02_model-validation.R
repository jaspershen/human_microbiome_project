library(magrittr)
library(brms)


# Markdown syntax
masstools::setwd_project()
setwd("data_analysis/brms/")
load("model-output.RData")

m <- brm
rm(brm)

# Change `days` back to numeric
m[["data"]]$days <- 
  m[["data"]]$days %>%
  as.character() %>%
  as.numeric()

# Model evaluation based on [the brms vignette https://cran.r-project.org/web/packages/brms/vignettes/brms_nonlinear.html]

summary(m)

# Plot population-level effects
plot(m, pars = "^b")

plot(marginal_effects(m), points = TRUE)

pp_check(m)

loo(m)


## Model evaluation based on [Mi et al., PLoS One, 2015. Goodness-of-Fit Tests and Model Diagnostics for Negative Binomial Regression of RNA Sequencing Data https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4365073/]
