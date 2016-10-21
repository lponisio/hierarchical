## library(devtools)
## install_github("nimble-dev/nimble",
##                ref = "devel",
##                subdir = "packages/nimble")

library(nimble)
library(igraph)
library(reshape)
source("../all/plotting.R")
source("../all/runNimble.R")
source('../all/samplers/sampler_RW_shift.R')
source("src/reformatData.R")
source("src/multispeciesOcc.R")

save.dir <-  "../../../saved/multiSpp-singleSea/saved"

survey.data <- read.csv("data/occupancy_data.csv")
species.groups <- read.csv("data/species_groups.csv")
survey.dates <- read.csv("data/survey_dates.csv")
habitat <- read.csv("data/habitat.csv")

## mcmc settings
scale <- 5e3
burnin <- 1e2*scale
niter <- (1e3)*scale

