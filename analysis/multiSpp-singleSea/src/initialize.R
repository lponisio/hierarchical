library(devtools)
## install_github("nimble-dev/nimble",
##                ref = "devel",
##                subdir = "packages/nimble")

library(nimble)
library(igraph)
library(reshape)
source("../all/plotting.R")
source("../all/runNimble.R")
source("../all/comparMCMCs_withMonitors.R")
source("../all/samplers/sampler_crossLevel_new.R")
source("../all/samplers/leaveOutSampler.R")
source("src/reformatData.R")
source("src/multispeciesOcc.R")

dir.create(file.path("../../../occupancy_saved/saved/multiSpp-singleSea/saved"),
           showWarnings = FALSE)
save.dir <-  "../../../occupancy_saved/saved/multiSpp-singleSea/saved"

survey.data <- read.csv("data/occupancy_data.csv")
species.groups <- read.csv("data/species_groups.csv")
survey.dates <- read.csv("data/survey_dates.csv")
habitat <- read.csv("data/habitat.csv")

## mcmc settings
scale <- 1e3
burnin <- 1e2*scale
niter <- (1e3)*scale

monitors <- c('mu.a1','mu.a2','mu.a3','mu.a4',
              'sigma.a1','sigma.a2','sigma.a3','sigma.a4',
              'cato.occ.mean', 'fcw.occ.mean',
              'cato.det.mean', 'fcw.det.mean',
              'sigma.ucato', 'sigma.vcato', 'sigma.ufcw',
              'sigma.vfcw', 'mu.b1', 'mu.b2', 'sigma.b1',
              'sigma.b2', 'u.cato', 'u.fcw',
              'a1', 'a2', 'a3', 'a4', 'v.cato', 'v.fcw', 'b1', 'b2')
