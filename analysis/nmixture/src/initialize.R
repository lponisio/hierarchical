library(devtools)
## install_github("nimble-dev/nimble",
##                ref = "devel",
##                subdir = "packages/nimble")

library(nimble)
library(igraph)
source("../all/plotting.R")
source("../all/runNimble.R")
source("../all/comparMCMCs_withMonitors.R")
source("../all/samplers/sampler_crossLevel_new.R")
source("../all/samplers/leaveOutSampler.R")

dir.create(file.path("../../../occupancy_saved/saved/nmixture/saved"),
           showWarnings = FALSE)
save.dir <-  "../../../occupancy_saved/saved/nmixture/saved"

## mcmc settings
scale <- 1e2
burnin <- 1e2*scale
niter <- (1e3)*scale

