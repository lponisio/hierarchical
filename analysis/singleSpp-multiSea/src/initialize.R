library(devtools)
## install_github("nimble-dev/nimble",
##                ref = "devel",
##                subdir = "packages/nimble")


## install_github("nimble-dev/nimble",
##                ref = "faster-vector-passing",
##                subdir = "packages/nimble")



library(nimble)
library(igraph)
source("src/dynamicOcc.R")
source("src/setup.R")
source("../all/plotting.R")
source("../all/samplers/sampler_crossLevel_new.R")
source("../all/samplers/leaveOutSampler.R")
source("../all/runNimble.R")

dir.create(file.path("../../../occupancy_saved/saved/singleSpp-multiSea/saved"),
           showWarnings = FALSE)
save.dir <-  "../../../occupancy_saved/saved/singleSpp-multiSea/saved"

## MCMC settings
scale <- 1e1
burnin <- 1e1*scale
niter <- (1e3)*scale
