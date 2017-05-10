library(devtools)
install_github("nimble-dev/nimble",
               ref = "devel",
               subdir = "packages/nimble")

library(nimble)
library(igraph)
source("src/dynamicOcc.R")
source("src/dataGen.R")
source("../all/plotting.R")
source("../all/samplers/sampler_crossLevel_new.R")
source("../all/samplers/leaveOutSampler.R")
source("../all/runNimble.R")

save.dir <- "../../../saved/singleSpp-multiSea/saved"

## MCMC settings
scale <- 1e1
burnin <- 1e1*scale
niter <- (1e3)*scale
