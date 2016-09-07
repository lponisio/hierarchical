## library(devtools)
## install_github("nimble-dev/nimble",
##                ref = "dlevel",
##                subdir = "packages/nimble")

library(nimble)
library(igraph)
source("src/dynamicOcc.R")
source("src/dataGen.R")
source("../all/plotting.R")
source("../all/runNimble.R")

save.dir <- "../../../saved/singleSpp-multiSea/saved"

## samplers
source("../all/samplers/sampler_z.R")
source("../all/samplers/sampler_reflective.R")

## MCMC settings
scale <- 1e1
burnin <- 1e1*scale
niter <- (1e3)*scale
