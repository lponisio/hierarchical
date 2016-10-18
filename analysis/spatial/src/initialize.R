library(devtools)
install_github("nimble-dev/nimble",
               ref = "devel",
               subdir = "packages/nimble")

library(nimble)
library(igraph)
library(raster)

source("../all/plotting.R")
source("../all/runNimble.R")
source("src/dataGen.R")

save.dir <-  "../../../saved/spatial/saved"

## MCMC settings
scale <- 1e2
burnin <- 1e1*scale
niter <- (1e3)*scale
