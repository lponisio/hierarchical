library(devtools)
## install_github("nimble-dev/nimble",
##                ref = "devel",
##                subdir = "packages/nimble")

library(nimble)
library(igraph)
source("src/dNmixture.R")
source("../all/plotting.R")
## source("src/plotting.R")

dir.create(file.path("../../../occupancy_saved/saved/nmixture/saved"),
           showWarnings = FALSE)
save.dir <-  "../../../occupancy_saved/saved/nmixture/saved"

## mcmc settings
scale <- 1
burnin <- 1e2*scale
niter <- (1e3)*scale

