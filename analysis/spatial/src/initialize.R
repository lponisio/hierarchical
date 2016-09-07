## library(devtools)
## install_github("nimble-dev/nimble",
##                ref = "devel",
##                subdir = "packages/nimble")

library(nimble)
library(igraph)
library(raster)

source("../all/plotting.R")
source("../all/runNimble.R")
source("src/dataGen.R")

## samplers
source("../all/samplers/sampler_z.R")
source("../all/samplers/sampler_reflective.R")

save.dir <-  "../../../saved/spatial/saved"
