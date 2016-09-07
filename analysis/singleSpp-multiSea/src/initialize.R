## library(devtools)
## install_github("nimble-dev/nimble",
##                ref = "dlevel",
##                subdir = "packages/nimble")

library(nimble)
library(mcmcplots)
library(igraph)
source("src/dynamicOcc.R")
source("src/dataGen.R")
source("../all/plotting.R")
source("../all/runNimble.R")

save.dir <- "../../../saved/singleSpp-multiSea/saved"

## samplers
source("../all/samplers/sampler_z.R")
source("../all/samplers/sampler_reflective.R")

## *********************************************************************
##  Dynamic (multi-season) site-occupancy models
## *********************************************************************
set.seed(4)
data <- genDynamicOccData()

## data zs with 0s set to NAs
zs <- apply(data$y, c(1, 3), max)
zs[zs == 0] <- NA

## initial condiations, NAs where 1s are in z, and 1s are where NA
zinits <- zs
zinits[zinits == 1] <- 2
zinits[is.na(zinits)] <- 1
zinits[zinits == 2] <- NA
inits <- list(z = zinits)

## model data
model.data <- list(y = data$y, z = zs)

## constants
constants <- list(nsite = dim(data$y)[1],
                  nrep = dim(data$y)[2],
                  nyear = dim(data$y)[3])

## parameters to monitor
monitors <- c("psi",
              "phi",
              "gamma",
              "p",
              "n.occ",
              "growthr",
              "turnover")


## MCMC settings
scale <- 1e1
burnin <- 1e1*scale
niter <- (1e3)*scale
