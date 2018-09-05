rm(list=ls())
setwd("~/Dropbox/occupancy")
setwd("analysis/nmixture")
source('src/initialize.R')

latent <- TRUE
hyper.param <- TRUE
source("src/setup.R")
source('src/models.R')

input1 <- c(code=nmixture,
            model.input)


## *********************************************************************
## original: vanilla nimble and JAGS
## *********************************************************************

nmixture.orig <- compareMCMCs(input1,
                           MCMCs=c('jags', 'nimble'),
                           niter=niter,
                           burnin = burnin,
                           summary=FALSE,
                           check=FALSE)



latent <- TRUE
hyper.param <- FALSE
source("src/setup.R")
source('src/models.R')

input1 <- c(code=nmixture,
            model.input)


nmixture.orig <- compareMCMCs(input1,
                           MCMCs=c('jags', 'nimble'),
                           niter=niter,
                           burnin = burnin,
                           summary=FALSE,
                           check=FALSE)


