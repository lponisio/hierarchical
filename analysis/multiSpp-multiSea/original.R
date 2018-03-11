## ************************************************************
## setwd('~/Dropbox/occupancy')
rm(list=ls())
setwd('analysis/multiSpp-multiSea')

source("src/initialize.R")
load("data/5-0-350.Rdata")
source("src/complete_allInt.R")

input1 <- c(code=ms.ms.occ,
            model.input)

## *****************************************************************
## vanilla nimble/jags
## *****************************************************************

ms.ms.nimble <- compareMCMCs_withMonitors(input1,
                                          MCMCs=c('nimble', 'jags'),
                                          niter=niter,
                                          burnin = burnin,
                                          summary=FALSE,
                                          check=FALSE,
                                          monitors=model.input$monitors)

save(ms.ms.orig, file=file.path(save.dir, 'orig.Rdata'))

## *********************************************************************
## sampler only a subset of latent states
## *********************************************************************

MCMCdefs.subsamp <- list('nimbleSubsamp' = quote({
    customSpec <- configureMCMC(Rmodel)
    customSpec$removeSamplers('Z')
    customSpec$addSampler('Z', type = 'sampler_latentSub',
                          control = list(leaveOutProportion = 0.05,
                                         control = list()))
    customSpec
}))

ms.ms.subsamp <- compareMCMCs_withMonitors(input1,
                              MCMCs=c('nimbleSubsamp'),
                              MCMCdefs = MCMCdefs.subsamp,
                              niter= niter,
                              burnin = burnin,
                              monitors=model.input$monitors,
                              summary=FALSE,
                              check=FALSE)

save(ms.ms.subsamp, file=file.path(save.dir, 'subsamp.Rdata'))


## *********************************************************************
## cross level sampler
## *********************************************************************
## sample latent and top-level parameters jointly

