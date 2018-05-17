## setwd("~/Dropbox/occupancy/")
rm(list=ls())
setwd("analysis/singleSpp-multiSea")
source('src/initialize.R')
set.seed(444)
data <- genDynamicOccData()
model.input <- prepModDataOcc(data, include.zs=FALSE)

## *********************************************************************
##  Multi-season occupancy model: remove latent states using
##  user-defined NIMBLE function dDynamicOcc
##  *********************************************************************

## Specify model in NIMBLE
ss.ms.occ <- nimbleCode({
  ##  priors
  psi1 ~ dunif(0, 1)

  psi[1] <- psi1
  for(k in 1:(nyear-1)){
    phi[k] ~ dunif(0, 1)
    gamma[k] ~ dunif(0, 1)
    p[k] ~ dunif(0, 1)
  }
  p[nyear] ~ dunif(0, 1)

  ## Ecological submodel: Define state conditional on parameters
  for(i in 1:nsite) {
    ## removes the z's and muZ's from the model and compute
    ## the probability of all reps over all years for one site.
    y[i, 1:nrep, 1:nyear] ~ dDynamicOccupancy(nrep,
                                              psi1,
                                              phi[1:(nyear-1)],
                                              gamma[1:(nyear-1)],
                                              p[1:nyear])
  }
})

input1 <- c(code=ss.ms.occ,
            model.input)

## *********************************************************************
## vanilla nimble and slice samplers
## *********************************************************************

ss.ms.filter <- compareMCMCs(input1,
                           MCMCs=c('nimble'),
                           niter=niter,
                           burnin = burnin,
                           summary=FALSE,
                           check=FALSE)

save(ss.ms.filter, file=file.path(save.dir, "filter.Rdata"))

## *********************************************************************
## phi and gam parameters together
## *********************************************************************

MCMCdefs.blocking <- list('blocking' = quote({
    customSpec <- configureMCMC(Rmodel)
        customSpec$removeSamplers(Rmodel$getNodeNames(includeData = FALSE), print=FALSE)
        customSpec$addSampler(target = Rmodel$getNodeNames(includeData = FALSE),
                              type = "AF_slice")
    customSpec
}))

ss.ms.filter.blocking <- compareMCMCs(input1,
                           MCMCs=c('blocking'),
                           MCMCdefs = MCMCdefs.blocking,
                           niter=niter,
                           burnin = burnin,
                           summary=FALSE,
                           check=FALSE)

save(ss.ms.filter.blocking, file=file.path(save.dir, "filter_blocking.Rdata"))
