rm(list=ls())

setwd("~/Dropbox/nimble/occupancy/analysis/singleSpp-multiSea")
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
                           MCMCs=c('nimble', 'nimble_slice'),
                           niter=niter,
                           burnin = burnin,
                           summary=FALSE,
                           check=FALSE)

save(ss.ms.filter, file=file.path(save.dir, "filter.Rdata"))


## *********************************************************************
## block together phi and gamma
## *********************************************************************

## MCMCdefs.opt5 <- list('nimbleOpt5' = quote({
##   customSpec <- configureMCMC(Rmodel)
##   ## find node names for random effects
##   parms.phi <- Rmodel$getNodeNames(includeData = FALSE)[grepl("^phi",
##                                      Rmodel$getNodeNames(includeData = FALSE))]
##   parms.gam <- Rmodel$getNodeNames(includeData = FALSE)[grepl("^gamma",
##                                      Rmodel$getNodeNames(includeData =
##                                                          FALSE))]
##   phi.gam <- cbind(parms.phi, parms.gam)[-1,]
##   for(i in 1:nrow(phi.gam)){
##     customSpec$removeSamplers(phi.gam[i,], print=FALSE)
##     customSpec$addSampler(target = phi.gam[i,],
##                           type = "RW_block")
##   }
##   customSpec
## }))

## ss.ms.opt5 <- compareMCMCs(input1,
##                            MCMCs=c('nimbleOpt5'),
##                            MCMCdefs = MCMCdefs.opt5,
##                            niter=niter,
##                            burnin = burnin,
##                            summary=FALSE,
##                            check=FALSE)

## save(ss.ms.opt5, file=file.path(save.dir, "opt5.Rdata"))


