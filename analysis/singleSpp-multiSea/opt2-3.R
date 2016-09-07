rm(list=ls())
setwd("~/Dropbox/nimble/occupancy/analysis/singleSpp-multiSea")
source('src/initialize.R')

data <- genDynamicOccData()
model.input <- prepModDataOcc(data)

## *********************************************************************
##  Multi-season occupancy model: custom z sampler
## *********************************************************************

ss.ms.occ <- nimbleCode({
  ## Specify priors
  psi1 ~ dunif(0, 1)

  for(k in 1:(nyear-1)){
    phi[k] ~ dunif(0, 1)
    gamma[k] ~ dunif(0, 1)
    p[k] ~ dunif(0, 1)
  }
  p[nyear] ~ dunif(0, 1)

  ## Ecological submodel: Define state conditional on parameters
  for (i in 1:nsite){
    z[i,1] ~ dbern(psi1)
    for (k in 2:nyear){
      muZ[i,k]<- z[i,k-1]*phi[k-1] + (1-z[i,k-1])*gamma[k-1]
      z[i,k] ~ dbern(muZ[i,k])
    }
  }

  ## Observation model
  for (i in 1:nsite){
    for (j in 1:nrep){
      for (k in 1:nyear){
        muy[i,j,k] <- z[i,k]*p[k]
        y[i,j,k] ~ dbern(muy[i,j,k])
      }
    }
  }
  psi[1] <- psi1
})

input1 <- c(code=ss.ms.occ,
           model.input)

## *********************************************************************
## opt 2: add custom z sampler and slice on uniform(0,1) nodes
## *********************************************************************

MCMCdefs.opt2 <- list('nimbleOpt2' = quote({
  customSpec <- configureMCMC(Rmodel)
  customSpec$removeSamplers('phi', print=FALSE)
  customSpec$removeSamplers('gamma', print=FALSE)
  customSpec$removeSamplers('p', print=FALSE)
  customSpec$removeSamplers('psi1', print=FALSE)
  ## happens to be all top nodes
  zeroOneNodes <- Rmodel$getNodeNames(topOnly = TRUE)
  for(zon in zeroOneNodes) customSpec$addSampler(target = zon,
                                                 type =
                                                   "slice",
                                                 print=FALSE)
  customSpec
}))

## *********************************************************************
## run with compareMCMCs

ss.ms.opt2 <- compareMCMCs(input1,
                           MCMCs=c('nimbleOpt2'),
                           MCMCdefs = MCMCdefs.opt2,
                           niter= niter,
                           burnin = burnin,
                           summary=FALSE,
                           check=FALSE)

save(ss.ms.opt2, file=file.path(save.dir, "opt2.Rdata"))


## *********************************************************************
## opt3:  add custom z sampler and reflective samplers
## *********************************************************************

MCMCdefs.opt3 <- list('nimbleOpt3' = quote({
  customSpec <- configureMCMC(Rmodel)
  customSpec$removeSamplers('phi', print=FALSE)
  customSpec$removeSamplers('gamma', print=FALSE)
  customSpec$removeSamplers('p', print=FALSE)
  customSpec$removeSamplers('psi1', print=FALSE)
  ## happens to be all top nodes
  zeroOneNodes <- Rmodel$getNodeNames(topOnly = TRUE)
  for(zon in zeroOneNodes) customSpec$addSampler(target = zon,
                                                 type =
                                                   sampler_RW_reflect,
                                                 print=FALSE)
  customSpec
}))


## *********************************************************************
## run with compareMCMCs

ss.ms.opt3 <- compareMCMCs(input1,
                           MCMCs=c('nimbleOpt3'),
                           MCMCdefs = MCMCdefs.opt3,
                           niter= niter,
                           burnin = burnin,
                           summary=FALSE,
                           check=FALSE)

save(ss.ms.opt3, file=file.path(save.dir, "opt3.Rdata"))

