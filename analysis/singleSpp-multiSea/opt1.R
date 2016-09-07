rm(list=ls())
setwd("~/Dropbox/nimble/occupancy/analysis/singleSpp-multiSea")
source('src/initialize.R')

data <- genDynamicOccData()
model.input <- prepModDataOcc(data)

## *********************************************************************
##  Multi-season occupancy model: option 1: user defined distribution
##  that computes all dbern(psi)s in one loop
##  *********************************************************************

## specify model in NIMBLE
ss.ms.occ <- nimbleCode({
  ## Specify priors
  psi1 ~ dunif(0, 1)

  for(k in 1:(nyear-1)){
    phi[k] ~ dunif(0, 1)
    gamma[k] ~ dunif(0, 1)
    p[k] ~ dunif(0, 1)
  }
  p[nyear] ~ dunif(0, 1)

  for (i in 1:nsite){
    z[i,1] ~ dbern(psi1)
    for (k in 2:nyear){
      muZ[i,k]<- z[i,k-1]*phi[k-1] + (1-z[i,k-1])*gamma[k-1]
      z[i,k] ~ dbern(muZ[i,k])
    }
  }

  for (i in 1:nsite){
      for (k in 1:nyear){
          muy[i,k] <- z[i,k]*p[k]
          y[i,1:nrep,k] ~ dbern_vec(muy[i,k], nrep)
    }
  }
  psi[1] <- psi1
})

## *********************************************************************
## run with compareMCMCs

input1 <- c(model.input,
            code=ss.ms.occ)

ss.ms.opt1 <- compareMCMCs(input1,
                                MCMCs=c('nimble'),
                                niter=niter,
                                burnin = burnin,
                                summary=FALSE,
                           check=FALSE)

save(ss.ms.opt1, file=file.path(save.dir, "opt1.Rdata"))
