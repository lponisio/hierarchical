rm(list=ls())
setwd("~/Dropbox/nimble-dev/occupancy/singleSpp-multiSea")
source('src/initialize.R')

## *********************************************************************
##  Multi-season occupancy model: option 3 remove latent states using
##  user-defined NIMBLE function
##  *********************************************************************

## Specify model in NIMBLE
ss.ms.occ <- nimbleCode({
  ##  priors
  psi1 ~ dunif(0, 1)

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

  ## Derived parameters: Sample and population occupancy, growth rate
  ## and turnover
  psi[1] <- psi1
  n.occ[1]<- sum(z[1:nsite,1])
  for (k in 2:nyear){
    psi[k] <- psi[k-1]*phi[k-1] + (1-psi[k-1])*gamma[k-1]
    n.occ[k] <- sum(z[1:nsite,k])
    growthr[k-1] <- psi[k]/psi[k-1]
    turnover[k-1] <- (1 - psi[k-1]) * gamma[k-1]/psi[k]
  }
})

## *********************************************************************
## run with compareMCMCs

input1 <- list(code=ss.ms.occ,
               constants=constants,
               data=model.data,
               inits=inits)


ss.ms.opt3 <- compareMCMCs(input1,
                           MCMCs=c('nimble', 'autoBlock'),
                           niter=niter,
                           burnin = burnin,
                           summary=FALSE,
                           check=FALSE)

save(ss.ms.opt3, file="saved/opt3.Rdata")
