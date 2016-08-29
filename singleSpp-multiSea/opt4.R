rm(list=ls())
setwd("~/Dropbox/occupancy-nimble/singleSpp-multiSea")
source('src/initialize.R')

## *********************************************************************
##  Multi-season occupancy model: option 4 remove latent states using
##  user-defined NIMBLE function, and add a block sampler to phi[i-1]
##  and gamma[k-1]
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
## ## add block samplers
MCMCdefs.opt4 <- list('nimbleOpt4' = quote({
  customSpec <- configureMCMC(Rmodel)
  ## find node names for random effects
  sp.parms.phi <- Rmodel$getNodeNames(includeData = FALSE))[grepl("\\<phi\\>",
                                    Rmodel$getNodeNames(includeData = FALSE))]
  sp.parms.gam <- Rmodel$getNodeNames(includeData = FALSE))[grepl("[[gamma]]",
                                     Rmodel$getNodeNames(includeData = FALSE))]
  customSpec$removeSamplers(c(sp.parms.phi, sp.parms.gam),
                            print=FALSE)

  customSpec$addSampler(target = sp.parms.a,
                        type = "RW_block", log=TRUE)
  customSpec$addSampler(target = sp.parms.b,
                        type = "RW_block", log=TRUE)
  customSpec
}))


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
