rm(list=ls())
gctorture()
setwd("~/Dropbox/nimble/occupancy/analysis/singleSpp-multiSea")
source('src/initialize.R')

## *********************************************************************
##  Multi-season occupancy model: option 4-5 remove latent states using
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
## opt 4 run with compareMCMCs
## *********************************************************************

input1 <- list(code=ss.ms.occ,
               constants=constants,
               data=model.data,
               inits=inits)

monitors <- c("phi", "gamma", "psi")

ss.ms.opt4 <- compareMCMCs(input1,
                           MCMCs=c('nimble', 'autoBlock'),
                           niter=niter,
                           burnin = burnin,
                           summary=FALSE,
                           check=FALSE)

save(ss.ms.opt4, file=file.path(save.dir, "opt4.Rdata"))


## *********************************************************************
## opt 5 add block samplers
## *********************************************************************
MCMCdefs.opt5 <- list('nimbleOpt5' = quote({
  customSpec <- configureMCMC(Rmodel)
  ## find node names for random effects
  parms.phi <- Rmodel$getNodeNames(includeData = FALSE)[grepl("^phi",
                                     Rmodel$getNodeNames(includeData = FALSE))]
  parms.gam <- Rmodel$getNodeNames(includeData = FALSE)[grepl("^gamma",
                                     Rmodel$getNodeNames(includeData =
                                                         FALSE))]
  phi.gam <- cbind(parms.phi, parms.gam)[-1,]
  for(i in 1:nrow(phi.gam)){
    customSpec$removeSamplers(phi.gam[i,], print=FALSE)
    customSpec$addSampler(target = phi.gam[i,],
                          type = "RW_block")
  }
  customSpec
}))


## *********************************************************************
## run with compareMCMCs
## *********************************************************************

input1 <- list(code=ss.ms.occ,
               constants=constants,
               data=model.data,
               inits=inits)


ss.ms.opt5 <- compareMCMCs(input1,
                           MCMCs=c('nimbleOpt5'),
                           MCMCdefs = MCMCdefs.opt5,
                           niter=niter,
                           burnin = burnin,
                           summary=FALSE,
                           check=FALSE)

save(ss.ms.opt5, file=file.path(save.dir, "opt5.Rdata"))


## *********************************************************************
options(error=recover)
## build model
R.model <- nimbleModel(code=ss.ms.occ,
                       constants=constants,
                       data=model.data,
                       inits=inits,
                       check=FALSE)
message('R model created')


## configure and build mcmc
mcmc.spec <- configureMCMC(R.model,
                           print=FALSE,
                           monitors = monitors,
                           thin=1)
mcmc <- buildMCMC(mcmc.spec)
message('MCMC built')

## compile model in C++
C.model <- compileNimble(R.model)
C.mcmc <- compileNimble(mcmc, project = R.model)
message('NIMBLE model compiled')

source('../cppp/src/calcCPPP.R', chdir = TRUE)
options(mc.cores=6)

generateCPPP(R.model,
             C.model,
             C.mcmc,
             mcmc,
             dataName = 'y',
             paramNames = monitors, 
             MCMCIter = 1000, 
             NSamp = 1000,
             NPDist = 100,
             thin = 1)


list(code=ss.ms.occ,
               constants=constants,
               data=model.data,
               inits=inits)
