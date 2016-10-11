rm(list=ls())

setwd("~/Dropbox/nimble/occupancy/analysis/singleSpp-multiSea")
source('src/initialize.R')
set.seed(444)
data <- genDynamicOccData()
model.input <- prepModDataOcc(data, include.zs=FALSE)

## *********************************************************************
##  Multi-season occupancy model: option 4-5 remove latent states using
##  user-defined NIMBLE function
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
## opt 4 run with compareMCMCs
## *********************************************************************

ss.ms.opt4 <- compareMCMCs(input1,
                           MCMCs=c('nimble', 'autoBlock', 'nimble_slice'),
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

ss.ms.opt5 <- compareMCMCs(input1,
                           MCMCs=c('nimbleOpt5'),
                           MCMCdefs = MCMCdefs.opt5,
                           niter=niter,
                           burnin = burnin,
                           summary=FALSE,
                           check=FALSE)

save(ss.ms.opt5, file=file.path(save.dir, "opt5.Rdata"))


## *********************************************************************
occ.R.model <- nimbleModel(code=ss.ms.occ,
                           constants=input1$constants,
                           data=input1$data,
                           inits=input1$inits,
                           check=FALSE)

occ.mcmc <- buildMCMC(occ.R.model)
occ.C.model <- compileNimble(occ.R.model)
occ.C.mcmc <- compileNimble(occ.mcmc, project = occ.R.model)
occ.C.mcmc$run(niter)

source('../cppp/src/calcCPPP.R', chdir = TRUE)
options(mc.cores=1)

test.opt4 <- generateCPPP(occ.R.model,
                          occ.C.model,
                          occ.C.mcmc,
                          occ.mcmc,
                          dataName = 'y',
                          paramNames = input1$monitors, 
                          MCMCIter = niter, 
                          NSamp = 10^3,
                          NPDist = 10^3,
                          burnInProportion = 0.10,
                          thin = 1,
                          averageParams = TRUE,
                          discFuncGenerator=likeDiscFuncGenerator)

save(test.opt4, file=file.path(save.dir, "ssms_noz_CPPP.Rdata"))

