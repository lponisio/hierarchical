rm(list=ls())
setwd("~/Dropbox/nimble/occupancy/analysis/singleSpp-multiSea")

source('src/initialize.R')
set.seed(444)
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
                                                 type ="slice",
                                                 print=FALSE)
  customSpec$removeSamplers('z')
  customSpec$addSampler('z', type = "sampler_latentSub",
                        control = list(leaveOutProportion = 0.7,
                          control = list()))
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

save(ss.ms.opt2, file=file.path(save.dir, "opt2_saddness.Rdata"))


## *********************************************************************
## model assessment
## *********************************************************************
## build model

## occ.R.model <- nimbleModel(code=ss.ms.occ,
##                            constants=input1$constants,
##                            data=input1$data,
##                            inits=input1$inits,
##                            check=FALSE)



## occ.mcmc <- buildMCMC(occ.R.model)
## occ.C.model <- compileNimble(occ.R.model)
## occ.C.mcmc <- compileNimble(occ.mcmc, project = occ.R.model)
## occ.C.mcmc$run(niter)

## as.matrix(occ.C.mcmc$mvSamples)

## source('../cppp/src/calcCPPP.R', chdir = TRUE)
## options(mc.cores=6)

## test.opt2 <- generateCPPP(occ.R.model,
##                           occ.C.model,
##                           occ.C.mcmc,
##                           occ.mcmc,
##                           dataName = 'y',
##                           paramNames = input1$monitors, 
##                           NSamp = 10^3,
##                           NPDist = 10^3,
##                           burnInProp = 0.10,
##                           thin = 1,
##                           averageParams = TRUE,
##                           discFuncGenerator=likeDiscFuncGenerator)

## save(test.opt2, file=file.path(save.dir, "ssms_CPPP.Rdata"))
