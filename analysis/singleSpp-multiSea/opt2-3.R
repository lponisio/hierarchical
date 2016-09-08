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
  
  psi[1] <- psi1
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
## model assessment
## *********************************************************************
## build model
R.model <- nimbleModel(code=ss.ms.occ,
                       constants=input1$constants,
                       data=input1$data,
                       inits=input1$inits,
                       check=FALSE)
message('R model created')

customSpec <- configureMCMC(R.model,
                            monitors=input1$monitors)
customSpec$removeSamplers('phi', print=FALSE)
customSpec$removeSamplers('gamma', print=FALSE)
customSpec$removeSamplers('p', print=FALSE)
customSpec$removeSamplers('psi1', print=FALSE)
## happens to be all top nodes
zeroOneNodes <- R.model$getNodeNames(topOnly = TRUE)
for(zon in zeroOneNodes) customSpec$addSampler(target = zon,
                                                 type =
                                                   "slice",
                                                 print=FALSE)

mcmc <- buildMCMC(customSpec)
message('MCMC built')

## compile model in C++
C.model <- compileNimble(R.model)
C.mcmc <- compileNimble(mcmc, project = R.model)
message('NIMBLE model compiled')

source('../cppp/src/calcCPPP.R', chdir = TRUE)
options(mc.cores=2)

test.opt2 <- generateCPPP(R.model,
             C.model,
             C.mcmc,
             mcmc,
             dataName = 'y',
             paramNames = input1$monitors, 
             MCMCIter = 1000, 
             NSamp = 100,
             NPDist = 10,
             thin = 1)


list(code=ss.ms.occ,
               constants=constants,
               data=model.data,
               inits=inits)


