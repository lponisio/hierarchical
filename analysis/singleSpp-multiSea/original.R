## setwd("~/Dropbox/occupancy")
rm(list=ls())
setwd("analysis/singleSpp-multiSea")
source('src/initialize.R')

set.seed(444)
data <- genDynamicOccData()
model.input <- prepModDataOcc(data)


## *********************************************************************
##  original multi-season occupancy : JAGS & NIMBLE
## *********************************************************************

ss.ms.occ <- nimbleCode({
  ## Specify priors
  psi1 ~ dunif(0, 1)

  for(year in 1:(nyear-1)){
    phi[year] ~ dunif(0, 1)
    gamma[year] ~ dunif(0, 1)
    p[year] ~ dunif(0, 1)
  }
  p[nyear] ~ dunif(0, 1)

  ## Ecological submodel: Define state conditional on parameters
  for (site in 1:nsite){
    z[site,1] ~ dbern(psi1)
    for (year in 2:nyear){
        muZ[site,year]<- z[site,year-1]*phi[year-1] +
            (1-z[site,year-1])*gamma[year-1]
      z[site,year] ~ dbern(muZ[site,year])
    }
  }

  ## Observation model
  for (site in 1:nsite){
    for (rep in 1:nrep){
      for (year in 1:nyear){
        muy[site,rep,year] <- z[site,year]*p[year]
        y[site,rep,year] ~ dbern(muy[site,rep,year])
      }
    }
  }

})

input1 <- c(code=ss.ms.occ,
            model.input)


## *********************************************************************
## original: vanilla nimble and JAGS
## *********************************************************************

ss.ms.orig <- compareMCMCs(input1,
                           MCMCs=c('jags', 'nimble'),
                           niter=niter,
                           burnin = burnin,
                           summary=FALSE,
                           check=FALSE)

save(ss.ms.orig, file=file.path(save.dir, "orig.Rdata"))


## *********************************************************************
## sampler only a subset of latent states
## *********************************************************************

MCMCdefs.subsamp <- list('nimbleSubsamp' = quote({
    customSpec <- configureMCMC(Rmodel)
    customSpec$removeSamplers('z')
    customSpec$addSampler('z', type = 'sampler_latentSub',
                          control = list(leaveOutProportion = 0.6,
                                         control = list()))
    customSpec
}))

ss.ms.subsamp <- compareMCMCs(input1,
                           MCMCs=c('nimbleSubsamp'),
                           MCMCdefs = MCMCdefs.subsamp,
                           niter= niter,
                           burnin = burnin,
                           summary=FALSE,
                           check=FALSE)

save(ss.ms.subsamp, file=file.path(save.dir, 'subsamp.Rdata'))



## *********************************************************************
## block together phi and gamma
## *********************************************************************

MCMCdefs.blocking <- list('blocking' = quote({
    customSpec <- configureMCMC(Rmodel)
    ## find node names for random effects
    for(i in 1:nrow(phi.gam)){
        customSpec$removeSamplers(phi.gam[i,], print=FALSE)
        customSpec$addSampler(target = phi.gam[i,],
                              type = "AF_slice")
    }
    customSpec
}))

ss.ms.blocking <- compareMCMCs(input1,
                           MCMCs=c('blocking'),
                           MCMCdefs = MCMCdefs.blocking,
                           niter=niter,
                           burnin = burnin,
                           summary=FALSE,
                           check=FALSE)

save(ss.ms.blocking, file=file.path(save.dir, "blocking.Rdata"))


