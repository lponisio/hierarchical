## setwd("~/Dropbox/occupancy")
rm(list=ls())
setwd("analysis/singleSpp-multiSea")
source('src/initialize.R')

set.seed(444)
sim.input <- genDynamicOccData()
model.input <- prepModDataOcc(sim.input)

## *********************************************************************
##  original multi-season occupancy : JAGS & NIMBLE
## *********************************************************************

ss.ms.occ <- nimbleCode({
    ## Specify priors
    psi1 ~ dunif(0, 1)

    mu.p     ~ dnorm(0,0.001)
    sigma.p     ~ dunif(0,100)
    tau.p <- 1/(sigma.p*sigma.p)

    mu.phi  ~ dnorm(0,0.001)
    mu.gamma  ~ dnorm(0,0.001)
    sigma.phi ~ dunif(0,100)
    sigma.gamma ~ dunif(0,100)
    tau.phi <-  1/(sigma.phi*sigma.phi)
    tau.gamma <-  1/(sigma.gamma*sigma.gamma)

    for(year in 1:(nyear -1)) {
        p[year]   ~ dnorm(mu.p,     tau.p)
        phi[year] ~ dnorm(mu.phi, tau.phi)
        gamma[year] ~ dnorm(mu.gamma, tau.gamma)
    }
    p[nyear]      ~ dnorm(mu.p,     tau.p)


    ## Ecological submodel: Define state conditional on parameters
    for (site in 1:nsite){
        z[site,1] ~ dbern(psi1)
        for (year in 2:nyear){
            logit(muZ[site,year]) <- z[site,year-1]*phi[year-1] +
                (1-z[site,year-1])*gamma[year-1]
            z[site,year] ~ dbern(muZ[site,year])
        }
    }

    ## Observation model
    for (site in 1:nsite){
        for (rep in 1:nrep){
            for (year in 1:nyear){
                logit(muy[site,rep,year]) <- z[site,year]*p[year]
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
                           check=TRUE)

save(ss.ms.orig, file=file.path(save.dir, "orig.Rdata"))

## *********************************************************************
## block together phi and gamma, af slice
## *********************************************************************

MCMCdefs.AFSS.block <- list('AFSS_block' = quote({
    customSpec <- configureMCMC(Rmodel)
    parms.phi <- Rmodel$getNodeNames(
                            includeData = FALSE)[grepl("^phi",
                                                       Rmodel$getNodeNames(includeData = FALSE))]
    parms.gam <- Rmodel$getNodeNames(
                            includeData = FALSE)[grepl("^gamma",
                                                       Rmodel$getNodeNames(includeData = FALSE))]
    phi.gam <- cbind(parms.phi, parms.gam)[-1,]
    ## find node names for random effects
    for(i in 1:nrow(phi.gam)){
        customSpec$removeSamplers(phi.gam[i,], print=FALSE)
        customSpec$addSampler(target = phi.gam[i,],
                              type = "AF_slice")
    }
    customSpec
}))

ss.ms.AFSblocking <- compareMCMCs(input1,
                                     MCMCs=c('AFSS_block'),
                                     MCMCdefs = MCMCdefs.AFSS.block,
                                     niter=niter,
                                     burnin = burnin,
                                     summary=FALSE,
                                     check=FALSE)

save(ss.ms.AFSblocking, file=file.path(save.dir, "AFSS_block.Rdata"))


## *********************************************************************
## block together phi and gamma, rw block
## *********************************************************************

MCMCdefs.RW.block <- list('RW_block' = quote({
    customSpec <- configureMCMC(Rmodel)
    parms.phi <- Rmodel$getNodeNames(
                            includeData = FALSE)[grepl("^phi",
                                                       Rmodel$getNodeNames(includeData = FALSE))]
    parms.gam <- Rmodel$getNodeNames(
                            includeData = FALSE)[grepl("^gamma",
                                                       Rmodel$getNodeNames(includeData = FALSE))]
    phi.gam <- cbind(parms.phi, parms.gam)[-1,]
    ## find node names for random effects
    for(i in 1:nrow(phi.gam)){
        customSpec$removeSamplers(phi.gam[i,], print=FALSE)
        customSpec$addSampler(target = phi.gam[i,],
                              type = "RW_block")
    }
    customSpec
}))

ss.ms.RWblocking <- compareMCMCs(input1,
                                     MCMCs=c('RW_block'),
                                     MCMCdefs = MCMCdefs.RW.block,
                                     niter=niter,
                                     burnin = burnin,
                                     summary=FALSE,
                                     check=FALSE)

save(ss.ms.RWblocking, file=file.path(save.dir, "RW_block.Rdata"))



## *********************************************************************
## sampler only a subset of latent states
## *********************************************************************

## MCMCdefs.subsamp <- list('nimbleSubsamp' = quote({
##     customSpec <- configureMCMC(Rmodel)
##     customSpec$removeSamplers('z')
##     customSpec$addSampler('z', type = 'sampler_latentSub',
##                           control = list(leaveOutProportion = 0.6,
##                                          control = list()))
##     customSpec
## }))

## ss.ms.subsamp <- compareMCMCs(input1,
##                               MCMCs=c('nimbleSubsamp'),
##                               MCMCdefs = MCMCdefs.subsamp,
##                               niter= niter,
##                               burnin = burnin,
##                               summary=FALSE,
##                               check=FALSE)

## save(ss.ms.subsamp, file=file.path(save.dir, 'subsamp.Rdata'))


