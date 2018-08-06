## setwd("~/Dropbox/occupancy/")
rm(list=ls())
setwd("analysis/singleSpp-multiSea")
source('src/initialize.R')
set.seed(444)
data <- genDynamicOccData()
model.input <- prepModDataOcc(data, include.zs=FALSE)

## *********************************************************************
##  Multi-season occupancy model: remove latent states using
##  user-defined NIMBLE function dDynamicOcc
##  *********************************************************************

ss.ms.occ <- nimbleCode({
    ##  priors
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
    for(i in 1:nsite) {
        ## removes the z's and muZ's from the model and compute
        ## the probability of all reps over all years for one site.
        y[i, 1:nrep, 1:nyear] ~ dDynamicOccupancy(nrep,
                                                  psi1,
                                                  expit(phi[1:(nyear-1)]),
                                                  expit(gamma[1:(nyear-1)]),
                                                  expit(p[1:nyear]))
    }
})

input1 <- c(code=ss.ms.occ,
            model.input)

## *********************************************************************
## vanilla nimble and slice samplers
## *********************************************************************

ss.ms.filter <- compareMCMCs(input1,
                             MCMCs=c('nimble'),
                             niter=niter,
                             burnin = burnin,
                             summary=FALSE,
                             check=FALSE)

save(ss.ms.filter, file=file.path(save.dir, "filter.Rdata"))

## *********************************************************************
## block together phi and gamma
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

ss.ms.filer.blocking <- compareMCMCs(input1,
                                     MCMCs=c('AFSS_block'),
                                     MCMCdefs = MCMCdefs.AFSS.block,
                                     niter=niter,
                                     burnin = burnin,
                                     summary=FALSE,
                                     check=FALSE)

save(ss.ms.filter.blocking, file=file.path(save.dir, "filter_AFSS_block.Rdata"))
