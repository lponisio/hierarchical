rm(list=ls())
setwd('~/Dropbox/nimble/occupancy/analysis/multiSpp-singleSea')

source('src/initialize.R')
## don't agument data
n.zeroes <- 0
model.input <- prepMutiSpData(survey.data,
                              survey.dates,
                              species.groups,
                              habitat,
                              n.zeros,
                              remove.zs=FALSE)

load('~/Dropbox/nimble/saved/multiSpp-singleSea/saved/orig.Rdata')

model.input$inits <- c(model.input$inits, ms.ss.orig[[1]]$summary["nimble", "mean",])

## *********************************************************************
## multi-species site-occupancy models: original
## *********************************************************************

ms.ss.occ <- nimbleCode({
    ## Define prior distributions for community-level model parameters
    cato.occ.mean ~ dunif(0,1)
    mu.ucato <- log(cato.occ.mean) - log(1-cato.occ.mean)
    sigma.ucato ~ dunif(0, 100)
    tau.ucato <-  1/(sigma.ucato*sigma.ucato)

    fcw.occ.mean ~ dunif(0,1)
    mu.ufcw <- log(fcw.occ.mean) - log(1-fcw.occ.mean)
    sigma.ufcw ~ dunif(0, 100)
    tau.ufcw <-  1/(sigma.ufcw*sigma.ufcw)

    cato.det.mean ~ dunif(0,1)
    mu.vcato <- log(cato.det.mean) - log(1-cato.det.mean)

    fcw.det.mean ~ dunif(0,1)
    mu.vfcw <- log(fcw.det.mean) - log(1-fcw.det.mean)

    ## random effects

    mu.a1 ~ dnorm(0, 0.001)
    sigma.a1 ~ dunif(0, 100)
    tau.a1 <-  1/(sigma.a1*sigma.a1)

    mu.a2 ~ dnorm(0, 0.001)
    sigma.a2 ~ dunif(0, 100)
    tau.a2 <-  1/(sigma.a2*sigma.a2)

    mu.a3 ~ dnorm(0, 0.001)
    sigma.a3 ~ dunif(0, 100)
    tau.a3 <-  1/(sigma.a3*sigma.a3)

    mu.a4 ~ dnorm(0, 0.001)
    sigma.a4 ~ dunif(0, 100)
    tau.a4 <-  1/(sigma.a4*sigma.a4)

    sigma.vcato ~ dunif(0, 100)
    sigma.vfcw ~ dunif(0, 100)
    tau.vcato <-  1/(sigma.vcato*sigma.vcato)
    tau.vfcw <-  1/(sigma.vfcw*sigma.vfcw)

    mu.b1 ~ dnorm(0, 0.001)
    sigma.b1 ~ dunif(0, 100)
    tau.b1 <-  1/(sigma.b1*sigma.b1)

    mu.b2 ~ dnorm(0, 0.001)
    sigma.b2 ~ dunif(0, 100)
    tau.b2 <-  1/(sigma.b2*sigma.b2)

    for (i in 1:(num.species)) {
        ## Create priors for species i from the community level prior
        ## distributions

        u.cato[i] ~ dnorm(mu.ucato, tau.ucato)
        u.fcw[i] ~ dnorm(mu.ufcw, tau.ufcw)
        a1[i] ~ dnorm(mu.a1, tau.a1)
        a2[i] ~ dnorm(mu.a2, tau.a2)
        a3[i] ~ dnorm(mu.a3, tau.a3)
        a4[i] ~ dnorm(mu.a4, tau.a4)
        v.cato[i] ~ dnorm(mu.vcato, tau.vcato)
        v.fcw[i] ~ dnorm(mu.vfcw, tau.vfcw)
        b1[i] ~ dnorm(mu.b1, tau.b1)
        b2[i] ~ dnorm(mu.b2, tau.b2)

        ## Create a loop to estimate the Z matrix (true occurrence for
        ## species i at point j).
        for (j in 1:num.points) {
            logit(psi[j,i]) <- u.cato[i]*(1-habitat.ind[j]) +
                u.fcw[i]*habitat.ind[j] +
                a1[i]*ufc.linear[j] +
                a2[i]*ufc.quadratic[j] +
                a3[i]*ba.linear[j] +
                a4[i]*ba.quadratic[j]
            mu.psi[j,i] <- psi[j,i]
            Z[j,i] ~ dbern(mu.psi[j,i])

            ## Create a loop to estimate detection for species i at point k
            ## during sampling period k.
            for (k in 1:num.reps[j]) {
                logit(p[j,k,i]) <-  v.cato[i]*(1-habitat.ind[j]) +
                    v.fcw[i]*habitat.ind[j] +
                    b1[i]*date.linear[j,k] +
                    b2[i]*date.quadratic[j,k]
                mu.p[j,k,i] <- p[j,k,i]*Z[j,i]
                X[j,k,i] ~ dbern(mu.p[j,k,i])
            }
        }
    }
})



input1 <- c(code=ms.ss.occ, model.input)

## *********************************************************************
## original model with nimble

ms.ss.orig <- compareMCMCs(input1,
                           MCMCs=c('jags', 'nimble'),
                           niter=niter,
                           burnin = burnin,
                           summary=FALSE,
                           check=FALSE)

save(ms.ss.orig, file=file.path(save.dir, "orig.Rdata"))

## *********************************************************************
## sampler only a subset of latent states
## *********************************************************************

MCMCdefs.subsamp <- list('nimbleSubsamp' = quote({
    customSpec <- configureMCMC(Rmodel)
    customSpec$removeSamplers('Z')
    customSpec$addSampler('Z', type = 'sampler_latentSub',
                          control = list(leaveOutProportion = 0.05,
                                         control = list()))
    customSpec
}))

ms.ss.subsamp <- compareMCMCs(input1,
                              MCMCs=c('nimbleSubsamp'),
                              MCMCdefs = MCMCdefs.subsamp,
                              niter= niter,
                              burnin = burnin,
                              summary=FALSE,
                              check=FALSE)

save(ms.ss.subsamp, file=file.path(save.dir, 'subsamp.Rdata'))


## *********************************************************************
## cross level sampler
## *********************************************************************
## sample latent and top-level parameters jointly


MCMCdefs.crosslevel <- list('nimbleCrosslevel' = quote({
    customSpec <- configureMCMC(Rmodel)
    customSpec$removeSamplers('Z')
    ## find node names of each species for random effects
    base.names <- c("a1", "a2", "a3", "a4", "b1", "b2", "u.cato",
                    "u.fcw", "v.cato", "v.fcw" )

    customSpec$removeSamplers(base.names)
    customSpec$addSampler(target = base.names,
                          type ='sampler_crossLevelBinary',
                          print=FALSE)
    customSpec
}))

ms.ss.crosslevel <- compareMCMCs(input1,
                                 MCMCs=c('nimbleCrosslevel'),
                                 MCMCdefs = MCMCdefs.crosslevel,
                                 niter= niter,
                                 burnin = burnin,
                                 summary=FALSE,
                                 check=FALSE)

save(ms.ss.crosslevel, file=file.path(save.dir, 'crosslevel.Rdata'))




