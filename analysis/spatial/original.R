rm(list=ls())
setwd("occupancy/analysis/spatial")
source('src/initialize.R')
library('rjags')
load.module("msm")
load.module("glm")

set.seed(444)
dats <- genSpatialOccData(alpha=0.6, sigma=5, delta=0.1)
model.input <- prepModData(dats$data, dats$y, dats$distance,
                           nsite=250, inits=dats$inits)

mexp <- nimbleFunction(
    run = function(A = double(2)){
        returnType(double(2))
        outMat <- exp(A)
        return(outMat)
    })

sp.mod <- nimbleCode({
    ## priors

    ## logSigma ~ dnorm(0, 0.001)
    ## sigma <- exp(logSigma)
    ## logDelta ~ dnorm(0, 0.001)
    ## delta <- exp(logDelta)

    delta ~ dunif(0, 100)
    sigma ~ dunif(0, 100)


    p ~ dunif(0, 1)
    alpha ~ dnorm(0, 0.001)

    rho[1:nsite] ~ dmnorm(zeros[1:nsite],
                          D.tau[1:nsite, 1:nsite])


    ## Likelihood
    ## Ecological model for true occurrence
    for (site in 1:nsite) {
        z[site] ~ dbern(psi[site])
        logit(psi[site]) <- alpha + rho[site]
        p.eff[site] <- z[site] * p

        ## Observation model for replicated detection/nondetection
        ## observations
        for (rep in 1:nreps) {
            y[site, rep] ~ dbern(p.eff[site])
        }
    }

    ## create covariance matrix based on distances (must be 1/cov for
    ## JAGS)

    ## mexp is jags's version fo matrix exponentiation, very sensitive
    prep.cov[1:nsite, 1:nsite] <- 1/(exp(delta)*mexp(D[1:nsite, 1:nsite]))
    D.cov[1:nsite, 1:nsite] <- (0.95*prep.cov[1:nsite, 1:nsite] +
                                0.05*DI[1:nsite, 1:nsite])

    D.tau[1:nsite, 1:nsite] <- (1/(sigma))*inverse(D.cov[1:nsite, 1:nsite])
})


input1 <- c(code=sp.mod,
            model.input)

## *********************************************************************
## opt 1:vanilla nimble and auto block
## *********************************************************************

sp.orig <- compareMCMCs(input1,
                        MCMCs=c("jags"),
                        niter=niter,
                        burnin = burnin,
                        summary=FALSE,
                        check=FALSE)

save(sp.orig, file=file.path(save.dir, "orig.Rdata"))

