rm(list=ls())
##setwd("~/Dropbox/occupancy-nimble/spatial")
source('src/initialize.R')

sp.mod <- nimbleCode({
  ## priors
  delta ~ dunif(0, 10)
  sigma ~ dunif(0, 10)
##  psi ~ dunif(0, 1)
  p ~ dunif(0, 1)
  alpha ~ dnorm(0, 0.001)
  b1 ~ dnorm(0, 0.001)

  ## Likelihood
  ## Ecological model for true occurrence
  for (i in 1:nsite) {
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- alpha + b1*elev[i] + rho[i]
    p.eff[i] <- z[i] * p

    ## Observation model for replicated detection/nondetection
    ## observations
    for (j in 1:nreps) {
      y[i,j] ~ dbern(p.eff[i])
    }
  }

  rho[1:nsite] ~ dmnorm(zeros[1:nsite],
                        cov = D.cov[1:nsite, 1:nsite]) 
  ## derived quantities
  ## turning the distance matrix to covariance matrix
  D.cov[1:nsite, 1:nsite] <- (sigma^2)*exp(-delta*D[1:nsite, 1:nsite])
  
})

input1 <- list(code=sp.mod,
               constants=constants,
               data=model.data,
               inits=inits)

sp.mod.1 <- nimbleModel(sp.mod, constants = constants, data = model.data, inits = inits, debug = TRUE)

sp.mod.opt1 <- compareMCMCs(input1,
                            MCMCs=c('nimble'),
                            niter=niter,
                            burnin = burnin,
                            summary=FALSE,
                            check=FALSE)

## save(sp.mod.opt1, file="saved/opt1.Rdata")

## checkChains(sp.mod.opt1[[1]]$samples,
##             f.path = "figures/chains/%s.pdf")
