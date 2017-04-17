rm(list=ls())
setwd("~/Dropbox/nimble/occupancy/analysis/spatial")
source('src/initialize.R')

set.seed(444)
dats <- genSpatialOccData()
model.input <- prepModData(dats$data, dats$y, dats$distance,
                           nsite=250, inits=dats$inits)

sp.mod <- nimbleCode({
  ## priors
  ## delta ~ dunif(0.1, 10)
  logSigma ~ dnorm(0, 1)
  sigma <- exp(logSigma)
  logDelta ~ dnorm(0, 1)
  delta <- exp(logDelta)
  ## sigma ~ dunif(0.1, 100)
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

  ## create covariance matrix based on distances
  prep.cov[1:nsite, 1:nsite] <- exp(-delta*D[1:nsite, 1:nsite])

  D.cov[1:nsite, 1:nsite] <- (sigma^2)*
    (0.95*prep.cov[1:nsite, 1:nsite] + 0.05*DI[1:nsite, 1:nsite])

})

input1 <- c(code=sp.mod,
            model.input)

## *********************************************************************
## opt 1:vanilla nimble and auto block
## *********************************************************************

sp.opt1 <- compareMCMCs(input1,
                        MCMCs=c("nimble"),
                        niter=niter,
                        burnin = burnin,
                        summary=FALSE,
                        check=FALSE)

save(sp.opt1, file=file.path(save.dir, "opt1.Rdata"))

## *********************************************************************
## opt 2: add custom z sampler and slice on uniform(0,1) nodes
## *********************************************************************


MCMCdefs.opt2 <- list('nimbleOpt2' = quote({
  customSpec <- configureMCMC(Rmodel)
  ## slice sampler for 0,1 varaibles
  customSpec$removeSamplers('p', print=FALSE)
  customSpec$addSampler(target = 'p',
                        type =
                        "slice",
                        print=FALSE)
  ## multivariate slice sampler
  customSpec$removeSamplers('rho', print=FALSE)
  customSpec$addSampler(target = 'rho',
                                                 type =
                                                   "AFSS_to_RM_block",
   control = list(AF_sliceControl =  list(sliceWidths = rep(1, length(pumpParams)), ## set initial slice widths
                                                       factorBurnIn = 10000,  ## number of iterations until factors stop adapting
                                                       factorAdaptInterval = 1000,  ## factors will be adapted every 1000 iterations
                                                       sliceBurnIn = 1000,
                                                       sliceMaxSteps = 100),
                               RWcontrol = list(propCov = diag(length(pumpParams)), scale = 1,
                                                adaptInterval = 500, adaptScaleOnly = F,
                                                adaptive = T),
                               nAFSSIters = 8000,
                               essThreshold = .5,
                               numESSAdaptations = 50,
                               timeSwitch = TRUE))
  customSpec
}))

## ## *********************************************************************
## ## run with compareMCMCs

sp.opt2 <- compareMCMCs(input1,
                        MCMCs=c('nimbleOpt2'),
                        MCMCdefs = MCMCdefs.opt2,
                        niter= niter,
                        burnin = burnin,
                        summary=FALSE,
                        check=FALSE)

save(sp.opt2, file=file.path(save.dir, "opt2.Rdata"))
