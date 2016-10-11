rm(list=ls())
setwd("~/Dropbox/nimble/occupancy/analysis/spatial")
source('src/initialize.R')
library(rjags)
load.module("msm")

set.seed(444)
dats <- genSpatialOccData()
model.input <- prepModData(dats$data, dats$y, dats$distance,
                           nsite=50)

mexp <- nimbleFunction(
  run = function(A = double(2)){
    returnType(double(2))
    outMat <- exp(A)
    return(outMat)
  })

sp.mod <- nimbleCode({
  ## priors
  delta ~ dunif(0.01, 10)
  sigma ~ dunif(0, 10)
  p ~ dunif(0, 1)
  alpha ~ dnorm(0, 0.001)
  b1 ~ dnorm(0, 0.001)

  rho[1:nsite] ~ dmnorm(zeros[1:nsite],
                        D.tau[1:nsite, 1:nsite])

  
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

  ## create covariance matrix based on distances (must be 1/cov for
  ## JAGS)

  ## mexp is jags's version fo matrix exponentiation, very sensitive 
  ## temp.cov[1:nsite, 1:nsite] <- -delta*D[1:nsite, 1:nsite]
  D.cov[1:nsite, 1:nsite]  <- (sigma^2)*mexp(-delta*D[1:nsite, 1:nsite])

  ## for(i in 1:nsite){
  ##   for(j in 1:nsite){
  ##     temp.cov[i, j] <- -delta*D[i, j]
  ##     D.cov[i, j]  <- (sigma^2)* exp(temp.cov[i, j])
  ##   }
  ## }
  
  D.tau[1:nsite, 1:nsite] <- inverse(D.cov[1:nsite, 1:nsite])
  
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
