rm(list=ls())
setwd("~/Dropbox/nimble/occupancy/analysis/spatial")
source('src/initialize.R')
library(rjags)
load.module("msm")

set.seed(444)
dats <- genSpatialOccData()
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

  logSigma ~ dnorm(0, 0.001)
  sigma <- exp(logSigma)
  logDelta ~ dnorm(0, 0.001)
  delta <- exp(logDelta)

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
  prep.cov[1:nsite, 1:nsite] <- 1/(exp(delta)*mexp(D[1:nsite, 1:nsite]))
  D.cov[1:nsite, 1:nsite] <- (sigma^2)*(0.95*prep.cov[1:nsite, 1:nsite] + 0.05*DI[1:nsite, 1:nsite])

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


## *********************************************************************
## model in jags
## *********************************************************************

library('R2jags')

sp.mod.jags <- function(dd,
                        ni, nt, nb, nc) {

    sink('model.jags')
    cat('model{

  ## priors
  delta ~ dunif(0.1, 10)
  sigma ~ dunif(0.1, 100)
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

    prep.cov[1:nsite, 1:nsite] <- 1/(exp(delta)*mexp(D[1:nsite, 1:nsite]))
  D.cov[1:nsite, 1:nsite] <- (sigma^2)*(0.95*prep.cov[1:nsite, 1:nsite] + 0.05*DI[1:nsite, 1:nsite])


  D.tau[1:nsite, 1:nsite] <- inverse(D.cov[1:nsite, 1:nsite])

  }',fill = TRUE)
    sink()
    jags(dd$data, dd$inits, dd$params,'model.jags',n.chains=nc,
         n.thin=nt,n.iter=ni,n.burnin=nb,working.directory=NULL)
}


analyse.jags <- function(d, ni, nt, nb, nc) {
    names(d)[names(d) == "monitors"] <-  "params"
    d$data <-  c(d$constants, d$data)
    d$constants <- NULL
    these.inits <- d$inits

    my.inits <- function() {
        these.inits
    }

    dd <- list(data=d$data, inits=my.inits, params=model.input$monitors)
    res <- list(data=d$data, bugs=sp.mod.jags(dd, ni=ni, nt=nt, nb=nb, nc=nc))
    summary <- list(data=d, bugs=res$bugs$BUGSoutput$summary)
    return(summary)
}


scale <- 1e3
save.dir <- '~/Dropbox/nimble/saved'


res <- analyse.jags(model.input,
                     ni=(1e3+1e1)*scale,
                     nt=scale,
                     nb=1e1*scale,
                     nc=3)


## ## cols <- c('mean', 'sd', '2.5%', '97.5%', 'Rhat', 'n.eff')
## ## res$bugs[,cols]
