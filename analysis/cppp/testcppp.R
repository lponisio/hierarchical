rm(list=ls())
library(nimble)
library(parallel)
options(mc.cores=2)
nthin <- 2
source('~/Dropbox/occupancy-nimble/cppp/src/calcCPPP.R')


pumpCode <- nimbleCode({
  for (i in 1:N){
    theta[i] ~ dgamma(alpha,beta)
    lambda[i] <- theta[i]*t[i]
    x[i] ~ dpois(lambda[i])
  }
  alpha ~ dexp(1.0)
  beta ~ dgamma(0.1,1.0)
})

pumpConsts <- list(N = 10,
                   t = c(94.3, 15.7, 62.9, 126, 5.24,
                     31.4, 1.05, 1.05, 2.1, 10.5))
pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22))
pumpInits <- list(alpha = 1, beta = 1,
                  theta = rep(0.1, pumpConsts$N))



## build model
R.model <- nimbleModel(code=pumpCode,
                       constants=pumpConsts,
                       data=pumpData,
                       inits=pumpInits,
                       check=FALSE)
message('R model created')

## configure and build mcmc
mcmc.spec <- configureMCMC(R.model,
                           print=FALSE,
                           monitors = c("alpha", "beta"),
                           thin=nthin)
mcmc <- buildMCMC(mcmc.spec)
message('MCMC built')

## compile model in C++
D.model <- compileNimble(R.model)
D.mcmc <- compileNimble(mcmc, project = R.model)
message('NIMBLE model compiled')
set.seed(4)
output <- generateCPPP(R.model,
             D.model,
             D.mcmc,
             mcmc,
             dataNames = 'x',
             paramNames = c('alpha','beta'), 
             MCMCIter = 3000, 
             NSamp = 3000,
             NPDist = 500,
             burnInProportion = .1,
             thin = nthin)
