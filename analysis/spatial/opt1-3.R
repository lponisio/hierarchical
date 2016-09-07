rm(list=ls())
gctorture()
setwd("~/Dropbox/nimble/occupancy/analysis/spatial")
source('src/initialize.R')

sp.mod <- nimbleCode({
  ## priors
  delta ~ dunif(0, 1)
  sigma ~ dunif(0, 10)
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

## *********************************************************************
## opt 1:vanilla nimble and auto block
## *********************************************************************

sp.opt1 <- compareMCMCs(input1,
                            MCMCs=c("nimble", "autoBlock"),
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
  ## identify samplers to replace
  znodes <- Rmodel$expandNodeNames('z')
  znodes <- znodes[!Rmodel$isData(znodes)]
  ## remove samplers
  customSpec$removeSamplers(znodes, print=FALSE)
  ## add custom samples
  for(znode in znodes) customSpec$addSampler(target = znode,
                                             type = custom_z_sampler,
                                             print=FALSE)
  ## slice sampler for 0,1 varaibles
  customSpec$removeSamplers('p', print=FALSE)
  customSpec$addSampler(target = 'p',
                                                 type =
                                                   "slice",
                                                 print=FALSE)
  ## multivariate normal sampler
  ## customSpec$removeSamplers('rho', print=FALSE)
  ## customSpec$addSampler(target = 'rho',
  ##                                                type =
  ##                                                  "ess",
  ##                                                print=FALSE)
  customSpec
}))

## *********************************************************************
## run with compareMCMCs

sp.opt2 <- compareMCMCs(input1,
                           MCMCs=c('nimbleOpt2'),
                           MCMCdefs = MCMCdefs.opt2,
                           niter= niter,
                           burnin = burnin,
                           summary=FALSE,
                           check=FALSE)

save(sp.opt2, file=file.path(save.dir, "opt2.Rdata"))

## *********************************************************************
## opt3:  add custom z sampler and reflective samplers
## *********************************************************************

MCMCdefs.opt3 <- list('nimbleOpt3' = quote({
  customSpec <- configureMCMC(Rmodel)
  ## identify samplers to replace
  znodes <- Rmodel$expandNodeNames('z')
  znodes <- znodes[!Rmodel$isData(znodes)]
  ## remove samplers
  customSpec$removeSamplers(znodes, print=FALSE)
  ## add custom samples
  for(znode in znodes) customSpec$addSampler(target = znode,
                                             type = custom_z_sampler,
                                             print=FALSE)

  ## reflective for 0,1 varaibles
  customSpec$removeSamplers('p', print=FALSE)
  customSpec$addSampler(target = 'p',
                                                 type =
                                                   "sampler_RW_reflect",
                                                 print=FALSE)
  ## multivariate normal sampler
  ## customSpec$removeSamplers('rho', print=FALSE)
  ## customSpec$addSampler(target = 'rho',
  ##                                                type =
  ##                                                  "ess",
  ##                                                print=FALSE)
  customSpec
}))


## *********************************************************************
## run with compareMCMCs

sp.opt3 <- compareMCMCs(input1,
                           MCMCs=c('nimbleOpt3'),
                           MCMCdefs = MCMCdefs.opt3,
                           niter= niter,
                           burnin = burnin,
                           summary=FALSE,
                           check=FALSE)

save(sp.opt3, file=file.path(save.dir, "opt3.Rdata"))


## *********************************************************************
## opt:  automated factor slice sampler
## *********************************************************************
