rm(list=ls())
setwd('~/Dropbox/nimble/occupancy/analysis/multiSpp-singleSea')

source('src/initialize.R')
## don't agument data
n.zeroes <- 0
model.input <- prepMutiSpData(survey.data,
                              survey.dates,
                              species.groups,
                              habitat,
                              n.zeros)


## *********************************************************************
## multi-species site-occupancy models: original
## *********************************************************************

ms.ss.occ <- nimbleCode({
  ## Define prior distributions for community-level model parameters
  cato.occ.mean ~ dunif(0,1)
  mu.ucato <- log(cato.occ.mean) - log(1-cato.occ.mean)

  fcw.occ.mean ~ dunif(0,1)
  mu.ufcw <- log(fcw.occ.mean) - log(1-fcw.occ.mean)

  cato.det.mean ~ dunif(0,1)
  mu.vcato <- log(cato.det.mean) - log(1-cato.det.mean)

  fcw.det.mean ~ dunif(0,1)
  mu.vfcw <- log(fcw.det.mean) - log(1-fcw.det.mean)

  ## random effects
  sigma.ucato ~ dunif(0, 100)
  tau.ucato <-  1/(sigma.ucato*sigma.ucato)

  sigma.ufcw ~ dunif(0, 100)
  tau.ufcw <-  1/(sigma.ufcw*sigma.ufcw)

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
      ## Occurence model: u.cato and u.fcw are the occurrence
      ## probabilities (on the logit scale) for species i at points in
      ## the CATO and FCW study area, respectively, for average values
      ## of UFC and BA. The coefficients for the four 'a' terms are
      ## the linear and squared effects of understory foliage and tree
      ## basal area on species i.
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
        ## Detection model for the observed data X: v.cato and v.fcw
        ## are the detection probabilities (on the logit scale) for
        ## species i at points j during sampling periods k in the CATO
        ## and FCW study area, respectively, for for linear and
        ## squared terms of Julian dates.
        logit(p[j,k,i]) <-  v.cato[i]*(1-habitat.ind[j]) +
          v.fcw[i]*habitat.ind[j] +
          b1[i]*date.linear[j,k] +
          b2[i]*date.quadratic[j,k]

        mu.p[j,k,i] <- p[j,k,i]*Z[j,i]
        X[j,k,i] ~ dbern(mu.p[j,k,i])
      }
    }
  }

  ## Derived quantities:
  ## Create a loop to determine point level
  ## richness estimates for the whole community
  ## and for subsets or assemblages of interest.
  ## for(j in 1:num.points){
  ##   N.site[j]<- sum(mu.psi[j,1:(num.species)])
  ##   N.ground[j]<- inprod(Z[j,1:num.species],ground[1:num.species])
  ##   N.mid[j]<- inprod(Z[j,1:num.species],mid[1:num.species])
  ## }
})

 
model.input$data[["onesRow"]] <- NULL
model.input$constants[["max.num.reps"]] <- NULL

input1 <- c(code=ms.ss.occ, model.input)


occ.R.model <- nimbleModel(code=ms.ss.occ,
                           constants=input1$constants,
                           data=input1$data,
                           inits=input1$inits,
                           check=FALSE)



## *********************************************************************
## original model with nimble

ms.ss.orig <- compareMCMCs(input1,
                           MCMCs=c('jags'),
                           niter=niter,
                           burnin = burnin,
                           summary=FALSE,
                           check=FALSE)

save(ms.ss.orig, file=file.path(save.dir, "orig.Rdata"))
