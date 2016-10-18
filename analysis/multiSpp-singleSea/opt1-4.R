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
## multi-species site-occupancy models: vectorized with custom
## function to remove zs
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
  sigma.ufcw ~ dunif(0, 100)
  
  mu.a1 ~ dnorm(0, 0.001)
  sigma.a1 ~ dunif(0, 100)
  
  mu.a2 ~ dnorm(0, 0.001)
  sigma.a2 ~ dunif(0, 100)
  
  mu.a3 ~ dnorm(0, 0.001)
  sigma.a3 ~ dunif(0, 100)
  
  mu.a4 ~ dnorm(0, 0.001)
  sigma.a4 ~ dunif(0, 100)

  sigma.vcato ~ dunif(0, 100)
  sigma.vfcw ~ dunif(0, 100)
  
  mu.b1 ~ dnorm(0, 0.001)
  sigma.b1 ~ dunif(0, 100)
  mu.b2 ~ dnorm(0, 0.001)
  sigma.b2 ~ dunif(0, 100)


  for (i in 1:(num.species)) {
    ## Create priors for species i from the community level prior
    ## distributions

    u.cato[i] ~ dnorm(mu.ucato, sd=sigma.ucato)
    u.fcw[i] ~ dnorm(mu.ufcw, sd=sigma.ufcw)
    a1[i] ~ dnorm(mu.a1, sd=sigma.a1)
    a2[i] ~ dnorm(mu.a2, sd=sigma.a2)
    a3[i] ~ dnorm(mu.a3, sd=sigma.a3)
    a4[i] ~ dnorm(mu.a4, sd=sigma.a4)

    v.cato[i] ~ dnorm(mu.vcato, sd=sigma.vcato)
    v.fcw[i] ~ dnorm(mu.vfcw, sd=sigma.vfcw)
    b1[i] ~ dnorm(mu.b1, sd=sigma.b1)
    b2[i] ~ dnorm(mu.b2, sd=sigma.b2)

    ## vectorize the calculation of psi.
    logit(psi[1:num.points,i]) <-
      u.cato[i]*(1-habitat.ind[1:num.points]) +
        u.fcw[i]*habitat.ind[1:num.points] +
          a1[i]*ufc.linear[1:num.points] +
            a2[i]*ufc.quadratic[1:num.points] +
              a3[i]*ba.linear[1:num.points] +
                a4[i]*ba.quadratic[1:num.points]
    ## vectorized calculation
    mu.psi[1:num.points,i] <- psi[1:num.points, i]

    ## For our purpose a better way to write this way is to
    ## not worry that some elements of date.linear and date.quadratic
    ## aren't used, since the benefit of vectorizing the computation
    ## should be much greater than the cost of a few extra elements
    logit(p[1:num.points, 1:max.num.reps, i]) <-
      (v.cato[i]*(1-habitat.ind[1:num.points]) +
       v.fcw[i]*habitat.ind[1:num.points]) %*%
         asRow(onesRow[1, 1:max.num.reps])+
           b1[i]*date.linear[1:num.points,1:max.num.reps] +
             b2[i]*date.quadratic[1:num.points,1:max.num.reps]

    ## user defined distribution to combine the bernoulli occupancy
    ## and detection events.  We can also make this is a single
    ## compuation for the entire matrix of locations-x-visits, for
    ## each species (i)
    
    X[1:num.points, 1:max.num.reps, i] ~ dBernDetectionMatrix(
      occProb = mu.psi[1:num.points,i],
      detectionProb = p[1:num.points, 1:max.num.reps,i],
      numReps = num.reps[1:num.points])
  }
  ## Derived quantities:
  ## for(j in 1:num.points){
  ##   N.site[j]<- sum(mu.psi[j,1:(num.species)])
  ##   N.ground[j]<- sum(mu.psi[j,1:num.species] * ground[1:num.species])
  ##   N.mid[j]<- sum(mu.psi[j,1:num.species] * mid[1:num.species])
  ## }
})


input1 <- c(code=ms.ss.occ, model.input)


## *********************************************************************
## option 1 (vanilla NIMBLE) using compare MCMC
## *********************************************************************

ms.ss.opt1 <- compareMCMCs(input1,
                           MCMCs=c('nimble', 'autoBlock'),
                           niter=niter,
                           burnin = burnin,
                           summary=FALSE,
                           check=FALSE)

save(ms.ss.opt1, file=file.path(save.dir, "opt1.Rdata"))


## *********************************************************************
## option 2: block sampler for species random effects for each species
## *********************************************************************

## remove the samples, add block samplers
MCMCdefs.opt2 <- list('nimbleOpt2' = quote({
  customSpec <- configureMCMC(Rmodel)
  ## find node names of each species for random effects
  base.names <- c("a1", "a2", "a3", "a4", "b1", "b2", "u.cato",
                  "u.fcw", "v.cato", "v.fcw" )
  exp.names.list <- list()
  for(bn in base.names){
    exp.names.list[[bn]] <- Rmodel$expandNodeNames(bn)
  }
  for(i in 1:length(exp.names.list[[1]])){
    blocknames <- unlist(lapply(exp.names.list, function(x) x[i]))
    customSpec$removeSamplers(blocknames, print=FALSE)
    customSpec$addSampler(target = blocknames,
                          type = "RW_block")
  }
  customSpec
}))

## run the model
ms.ss.opt2 <- compareMCMCs(input1,
                           MCMCs=c('nimbleOpt2'),
                           MCMCdefs = MCMCdefs.opt2,
                           niter= niter,
                           burnin = burnin,
                           summary=FALSE,
                           check=FALSE)

save(ms.ss.opt2, file=file.path(save.dir, "opt2.Rdata"))

## *********************************************************************
## ## option 3: block sampler for species random effects for each
## "type"
## *********************************************************************

## remove the samples, add block samplers
MCMCdefs.opt3 <- list('nimbleOpt3' = quote({
  customSpec <- configureMCMC(Rmodel)
  ## find node names for random effects
  sp.parms.a <- Rmodel$getNodeNames(includeData = FALSE)[grepl("^a",
                                      Rmodel$getNodeNames(includeData = FALSE))]
  sp.parms.b <- Rmodel$getNodeNames(includeData = FALSE)[grepl("^b",
                                      Rmodel$getNodeNames(includeData = FALSE))]
  sp.parms.u <- Rmodel$getNodeNames(includeData = FALSE)[grepl("^u",
                                      Rmodel$getNodeNames(includeData = FALSE))]
  sp.parms.v <- Rmodel$getNodeNames(includeData = FALSE)[grepl("^v",
                                      Rmodel$getNodeNames(includeData = FALSE))]
  customSpec$removeSamplers(c(sp.parms.a, sp.parms.b, sp.parms.u,
                              sp.parms.v),
                            print=FALSE)

  customSpec$addSampler(target = sp.parms.a,
                        type = "RW_block")
  customSpec$addSampler(target = sp.parms.b,
                        type = "RW_block")
  customSpec$addSampler(target = sp.parms.u,
                        type = "RW_block")
  customSpec$addSampler(target = sp.parms.v,
                        type = "RW_block")
  customSpec
}))

## run the model
ms.ss.opt3 <- compareMCMCs(input1,
                           MCMCs=c('nimbleOpt3'),
                           MCMCdefs = MCMCdefs.opt3,
                           niter= niter,
                           burnin = burnin,
                           summary=FALSE,
                           check=FALSE)

save(ms.ss.opt3, file=file.path(save.dir, "opt3.Rdata"))

## *********************************************************************
## ## option 4: sigma sampler
## *********************************************************************

## remove the samples, add block samplers
MCMCdefs.opt4 <- list('nimbleOpt4' = quote({
  customSpec <- configureMCMC(Rmodel)
  ## find node names for random effects
  sigma.parms <- Rmodel$getNodeNames(includeData = FALSE)[grepl("^sigma",
                                       Rmodel$getNodeNames(includeData = FALSE))]
  mu.parms <- paste("mu",
                    sapply(strsplit(sigma.parms, "[.]"), function(x) x[2]),
                    sep=".")
  customSpec$removeSamplers(sigma.parms,
                            print=FALSE)
  
  for(i in 1:length(sigma.parms)){
    customSpec$addSampler(target = sigma.parms[i],
                          type = "sampler_RW_log_shift",
                          control = list(shiftNodes = mu.parms[i]))
    
  }
  customSpec
}))

## run the model
ms.ss.opt4 <- compareMCMCs(input1,
                           MCMCs=c('nimbleOpt4'),
                           MCMCdefs = MCMCdefs.opt4,
                           niter= niter,
                           burnin = burnin,
                           summary=FALSE,
                           check=FALSE)

save(ms.ss.opt4, file=file.path(save.dir, "opt4.Rdata"))


## *********************************************************************
## ## CPPP
## *********************************************************************

occ.R.model <- nimbleModel(code=ms.ss.occ,
                           constants=input1$constants,
                           data=input1$data,
                           inits=input1$inits,
                           check=FALSE)

occ.mcmc <- buildMCMC(occ.R.model)
occ.C.model <- compileNimble(occ.R.model)
occ.C.mcmc <- compileNimble(occ.mcmc, project = occ.R.model)
occ.C.mcmc$run(niter)

source('../cppp/src/calcCPPP.R', chdir = TRUE)
options(mc.cores=8)

test.opt2 <- generateCPPP(occ.R.model,
                          occ.C.model,
                          occ.C.mcmc,
                          occ.mcmc,
                          dataName = 'X',
                          paramNames = input1$monitors, 
                          MCMCIter = niter, 
                          NSamp = 10^3,
                          NPDist = 10^3,
                          burnInProportion = 0.10,
                          thin = 1,
                          averageParams = TRUE,
                          discFuncGenerator=likeDiscFuncGenerator,
                          returnChains=TRUE)

save(test.opt2, file=file.path(save.dir, "ms_ss_noz_CPPP.Rdata"))


## *********************************************************************
## ## cross validation
## *********************************************************************
## if not already run above
occ.R.model <- nimbleModel(code=ms.ss.occ,
                           constants=input1$constants,
                           data=input1$data,
                           inits=input1$inits,
                           check=FALSE)

options(mc.cores=1)

source('../crossValidation/crossValidationFunction.R')

output <- crossValidateOne(model=occ.R.model,
                           dataNames= "X",
                           MCMCIter= niter,
                           burnIn=niter*0.1,
                           thin=1,
                           leaveOutIndex=3,
                           MCMCdefs=NULL)

## *********************************************************************

ms.ss.occ.simp <- nimbleCode({

  ## random effects  
  a1 ~ dnorm(0, 0.001)  
  a2 ~ dnorm(0, 0.001)
  a3 ~ dnorm(0, 0.001)
  a4 ~ dnorm(0, 0.001)
  b1 ~ dnorm(0, 0.001)
  b2 ~ dnorm(0, 0.001)


  for (i in 1:(num.species)) {
    ## Create priors for species i from the community level prior
    ## distributions

    ## vectorize the calculation of psi.
    logit(psi[1:num.points]) <-
      u.cato*(1-habitat.ind[1:num.points]) +
        u.fcw*habitat.ind[1:num.points] +
          a1*ufc.linear[1:num.points] +
            a2*ufc.quadratic[1:num.points] +
              a3*ba.linear[1:num.points] +
                a4*ba.quadratic[1:num.points]
    ## vectorized calculation
    mu.psi[1:num.points] <- psi[1:num.points]

    ## For our purpose a better way to write this way is to
    ## not worry that some elements of date.linear and date.quadratic
    ## aren't used, since the benefit of vectorizing the computation
    ## should be much greater than the cost of a few extra elements
    logit(p[1:num.points, 1:max.num.reps]) <-
      (v.cato*(1-habitat.ind[1:num.points]) +
       v.fcw*habitat.ind[1:num.points]) %*%
         asRow(onesRow[1, 1:max.num.reps])+
           b1*date.linear[1:num.points,1:max.num.reps] +
             b2*date.quadratic[1:num.points,1:max.num.reps]

    ## user defined distribution to combine the bernoulli occupancy
    ## and detection events.  We can also make this is a single
    ## compuation for the entire matrix of locations-x-visits, for
    ## each species (i)
    
    X[1:num.points, 1:max.num.reps, i] ~ dBernDetectionMatrix(
      occProb = mu.psi[1:num.points],
      detectionProb = p[1:num.points, 1:max.num.reps],
      numReps = num.reps[1:num.points])
  }
  ## Derived quantities:
  ## for(j in 1:num.points){
  ##   N.site[j]<- sum(mu.psi[j,1:(num.species)])
  ##   N.ground[j]<- sum(mu.psi[j,1:num.species] * ground[1:num.species])
  ##   N.mid[j]<- sum(mu.psi[j,1:num.species] * mid[1:num.species])
  ## }
})



occ.R.model.simp <- nimbleModel(code=ms.ss.occ.simp,
                           constants=input1$constants,
                           data=input1$data,
                           inits=input1$inits,
                           check=FALSE)

occ.mcmc <- buildMCMC(occ.R.model.simp)
occ.C.model.simp <- compileNimble(occ.R.model.simp)
occ.C.mcmc <- compileNimble(occ.mcmc, project = occ.R.model.simp)


output.simp <- crossValidateOne(model=occ.R.model.simp,
                           dataNames= "X",
                           MCMCIter= 10^4,
                           burnIn=10^2,
                           thin=2,
                           leaveOutIndex=3,
                           MCMCdefs=NULL)
