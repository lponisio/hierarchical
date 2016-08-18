rm(list=ls())
setwd('~/Dropbox/nimble-dev/occupancy/multiSpp-singleSea')

## don't agument data
n_zeroes <- 0
source('src/initialize.R')

## *********************************************************************
## multi-species site-occupancy models: vectorized with custom
## function to remove zs
## *********************************************************************

## Since calculations with date_linear and date_quadratic are now
## vectorized, we'll set the NAs to 0
model_data$date_linear[is.na(model_data$date_linear)] <- 0
model_data$date_quadratic[is.na(model_data$date_quadratic)] <- 0

ms.ss.occ <- nimbleCode({
  ## Define prior distributions for community-level model parameters
  cato_occ_mean ~ dunif(0,1)
  mu_ucato <- log(cato_occ_mean) - log(1-cato_occ_mean)

  fcw_occ_mean ~ dunif(0,1)
  mu_ufcw <- log(fcw_occ_mean) - log(1-fcw_occ_mean)

  cato_det_mean ~ dunif(0,1)
  mu_vcato <- log(cato_det_mean) - log(1-cato_det_mean)

  fcw_det_mean ~ dunif(0,1)
  mu_vfcw <- log(fcw_det_mean) - log(1-fcw_det_mean)

  ## random effects
  sigma_ucato ~ dunif(0, 100)
  sigma_ufcw ~ dunif(0, 100)
  
  mu_a1 ~ dnorm(0, 0.001)
  sigma_a1 ~ dunif(0, 100)
  
  mu_a2 ~ dnorm(0, 0.001)
  sigma_a2 ~ dunif(0, 100)
  
  mu_a3 ~ dnorm(0, 0.001)
  sigma_a3 ~ dunif(0, 100)
  
  mu_a4 ~ dnorm(0, 0.001)
  sigma_a4 ~ dunif(0, 100)

  sigma_vcato ~ dunif(0, 100)
  sigma_vfcw ~ dunif(0, 100)
  
  mu_b1 ~ dnorm(0, 0.001)
  sigma_b1 ~ dunif(0, 100)
  mu_b2 ~ dnorm(0, 0.001)
  sigma_b2 ~ dunif(0, 100)


  for (i in 1:(num_species)) {
    ## Create priors for species i from the community level prior
    ## distributions

    u_cato[i] ~ dnorm(mu_ucato, sd=sigma_ucato)
    u_fcw[i] ~ dnorm(mu_ufcw, sd=sigma_ufcw)
    a1[i] ~ dnorm(mu_a1, sd=sigma_a1)
    a2[i] ~ dnorm(mu_a2, sd=sigma_a2)
    a3[i] ~ dnorm(mu_a3, sd=sigma_a3)
    a4[i] ~ dnorm(mu_a4, sd=sigma_a4)

    v_cato[i] ~ dnorm(mu_vcato, sd=sigma_vcato)
    v_fcw[i] ~ dnorm(mu_vfcw, sd=sigma_vfcw)
    b1[i] ~ dnorm(mu_b1, sd=sigma_b1)
    b2[i] ~ dnorm(mu_b2, sd=sigma_b2)

    ## Major change: We can vectorize the calculation of psi.  This
    ## means we don't need the for loop.  It reduces the number of
    ## nodes in the model by a factor of num_points
    logit(psi[1:num_points,i]) <-
      u_cato[i]*(1-habitat_ind[1:num_points]) +
      u_fcw[i]*habitat_ind[1:num_points] +
      a1[i]*ufc_linear[1:num_points] +
      a2[i]*ufc_quadratic[1:num_points] +
      a3[i]*ba_linear[1:num_points] +
      a4[i]*ba_quadratic[1:num_points]
    ## We can also vectorize this
    mu_psi[1:num_points,i] <- psi[1:num_points, i]

    ## and this: For our purpose a better way to write this way is to
    ## not worry that some elements of date_linear and date_quadratic
    ## aren't used, since the benefit of vectorizing the computation
    ## should be much greater than the cost of a few extra elements
    logit(p[1:num_points, 1:max_num_reps, i]) <-
      (v_cato[i]*(1-habitat_ind[1:num_points]) +
       v_fcw[i]*habitat_ind[1:num_points]) %*%
      asRow(onesRow[1, 1:max_num_reps])+
      b1[i]*date_linear[1:num_points,1:max_num_reps] +
      b2[i]*date_quadratic[1:num_points,1:max_num_reps]

    ## This is the biggest change: We can write our own distribution
    ## to combine the bernoulli occupancy and detection events.  We
    ## can also make this is a single compuation for the entire matrix
    ## of locations-x-visits, for each species (i) The code to define
    ## dBernDetectionMatrix is below

    X[1:num_points, 1:max_num_reps, i] ~ dBernDetectionMatrix(
      occProb = mu_psi[1:num_points,i],
      detectionProb = p[1:num_points, 1:max_num_reps,i],
      num_reps = num_reps[1:num_points])
  }
  ## Derived quantities:
  ## since we don't have Z's any more, I'm going to define these
  ## derived quantities differently, as expected values.  One could
  ## simply generate Z's at this predictive stage to capture that
  ## additional variation or figure out derived quantities as wanted
  ## for(j in 1:num_points){
  ##   N_site[j]<- sum(mu_psi[j,1:(num_species)])
  ##   N_ground[j]<- sum(mu_psi[j,1:num_species] * ground[1:num_species])
  ##   N_mid[j]<- sum(mu_psi[j,1:num_species] * mid[1:num_species])
  ## }
})


## zs are removed from these models
model_data[["Z"]] <- NULL
model_data[["ground"]] <- NULL
model_data[["mid"]] <- NULL
inits[["Z"]] <- NULL

input1 <- list(code=ms.ss.occ,
               constants=constants,
               data=model_data,
               inits=inits)


## *********************************************************************
## option 1 (vanilla NIMBLE) using compare MCMC

ms.ss.opt1 <- compareMCMCs(input1,
                           MCMCs=c('nimble', 'autoBlock'),
                           niter=niter,
                           burnin = burnin,
                           summary=FALSE,
                           check=FALSE)

save(ms.ss.opt1, file="saved/opt1.Rdata")


## *********************************************************************
## option 2: block sampler for species random effects for each "type"

## remove the samples, add block samplers
MCMCdefs.opt2 <- list('nimbleOpt2' = quote({
  customSpec <- configureMCMC(Rmodel)
  ## find node names for random effects
  sp.parms.a <- Rmodel$getNodeNames()[grepl("^a",
                                    Rmodel$getNodeNames(includeData = FALSE))]
  sp.parms.b <- Rmodel$getNodeNames()[grepl("^b",
                                     Rmodel$getNodeNames(includeData = FALSE))]
  sp.parms.u <- Rmodel$getNodeNames()[grepl("^u",
                                     Rmodel$getNodeNames(includeData = FALSE))]
  sp.parms.v <- Rmodel$getNodeNames()[grepl("^v",
                                     Rmodel$getNodeNames(includeData = FALSE))]
  customSpec$removeSamplers(c(sp.parms.a, sp.parms.b, sp.parms.u,
                              sp.parms.v),
                            print=FALSE)

  customSpec$addSampler(target = sp.parms.a,
                        type = "RW_block", log=TRUE)
  customSpec$addSampler(target = sp.parms.b,
                        type = "RW_block", log=TRUE)
  customSpec$addSampler(target = sp.parms.u,
                        type = "RW_block", log=TRUE)
  customSpec$addSampler(target = sp.parms.v,
                        type = "RW_block", log=TRUE)
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

save(ms.ss.opt2, file="saved/opt2.Rdata")

## *********************************************************************
## option 3: block sampler for species random effects for each species

## remove the samples, add block samplers
MCMCdefs.opt3 <- list('nimbleOpt3' = quote({
  customSpec <- configureMCMC(Rmodel)
  ## find node names of each species for random effects
  base.names <- c("a1", "a2", "a3", "a4", "b1", "b2", "u_cato",
                  "u_fcw", "v_cato", "v_fcw" )
  exp.names.list <- list()
  for(bn in base.names){
    exp.names.list[[bn]] <- Rmodel$expandNodeNames(bn)
  }
  for(i in 1:length(exp.names.list[[1]])){
    blocknames <- unlist(lapply(exp.names.list, function(x) x[i]))
    customSpec$removeSamplers(blocknames, print=FALSE)
    customSpec$addSampler(target = blocknames, type = "RW_block", log=TRUE)
  }
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

save(ms.ss.opt3, file="saved/opt3.Rdata")
