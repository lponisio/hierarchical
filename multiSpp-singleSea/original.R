rm(list=ls())
setwd('~/Dropbox/occupancy-nimble/multiSpp-singleSea')

## don't augment data
n_zeroes <- 0
source('src/initialize.R')

## *********************************************************************
## multi-species site-occupancy models: original
## *********************************************************************

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
  tau_ucato <-  1/(sigma_ucato*sigma_ucato)

  sigma_ufcw ~ dunif(0, 100)
  tau_ufcw <-  1/(sigma_ufcw*sigma_ufcw)

  mu_a1 ~ dnorm(0, 0.001)
  sigma_a1 ~ dunif(0, 100)
  tau_a1 <-  1/(sigma_a1*sigma_a1)

  mu_a2 ~ dnorm(0, 0.001)
  sigma_a2 ~ dunif(0, 100)
  tau_a2 <-  1/(sigma_a2*sigma_a2)

  mu_a3 ~ dnorm(0, 0.001)
  sigma_a3 ~ dunif(0, 100)
  tau_a3 <-  1/(sigma_a3*sigma_a3)

  mu_a4 ~ dnorm(0, 0.001)
  sigma_a4 ~ dunif(0, 100)
  tau_a4 <-  1/(sigma_a4*sigma_a4)

  sigma_vcato ~ dunif(0, 100)
  sigma_vfcw ~ dunif(0, 100)
  tau_vcato <-  1/(sigma_vcato*sigma_vcato)
  tau_vfcw <-  1/(sigma_vfcw*sigma_vfcw)

  mu_b1 ~ dnorm(0, 0.001)
  sigma_b1 ~ dunif(0, 100)
  tau_b1 <-  1/(sigma_b1*sigma_b1)

  mu_b2 ~ dnorm(0, 0.001)
  sigma_b2 ~ dunif(0, 100)
  tau_b2 <-  1/(sigma_b2*sigma_b2)

  for (i in 1:(num_species)) {
    ## Create priors for species i from the community level prior
    ## distributions

    u_cato[i] ~ dnorm(mu_ucato, tau_ucato)
    u_fcw[i] ~ dnorm(mu_ufcw, tau_ufcw)
    a1[i] ~ dnorm(mu_a1, tau_a1)
    a2[i] ~ dnorm(mu_a2, tau_a2)
    a3[i] ~ dnorm(mu_a3, tau_a3)
    a4[i] ~ dnorm(mu_a4, tau_a4)

    v_cato[i] ~ dnorm(mu_vcato, tau_vcato)
    v_fcw[i] ~ dnorm(mu_vfcw, tau_vfcw)
    b1[i] ~ dnorm(mu_b1, tau_b1)
    b2[i] ~ dnorm(mu_b2, tau_b2)

    ## Create a loop to estimate the Z matrix (true occurrence for
    ## species i at point j).
    for (j in 1:num_points) {
      ## Occurence model: u_cato and u_fcw are the occurrence
      ## probabilities (on the logit scale) for species i at points in
      ## the CATO and FCW study area, respectively, for average values
      ## of UFC and BA. The coefficients for the four 'a' terms are
      ## the linear and squared effects of understory foliage and tree
      ## basal area on species i.
      logit(psi[j,i]) <- u_cato[i]*(1-habitat_ind[j]) +
        u_fcw[i]*habitat_ind[j] +
        a1[i]*ufc_linear[j] +
        a2[i]*ufc_quadratic[j] +
        a3[i]*ba_linear[j] +
        a4[i]*ba_quadratic[j]

      mu_psi[j,i] <- psi[j,i]
      Z[j,i] ~ dbern(mu_psi[j,i])

      ## Create a loop to estimate detection for species i at point k
      ## during sampling period k.
      for (k in 1:num_reps[j]) {
        ## Detection model for the observed data X: v_cato and v_fcw
        ## are the detection probabilities (on the logit scale) for
        ## species i at points j during sampling periods k in the CATO
        ## and FCW study area, respectively, for for linear and
        ## squared terms of Julian dates.
        logit(p[j,k,i]) <-  v_cato[i]*(1-habitat_ind[j]) +
          v_fcw[i]*habitat_ind[j] +
          b1[i]*date_linear[j,k] +
          b2[i]*date_quadratic[j,k]

        mu_p[j,k,i] <- p[j,k,i]*Z[j,i]
        X[j,k,i] ~ dbern(mu_p[j,k,i])
      }
    }
  }

  ## Derived quantities:
  ## Create a loop to determine point level
  ## richness estimates for the whole community
  ## and for subsets or assemblages of interest.
  ## for(j in 1:num_points){
  ##   N_site[j]<- sum(mu_psi[j,1:(num_species)])
  ##   N_ground[j]<- inprod(Z[j,1:num_species],ground[1:num_species])
  ##   N_mid[j]<- inprod(Z[j,1:num_species],mid[1:num_species])
  ## }
})


model_data[["onesRow"]] <- NULL
constants[["max_num_reps"]] <- NULL
model_data[["ground"]] <- NULL
model_data[["mid"]] <- NULL

## *********************************************************************
## original model with nimble

input1 <- list(code=ms.ss.occ,
               constants=constants,
               data=model_data,
               inits=inits)

ms.ss.orig <- compareMCMCs(input1,
                           MCMCs=c('jags','nimble'),
                           niter=niter,
                           burnin = burnin,
                           summary=FALSE,
                           check=FALSE)

save(ms.ss.orig, file="saved/orig.Rdata")
