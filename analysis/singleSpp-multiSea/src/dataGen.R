logit <- function(x) {
  log(x/(1 - x))
}

antilogit <- function(x) {
  exp(x)/(1 + exp(x))
}

##  Generation and analysis of simulated data for multi season
##  occupancy model
genDynamicOccData <- function(R = 250,
                              J = 3, K = 10,
                              psi1 = 0.4,
                              range.p = c(0.2, 0.4),
                              range.phi = c(0.6, 0.8),
                              range.gamma = c(0, 0.1)) {
  ## Function to simulate detection/nondetection data for dynamic site-occ model
  ## Annual variation in probabilities of patch survival, colonization and
  ## detection is specified by the bounds of a uniform distribution.

  ## Function arguments:
  ## R - Number of sites
  ## J - Number of replicate surveys
  ## K - Number of years
  ## psi1 - occupancy probability in first year
  ## range.p - bounds of uniform distribution from which annual p drawn
  ## range.psi and range.gamma - same for survival and colonization probability

  ## Set up some required arrays
  site <- 1:R					## Sites
  year <- 1:K					## Years
  psi <- rep(NA, K)				## Occupancy probability
  muZ <- z <- array(dim = c(R, K))	## Expected and realized occurrence
  y <- array(NA, dim = c(R, J, K))	## Detection histories

  ## Determine initial occupancy and demographic parameters
  psi[1] <- psi1				## Initial occupancy probability
  p <- runif(n = K, min = range.p[1], max = range.p[2])
  phi <- runif(n = K-1, min = range.phi[1], max = range.phi[2])
  gamma <- runif(n = K-1, min = range.gamma[1], max = range.gamma[2])

  ## Generate latent states of occurrence
  ## First year
  z[,1] <- rbinom(R, 1, psi[1])		## Initial occupancy state
  ## Later years
  for(i in 1:R){				## Loop over sites
    for(k in 2:K){				## Loop over years
      muZ[k] <- z[i, k-1]*phi[k-1] + (1-z[i, k-1])*gamma[k-1] ## Prob for occ.
      z[i,k] <- rbinom(1, 1, muZ[k])
    }
  }

  ## Generate detection/nondetection data
  for(i in 1:R){
    for(k in 1:K){
      prob <- z[i,k] * p[k]
      for(j in 1:J){
        y[i,j,k] <- rbinom(1, 1, prob)
      }
    }
  }

  ## Compute annual population occupancy
  for (k in 2:K){
    psi[k] <- psi[k-1]*phi[k-1] + (1-psi[k-1])*gamma[k-1]
  }

  ## Plot apparent occupancy
  psi.app <- apply(apply(y, c(1,3), max), 2, mean)

  return(list(R = R, J = J, K = K,
              psi = psi, psi.app = psi.app, z = z,
              phi = phi, gamma = gamma, p = p, y = y))
}

## http://www.r-bloggers.com/multi-species-dynamic-occupancy-model-with-r-and-jags/

genMultiSpecOccData <- function(p_beta = 0.7,
                                sdbeta = 2,
                                p_rho = 0.8,
                                sdrho = 1,
                                nsite = 150,
                                nspec = 6,
                                nyear = 1,
                                nrep = 5,
                                p_p = 0.7,
                                sdp = 1.5) {
  ## community level hyperparameters
  mubeta <- logit(p_beta)
  murho <- logit(p_rho)

  ## species specific random effects
  beta <- rnorm(nspec, mubeta, sdbeta)
  rho <- rnorm(nspec, murho, sdrho)

  mup <- logit(p_p)
  lp <- rnorm(nspec, mup, sdp)
  p <- antilogit(lp)

  ## initial occupancy states
  rho0 <- runif(nspec, 0, 1)
  z0 <- array(dim = c(nsite, nspec))
  for (i in 1:nspec) {
    z0[, i] <- rbinom(nsite, 1, rho0[i])
  }

  ## subsequent occupancy
  z <- array(dim = c(nsite, nspec, nyear))
  lpsi <- array(dim = c(nsite, nspec, nyear))
  psi <- array(dim = c(nsite, nspec, nyear))
  for (j in 1:nsite) {
    for (i in 1:nspec) {
      for (t in 1:nyear) {
        if (t == 1) {
          lpsi[j, i, t] <- beta[i] + rho[i] * z0[j, i]
          psi[j, i, t] <- antilogit(lpsi[j, i, t])
          z[j, i, t] <- rbinom(1, 1, psi[j, i, t])
        } else {
          lpsi[j, i, t] <- beta[i] + rho[i] * z[j, i, t - 1]
          psi[j, i, t] <- antilogit(lpsi[j, i, t])
          z[j, i, t] <- rbinom(1, 1, psi[j, i, t])
        }
      }
    }
  }
  x <- array(dim = c(nsite, nspec, nyear, nrep))
  for (j in 1:nsite) {
    for (i in 1:nspec) {
      for (t in 1:nyear) {
        for (k in 1:nrep) {
          x[j, i, t, k] <- rbinom(1, 1, p[i] * z[j, i, t])
        }
      }
    }
  }
  return(x)
}
