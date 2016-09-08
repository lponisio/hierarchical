library(nimble)
library(igraph)

distributionsInputList <- nimble:::distributionsInputList
## OR
## library(devtools)
## install_github("lponisio/nimble",
##                ref = "userDistWorkArnd",
##                subdir = "packages/nimble")

## *********************************************************************
## user defined functions to remove latent states
## *********************************************************************

## DynamicOccupancy removes the z's and muZ's from the model and computes
## the probability of all reps over all years for one site.
dDynamicOccupancy <- nimbleFunction(
  ## I've checked that this runs and compiles, but I haven't tested if
  ## I got the logic right!
  run = function(x = double(2),
    nrep = integer(),
    psi1 = double(),
    phi = double(1),
    gamma = double(1),
    p = double(1),
    log = integer(0, default = 0)) {
    prob1 <- psi1 * p[1]
    numObs <- sum(x[,1]) ## do I have the right orientation?
    ## prob of the occupied sites out of the total sites given p
    ProbOccAndCount <- psi1 * dbinom(numObs, size = nrep, p = p[1], log = 0)
    ## prob of the empty sites
    ProbUnoccAndCount <- (1-psi1) * (numObs == 0)
    ## probably of the observed states
    ProbCount <- ProbOccAndCount + ProbUnoccAndCount
    ProbOccGivenCount <- ProbOccAndCount / ProbCount
    ProbOccNextTime <- ProbOccGivenCount * phi[1] +
      (1-ProbOccGivenCount) * gamma[1]
    ll <- log(ProbCount)
    nyears <- dim(x)[2]
    for(t in 2:nyears) {
      numObs <- sum(x[,t])
      ProbOccAndCount <- ProbOccNextTime *
        dbinom(numObs, size = nrep, p = p[t], log = 0)
      ProbUnoccAndCount <- (1-ProbOccNextTime) * (numObs == 0)
      ProbCount <- ProbOccAndCount + ProbUnoccAndCount
      ProbOccGivenCount <- ProbOccAndCount / ProbCount
      ll <- ll + log(ProbCount)
      if(t < nyears) ProbOccNextTime <- ProbOccGivenCount * phi[t] +
        (1-ProbOccGivenCount) * gamma[t]
    }
    if(log) return(ll)
    else return(exp(ll))
    returnType(double())
  }
  )

rDynamicOccupancy <- nimbleFunction(
  run = function(n = integer(),
    nrep = integer(),
    psi1 = double(),
    phi = double(1),
    gamma = double(1),
    p = double(1),
    log = integer(0, default = 0)) {
    nyear <- length(p)
    ans <- double(2)
    setSize(ans, nrep, nyear)
    ## could populate ans here, but I'm just doing this as a placeholder
    returnType(double(2))
    return(ans)
  }
  )

registerDistributions(list(
  dDynamicOccupancy = list(
    BUGSdist = "dDynamicOccupancy(nrep, psi1, phi, gamma, p)",
    Rdist = "dDynamicOccupancy(nrep, psi1, phi, gamma, p)",
    types = c('value = double(2)',
      'nrep = integer(0)',
      'psi1 = double()',
      'phi = double(1)',
      'gamma = double(1)',
      'p = double(1)'))
  ))

## *********************************************************************
## data prep
## *********************************************************************

logit <- function(x) {
  log(x/(1 - x))
}

expit <- function(x) {
  exp(x)/(1 + exp(x))
}

## Function to simulate detection/nondetection data for dynamic site-occ model
## Annual variation in probabilities of patch survival, colonization and
## detection is specified by the bounds of a uniform distribution.
genDynamicOccData <- function(R = 250,
                              J = 3, K = 10,
                              psi1 = 0.4,
                              range.p = c(0.2, 0.4),
                              range.phi = c(0.6, 0.8),
                              range.gamma = c(0, 0.1)) {
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


## prep data for nimble model
prepModDataOcc <- function(data,
                           monitors = c("psi",
                             "phi",
                             "gamma",
                             "p"),
                           include.zs=TRUE){
  ## data zs with 0s set to NAs
  zs <- apply(data$y, c(1, 3), max)
  zs[zs == 0] <- NA

  ## initial condiations, NAs where 1s are in z, and 1s are where NA
  zinits <- zs
  zinits[zinits == 1] <- 2
  zinits[is.na(zinits)] <- 1
  zinits[zinits == 2] <- NA
  inits <- list(z = zinits)

  ## constants
  constants <- list(nsite = dim(data$y)[1],
                    nrep = dim(data$y)[2],
                    nyear = dim(data$y)[3])
  if(include.zs){
    model.data <- list(y = data$y, z = zs)
    inits <- list(z = zinits)
  } else{
    model.data <- list(y = data$y)
    inits <- list()
  }
  model.input <- list(data=model.data,
                      monitors=monitors,
                      constants=constants,
                      inits=inits)
  return(model.input)
}

data <- genDynamicOccData()
model.input <- prepModDataOcc(data, include.zs=FALSE)


## MCMC settings
scale <- 1e1
burnin <- 1e1*scale
niter <- (1e3)*scale

## *********************************************************************
##  Multi-season occupancy model: option 4-5 remove latent states using
##  user-defined NIMBLE function
##  *********************************************************************

## Specify model in NIMBLE
ss.ms.occ <- nimbleCode({
  ##  priors
  psi1 ~ dunif(0, 1)

  psi[1] <- psi1
  for(k in 1:(nyear-1)){
    phi[k] ~ dunif(0, 1)
    gamma[k] ~ dunif(0, 1)
    p[k] ~ dunif(0, 1)
  }
  p[nyear] ~ dunif(0, 1)

  ## Ecological submodel: Define state conditional on parameters
  for(i in 1:nsite) {
    ## removes the z's and muZ's from the model and compute
    ## the probability of all reps over all years for one site.
    y[i, 1:nrep, 1:nyear] ~ dDynamicOccupancy(nrep,
                                              psi1,
                                              phi[1:(nyear-1)],
                                              gamma[1:(nyear-1)],
                                              p[1:nyear])
  }
})

## *********************************************************************
## opt 4 run with compareMCMCs
## *********************************************************************

input1 <- c(code=ss.ms.occ,
            model.input)

ss.ms.opt4 <- compareMCMCs(input1,
                           MCMCs=c('nimble', 'autoBlock', 'nimble_slice'),
                           niter=niter,
                           burnin = burnin,
                           summary=FALSE,
                           check=FALSE)
