rm(list=ls())
library(nimble)
library(igraph)
library(raster)

## MCMC settings
scale <- 1e2
burnin <- 1e1*scale
niter <- (1e3)*scale

## *********************************************************************
## simulate data
## *********************************************************************
expit <- function(x) 1/(1 + exp(-x))

## random multivariate normal draw
rmvn <- function(n, mu = 0, V = matrix(1)) {
  p <- length(mu)
  if (any(is.na(match(dim(V), p)))) 
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n * p), ncol = p) %*% D + rep(mu, rep(n, p)))
}

## simulate detection with covariate and spatial random effect data
## set
genSpatialOccData <- function(ngrid = 50, ## number of spatial grid cells
                              nreps = 10, ## number of repeated
                              ## samples at sites
                              alpha = 1, ## intercept 
                              beta1 = 2, ## covariate coefficient
                              p = 0.6, ## prob detection
                              sigma = 0.5, 
                              delta = 0.5){ ## decay parameter

  ## Set up a square lattice region
  simgrid <- expand.grid(1:ngrid, 1:ngrid)
  n <- nrow(simgrid)

  ## Set up distance matrix
  distance <- as.matrix(dist(simgrid))

  ## Generate spatial random effect
  X <- rmvn(1, rep(0, n),  (sigma^2)*exp(-delta * distance))
  Xraster <- rasterFromXYZ(cbind(simgrid[, 1:2] - 0.5, X))

  ## simulate "elevation" data
  elev <- raster(matrix(rnorm(n), ngrid, ngrid),
                 xmn = 0, xmx = ngrid,
                 ymn = 0, ymx = ngrid)
  elev <- scale(elev)

  ## calculate probabilities of occurrence
  psi <- expit(alpha + beta1 * raster::values(elev) +
               raster::values(Xraster))

  ## Latent occurrence state
  z <- rbinom(n = n, size = 1, prob = psi) 
  z <- rasterFromXYZ(cbind(coordinates(elev), z))

  coords <- coordinates(z)
  fulldata <- data.frame(coords,
                         elevation = extract(elev, coords),
                         Xvar = extract(Xraster, coords),
                         z = extract(z, coords))

  ## Observation process: Sample detection/nondetection observations
  ## from a Bernoulli(with p) if z=1
  y <- matrix(NA, nrow = n, ncol = nreps)

  for (j in 1:nreps){
    y[,j] <- rbinom(n = n, size = 1, prob = fulldata$z * p)
  }

  return(list(data=fulldata,
              y=y,
              distance=distance))
}


## preps data for nimble
prepModData <- function(fulldata, y, distance, nsite,
                        monitors=c("delta", "sigma", "psi",
                          "p", "alpha", "b1")){
  ## subsample at "sites" (create a grid of sites to avoid any that
  ## are too close to eachother)
  sites <- round(seq(from=1, to=nrow(fulldata), length=nsite))
  y <- y[sites,]
  
  ## data zs with 0s set to NAs
  zs <- apply(y, 1, max)
  zs[zs == 0] <- NA

  ## initial conditions, NAs where 1s are in y, and 1s are where NA
  zinits <- zs
  zinits[zinits == 1] <- 2
  zinits[is.na(zinits)] <- 1
  zinits[zinits == 2] <- NA
  inits <- list(z = zinits)

  model.data <- list(D = distance[sites, sites],
                     y = y,
                     z = zs,
                     zeros=rep(0, nsite),
                     elev = fulldata$elevation[sites])
  
  model.input <- list(constants=list(nsite = nsite,
                        nreps=ncol(y)),
                      data=model.data,
                      monitors=monitors,
                      inits=inits)
  
  return(model.input)
}


## simulate data
set.seed(444)
dats <- genSpatialOccData()
model.input <- prepModData(dats$data, dats$y, dats$distance,
                           nsite=50)

## *********************************************************************
## distance decay model for nimble
## *********************************************************************
sp.mod <- nimbleCode({
  ## priors
  delta ~ dunif(0.01, 10)
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

  ## create covariance matrix based on distances
  D.cov[1:nsite, 1:nsite] <- (sigma^2)*exp(-delta*D[1:nsite, 1:nsite])
  
})

input1 <- c(code=sp.mod,
            model.input)

## *********************************************************************
## opt 1:vanilla nimble and auto block, runs but chains look terrible
## in various ways. Some do not converge, delta in particular seems to
## have nothing informing it (posterior looks like prior)
## *********************************************************************

sp.opt1 <- compareMCMCs(input1,
                        MCMCs=c("nimble"),
                        niter=niter,
                        burnin = burnin,
                        summary=FALSE,
                        check=FALSE)


## *********************************************************************
## jags model
## *********************************************************************

library(rjags)
## needed for mexp function
load.module("msm")

## needed for to trick nimble to read R model with mexp
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
                        D.tau[1:nsite, 1:nsite])

  ## create covariance matrix based on distances (must be 1/cov for
  ## JAGS)

  ## mexp is jags's version fo matrix exponentiation, for some reason
  ## the below will not work
  
  ## D.cov[1:nsite, 1:nsite]  <- (sigma^2)*mexp(-delta*D[1:nsite, 1:nsite])

  ## terribly inefficeint but I cannot get mexp to work properly.
  for(i in 1:nsite){
    for(j in 1:nsite){
      temp.cov[i, j] <- -delta*D[i, j]
      D.cov[i, j]  <- (sigma^2)* exp(temp.cov[i, j])
    }
  }
  
  D.tau[1:nsite, 1:nsite] <- inverse(D.cov[1:nsite, 1:nsite])
  
})


input1 <- c(code=sp.mod,
            model.input)


## *********************************************************************
## jags
## *********************************************************************

sp.orig <- compareMCMCs(input1,
                        MCMCs=c("jags"),
                        niter=niter,
                        burnin = burnin,
                        summary=FALSE,
                        check=FALSE)

