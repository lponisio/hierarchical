rm(list=ls())
library(raster)
library('rjags')
library('R2jags')
load.module("msm")

## *********************************************************************
## simulate data
## *********************************************************************
expit <- function(x) 1/(1 + exp(-x))

rmvn <- function(n, mu = 0, V = matrix(1)) {
  p <- length(mu)
  if (any(is.na(match(dim(V), p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n * p), ncol = p) %*% D + rep(mu, rep(n, p)))
}


genSpatialOccData <- function(ngrid = 50,
                              nreps = 50,
                              alpha = 1,
                              b1 = 0.5,
                              p = 0.8,
                              sigma = 8,
                              delta = 0.1){
  ## Set up a square lattice region
  simgrid <- expand.grid(1:ngrid, 1:ngrid)
  n <- nrow(simgrid)

  ## Set up distance matrix
  dist.mat <- as.matrix(dist(simgrid))

  prep.cor.mat <- exp(-delta * dist.mat)
  cor.mat <- sigma^2*(0.95*prep.cor.mat +
                      0.05*diag(nrow(prep.cor.mat)))
  ## Generate spatial random effect
  X <- rmvn(1, rep(0, n), cor.mat)

  Xraster <- rasterFromXYZ(cbind(simgrid[, 1:2] - 0.5, X))

  ## simulate elevation data
  elev <- raster(matrix(rnorm(n), ngrid, ngrid),
                 xmn = 0, xmx = ngrid,
                 ymn = 0, ymx = ngrid)
  elev <- scale(elev)

  ## calculate probabilities of occurrence
  psi <- expit(alpha + b1 * raster::values(elev) +
               raster::values(Xraster))


  ## Latent occurrence state
  z <- rbinom(n = n, size = 1, prob = psi)
  z <- rasterFromXYZ(cbind(coordinates(elev), z))
  quartz()
  plot(z)

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
              distance=dist.mat,
         inits=list(delta=delta,
           sigma=sigma,
           alpha=alpha,
           b1=b1,
           p=p)))
}


## preps data for nimble
prepModData <- function(fulldata, y, dist.mat, nsite, inits,
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
  inits$z <- zinits

  model.data <- list(D = dist.mat[sites, sites],
                     y = y,
                     z = zs,
                     zeros=rep(0, nsite),
                     elev = fulldata$elevation[sites],
                     DI = diag(nsite))

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
                           nsite=250, inits=dats$inits)

## *********************************************************************
## jags model
## *********************************************************************

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

  dist.mat[1:nsite, 1:nsite] <- -delta*D[1:nsite, 1:nsite]
  prep.cov[1:nsite, 1:nsite] <- mexp(dist.mat[1:nsite, 1:nsite])
  D.cov[1:nsite, 1:nsite] <- (sigma^2)*(0.95*prep.cov[1:nsite, 1:nsite] + 0.05*DI[1:nsite, 1:nsite])

 ##   for(i in 1:nsite){
 ##   for(j in 1:nsite){
 ##     prep.cov[i, j]  <- exp(-delta*D[i, j])
 ##     D.cov[i, j] <- (sigma^2)*(0.95*prep.cov[i, j] + 0.05*DI[i, j])
 ##   }
 ## }



  D.tau[1:nsite, 1:nsite] <- (D.cov[1:nsite, 1:nsite])

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
res <- analyse.jags(model.input,
                     ni=(1e3+1e1)*scale,
                     nt=scale,
                     nb=1e1*scale,
                     nc=3)


