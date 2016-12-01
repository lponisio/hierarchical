
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
                              beta1 = 6,
                              p = 0.8,
                              sigma = 0.5,
                              delta = 0.5){
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
  quartz()
  plot(Xraster)
  ## simulate elevation data
  elev <- raster(matrix(rnorm(n), ngrid, ngrid),
                 xmn = 0, xmx = ngrid,
                 ymn = 0, ymx = ngrid)
  elev <- scale(elev)

  ## calculate probabilities of occurrence
  psi <- expit(alpha + beta1 * raster::values(elev) +
               raster::values(Xraster))
  
  ## quartz()
  ## rasterFromXYZ(cbind(coordinates(elev), psi))

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
              distance=dist.mat))
}


## preps data for nimble
prepModData <- function(fulldata, y, dist.mat, nsite,
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
  inits <- list(z = zinits,
                sigma=0.1,
                delta=0.1)

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


