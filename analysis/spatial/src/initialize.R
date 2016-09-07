## library(devtools)
## install_github("nimble-dev/nimble",
##                ref = "devel",
##                subdir = "packages/nimble")


library(nimble)
library(igraph)
library(raster)

source("../all/plotting.R")
source("../all/runNimble.R")

save.dir <-  "../../../saved/spatial/saved"

## samplers
source("../all/samplers/sampler_z.R")
source("../all/samplers/sampler_reflective.R")

## parameters
ngrid <- 25
nsite <- 20
nreps <- 5

alpha <- 0.1
beta1 <- 2
p <- 0.6
sigma <- 0.5
delta <- 0.05

set.seed(44)

rmvn <- function(n, mu = 0, V = matrix(1)) {
  p <- length(mu)
  if (any(is.na(match(dim(V), p)))) 
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n * p), ncol = p) %*% D + rep(mu, rep(n, p)))
}

## Set up a square lattice region
simgrid <- expand.grid(1:ngrid, 1:ngrid)
n <- nrow(simgrid)

## Set up distance matrix
distance <- as.matrix(dist(simgrid))

## Generate spatial random variable
X <- rmvn(1, rep(0, n),  (sigma^2)*exp(-delta * distance))
Xraster <- rasterFromXYZ(cbind(simgrid[, 1:2] - 0.5, X))

## simulate elevation data
elev <- raster(matrix(rnorm(n), ngrid, ngrid),
               xmn = 0, xmx = ngrid,
               ymn = 0, ymx = ngrid)
elev <- scale(elev)

## calculate probabilities of occurrence
expit <- function(x) 1/(1 + exp(-x))
psi <- expit(alpha + beta1 * raster::values(elev) +  raster::values(Xraster))

## Latent occurrence state
z <- rbinom(n = n, size = 1, prob = psi) 
z <- rasterFromXYZ(cbind(coordinates(elev), z))

coords <- coordinates(z)
fulldata <- data.frame(coords,
                       elevation = extract(elev, coords),
                       Xvar = extract(Xraster, coords),
                       z = extract(z, coords))

## subsample at "sites"
sites <- sample(1:n, nsite)

## Observation process: Sample detection/nondetection observations
## from a Bernoulli(with p) if z=1
y <- matrix(NA, nrow = n, ncol = nreps)

for (j in 1:nreps){
  y[,j] <- rbinom(n = n, size = 1, prob = fulldata$z * p)
}

## NA for "sites" that are not sampled
ysample <- y
ysample[-sites, ] <- NA
## drop ys that are not sampled for this analysis
y <- ysample[sites,]

## data zs with 0s set to NAs
zs <- apply(y, 1, max)
zs[zs == 0] <- NA

## initial conditions, NAs where 1s are in y, and 1s are where NA
zinits <- zs
zinits[zinits == 1] <- 2
zinits[is.na(zinits)] <- 1
zinits[zinits == 2] <- NA
inits <- list(z = zinits)

## model data
model.data <- list(D = distance[sites, sites],
                   y = y,
                   z = zs,
                   zeros=rep(0, nsite),
                   elev = fulldata$elevation[sites])
## constants
constants <- list(nsite = nsite, nreps=nreps)

## parameters to monitor
monitors <- c("delta", "sigma", "psi", "p", "alpha", "b1")


## MCMC settings
scale <- 5e1^2
burnin <- 1e1*scale
niter <- (1e3)*scale

