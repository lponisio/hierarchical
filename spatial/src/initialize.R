## library(devtools)
## install_github("nimble-dev/nimble",
##                ref = "devel",
##                subdir = "packages/nimble")

library(nimble)
library(igraph)
library(raster)

expit <- function(x) 1/(1+exp(-x))
logit <- function(x) log(x/(1-x))

## parameters
ngrid <- 25
nsite <- 20
nreps <- 5

alpha <- 0.1
beta1 <- 2
p <- 0.6
psi <- 0.7
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

## calculate probabilities of occurrence
psi <- expit(alpha + beta1 * values(elev) +  values(Xraster))

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

## model data
model.data <- list(D = distance[sites, sites],
                   y = ysample[sites,],
                   zeros=rep(0, nsite),
                   elev = fulldata$elevation[sites])
## constants
constants <- list(nsite = nsite, nreps=nreps)

## parameters to monitor
monitors <- c("delta", "sigma", "psi", "p", "alpha", "b1")

## inits
inits <- list()

## MCMC settings
scale <- 1e1
burnin <- 1e1*scale
niter <- (1e3)*scale

