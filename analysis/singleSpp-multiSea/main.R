rm(list=ls())
setwd("occupancy/analysis/singleSpp-multiSea")

## vanilla nimble and jags, subsample latent states, cross level
## sampler
source('original.R')

## custom function for latent state
source('filtering.R')

