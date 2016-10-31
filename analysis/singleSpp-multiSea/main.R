rm(list=ls())
setwd("~//Dropbox//nimble//occupancy//analysis//singleSpp-multiSea")

## original model
source('original.R')


## custom samplers for zs and reflective sampler, or slice sampler
source('opt2-3.R')

## custom function for latent state,  block sampler of phi[i-1],
## gamma[i-1]
source('opt4-5.R')

