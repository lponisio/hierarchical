rm(list=ls())
setwd("~/Dropbox/nimble-dev/occupancy/singleSpp-multiSea")

## original model
source('original.R')

## vectorized Bernoulli calls
source('opt1.R')

## custom samplers for zs and reflective sampler
source('opt2.R')

## custom function for latent state
source('opt3.R')
