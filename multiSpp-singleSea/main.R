rm(list=ls())
setwd('~/Dropbox/occupancy-nimble/multiSpp-singleSea')

## original model
source('original.R')

## remove Zs and add block samplers to species random effects
source('opt1-4.R')

## comparisons
source('compare.R')
