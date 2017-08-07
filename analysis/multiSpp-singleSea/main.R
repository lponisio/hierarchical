rm(list=ls())
setwd('occupancy/analysis/multiSpp-singleSea')

## original model
source('original.R')

## remove Zs and add block samplers to species random effects
source('filtering.R')

## comparisons
source('compare.R')
