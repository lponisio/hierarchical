rm(list=ls())
setwd("~/Dropbox/nimble/occupancy/analysis/spatial")

## original model
source('original.R')

## custom samplers for zs and reflective sampler, or slice sampler
source('opt1-3.R')

