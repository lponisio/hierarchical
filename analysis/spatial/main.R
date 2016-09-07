rm(list=ls())
setwd("~/Dropbox/occupancy-nimble/spatial")

## original model
source('original.R')

## custom samplers for zs and reflective sampler, or slice sampler
source('opt1-3.R')

