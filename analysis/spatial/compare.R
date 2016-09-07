rm(list=ls())
setwd('~/Dropbox/occupancy-nimble/spatial')
source('src/initialize.R')

## original model jags and nimble
load(file="saved/orig.Rdata")

## vanilla nimble and auto block
load(file="saved/opt1.Rdata")

## custom sampler for zs, slice for other parms
load(file="saved/opt2.Rdata")

## custom sampler for zs, reflective sampler for other parms
load(file="saved/opt3.Rdata")

## rename results
sp.opt2[[1]] <- rename_MCMC_comparison_method('nimbleOpt2',
                                                 'slice',
                                                 comparison=sp.opt2[[1]])
sp.opt3[[1]] <- rename_MCMC_comparison_method('nimbleOpt3',
                                                 'reflective',
                                                 comparison=sp.opt3[[1]])
## compare mcmcs
sp.occ.all <- combine_MCMC_comparison_results(sp.orig[[1]],
                                                 sp.opt1[[1]],
                                                 sp.opt2[[1]],
                                                 sp.opt3[[1]],
                                                 name = "sp" )

make_MCMC_comparison_pages(sp.occ.all, dir="figures/comparisons")
