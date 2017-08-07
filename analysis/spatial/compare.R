rm(list=ls())
library(nimble)
setwd("occupancy/analysis/spatial")
source('../all/plotting.R')
source('src/initialize.R')


## original model jags and nimble
load(file=file.path(save.dir, "orig.Rdata"))

## ## vanilla nimble and auto block
load(file=file.path(save.dir, "nimble.Rdata"))

## ## custom sampler for zs, slice for other parms
load(file=file.path(save.dir, "AFSlice.Rdata"))

## ## rename results
sp.AFSlice[[1]] <- rename_MCMC_comparison_method('nimbleAFSlice',
                                              'AF slice',
                                              comparison=sp.AFSlice[[1]])
## compare mcmcs
sp.occ.all <- combine_MCMC_comparison_results(sp.orig[[1]],
                                              sp.nimble[[1]],
                                              sp.AFSlice[[1]],
                                              name = "sp" )
make_MCMC_comparison_pages(sp.occ.all,
                           dir=file.path(save.dir,
                                         "../figures/comparisons"))

checkChains(sp.occ.all[[1]]$samples,
            f.path = file.path(save.dir,
                               "../figures/chains/%s.pdf"))


## $logDelta
## [1] -2.302585

## $logSigma
## [1] 2.079442

## $alpha
## [1] 0.5

## $p
## [1] 0.8
