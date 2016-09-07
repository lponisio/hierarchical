rm(list=ls())
setwd("~/Dropbox/nimble/occupancy/analysis/spatial")
save.dir <-  "../../../saved/spatial/saved"

## original model jags and nimble
load(file=file.path(save.dir, "orig.Rdata"))

## vanilla nimble and auto block
load(file=file.path(save.dir, "opt1.Rdata"))

## custom sampler for zs, slice for other parms
load(file=file.path(save.dir, "opt2.Rdata"))

## custom sampler for zs, reflective sampler for other parms
load(file=file.path(save.dir, "opt3.Rdata"))

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

make_MCMC_comparison_pages(sp.occ.all,
                           dir=file.path(save.dir, "figures/comparisons"))
