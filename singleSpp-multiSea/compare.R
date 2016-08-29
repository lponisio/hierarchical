rm(list=ls())
setwd('~/Dropbox/occupancy-nimble/singleSpp-multiSea')
source('src/initialize.R')

## original model jags and nimble
load(file="saved/orig.Rdata")

## vectorized bernuli calls
load(file="saved/opt1.Rdata")

## custom sampler for zs, slice for other parms
load(file="saved/opt2.Rdata")

## custom sampler for zs, reflective sampler for other parms
load(file="saved/opt3.Rdata")

## costum function for latent state
load(file="saved/opt4.Rdata")

## costum function for latent state + block samplers on phi[i-1],
## gamma[i-1]
load(file="saved/opt5.Rdata") 

## rename results
ss.ms.opt1[[1]] <- rename_MCMC_comparison_method('nimble', 'vectorized',
                                                 comparison=ss.ms.opt1[[1]])
ss.ms.opt2[[1]] <- rename_MCMC_comparison_method('nimbleOpt2',
                                                 'slice',
                                                 comparison=ss.ms.opt2[[1]])
ss.ms.opt3[[1]] <- rename_MCMC_comparison_method('nimbleOpt3',
                                                 'reflective',
                                                 comparison=ss.ms.opt3[[1]])
ss.ms.opt4[[1]] <- rename_MCMC_comparison_method(c('nimble', 'autoBlock'),
                                                 c('no z',
                                                   'no z, autoBlock'),
                                                 comparison=ss.ms.opt4[[1]])
ss.ms.opt5[[1]] <- rename_MCMC_comparison_method('nimbleOpt5',
                                                 'block phi gam',
                                                 comparison=ss.ms.opt5[[1]])
## compare mcmcs
ss.ms.occ.all <- combine_MCMC_comparison_results(ss.ms.orig[[1]],
                                                 ss.ms.opt1[[1]],
                                                 ss.ms.opt2[[1]],
                                                 ss.ms.opt3[[1]],
                                                 ss.ms.opt4[[1]],
                                                 ss.ms.opt5[[1]],
                                                 name = "ss.ms" )

make_MCMC_comparison_pages(ss.ms.occ.all, dir="figures/comparisons")
