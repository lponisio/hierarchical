rm(list=ls())
library(nimble)
setwd("~/Dropbox/nimble/occupancy/analysis/singleSpp-multiSea")
source('../all/plotting.R')
save.dir <- "../../../saved/singleSpp-multiSea/saved"

## original model jags and nimble
load(file=file.path(save.dir, "orig.Rdata"))

## custom sampler for zs, slice for other parms
load(file=file.path(save.dir, "opt2.Rdata"))

## custom function for latent state
load(file=file.path(save.dir, "opt4.Rdata"))

## costum function for latent state + block samplers on phi[i-1],
## gamma[i-1]
load(file=file.path(save.dir, "opt5.Rdata")) 

## rename results

ss.ms.orig[[1]] <- rename_MCMC_comparison_method(c('nimble', 'jags'),
                                                 c('NIMBLE-latent',
                                                   'JAGS-latent'),
                                                 comparison=ss.ms.orig[[1]])


ss.ms.opt2[[1]] <- rename_MCMC_comparison_method('nimbleOpt2',
                                                 'slice',
                                                 comparison=ss.ms.opt2[[1]])

ss.ms.opt4[[1]] <- rename_MCMC_comparison_method(c('nimble',
                                                   'autoBlock',
                                                   'nimble_slice'),
                                                 c('filter',
                                                   'filter + autoblock',
                                                   'filter + slice'),
                                                 comparison=ss.ms.opt4[[1]])

## ss.ms.opt5[[1]] <- rename_MCMC_comparison_method('nimbleOpt5',
##                                                  'block phi gam',
##                                                  comparison=ss.ms.opt5[[1]])
## compare mcmcs
ss.ms.occ.all <- combine_MCMC_comparison_results(ss.ms.orig[[1]],
                                                 ## ss.ms.opt2[[1]],
                                                 ss.ms.opt4[[1]],
                                                 ## ss.ms.opt5[[1]],
                                                 name = "ss.ms" )

make_MCMC_comparison_pages(ss.ms.occ.all,
                           dir=file.path(save.dir, "../figures/comparisons"))


checkChains(ss.ms.occ.all[[1]]$samples,
            f.path = file.path(save.dir,
              "../figures/chains/%s.pdf")
            )

## ****************************************
## custom figs
## ****************************************

by.param <- apply(ss.ms.occ.all[[1]]$samples, c(1,2), effectiveSize)
by.config <- ss.ms.occ.all[[1]]$efficiency

plotEffSize(by.param, by.config, f.path= file.path(save.dir,
              "../figures/comparisons/%s.pdf"), name="singleSpp-MultiSea",
            at=9)
