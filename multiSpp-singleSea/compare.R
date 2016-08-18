rm(list=ls())
library(nimble)
setwd('~/Dropbox/nimble-dev/occupancy/multiSpp-singleSea')
source('src/plotting.R')

## original model jags and nimble
load(file="saved/orig.Rdata")

## vectorized, likelihood for latent state, derived quantity
## calculation
load(file="saved/opt1.Rdata")

## option 1 + custom block sampler on species random effect for each
## random effect type
load(file="saved/opt2.Rdata")

## option 1 + custom block sampler on species random effect for each
## species
load(file="saved/opt3.Rdata")

ms.ss.opt1[[1]] <- rename_MCMC_comparison_method('nimble', 'remove_z',
                                                 comparison=ms.ss.opt1[[1]])
ms.ss.opt2[[1]] <- rename_MCMC_comparison_method('nimbleOpt2', 'block_1',
                                                 comparison=ms.ss.opt2[[1]])
ms.ss.opt3[[1]] <- rename_MCMC_comparison_method('nimbleOpt3',
                                                 'block_2',
                                                 comparison=ms.ss.opt3[[1]])

## compare mcmcs
ms.ss.occ.all <- combine_MCMC_comparison_results(ms.ss.orig[[1]],
                                                 ms.ss.opt1[[1]],
                                                 ms.ss.opt2[[1]],
                                                 ms.ss.opt3[[1]],
                                                 name = "ms.ss" )

make_MCMC_comparison_pages(ms.ss.occ.all, dir="figures/comparisons")


## look at samples

checkChains(ms.ss.occ.all[[1]]$samples,
            f.path = "figures/comparisons/chains/%s.pdf")
