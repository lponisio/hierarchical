rm(list=ls())
## gctorture()

library(nimble)
setwd('~/Dropbox/nimble/occupancy/analysis/multiSpp-singleSea')
source('../all/plotting.R')
save.dir <-  "../../../saved/multiSpp-singleSea/saved"

## original model jags and nimble
load(file=file.path(save.dir, "orig.Rdata"))

## vectorized, likelihood for latent state, derived quantity
## calculation
load(file=file.path(save.dir, "opt1.Rdata"))

## option 1 + custom block sampler on species random effect for each
## species
load(file=file.path(save.dir, "opt2.Rdata"))

## ## option 1 + custom block sampler on species random effect for each
## ## random effect type
## load(file=file.path(save.dir, "opt3.Rdata")

## option 1 + sigma sampler on random effects
load(file=file.path(save.dir, "opt4.Rdata"))

ms.ss.opt1[[1]] <- rename_MCMC_comparison_method('nimble', 'remove_z',
                                                 comparison=ms.ss.opt1[[1]])
ms.ss.opt2[[1]] <- rename_MCMC_comparison_method('nimbleOpt2', 'block_1',
                                                 comparison=ms.ss.opt2[[1]])
## ms.ss.opt3[[1]] <- rename_MCMC_comparison_method('nimbleOpt3',
##                                                  'block_2',
##                                                  comparison=ms.ss.opt3[[1]])
ms.ss.opt4[[1]] <- rename_MCMC_comparison_method('nimbleOpt4',
                                                 'sigma sampler',
                                                 comparison=ms.ss.opt4[[1]])

## compare mcmcs
ms.ss.occ.all <- combine_MCMC_comparison_results(ms.ss.orig[[1]],
                                                 ms.ss.opt1[[1]],
                                                 ms.ss.opt2[[1]],
                                                 ## ms.ss.opt3[[1]],
                                                 ms.ss.opt4[[1]],
                                                 name = "ms.ss" )

make_MCMC_comparison_pages(ms.ss.occ.all,
                           dir=file.path(save.dir, "../figures/comparisons"))


## look at samples

checkChains(ms.ss.occ.all[[1]]$samples,
            f.path = file.path(save.dir,
            "../figures/chains/%s.pdf")
)
