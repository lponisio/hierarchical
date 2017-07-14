rm(list=ls())
library(nimble)
library(coda)

setwd('~/Dropbox/nimble/occupancy/analysis/multiSpp-singleSea')
source('../all/plotting.R')
save.dir <-  "../../../saved/multiSpp-singleSea/saved"

## original model jags and nimble
load(file=file.path(save.dir, "orig.Rdata"))

## vectorized, likelihood for latent state, derived quantity
## calculation
load(file=file.path(save.dir, "filter.Rdata"))

## cross level sampler
load(file=file.path(save.dir, "crosslevel.Rdata"))

## subsampling
load(file=file.path(save.dir, "subsamp.Rdata"))

## filter + custom block sampler on species random effect for each
## species
## load(file=file.path(save.dir, "opt2.Rdata"))

## filter + custom block sampler on species random effect for each
## random effect type
## load(file=file.path(save.dir, "opt3.Rdata")

## filter + sigma sampler on random effects
## load(file=file.path(save.dir, "opt4.Rdata"))



ms.ss.orig[[1]] <- rename_MCMC_comparison_method(c('nimble', 'jags'),
                                                 c('NIMBLE-latent',
                                                   'JAGS-latent'),
                                                 comparison=ms.ss.orig[[1]])

ms.ss.filter[[1]] <- rename_MCMC_comparison_method(c('nimble'),
                                                 c('NIMBLE-filter'),
                                                 comparison=ms.ss.filter[[1]])

ms.ss.crosslevel[[1]] <- rename_MCMC_comparison_method('nimbleCrosslevel',
                                                       'NIMBLE-cross-level',
                                                       comparison=ms.ss.crosslevel[[1]])

ms.ss.subsamp[[1]] <- rename_MCMC_comparison_method('nimbleSubsamp',
                                                    'NIMBLE-subsample',
                                                    comparison=ms.ss.subsamp[[1]])


## ms.ss.opt2[[1]] <- rename_MCMC_comparison_method('nimbleOpt2',
##                                                  'filter + sp. block',
##                                                  comparison=ms.ss.opt2[[1]])
## ms.ss.opt3[[1]] <- rename_MCMC_comparison_method('nimbleOpt3',
##                                                  'block_2',
##                                                  comparison=ms.ss.opt3[[1]])
## ms.ss.opt4[[1]] <- rename_MCMC_comparison_method('nimbleOpt4',
##                                                  'sigma sampler',
##                                                  comparison=ms.ss.opt4[[1]])

## compare mcmcs
ms.ss.occ.all <- combine_MCMC_comparison_results(ms.ss.orig[[1]],
                                                 ms.ss.filter[[1]],
                                                 ms.ss.crosslevel[[1]],
                                                 ms.ss.subsamp[[1]],
                                                 name = "ms.ss")

make_MCMC_comparison_pages(ms.ss.occ.all,
                           dir=file.path(save.dir, "../figures/comparisons"))


## look at samples

checkChains(ms.ss.occ.all$ms.ss$samples,
            f.path = file.path(save.dir,
            "../figures/chains/%s.pdf")
)


## ****************************************
## custom figs
## ****************************************

by.param <- apply(ms.ss.occ.all$ms.ss$samples, c(1,2), effectiveSize)/
  ms.ss.occ.all$ms.ss$timing
by.config <- ms.ss.occ.all$ms.ss$efficiency

source('../all/plotting.R')
plotEffSize(by.config, by.param, f.path= file.path(save.dir,
              "../figures/comparisons/%s%s.pdf"), "MultiSpp-SingleSea",
            at=0.4, adj1=0.03, adj2=0.1)

(ms.ss.occ.all$ms.ss$efficiency$mean["NIMBLE-latent"] - ms.ss.occ.all$ms.ss$efficiency$mean["JAGS-latent"])/max(ms.ss.occ.all$ms.ss$efficiency$mean["NIMBLE-latent"])

(ms.ss.occ.all$ms.ss$efficiency$mean["NIMBLE-filter"] - ms.ss.occ.all$ms.ss$efficiency$mean["JAGS-latent"])/max(ms.ss.occ.all$ms.ss$efficiency$mean["NIMBLE-filter"])

