rm(list=ls())
library(nimble)
library(coda)
setwd("~/Dropbox/nimble/occupancy/analysis/singleSpp-multiSea")
source('../all/plotting.R')
save.dir <- "../../../saved/singleSpp-multiSea/saved"

## original model jags and nimble
load(file=file.path(save.dir, "orig.Rdata"))

## subsample latent states
load(file=file.path(save.dir, "subsamp.Rdata"))

## cross level sampler
load(file=file.path(save.dir, "crosslevel.Rdata"))

## filter over latent states
load(file=file.path(save.dir, "filter.Rdata"))


## rename results

ss.ms.orig[[1]] <- rename_MCMC_comparison_method(c('nimble', 'jags'),
                                                 c('NIMBLE-latent',
                                                   'JAGS-latent'),
                                                 comparison=ss.ms.orig[[1]])


ss.ms.crosslevel[[1]] <- rename_MCMC_comparison_method('nimbleCrosslevel',
                                                       'Cross-level',
                                                       comparison=ss.ms.crosslevel[[1]])

ss.ms.subsamp[[1]] <- rename_MCMC_comparison_method('nimbleSubsamp',
                                                    'Subsample latent',
                                                    comparison=ss.ms.subsamp[[1]])

ss.ms.filter[[1]] <- rename_MCMC_comparison_method('nimble',
                                                   'Filter',
                                                   comparison=ss.ms.filter[[1]])

## compare mcmcs
ss.ms.occ.all <- combine_MCMC_comparison_results(ss.ms.orig[[1]],
                                                 ss.ms.crosslevel[[1]],
                                                 ss.ms.subsamp[[1]],
                                                 ss.ms.filter[[1]],
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

by.param <- apply(ss.ms.occ.all[[1]]$samples, c(1,2), effectiveSize)/
    ss.ms.occ.all[[1]]$timing
by.config <- ss.ms.occ.all[[1]]$efficiency


source('../all/plotting.R')
plotEffSize(by.config, by.param, f.path= file.path(save.dir,
                                                   "../figures/comparisons/%s%s.pdf"),
            "SingleSpp-MultiSea",
            at=9, adj1=1, adj2=0.3, widths=c(4.5, 8.5))
