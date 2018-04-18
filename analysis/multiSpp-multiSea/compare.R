## setwd('~/Dropbox/occupancy')
rm(list=ls())
library(nimble)
library(coda)
setwd('analysis/multiSpp-multiSea')
source('../all/plotting.R')
source('src/initialize.R')

## original model in nimble
load(file=file.path(save.dir, "jags.Rdata"))

## original model in nimble
load(file=file.path(save.dir, "nimble.Rdata"))

## vectorized, likelihood for latent state, derived quantity
## calculation
load(file=file.path(save.dir, "filtering.Rdata"))

## subsampling
load(file=file.path(save.dir, "subsamp.Rdata"))

ms.ms.orig.nim[[1]] <- rename_MCMC_comparison_method(c('nimble'),
                                                 c('NIMBLE-latent' ),
                                              comparison=ms.ms.orig.nim[[1]])


ms.ms.orig[[1]] <- rename_MCMC_comparison_method(c('jags'),
                                                 c('JAGS-latent' ),
                                              comparison=ms.ms.orig[[1]])

ms.ms.filter <- ms.ms.nimble
ms.ms.filter[[1]] <- rename_MCMC_comparison_method(c('nimble'),
                                                 c('NIMBLE-filter'),
                                              comparison=ms.ms.filter[[1]])

## ms.ms.crosslevel[[1]] <- rename_MCMC_comparison_method('nimbleCrosslevel',
##                                                        'NIMBLE-cross-level',
##                                               comparison=ms.ms.crosslevel[[1]])

## ms.ms.subsamp[[1]] <- rename_MCMC_comparison_method('nimbleSubsamp',
##                                                     'NIMBLE-subsample',
##                                               comparison=ms.ms.subsamp[[1]])

## compare mcmcs
ms.ms.occ.all <- combine_MCMC_comparison_results(ms.ms.orig[[1]],
                                                 ms.ms.orig.nim[[1]],
                                                 ms.ms.filter[[1]],
                                                 ## ms.ms.crosslevel[[1]],
                                                 ## ms.ms.subsamp[[1]],
                                                 name = "ms.ms")

save(ms.ms.occ.all, file=file.path(save.dir, "combined.Rdata"))

make_MCMC_comparison_pages(ms.ms.occ.all,
                           dir=file.path(save.dir, "../figures/comparisons"))


## look at samples

checkChains(ms.ms.occ.all$ms.ms$samples,
            f.path = file.path(save.dir,
                               "../figures/chains/%s.pdf"))

## ****************************************
## custom figs
## ****************************************

by.param <- apply(ms.ms.occ.all$ms.ms$samples, c(1,2), effectiveSize)/
  ms.ms.occ.all$ms.ms$timing
by.config <- ms.ms.occ.all$ms.ms$efficiency

source('../all/plotting.R')
plotEffSize(by.config, by.param, f.path= file.path(save.dir,
              "../figures/comparisons/%s%s.pdf"), "MultiSpp-multiSea",
            at=0.3, adj1=0.03, adj2=0.1)

(ms.ms.occ.all$ms.ms$efficiency$mean["NIMBLE-latent"] -
 ms.ms.occ.all$ms.ms$efficiency$mean["JAGS-latent"])/
    max(ms.ms.occ.all$ms.ms$efficiency$mean["NIMBLE-latent"])

(ms.ms.occ.all$ms.ms$efficiency$mean["NIMBLE-filter"] -
 ms.ms.occ.all$ms.ms$efficiency$mean["JAGS-latent"])/
    max(ms.ms.occ.all$ms.ms$efficiency$mean["NIMBLE-filter"])

