rm(list=ls())
library(nimble)
library(coda)
setwd("~/Dropbox/occupancy")
setwd("analysis/singleSpp-multiSea")

source('../all/plotting.R')
source('src/initialize.R')


## original model jags and nimble
load(file=file.path(save.dir, "orig.Rdata"))

## subsample latent states
load(file=file.path(save.dir, "subsamp.Rdata"))

## cross level sampler
## load(file=file.path(save.dir, "crosslevel.Rdata"))

## filter over latent states
load(file=file.path(save.dir, "filter.Rdata"))

## blocking with AFSS
load(file=file.path(save.dir, "blocking.Rdata"))


##  slice sampleers
load(file=file.path(save.dir, "slice.Rdata"))



## rename results
ss.ms.orig[[1]] <- rename_MCMC_comparison_method(c('nimble', 'jags'),
                                                 c('NIMBLE-latent',
                                                   'JAGS-latent'),
                                              comparison=ss.ms.orig[[1]])


## ss.ms.crosslevel[[1]] <- rename_MCMC_comparison_method('nimbleCrosslevel',
##                                                        'NIMBLE-cross-level',
##                                               comparison=ss.ms.crosslevel[[1]])

ss.ms.subsamp[[1]] <- rename_MCMC_comparison_method('nimbleSubsamp',
                                                    'NIMBLE-subsample',
                                              comparison=ss.ms.subsamp[[1]])

ss.ms.filter[[1]] <- rename_MCMC_comparison_method('nimble',
                                                   'NIMBLE-filter',
                                                   comparison=ss.ms.filter[[1]])

ss.ms.blocking[[1]] <- rename_MCMC_comparison_method('blocking',
                                                   'NIMBLE-blocking',
                                              comparison=ss.ms.blocking[[1]])


ss.ms.slice[[1]] <- rename_MCMC_comparison_method(c('nim_slice', 'nim_AFSS'),
                                                   c('NIMBLE-slice', 'NIMBLE-AFSS'),
                                              comparison=ss.ms.slice[[1]])


## compare mcmcs
ss.ms.occ.all <- combine_MCMC_comparison_results(ss.ms.orig[[1]],
                                                 ss.ms.blocking[[1]],
                                                 ss.ms.slice[[1]],
                                                 ## ss.ms.crosslevel[[1]],
                                                 ss.ms.subsamp[[1]],
                                                 ss.ms.filter[[1]],
                                                 name = "ss.ms" )

make_MCMC_comparison_pages(ss.ms.occ.all,
                   dir=file.path(save.dir, "../figures/comparisons"))


checkChains(ss.ms.occ.all[[1]]$samples,
            f.path = file.path(save.dir,
                               "../figures/chains/%s.pdf"))


## ****************************************
## custom figs
## ****************************************

by.param <- apply(ss.ms.occ.all$ss.ms$samples, c(1,2), effectiveSize)/
    ss.ms.occ.all$ss.ms$timing
by.config <- ss.ms.occ.all$ss.ms$efficiency


source('../all/plotting.R')
plotEffSize(by.config, by.param, f.path= file.path(save.dir,
                               "../figures/comparisons/%s%s.pdf"),
                                "SingleSpp-MultiSea",
            at=50, adj1=1, adj2=0.3, widths=c(4.5, 8.5))

## difference between default Jags and Nimble
(ss.ms.occ.all$ss.ms$efficiency$mean["JAGS-latent"] -
 ss.ms.occ.all$ss.ms$efficiency$mean["NIMBLE-latent"])/
    max(ss.ms.occ.all$ss.ms$efficiency$mean["JAGS-latent"])

## jags and filtering
(ss.ms.occ.all$ss.ms$efficiency$mean["JAGS-latent"] -
 ss.ms.occ.all$ss.ms$efficiency$mean["NIMBLE-filter"])/
    max(ss.ms.occ.all$ss.ms$efficiency$mean["JAGS-latent"])

## jags and subsampling
(ss.ms.occ.all$ss.ms$efficiency$mean["JAGS-latent"] -
 ss.ms.occ.all$ss.ms$efficiency$mean["NIMBLE-subsample"])/
    max(ss.ms.occ.all$ss.ms$efficiency$mean["JAGS-latent"])

## subsampling and nimble defaults
(ss.ms.occ.all$ss.ms$efficiency$mean["NIMBLE-subsample"] -
 ss.ms.occ.all$ss.ms$efficiency$mean["NIMBLE-latent"])/
    max(ss.ms.occ.all$ss.ms$efficiency$mean["NIMBLE-subsample"])

