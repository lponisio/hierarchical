## setwd("~/Dropbox/occupancy")
rm(list=ls())
library(nimble)
library(coda)
setwd("analysis/singleSpp-multiSea")

source('../all/plotting.R')
source('src/initialize.R')


## original model jags and nimble
load(file=file.path(save.dir, "orig.Rdata"))

## subsample latent states
## load(file=file.path(save.dir, "subsamp.Rdata"))

## filter over latent states
load(file=file.path(save.dir, "filter.Rdata"))

##  AFSS block
load(file=file.path(save.dir, "AFSS_block.Rdata"))

##  AFSS block and filtering
load(file=file.path(save.dir, "filter_AFSS_block.Rdata"))

##  RW block
load(file=file.path(save.dir, "RW_block.Rdata"))

##  RW block and filtering
load(file=file.path(save.dir, "filter_RW_block.Rdata"))

## slice
load(file=file.path(save.dir, "slice.Rdata"))

## ##  slice and filtering
## load(file=file.path(save.dir, "filter_slice.Rdata"))


## rename results
ss.ms.orig[[1]] <- rename_MCMC_comparison_method(c('nimble', 'jags'),
                                                 c('NIMBLE-latent',
                                                   'JAGS-latent'),
                                              comparison=ss.ms.orig[[1]])

ss.ms.filter[[1]] <- rename_MCMC_comparison_method('nimble',
                                                   'NIMBLE-rm-latent',
                                                   comparison=ss.ms.filter[[1]])

ss.ms.AFSblocking[[1]] <- rename_MCMC_comparison_method('AFSS_block',
                                                   'NIMBLE-AFSS-block',
                                              comparison=ss.ms.AFSblocking[[1]])

ss.ms.filter.AFSblocking[[1]] <- rename_MCMC_comparison_method('AFSS_block',
                                                   'NIMBLE-rm-latent-AFSS-block',
                                              comparison=ss.ms.filter.AFSblocking[[1]])

ss.ms.RWblocking[[1]] <- rename_MCMC_comparison_method('RW_block',
                                                   'NIMBLE-RW-block',
                                              comparison=ss.ms.RWblocking[[1]])

ss.ms.filter.RWblocking[[1]] <- rename_MCMC_comparison_method('RW_block',
                                                   'NIMBLE-rm-latent-RW-block',
                                              comparison=ss.ms.filter.RWblocking[[1]])


ss.ms.slice[[1]] <- rename_MCMC_comparison_method('slice_jags',
                                                   'JAGS-like-NIMBLE',
                                              comparison=ss.ms.slice[[1]])



## ss.ms.subsamp[[1]] <- rename_MCMC_comparison_method('nimbleSubsamp',
##                                                     'NIMBLE-subsample',
##                                               comparison=ss.ms.subsamp[[1]])


## compare mcmcs
ss.ms.occ.all <- combine_MCMC_comparison_results(ss.ms.orig[[1]],
                                                 ss.ms.AFSblocking[[1]],
                                                 ss.ms.RWblocking[[1]],
                                                 ss.ms.slice[[1]],
                                                 ss.ms.filter[[1]],
                                                 ## ss.ms.filter.AFSblocking[[1]],
                                                 ss.ms.filter.RWblocking[[1]],
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

