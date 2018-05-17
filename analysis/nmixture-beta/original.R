rm(list=ls())
setwd("~/Dropbox/occupancy")
setwd("analysis/nmixture-beta")
source('src/initialize.R')
source("src/setup.R")
model.input <- prepData(UIF.counts)

## This example is adapted for NIMBLE from the Gomez et al. 2017
## (https://doi.org/10.1111/2041-210X.12856) orginally from Yamaura et
## al. 2016 model
## ------------------------------------------------------------------------

ms.nmixture <- nimbleCode( {
    ## prior distributions on community level estimates

    ## mean value
    ## parameter related to abundance
    mu.a0 ~ dnorm(0,0.001)	## intercept
    ## parameter related to detectability
    mu.p0 ~ dnorm(0,0.001)

    ## standard deviation
    ## parameter related to abundance
    sigma.a0 ~ dunif(0,10)	## intercept

    ## parameter related to detectability
    sigma.p0 ~ dunif(0,10)

    ## for jags compatibility
    tau.a0 <- pow(sigma.a0,-2)
    tau.p0 <- pow(sigma.p0,-2)

        ## Estimate parameters fo each species
        for(i in 1:nspp){

            ## generating parameters of each species related to abundance
            a0[i] ~ dnorm(mu.a0, tau.a0)

            ## Estimating lambda from N which is a latent variable
            log(lambda[i]) <- a0[i]	## equation (4)

            ## generating parameters of each species related to detectability
            p0[i] ~ dnorm(mu.p0, tau.p0)

            ##values of p for each species (i) and each site (j)
            p[i] <- 1/(1+exp(-(p0[i])))	## equation (5)
            ## Loop along sites (point counts)
            for(j in 1:nsites){

                N[i,j] ~ dpois(lambda[i])	## latent abundance of each species in each year

                ## Loop along replicate counts
                for(t in 1:nvisits){

                    ## detection process model
                    counts[j,t,i] ~ dbin(p[i], N[i,j])	## detection process in each site

                }

            }

        }

    }
)


input1 <- c(code=ms.nmixture,
            model.input)


## *********************************************************************
## original: vanilla nimble and JAGS
## *********************************************************************

ms.nmixture.orig <- compareMCMCs(input1,
                           MCMCs=c('nimble'),
                           niter=niter,
                           burnin = burnin,
                           summary=FALSE,
                           check=FALSE)
