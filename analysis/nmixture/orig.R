rm(list=ls())
setwd("Dropbox/occupancy")
setwd("analysis/nmixture")
source('src/initialize.R')
source("src/setup.R")


## This example is adapted for NIMBLE from the AHM book by Jacob Levine and Perry de Valpine
## 6.11.1 Bayesian fitting of the basic ZIP N-mixture model
## ------------------------------------------------------------------------
## From section 6.9.2 - for data creation/organization:

## Specify model in BUGS language
## "ZIPNmix.txt"
nmixture <- nimbleCode( {

    ## Specify priors
    ## zero-inflation/suitability
    phi ~ dunif(0,1)          ## proportion of suitable sites (probability of being not a structural 0)
    theta <- 1-phi            ## zero-inflation (proportion of unsuitable)
    ltheta <- logit(theta)

    ## abundance
    beta0 ~ dnorm(0, 0.1)     ## log(lambda) intercept
    for(k in 1:7){            ## Regression params in lambda
        beta[k] ~ dnorm(0, 1)
    }
    tau.lam <- pow(sd.lam, -2)
    sd.lam ~ dunif(0, 2)      ## site heterogeneity in lambda

    ## detection
    for(j in 1:3){
        alpha0[j] <- logit(mean.p[j])
        mean.p[j] ~ dunif(0, 1) ## p intercept for occasions 1-3
    }
    for(k in 1:13){           ## Regression params in p
        alpha[k] ~ dnorm(0, 1)
    }
    tau.p.site <- pow(sd.p.site, -2)
    sd.p.site ~ dunif(0, 2)   ## site heterogeneity in p
    tau.p.survey <- pow(sd.p.survey, -2)
    sd.p.survey ~ dunif(0, 2) ## site-survey heterogeneity in p

    ## ZIP model for abundance
    for (i in 1:nsite){
        a[i] ~ dbern(phi)
        eps.lam[i] ~ dnorm(0, tau.lam)       ## Random site effects in log(abundance)
        loglam[i] <- beta0 + inprod(beta[1:7], lamDM[i, 1:7]) + eps.lam[i]
        loglam.lim[i] <- min(250, max(-250, loglam[i]))  ## Stabilize log
        lam[i] <- exp(loglam.lim[i])
        mu.poisson[i] <- a[i] * lam[i]
        N[i] ~ dpois(mu.poisson[i])
    }

    ## Measurement error model
    for (i in 1:nsite){
        eps.p.site[i] ~ dnorm(0, tau.p.site) ## Random site effects in logit(p)
        for (j in 1:nrep){
            y[i,j] ~ dbin(p[i,j], N[i])
            p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))
            lp.lim[i,j] <- min(250, max(-250, lp[i,j]))  ## Stabilize logit
            lp[i,j] <- alpha0[j] + alpha[1] * elev[i] + alpha[2] * elev2[i] +
                alpha[3] * date[i,j] + alpha[4] * date2[i,j] +
                alpha[5] * dur[i,j] + alpha[6] * dur2[i,j] +
                alpha[7] * elev[i] * date[i,j] + alpha[8] * elev2[i] * date[i,j] +
                alpha[9] * elev[i] * dur[i,j] + alpha[10] * elev[i] * dur2[i,j] +
                alpha[11] * elev2[i] * dur[i,j] + alpha[12] * date[i,j] * dur[i,j] +
                alpha[13] * date[i,j] * dur2[i,j] +
                eps.p.site[i] + eps.p.survey[i,j]
            eps.p.survey[i,j] ~ dnorm(0, tau.p.survey) ## Random site-survey effects
        }
    }
}
)



input1 <- c(code=nmixture,
            model.input)


## *********************************************************************
## original: vanilla nimble and JAGS
## *********************************************************************

nmixture.orig <- compareMCMCs(input1,
                           MCMCs=c('jags', 'nimble'),
                           niter=niter,
                           burnin = burnin,
                           summary=FALSE,
                           check=FALSE)





## Step 5: rewrite the model to group shared calculations
Section6p11_code_grouped <- nimbleCode( {

    ## Specify priors
    ## zero-inflation/suitability
    phi ~ dunif(0,1)          ## proportion of suitable sites (probability of being not a structural 0)
    theta <- 1-phi            ## zero-inflation (proportion of unsuitable)
    ltheta <- logit(theta)

    ## abundance
    beta0 ~ dnorm(0, 0.1)     ## log(lambda) intercept
    for(k in 1:7){            ## Regression params in lambda
        beta[k] ~ dnorm(0, 1)
    }
    tau.lam <- pow(sd.lam, -2)
    sd.lam ~ dunif(0, 2)      ## site heterogeneity in lambda

    ## detection
    for(j in 1:3){
        alpha0[j] <- logit(mean.p[j])
        mean.p[j] ~ dunif(0, 1) ## p intercept for occasions 1-3
    }
    for(k in 1:13){           ## Regression params in p
        alpha[k] ~ dnorm(0, 1)
    }
    tau.p.site <- pow(sd.p.site, -2)
    sd.p.site ~ dunif(0, 2)   ## site heterogeneity in p
    tau.p.survey <- pow(sd.p.survey, -2)
    sd.p.survey ~ dunif(0, 2) ## site-survey heterogeneity in p

    ## ZIP model for abundance
    ## new
    loglam.fixed[1:nsite] <- beta0 + (lamDM[1:nsite, 1:7] %*% beta[1:7])[1:nsite,1]
    for (i in 1:nsite){
        a[i] ~ dbern(phi)
        eps.lam[i] ~ dnorm(0, tau.lam)       ## Random site effects in log(abundance)

        ## original
        ## loglam[i] <- beta0 + inprod(beta[1:7], lamDM[i, 1:7]) + eps.lam[i] * hlam.on
        ## new
        loglam[i] <- loglam.fixed[i] + eps.lam[i] * hlam.on

        loglam.lim[i] <- min(250, max(-250, loglam[i]))  ## Stabilize log
        lam[i] <- exp(loglam.lim[i])
        mu.poisson[i] <- a[i] * lam[i]
        N[i] ~ dpois(mu.poisson[i])
    }

    ## Measurement error model
    ## new
    for(j in 1:nrep) {
        lp.fixed[1:nsite,j] <- alpha0[j] + alpha[1] * elev[1:nsite] + alpha[2] * elev2[1:nsite] +
            alpha[3] * date[1:nsite,j] + alpha[4] * date2[1:nsite,j] +
            alpha[5] * dur[1:nsite,j] + alpha[6] * dur2[1:nsite,j] +
            alpha[7] * elev[1:nsite] * date[1:nsite,j] + alpha[8] * elev2[1:nsite] * date[1:nsite,j] +
            alpha[9] * elev[1:nsite] * dur[1:nsite,j] + alpha[10] * elev[1:nsite] * dur2[1:nsite,j] +
            alpha[11] * elev2[1:nsite] * dur[1:nsite,j] + alpha[12] * date[1:nsite,j] * dur[1:nsite,j] +
            alpha[13] * date[1:nsite,j] * dur2[1:nsite,j]
    }

    for (i in 1:nsite){
        eps.p.site[i] ~ dnorm(0, tau.p.site) ## Random site effects in logit(p)
        for (j in 1:nrep){
            y[i,j] ~ dbin(p[i,j], N[i])
            p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))
            lp.lim[i,j] <- min(250, max(-250, lp[i,j]))  ## Stabilize logit
            ## original
            ## lp[i,j] <- alpha0[j] + alpha[1] * elev[i] + alpha[2] * elev2[i] +
            ##   alpha[3] * date[i,j] + alpha[4] * date2[i,j] +
            ##   alpha[5] * dur[i,j] + alpha[6] * dur2[i,j] +
            ##   alpha[7] * elev[i] * date[i,j] + alpha[8] * elev2[i] * date[i,j] +
            ##   alpha[9] * elev[i] * dur[i,j] + alpha[10] * elev[i] * dur2[i,j] +
            ##   alpha[11] * elev2[i] * dur[i,j] + alpha[12] * date[i,j] * dur[i,j] +
            ##   alpha[13] * date[i,j] * dur2[i,j] +
            ##   eps.p.site[i] * hp.site.on + eps.p.survey[i,j] * hp.survey.on
            ## new
            lp[i,j] <- lp.fixed[i, j] +
                eps.p.site[i] * hp.site.on + eps.p.survey[i,j] * hp.survey.on

            eps.p.survey[i,j] ~ dnorm(0, tau.p.survey) ## Random site-survey effects
        }
    }
    if(DO_POSTERIOR_PREDICTION) {
        ## Posterior predictive distributions of chi2 discrepancy
        for (i in 1:nsite) {
            for (j in 1:nrep) {
                y.sim[i,j] ~ dbin(p[i,j], N[i]) ## Create new data set under model
                e.count[i,j] <- N[i] * p[i,j]   ## Expected datum
                ## Chi-square discrepancy for the actual data
                chi2.actual[i,j] <- pow((y[i,j]-e.count[i,j]),2) / (e.count[i,j]+e)
                ## Chi-square discrepancy for the simulated ('perfect') data
                chi2.sim[i,j] <- pow((y.sim[i,j]-e.count[i,j]),2) / (e.count[i,j]+e)
                ## Add small value e to denominator to avoid division by zero
            }
        }
        ## Add up individual chi2 values for overall fit statistic
        fit.actual <- sum(chi2.actual[1:263, 1:3])  ## Fit statistic for actual data set
        fit.sim <- sum(chi2.sim[1:263, 1:3])        ## Fit statistic for a fitting model
        bpv <- step(fit.sim-fit.actual)    ## Bayesian p-value
        c.hat <- fit.actual/fit.sim        ## c-hat estimate

        ## Derived parameters: Total abundance at 263 sampled sites
        Ntotal263 <- sum(N[1:263])
    }
}
)

## m3 <- nimbleModel(Section6p11_code_grouped,
##                  constants = win.data1,
##                  inits = initsFull())

## MCMCconf3 <- configureMCMC(m3)
## ## Sampling effort is being wasted on the nodes that are "off" in the model:
## MCMCconf3$printSamplers("eps.lam")
## ## So let's remove those samplers
## MCMCconf3$removeSamplers("eps.lam")
## MCMCconf3$removeSamplers("eps.p.site")
## MCMCconf3$removeSamplers("eps.p.survey")
## MCMCconf3$printSamplers()
## MCMC3 <- buildMCMC(MCMCconf3)
## compiled3 <- compileNimble(m3, MCMC3)
## system.time(compiled3$MCMC$run(20)) ## get the NAs filled in before timing
## t200 <- system.time(compiled3$MCMC$run(200))
## t200 * 50000 / 200
## ## a little under 3 minutes

## Step 6: Try some different sampler configurations

## removeUnusedSamplers <- function(MCMCconf) {
##     MCMCconf$removeSamplers("eps.lam")
##     MCMCconf$removeSamplers("eps.p.site")
##     MCMCconf$removeSamplers("eps.p.survey")
##     MCMCconf$removeSamplers("sd.lam")
##     MCMCconf$removeSamplers("sd.p.site")
##     MCMCconf$removeSamplers("sd.p.survey")
##     MCMCconf
## }

## assignAFSS <- function(MCMCconf) {
##     MCMCconf$removeSamplers(c("beta0","beta"))
##     MCMCconf$addSampler(c("beta0","beta"), type = "AF_slice")
##     MCMCconf$removeSamplers(c("mean.p", "alpha"))
##     MCMCconf$addSampler(c("mean.p", "alpha"), type = "AF_slice")
##     MCMCconf
## }

## assignRWB <- function(MCMCconf) {
##     MCMCconf$removeSamplers(c("beta0","beta"))
##     for(i in 1:5)
##         MCMCconf$addSampler(c("beta0","beta"), type = "AF_slice")
##     MCMCconf$removeSamplers(c("mean.p", "alpha"))
##     for(i in 1:5)
##         MCMCconf$addSampler(c("mean.p", "alpha"), type = "AF_slice")
##     MCMCconf
## }


## ## m3 <- nimbleModel(Section6p11_code2,
## ##                  constants = win.data1,
## ##                  inits = initsFull())

## ## MCMCconf <- assignAFSS(removeUnusedSamplers(configureMCMC(m3, useConjugacy = FALSE)))
## ## MCMC <- buildMCMC(MCMCconf)
## ## MCMC$run(2)

## ## compiled <- compileNimble(m3, MCMC)
## ## compiled$MCMC$run(1000, reset = FALSE)


## MCMCdefs = list(
##     default = quote({
##         removeUnusedSamplers(configureMCMC(Rmodel))
##     }),
##     AFSS = quote({
##         assignAFSS(removeUnusedSamplers(configureMCMC(Rmodel)))
##     }),
##     RWB5 = quote({
##         assignRWB(removeUnusedSamplers(configureMCMC(Rmodel)))
##     })
## )

## nb = 1000
## ni = 10000

## ZIP_Nmix_comparisons <- compareMCMCs(
##     modelInfo = list(
##         code = Section6p11_code2,
##         data = win.data1,
##         inits = initsFull()
##         ),
##     monitors = params,
##     MCMCdefs = MCMCdefs,
##     MCMCs = c('default','AFSS', 'RWB5'),
##     summary = TRUE,
##     burnin = nb,
##     niter = ni
## )

## make_MCMC_comparison_pages(ZIP_Nmix_comparisons,
##                            "ZIP_Nmix_comparisons",
##                            modelNames = "ZIP_Nmix")

## browseURL(file.path("ZIP_Nmix_comparisons","ZIP_Nmix.html"))

Section6p11_code_jags_compatible <- nimbleCode( {

    ## Specify priors
    ## zero-inflation/suitability
    phi ~ dunif(0,1)          ## proportion of suitable sites (probability of being not a structural 0)
    theta <- 1-phi            ## zero-inflation (proportion of unsuitable)
    ltheta <- logit(theta)

    ## abundance
    beta0 ~ dnorm(0, 0.1)     ## log(lambda) intercept
    for(k in 1:7){            ## Regression params in lambda
        beta[k] ~ dnorm(0, 1)
    }
    ##    tau.lam <- pow(sd.lam, -2)
    ##    sd.lam ~ dunif(0, 2)      ## site heterogeneity in lambda

    ## detection
    for(j in 1:3){
        alpha0[j] <- logit(mean.p[j])
        mean.p[j] ~ dunif(0, 1) ## p intercept for occasions 1-3
    }
    for(k in 1:13){           ## Regression params in p
        alpha[k] ~ dnorm(0, 1)
    }
    ##    tau.p.site <- pow(sd.p.site, -2)
    ##    sd.p.site ~ dunif(0, 2)   ## site heterogeneity in p
    ##    tau.p.survey <- pow(sd.p.survey, -2)
    ##    sd.p.survey ~ dunif(0, 2) ## site-survey heterogeneity in p

    ## ZIP model for abundance
    for (i in 1:nsite){
        a[i] ~ dbern(phi)
        ##  eps.lam[i] ~ dnorm(0, tau.lam)       ## Random site effects in log(abundance)
        loglam[i] <- beta0 + inprod(beta[1:7], lamDM[i, 1:7]) ##+ eps.lam[i] * hlam.on
        loglam.lim[i] <- min(250, max(-250, loglam[i]))  ## Stabilize log
        lam[i] <- exp(loglam.lim[i])
        mu.poisson[i] <- a[i] * lam[i]
        N[i] ~ dpois(mu.poisson[i])
    }

    ## Measurement error model
    for (i in 1:nsite){
        ## eps.p.site[i] ~ dnorm(0, tau.p.site) ## Random site effects in logit(p)
        for (j in 1:nrep){
            y[i,j] ~ dbin(p[i,j], N[i])
            p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))
            lp.lim[i,j] <- min(250, max(-250, lp[i,j]))  ## Stabilize logit
            lp[i,j] <- alpha0[j] + alpha[1] * elev[i] + alpha[2] * elev2[i] +
                alpha[3] * date[i,j] + alpha[4] * date2[i,j] +
                alpha[5] * dur[i,j] + alpha[6] * dur2[i,j] +
                alpha[7] * elev[i] * date[i,j] + alpha[8] * elev2[i] * date[i,j] +
                alpha[9] * elev[i] * dur[i,j] + alpha[10] * elev[i] * dur2[i,j] +
                alpha[11] * elev2[i] * dur[i,j] + alpha[12] * date[i,j] * dur[i,j] +
                alpha[13] * date[i,j] * dur2[i,j] ##+
            ##eps.p.site[i] * hp.site.on + eps.p.survey[i,j] * hp.survey.on
            ##   eps.p.survey[i,j] ~ dnorm(0, tau.p.survey) ## Random site-survey effects
        }
    }
}
)
