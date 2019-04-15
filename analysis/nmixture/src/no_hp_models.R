## latent, no hyper param
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

    ## detection
    for(j in 1:3){
        alpha0[j] <- logit(mean.p[j])
        mean.p[j] ~ dunif(0, 1) ## p intercept for occasions 1-3
    }
    for(k in 1:13){           ## Regression params in p
        alpha[k] ~ dnorm(0, 1)
    }

    ## ZIP model for abundance
    for (i in 1:nsite){
        a[i] ~ dbern(phi)
        loglam[i] <- beta0 + inprod(beta[1:7], lamDM[i, 1:7])
        loglam.lim[i] <- min(250, max(-250, loglam[i]))  ## Stabilize log
        lam[i] <- exp(loglam.lim[i])
        mu.poisson[i] <- a[i] * lam[i]
        N[i] ~ dpois(mu.poisson[i])
    }

    ## Measurement error model
    for (i in 1:nsite){
        for (j in 1:nrep[i]){
            y[i,j] ~ dbin(p[i,j], N[i])
            p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))
            lp.lim[i,j] <- min(250, max(-250, lp[i,j]))  ## Stabilize logit
            lp[i,j] <- alpha0[j] + alpha[1] * elev[i] + alpha[2] * elev2[i] +
                alpha[3] * date[i,j] + alpha[4] * date2[i,j] +
                alpha[5] * dur[i,j] + alpha[6] * dur2[i,j] +
                alpha[7] * elev[i] * date[i,j] + alpha[8] * elev2[i] * date[i,j] +
                alpha[9] * elev[i] * dur[i,j] + alpha[10] * elev[i] * dur2[i,j] +
                alpha[11] * elev2[i] * dur[i,j] + alpha[12] * date[i,j] * dur[i,j] +
                alpha[13] * date[i,j] * dur2[i,j]
        }
    }
}
)


## no latent, no hyper param

nmixture <- nimbleCode( {
    lambdaRange[1] <- 0.01
    lambdaRange[2] <- 10000
    ## Specify priors
    ## zero-inflation/suitability
    phi ~ dunif(0,1)          ## proportion of suitable
    ## sites (probability of
    ## being not a structural 0)
    theta <- 1-phi            ## zero-inflation (proportion of unsuitable)
    ltheta <- logit(theta)

    ## abundance
    beta0 ~ dnorm(0, 0.1)     ## log(lambda) intercept
    for(k in 1:7){            ## Regression params in lambda
        beta[k] ~ dnorm(0, 1)
    }

    ## detection
    for(j in 1:3){
        alpha0[j] <- logit(mean.p[j])
        mean.p[j] ~ dunif(0, 1) ## p intercept for occasions 1-3
    }
    for(k in 1:13){           ## Regression params in p
        alpha[k] ~ dnorm(0, 1)
    }

    ## ZIP model for abundance
    for (i in 1:nsite){
        a[i] ~ dbern(phi)

        loglam[i] <- beta0 + inprod(beta[1:7], lamDM[i, 1:7])
        loglam.lim[i] <- min(250, max(-250, loglam[i]))  ## Stabilize log
        lam[i] <- exp(loglam.lim[i])
        mu.poisson[i] <- a[i] * lam[i]
    }



    ## Measurement error model
    for (i in 1:nsite){
        ## latent state integration
        y[i,1:nrep[i]] ~ dNmixtureRep(p[i,1:nrep[i]], lam[i], a[i], lambdaRange[1:2])
        for (j in 1:nrep[i]){
            p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))
            lp.lim[i,j] <- min(250, max(-250, lp[i,j]))  ## Stabilize logit
            lp[i,j] <- alpha0[j] + alpha[1] * elev[i] + alpha[2] * elev2[i] +
                alpha[3] * date[i,j] + alpha[4] * date2[i,j] +
                alpha[5] * dur[i,j] + alpha[6] * dur2[i,j] +
                alpha[7] * elev[i] * date[i,j] + alpha[8] * elev2[i] * date[i,j] +
                alpha[9] * elev[i] * dur[i,j] + alpha[10] * elev[i] * dur2[i,j] +
                alpha[11] * elev2[i] * dur[i,j] + alpha[12] * date[i,j] * dur[i,j] +
                alpha[13] * date[i,j] * dur2[i,j]
        }
    }
}
)
}


