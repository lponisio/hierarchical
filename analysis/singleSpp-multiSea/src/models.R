makeModel <- function(ffilter, hyper.param){
    if(!ffilter){
        if(hyper.param){
            ## ***************************************************************
            ## no filtering, including hyper param
            ## ***************************************************************
            ss.ms.occ <- nimbleCode({
                ## Specify priors
                psi1 ~ dunif(0, 1)
                mu.p     ~ dnorm(0,0.001)
                sigma.p     ~ dunif(0,100)
                tau.p <- 1/(sigma.p*sigma.p)
                mu.phi  ~ dnorm(0,0.001)
                mu.gamma  ~ dnorm(0,0.001)
                sigma.phi ~ dunif(0,100)
                sigma.gamma ~ dunif(0,100)
                tau.phi <-  1/(sigma.phi*sigma.phi)
                tau.gamma <-  1/(sigma.gamma*sigma.gamma)
                for(year in 1:(nyear -1)) {
                    p[year]  ~ dnorm(mu.p,     tau.p)
                    phi[year] ~ dnorm(mu.phi, tau.phi)
                    gamma[year] ~ dnorm(mu.gamma, tau.gamma)
                }
                p[nyear]      ~ dnorm(mu.p,     tau.p)

                ## Ecological submodel: Define state conditional on parameters
                for (site in 1:nsite){
                    z[site,1] ~ dbern(psi1)
                    for (year in 2:nyear){
                        logit(muZ[site,year]) <- z[site,year-1]*phi[year-1] +
                            (1-z[site,year-1])*gamma[year-1]
                        z[site,year] ~ dbern(muZ[site,year])
                    }
                }

                ## Observation model
                for (site in 1:nsite){
                    for (rep in 1:nrep){
                        for (year in 1:nyear){
                            logit(muy[site,rep,year]) <- z[site,year]*p[year]
                            y[site,rep,year] ~ dbern(muy[site,rep,year])
                        }
                    }
                }
            })

        } else{
            ## ***************************************************************
            ## no filtering, no hyper param
            ## ***************************************************************
            ss.ms.occ <- nimbleCode({
                ## Specify priors
                psi1 ~ dunif(0, 1)
                mu.p     ~ dnorm(0,0.001)
                mu.phi  ~ dnorm(0,0.001)
                mu.gamma  ~ dnorm(0,0.001)

                ## Ecological submodel: Define state conditional on parameters
                for (site in 1:nsite){
                    z[site,1] ~ dbern(psi1)
                    for (year in 2:nyear){
                        logit(muZ[site,year]) <- z[site,year-1]*mu.phi +
                            (1-z[site,year-1])*mu.gamma
                        z[site,year] ~ dbern(muZ[site,year])
                    }
                }
                ## Observation model
                for (site in 1:nsite){
                    for (rep in 1:nrep){
                        for (year in 1:nyear){
                            logit(muy[site,rep,year]) <- z[site,year]*mu.p
                            y[site,rep,year] ~ dbern(muy[site,rep,year])
                        }
                    }
                }
            })
        }
    }

    if(ffilter){
        if(hyper.param){
            ## ***************************************************************
            ## filtering and hyper param
            ## ***************************************************************
            ss.ms.occ <- nimbleCode({
                ##  priors
                psi1 ~ dunif(0, 1)
                mu.p     ~ dnorm(0,0.001)
                sigma.p     ~ dunif(0,100)
                tau.p <- 1/(sigma.p*sigma.p)
                mu.phi  ~ dnorm(0,0.001)
                mu.gamma  ~ dnorm(0,0.001)
                sigma.phi ~ dunif(0,100)
                sigma.gamma ~ dunif(0,100)
                tau.phi <-  1/(sigma.phi*sigma.phi)
                tau.gamma <-  1/(sigma.gamma*sigma.gamma)
                for(year in 1:(nyear -1)) {
                    p[year]   ~ dnorm(mu.p,     tau.p)
                    phi[year] ~ dnorm(mu.phi, tau.phi)
                    gamma[year] ~ dnorm(mu.gamma, tau.gamma)
                }
                p[nyear]      ~ dnorm(mu.p,     tau.p)

                ## Ecological submodel: Define state conditional on parameters
                for(i in 1:nsite) {
                    ## removes the z's and muZ's from the model and compute
                    ## the probability of all reps over all years for one site.
                    y[i, 1:nrep, 1:nyear] ~ dDynamicOccupancy(nrep,
                                                              psi1,
                                                              expit(phi[1:(nyear-1)]),
                                                              expit(gamma[1:(nyear-1)]),
                                                              expit(p[1:nyear]))
                }
            })
        }else{
            ## ***************************************************************
            ## filtering, no hyper param
            ## ***************************************************************
            ss.ms.occ <- nimbleCode({
                ##  priors
                psi1 ~ dunif(0, 1)
                mu.p     ~ dnorm(0,0.001)
                mu.phi  ~ dnorm(0,0.001)
                mu.gamma  ~ dnorm(0,0.001)
                ## Ecological submodel: Define state conditional on parameters
                for(i in 1:nsite) {
                    ## removes the z's and muZ's from the model and compute
                    ## the probability of all reps over all years for one site.
                    y[i, 1:nrep, 1:nyear] ~ dDynamicOccupancyNoYr(nrep,
                                                                  psi1,
                                                                  expit(mu.phi),
                                                                  expit(mu.gamma),
                                                                  expit(mu.p))
                }
            })
        }
    }
    return(ss.ms.occ)
}
