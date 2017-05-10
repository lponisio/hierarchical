library(nimble)

source('~//postdoc//afslice//afESS2.R')

pumpCode <- nimbleCode({
    for (i in 1:N){
        theta[i] ~ dgamma(alpha,beta)
        lambda[i] <- theta[i]*t[i]
        x[i] ~ dpois(lambda[i])
    }
    alpha ~ dexp(1.0)
    beta ~ dgamma(0.1,1.0)
})
pumpConsts <- list(N = 10,
                   t = c(94.3, 15.7, 62.9, 126, 5.24,
                         31.4, 1.05, 1.05, 2.1, 10.5))
pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22))

pumpInitial <- list(alpha = 0.8237675, beta = 1.2634812,  theta = rep(0.1, pumpConsts$N))

pumpParams <- c('alpha', 'beta')
pumpModel <- nimbleModel(code = pumpCode, name = 'pump', constants = pumpConsts,
                         data = pumpData, inits = pumpInitial)

compileNimble(pumpModel)
conf <- configureMCMC(pumpModel)
conf$removeSamplers(pumpParams)
conf$addSampler(pumpParams, type = 'AFSS_to_RW_block',
                control = list(AF_sliceControl =  list(sliceWidths = rep(1, length(pumpParams)),
                                                       factorBurnIn = 10000,
                                                       factorAdaptInterval = 1000,
                                                       sliceBurnIn = 1000,
                                                       sliceMaxSteps = 100),
                               RWcontrol = list(propCov = diag(length(pumpParams)), scale = 1,
                                                adaptInterval = 500, adaptScaleOnly = F,
                                                adaptive=T),
                               nAFSSIters
                               =
                               8000,
                               essThreshold = .5,
                               numESSAdaptations = 50,
                               timeSwitch = TRUE))
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc)
print(system.time(cmcmc$run(10000)))

