

dDynamicOccupancy <- nimbleFunction(
    ## DynamicOccupancy removes the z's and muZ's from the model and computes
    ## the probability of all reps over all years for one site.
    run = function(x = double(2),
                   nrep = double(),
                   psi1 = double(),
                   phi = double(1),
                   gamma = double(1),
                   p = double(1),
                   log = double(0, default = 0)) {

        numObs <- sum(x[,1])
        ## prob of the occupied sites out of the total sites given p
        ProbOccAndCount <- psi1 * dbinom(numObs, size = nrep, p = p[1], log = 0)
        ## prob of the empty sites
        ProbUnoccAndCount <- (1-psi1) * (numObs == 0)
        ## probably of the observed states
        ProbCount <- ProbOccAndCount + ProbUnoccAndCount
        ProbOccGivenCount <- ProbOccAndCount / ProbCount
        ProbOccNextTime <- ProbOccGivenCount * phi[1] +
            (1-ProbOccGivenCount) * gamma[1]
        ll <- log(ProbCount)
        nyears <- dim(x)[2]
        for(t in 2:nyears) {
            numObs <- sum(x[,t])
            ## similar change to this dbinom
            ProbOccAndCount <- ProbOccNextTime *
                dbinom(numObs, size = nrep, p = p[t], log = 0)
            ProbUnoccAndCount <- (1-ProbOccNextTime) * (numObs == 0)
            ProbCount <- ProbOccAndCount + ProbUnoccAndCount
            ProbOccGivenCount <- ProbOccAndCount / ProbCount
            ll <- ll + log(ProbCount)
            if(t < nyears) ProbOccNextTime <- ProbOccGivenCount * phi[t] +
                               (1-ProbOccGivenCount) * gamma[t]
        }
        if(log) return(ll)
        else return(exp(ll))
        returnType(double())
    }
)

## make matching changes here (at least to the input arguments)
## and check UserManual for non-needed "r" functions
rDynamicOccupancy <- nimbleFunction(
    run = function(n = double(),
                   nrep = double(),
                   psi1 = double(),
                   phi = double(1),
                   gamma = double(1),
                   p = double(1),
                   log = double(0, default = 0)) {
        nyear <- length(p)
        ans <- matrix(rbinom(nrep*nyear, 1, prob=p),
                      nrow=nrep, ncol=nyear)
        ## setSize(ans, nrep, nyear)
        returnType(double(2))
        return(ans)
    }
)

registerDistributions(list(
    dDynamicOccupancy = list(
        BUGSdist = "dDynamicOccupancy(nrep, psi1, phi, gamma, p)",
        Rdist = "dDynamicOccupancy(nrep, psi1, phi, gamma, p)",
        types = c('value = double(2)',
                  'nrep = double(0)',
                  'psi1 = double()',
                  'phi = double(1)',
                  'gamma = double(1)',
                  'p = double(1)'))
))



## ******************************************************************
## seperate function for when gamma, phi, p does not vary by year

dDynamicOccupancyNoYr <- nimbleFunction(
    ## DynamicOccupancy removes the z's and muZ's from the model and computes
    ## the probability of all reps over all years for one site.
    run = function(x = double(2),
                   nrep = double(),
                   psi1 = double(),
                   phi = double(),
                   gamma = double(),
                   p = double(),
                   log = double(0, default = 0)) {

        numObs <- sum(x[,1])
        ## prob of the occupied sites out of the total sites given p
        ProbOccAndCount <- psi1 * dbinom(numObs, size = nrep, p = p, log = 0)
        ## prob of the empty sites
        ProbUnoccAndCount <- (1-psi1) * (numObs == 0)
        ## probably of the observed states
        ProbCount <- ProbOccAndCount + ProbUnoccAndCount
        ProbOccGivenCount <- ProbOccAndCount / ProbCount
        ProbOccNextTime <- ProbOccGivenCount * phi +
            (1-ProbOccGivenCount) * gamma
        ll <- log(ProbCount)
        nyears <- dim(x)[2]
        for(t in 2:nyears) {
            numObs <- sum(x[,t])
            ## similar change to this dbinom
            ProbOccAndCount <- ProbOccNextTime *
                dbinom(numObs, size = nrep, p = p, log = 0)
            ProbUnoccAndCount <- (1-ProbOccNextTime) * (numObs == 0)
            ProbCount <- ProbOccAndCount + ProbUnoccAndCount
            ProbOccGivenCount <- ProbOccAndCount / ProbCount
            ll <- ll + log(ProbCount)
            if(t < nyears) ProbOccNextTime <- ProbOccGivenCount * phi +
                               (1-ProbOccGivenCount) * gamma
        }
        if(log) return(ll)
        else return(exp(ll))
        returnType(double())
    }
)

## make matching changes here (at least to the input arguments)
## and check UserManual for non-needed "r" functions
## rDynamicOccupancyNoYr <- nimbleFunction(
##     run = function(n = double(),
##                    nrep = double(),
##                    psi1 = double(),
##                    phi = double(),
##                    gamma = double(),
##                    p = double(),
##                    log = double(0, default = 0)) {
##         nyear <- length(p)
##         ans <- matrix(rbinom(nrep*nyear, 1, prob=p),
##                       nrow=nrep, ncol=nyear)
##         ## setSize(ans, nrep, nyear)
##         returnType(double(2))
##         return(ans)
##     }
## )

registerDistributions(list(
    dDynamicOccupancyNoYr = list(
        BUGSdist = "dDynamicOccupancyNoYr(nrep, psi1, phi, gamma, p)",
        Rdist = "dDynamicOccupancyNoYr(nrep, psi1, phi, gamma, p)",
        types = c('value = double(2)',
                  'nrep = double()',
                  'psi1 = double()',
                  'phi = double()',
                  'gamma = double()',
                  'p = double()'))
))



