
dNmixture <- nimbleFunction(
    run = function(x = double(0),
                   prob = double(),
                   lambda = double(),
                   notZero = double(),
                   log = integer(0, default = 0)) {
        if(is.na(x) | is.nan(x)) {
            if(log) return(-Inf)
            return (0)
        }
        if(notZero == 0) { ## It is structural zero
            if(x == 0) {
                if(log) return(0)
                else return(1)
            } else {
                if(log) return(-Inf)
                else return(0)
            }
        }
        ## Then all lambda < 0.01 will never be accepted.
        if(lambda < 0.01) {
            if(log) return(-Inf)
            else return(0)
        }
        ## It is not a structural zero
        obsProb <- 0
        if(x > 0){
            cumPoisProb <- ppois(x-1, lambda)
        } else {
            cumPoisProb <- 0
        }
        poissonTolerance <- 0.9999
        N <- x
        iterations <- 0
        while(cumPoisProb < poissonTolerance) {
            thisPoisProb <- dpois(N, lambda)
            thisObsProb <- thisPoisProb * dbinom(x, size = N, prob = prob)
            obsProb <- obsProb + thisObsProb
            cumPoisProb <- cumPoisProb + thisPoisProb
            N <- N + 1
            iterations <- iterations + 1
            if(iterations > 5000) {
                fractional <- (iterations - 5000)/1000
                if(abs(fractional - round(fractional)) < 0.00001)
                    print("On iteration ", iterations,
                          " with lambda = ",
                          lambda, " x = ", x, " cumPoisProb = ", cumPoisProb)
            }
        }
        if(log) return(log(obsProb))
        else return(obsProb)
        returnType(double(0))
    }
)


registerDistributions(list(
    dNmixture = list(
        BUGSdist = "dNmixture(prob, lambda, notZero)",
        Rdist = "dNmixture(prob, lambda, notZero)",
        types = c('value = double(0)',
                  'prob = double()',
                  'lambda = double()',
                  'notZero = double()'))
))

dNmixtureRep <- nimbleFunction(
    run = function(x = double(1),
                   prob = double(1),
                   lambda = double(),
                   notZero = double(),
                   lambdaRange = double(1),
                   log = integer(0, default = 0)) {
        ## Reject lambda < 0.01:
        if(lambda < lambdaRange[1] | lambda > lambdaRange[2]) {
            if(log) return(-Inf)
            else return(0)
        }
        if(is.na.vec(x) | is.nan.vec(x)) {
            if(log) return(-Inf)
            return (0)
        }
        if(notZero == 0) { ## It is structural zero
            if(all(x == 0)) {
                if(log) return(0)
                else return(1)
            } else {
                if(log) return(-Inf)
                else return(0)
            }
        }
        ## It is not a structural zero.
        ##
        ## For each x, the conditional distribution of (N - x | x) is pois(lambda * (1-p))
        ## We determine the lowest N and highest N at extreme quantiles and sum over those.
        minN <- min(x + qpois(0.00001, lambda * (1-prob)))
        maxN <- max(x + qpois(0.99999, lambda * (1-prob)))
        minN <- max( max(x), minN ) ## set minN to at least the largest x

        obsProb <- 0
        if(maxN > minN) { ## should normally be true, but check in case it isn't in some corner case.
        ##    print("counting from ", minN, " to ", maxN, " with lambda = ", lambda)
            for(N in minN:maxN) {
                thisObsProb <- dpois(N, lambda) * prod(dbinom(x, size = N, prob = prob))
                obsProb <- obsProb + thisObsProb
            }
        } else {
            ## return a potentially non-zero obsProb
            N <- max(x)
            obsProb <- dpois(N, lambda) * prod(dbinom(x, size = N, prob = prob))
        }
        if(log) return(log(obsProb))
        else return(obsProb)
        returnType(double(0))
    }
)



registerDistributions(list(
    dNmixtureRep = list(
        BUGSdist = "dNmixtureRep(prob, lambda, notZero, lambdaRange)",
        Rdist = "dNmixtureRep(prob, lambda, notZero, lambdaRange)",
        types = c('value = double(1)',
                  'prob = double(1)',
                  'lambda = double()',
                  'notZero = double()',
                  'lambdaRange = double(1)'))
))
