
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
                    print("On iteration ", iterations, " with lambda = ", lambda, " x = ", x, " cumPoisProb = ", cumPoisProb)
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
        ## It is not a structural zero
        obsProb <- 0
        N <- max(x)
        poissonTolerance <- 0.9999
        if(N > 0){
            cumPoisProb <- ppois(N-1, lambda)
            if(cumPoisProb >= poissonTolerance) {
                ## If there is virtually 0 probability that N is large enough to yield the observed counts,
                ## set obsProb as if N is known to be the largest count, so that a non-zero probability is returned.
                obsProb <- dpois(N, lambda) * prod(dbinom(x, size = N, prob = prob))
            }
        } else {
            cumPoisProb <- 0
        }
        iterations <- 0
        while(cumPoisProb < poissonTolerance) {
            thisPoisProb <- dpois(N, lambda)
            thisObsProb <- thisPoisProb * prod(dbinom(x, size = N, prob = prob))
            obsProb <- obsProb + thisObsProb
            cumPoisProb <- cumPoisProb + thisPoisProb
            N <- N + 1
            iterations <- iterations + 1
            if(iterations > 1000000) {
                fractional <- iterations/1000000
                if(abs(fractional - round(fractional)) < 0.0000001)
                    print("Warning: dNmixture is on N = ", N, " with lambda = ", lambda, " x = ", x, " cumPoisProb = ", cumPoisProb, ". Is it stuck forever?")
            }
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
