
dNmixture <- nimbleFunction(
    run = function(x = double(0),
                   prob = double(0),
                   lambda = double(0),
                   notZero = double(0),
                   log = integer(0, default = 0)) {
        if(is.finite(x) == 0){
            return(NA)
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
        while(cumPoisProb < poissonTolerance) {
            thisPoisProb <- dpois(N, lambda)
            thisObsProb <- thisPoisProb * dbinom(x, size = N, prob = prob)
            obsProb <- obsProb + thisObsProb
            cumPoisProb <- cumPoisProb + thisPoisProb
            N <- N + 1
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
                  'prob = double(0)',
                  'lambda = double(0)',
                  'notZero = double(0)'))
))



