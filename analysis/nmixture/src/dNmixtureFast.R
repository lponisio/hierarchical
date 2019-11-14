## For uncompiled, change is.na.vec to nimble:::is.na.vec.  just filed an issue on this

dNmixtureRepOrig <- nimbleFunction(
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


dNmixtureRepFast <- nimbleFunction(
  run = function(x = double(1),
                 prob = double(1),
                 lambda = double(),
                 notZero = double(),
                 lambdaRange = double(1),
                 log = integer(0, default = 0)) {
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
    minN <- min(x + qpois(0.00001, lambda * (1-prob)))
    maxN <- max(x + qpois(0.99999, lambda * (1-prob)))
    minN <- max( max(x), minN ) ## set minN to at least the largest x
    numReps <- length(x)
    logProb <- -Inf
    if(maxN > minN) {
      fac <- 1
      ## ff = lambda prod((1-p_i)) N^(numReps - 1) 
      ff <- lambda * prod( (1-prob) )
      numN <- maxN - minN + 1 - 1 ## remember: +1 for the count, but -1 because the summation should run from N = maxN to N = minN + 1
      for(i in 1:numN) {
        N <- maxN - i + 1
        fac <- 1 + fac * ff * prod(N/(N - x)) / N
      }
      logProb <- dpois(minN, lambda, log = TRUE) + sum(dbinom(x, size = minN, prob = prob, log = TRUE)) + log(fac)
    }
    if(log) return(logProb)
    else return(exp(logProb))
    returnType(double(0))
  }
)
