

dNmixtureRep <- nimbleFunction(
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
      logProb <- dpois(minN, lambda, log = TRUE) +
          sum(dbinom(x, size = minN, prob = prob, log = TRUE)) + log(fac)
    }
    if(log) return(logProb)
    else return(exp(logProb))
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
