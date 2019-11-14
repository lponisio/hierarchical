## a quick test

library(nimble)
source("dNmixtureFast.R")

## First case: I used N <- 20, and lambda = 24 in the dNmixture calls
## second case: I used N <- 100 and lambda = 150 inthe dNmixture calls
N <- 100
p <- c(.2, .4, .6, .8)
x <- rbinom(4, size = N, prob = p)

dNmixtureRepOrig(x = x, prob = p, lambda = 24, notZero = 1, lambdaRange = c(0.01, 100), log = TRUE)

dNmixtureRepFast(x = x, prob = p, lambda = 24, notZero = 1, lambdaRange = c(0.01, 100), log = TRUE)

timeOrig <- nimbleFunction(
  run = function(m = integer(),
                 x = double(1),
                 prob = double(1),
                 lambda = double(),
                 notZero = double(),
                 lambdaRange = double(1),
                 log = integer(0, default = 0)) {
    ans <- run.time({
      for(i in 1:m) {
        logProb <- dNmixtureRepOrig(x, prob, lambda, notZero, lambdaRange, log)
      }
    })
    return(ans)
    returnType(double())

  })

Corig <- compileNimble(dNmixtureRepOrig, timeOrig)

Corig$dNmixtureRepOrig(x = x, prob = p, lambda = 150, notZero = 1, lambdaRange = c(0.01, 200), log = TRUE)

Corig$timeOrig(m = 10000, x = x, prob = p, lambda = 150, notZero = 1, lambdaRange = c(0.01, 200), log = TRUE)




timeFast <- nimbleFunction(
  run = function(m = integer(),
                 x = double(1),
                 prob = double(1),
                 lambda = double(),
                 notZero = double(),
                 lambdaRange = double(1),
                 log = integer(0, default = 0)) {
    ans <- run.time({
      for(i in 1:m) {
        logProb <- dNmixtureRepFast(x, prob, lambda, notZero, lambdaRange, log)
      }
    })
    return(ans)
    returnType(double())

  })

Cfast <- compileNimble(dNmixtureRepFast, timeFast)

Cfast$dNmixtureRepFast(x = x, prob = p, lambda = 150, notZero = 1, lambdaRange = c(0.01, 200), log = TRUE)

Cfast$timeFast(m = 10000, x = x, prob = p, lambda = 150, notZero = 1, lambdaRange = c(0.01, 200), log = TRUE)
