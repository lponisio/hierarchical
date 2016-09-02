## *********************************************************************
## vectorize Bernoulli draws option 1
## *********************************************************************
## dbern_vec is a user-defined distribution that computes all the
## dbern(psi)'s in one loop

## accepts a scalar p

## we need to take length as an argument because for rbern_vec, there
## is other way to determine the length (it is not the intended use of
## n)

## and the dbern_vec needs to have the save parameter arguments as
## rbern_vec

dbern_vec <- nimbleFunction(
  run = function(x = double(1),
    p = double(),
    length = double(),
    log = integer(0, default = 0)) {
    if(length != dim(x)[1]) stop('In dbern_vec, length != length of x')

    ## For the scalar p case, one could also replace the following 6
    ## lines with

    ##sumx <- sum(x)
    ##ans <- log(p) * sumx + log(1-p) * (length-sumx)

    ans <- 0
    logp <- log(p)
    log1mp <- log(1-p)
    for(i in 1:length)
      if(x[i] == 0) ans <- ans + log1mp
      else ans <- ans + logp

    returnType(double())
    if(log) return(ans)
    else return(exp(ans))
  })

## note the n is really a dummy variable for future extensions.  we
## always treat it as 1.  It should not be used as the length needed
rbern_vec <- nimbleFunction(
  run = function(n = integer(), p = double(), length = double()) {
    ans <- double(1)
    setSize(ans, length)
    for(i in 1:length) ans[i] <- rbinom(1, size = 1, prob = p)
    returnType(double(1))
    return(ans)
  })

registerDistributions(list(dbern_vec = list(
                             BUGSdist = "dbern_vec(p, length)",
                             Rdist = "dbern_vec(p, length)",
                             range = c(0, 1),
                             types = c('value = double(1)',
                               'p = double(0)',
                               'length = double(0)'))
                           ))

## In this case, p provides information on the sizes of x, so we don't
## need separate arguments to provide that information

dbern_matrix <- nimbleFunction(
  run = function(x = double(2),
    p = double(2),
    log = integer(0, default = 0)) {
    nrow = dim(p)[1]
    ncol = dim(p)[2]
    if(nrow != dim(x)[1]){
      stop('In dbern_matrix, dim[1] != number of rows in x')
    }
    if(ncol != dim(x)[1]){
      stop('In dbern_matrix, dim[1] != number of cols in x')
    }
    ans <- 0
    for(i in 1:nrow)
      for(j in 1:ncol)
        if(x[i,j] == 0) ans <- ans + log(1-p[i,j])
        else ans <- ans + log(p[i,j])

    returnType(double())
    if(log) return(ans)
    else return(exp(ans))
  })

rbern_matrix <- nimbleFunction(
  run = function(n = integer(), p = double(2)) {
    ans = double(2)
    nrow = p[1]
    ncol = p[2]
    setSize(ans, nrow, ncol)
    for(i in 1:nrow)
      for(j in 1:ncol)
        ans[i,j] <- rbinom(1, size = 1, prob = p)
    returnType(double(2))
    return(ans)
  })

registerDistributions(list(dbern_matrix = list(
                             BUGSdist = "dbern_matrix(p)",
                             Rdist = "dbern_matrix(p)",
                             range = c(0, 1),
                             types = c('value = double(2)',
                               'p = double(2)'))
                           ))


## *********************************************************************
## custom sampler for zs, option 2
## *********************************************************************
custom_z_sampler <- nimbleFunction(
  contains = sampler_BASE,      ## required (class inheritance)
  setup = function(model, mvSaved, target, control) { ## required arguments
    calcNodes <- model$getDependencies(target)
    ## to be used for logProbs for values of 0 and 1 respectively
    logProbs <- c(0.0, 0.0)
  },
  run = function() {
    ## get the current value
    currentValue <- model[[target]]
    ## and the related sum of logProb's
    currentLogProb <- getLogProb(model, calcNodes)
    for(i in 1:2) {   ## check both possible values
      ## for the current value, record the current sum of logProb
      if(i - 1 == currentValue)
        logProbs[i] <<- currentLogProb
      else {
        ## for the other value, put it in the model
        model[[target]] <<- i - 1
        ## and calculate the sum of logProb's
        logProbs[i] <<- calculate(model, calcNodes)
      }
    }
    ## pick a random number between 0 and 1
    u <- runif(1, 0, 1)
    ## choose 0 with this probability, 1 otherwise
    if(u < exp(logProbs[1])/sum(exp(logProbs[1:2])))
      newValue <- 0
    else
      newValue <- 1

    if(newValue != currentValue) ## update the saved model states
      copy(from = model,
           to = mvSaved, row = 1,
           nodes = calcNodes, logProb = TRUE)
    else
      copy(from = mvSaved,
           to = model, row = 1,
           nodes = calcNodes, logProb = TRUE)
  },
  ## required for the MCMC system, but it doesn't have to do anything
  methods = list(
    reset = function () {}
    )
  )

## *********************************************************************
## relected custom z sampler, option 2b
## *********************************************************************

sampler_RW_reflect <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ###  control list extraction  ###
        adaptive      <- control$adaptive
        adaptInterval <- control$adaptInterval
        scale         <- control$scale
        ###  node list generation  ###
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        if(length(targetAsScalar) > 1)     stop('more than one target; cannot use RW sampler, try RW_block sampler')
        calcNodes  <- model$getDependencies(target)
        ###  numeric value generation  ###
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        scaleHistory          <- c(0, 0)
        acceptanceRateHistory <- c(0, 0)
        ## variables previously inside of nested functions:
        optimalAR <- 0.44
        gamma1    <- 0
        ## adapted from Chris' user-defined sampler training module
        dist <- model$getNodeDistribution(target)
        scalar <- getDistribution(dist)$types$value$nDim == 0
        if(!scalar) stop('Cannot use reflected sampler for non-scalar target')
        rg <- getDistribution(dist)$range
        reflectMin <- rg[1]
        reflectMax <- rg[2]
    },

    run = function() {
        propValue <- rnorm(1, mean = model[[target]], sd = scale)

        while(propValue < reflectMin | propValue > reflectMax) {
            if(propValue < reflectMin) propValue <- 2*reflectMin - propValue
            if(propValue > reflectMax) propValue <- 2*reflectMax - propValue
        }

     	model[[target]] <<- propValue
        logMHR <- calculateDiff(model, calcNodes)
        jump <- decide(logMHR)
        if(jump)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        else
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
        if(adaptive)     adaptiveProcedure(jump)
    },

    methods = list(

        adaptiveProcedure = function(jump = logical()) {
            timesRan <<- timesRan + 1
            if(jump)     timesAccepted <<- timesAccepted + 1
            if(timesRan %% adaptInterval == 0) {
                acceptanceRate <- timesAccepted / timesRan
                timesAdapted <<- timesAdapted + 1
                setSize(scaleHistory,          timesAdapted)
                setSize(acceptanceRateHistory, timesAdapted)
                scaleHistory[timesAdapted] <<- scale
                acceptanceRateHistory[timesAdapted] <<- acceptanceRate
                gamma1 <<- 1/((timesAdapted + 3)^0.8)
                gamma2 <- 10 * gamma1
                adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
                scale <<- scale * adaptFactor
                timesRan <<- 0
                timesAccepted <<- 0
            }
        },

        reset = function() {
            scale <<- scaleOriginal
            timesRan      <<- 0
            timesAccepted <<- 0
            timesAdapted  <<- 0
            scaleHistory          <<- scaleHistory          * 0
            acceptanceRateHistory <<- acceptanceRateHistory * 0
            gamma1 <<- 0
        }
    ), where = getLoadingNamespace()
)

## *********************************************************************
## remove Zs, option 3
## *********************************************************************
## DynamicOccupancy removes the z's and muZ's from the model and computes
## the probability of all reps over all years for one site.
dDynamicOccupancy <- nimbleFunction(
  ## I've checked that this runs and compiles, but I haven't tested if
  ## I got the logic right!
  run = function(x = double(2),
    nrep = integer(),
    psi1 = double(),
    phi = double(1),
    gamma = double(1),
    p = double(1),
    log = integer(0, default = 0)) {
    prob1 <- psi1 * p[1]
    numObs <- sum(x[,1]) ## do I have the right orientation?
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

rDynamicOccupancy <- nimbleFunction(
  run = function(n = integer(),
    nrep = integer(),
    psi1 = double(),
    phi = double(1),
    gamma = double(1),
    p = double(1),
    log = integer(0, default = 0)) {
    nyear <- length(p)
    ans <- double(2)
    setSize(ans, nrep, nyear)
    ## could populate ans here, but I'm just doing this as a placeholder
    returnType(double(2))
    return(ans)
  }
  )

registerDistributions(list(
  dDynamicOccupancy = list(
    BUGSdist = "dDynamicOccupancy(nrep, psi1, phi, gamma, p)",
    Rdist = "dDynamicOccupancy(nrep, psi1, phi, gamma, p)",
    types = c('value = double(2)',
      'nrep = integer(0)',
      'psi1 = double()',
      'phi = double(1)',
      'gamma = double(1)',
      'p = double(1)'))
  ))
