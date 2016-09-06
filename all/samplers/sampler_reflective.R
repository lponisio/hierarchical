
## *********************************************************************
## relected custom z sampler
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
