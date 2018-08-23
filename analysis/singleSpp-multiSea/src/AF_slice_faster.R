sampler_AF_slice_faster <- nimbleFunction(
    name = 'sampler_AF_slice_faster',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        widthVec            <- if(!is.null(control$sliceWidths))              control$sliceWidths              else 'oneVec'
        maxSteps            <- if(!is.null(control$sliceMaxSteps))            control$sliceMaxSteps            else 100
        adaptFactorMaxIter  <- if(!is.null(control$sliceAdaptFactorMaxIter))  control$sliceAdaptFactorMaxIter  else 15000
        adaptFactorInterval <- if(!is.null(control$sliceAdaptFactorInterval)) control$sliceAdaptFactorInterval else 1000
        adaptWidthMaxIter   <- if(!is.null(control$sliceAdaptWidthMaxIter))   control$sliceAdaptWidthMaxIter   else 512
        adaptWidthTolerance <- if(!is.null(control$sliceAdaptWidthTolerance)) control$sliceAdaptWidthTolerance else 0.1
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodes      <- model$getDependencies(target)
        iLastTarget <- max(which(calcNodes %in% target))
        calcNodesTarget <- calcNodes[1:iLastTarget]
        calcNodesAfterTarget <- calcNodes[-(1:iLastTarget)]
        ## calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
        ## numeric value generation
        d                  <- length(targetAsScalar)
        discrete           <- sapply(targetAsScalar, function(x) model$isDiscrete(x))
        anyDiscrete        <- any(discrete)
        gammaMatrix        <- diag(d)         # matrix of orthogonal bases
        if(is.character(widthVec) && widthVec == 'oneVec')   widthVec <- rep(1,d)
        widthVecOriginal   <- widthVec
        nExpansions        <- rep(0, d)       # number of expansions
        nContracts         <- rep(0, d)       # number of contractions
        adaptFactorMaxIterOriginal <- adaptFactorMaxIter
        factorCounter      <- 0               # number of iterations since last factor adaptation
        factorTimesAdapted <- 0               # number of times factors have adapted
        empirSamp          <- matrix(0, nrow=adaptFactorInterval, ncol=d)   # matrix of posterior samples
        empirCov           <- diag(d)
        allWidthsAdapted   <- 0               # indicates whether all widths have finished adapting
        widthCounter       <- 0               # number of iterations since last width adaptation
        adaptWidthMaxIterOriginal <- adaptWidthMaxIter
        adaptWidthInterval <- 1               # interval to adapt widths; doubles each time widths are adaptated
        widthIndicatorVec  <- rep(1, d)       # indicator of which widths are still adapting
        ## checks
        if(d <= 1)                         stop('AF_slice sampler must be used on at least two target nodes')
        if(class(widthVec) != 'numeric')   stop('sliceWidths must be a numeric vector')
        if(length(widthVec) != d)          stop('sliceWidths must have length = ', d)
    },
    run = function() {
        for(i in 1:d) {
            eigenVec <- gammaMatrix[, i]
            width <- widthVec[i]
            u <- getLogProb(model, calcNodes) - rexp(1, 1)   # generate (log)-auxiliary variable: exp(u) ~ uniform(0, exp(lp))
            x0 <- values(model, target)                      # create random interval (L,R), of width 'width', around current value of target
            Lbound <- -1.0 * runif(1, 0, 1) * width
            Rbound <- Lbound + width
            L <- x0 + Lbound * eigenVec
            R <- x0 + Rbound * eigenVec
            maxStepsL <- floor(runif(1, 0, 1) * maxSteps)    # randomly allot (maxSteps-1) into maxStepsL and maxStepsR
            maxStepsR <- maxSteps - 1 - maxStepsL
            lp <- setAndCalculateTarget(L)
            while(maxStepsL > 0 & !is.nan(lp) & lp >= u) {   # step L left until outside of slice (max maxStepsL steps)
                Lbound <- Lbound - width
                L <- x0 + Lbound * eigenVec
                lp <- setAndCalculateTarget(L)
                maxStepsL <- maxStepsL - 1
                nExpansions[i] <<- nExpansions[i] + 1
            }
            lp <- setAndCalculateTarget(R)
            while(maxStepsR > 0 & !is.nan(lp) & lp >= u) {   # step R right until outside of slice (max maxStepsR steps)
                Rbound <- Rbound + width
                R <- x0 + Rbound * eigenVec
                lp <- setAndCalculateTarget(R)
                maxStepsR <- maxStepsR - 1
                nExpansions[i] <<- nExpansions[i] + 1
            }
            prop <- Lbound + runif(1, 0, 1) * (Rbound - Lbound)
            x1 <- x0 + prop * eigenVec
            lp <- setAndCalculateTarget(x1)
            while(is.nan(lp) | lp < u) {   # must be is.nan()
                if(prop < 0) { Lbound <- prop }
                else         { Rbound <- prop }
                nContracts[i] <<- nContracts[i] + 1
                prop <- Lbound + runif(1, 0, 1) * (Rbound - Lbound)
                x1 <- x0 + prop * eigenVec
                lp <- setAndCalculateTarget(x1)
            }
        }
        nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        if(allWidthsAdapted == 0)   adaptWidths()
        if(adaptFactorMaxIter > 0)  adaptFactors()
    },
    methods = list(
        setAndCalculateTarget = function(targetValues = double(1)) {
            if(anyDiscrete == 1)
                for(i in 1:d)
                    if(discrete[i] == 1)   targetValues[i] <- floor(targetValues[i])            
            values(model, target) <<- targetValues
            lp <- model$calculate(calcNodesTarget)
            if(lp == -Inf) return(-Inf)
            lp <- lp + calculate(model, calcNodesAfterTarget)
            ## Following lines were intended to prevent bugs in dynamic index cases,
            ## but in other cases they violate topological ordering.
            ##            lp <- calculate(model, target)
            ##            if(lp == -Inf) return(-Inf) # deals with dynamic index out of bounds
            ##            lp <- lp + calculate(model, calcNodesNoSelf)
            returnType(double())
            return(lp)
        },
        adaptFactors = function() {
            adaptFactorMaxIter <<- adaptFactorMaxIter - 1
            factorCounter <<- factorCounter + 1
            empirSamp[factorCounter, 1:d] <<- values(model, target)
            if(factorCounter == adaptFactorInterval) {
                for(i in 1:d)   empirSamp[, i] <<- empirSamp[, i] - mean(empirSamp[, i])
                empirCov <<- (t(empirSamp) %*% empirSamp) / (adaptFactorInterval - 1)
                gammaMatrix <<- eigen(empirCov)$vectors  # replace old factors with new factors
                factorTimesAdapted <<- factorTimesAdapted + 1
                factorCounter      <<- 0
                nExpansions        <<- rep(0, d)
                nContracts         <<- rep(0, d)
                allWidthsAdapted   <<- 0
                widthCounter       <<- 0
                adaptWidthMaxIter  <<- adaptWidthMaxIterOriginal
                adaptWidthInterval <<- 1
                widthIndicatorVec  <<- rep(1, d)
            }
        },
        adaptWidths = function() {
            adaptWidthMaxIter <<- adaptWidthMaxIter - 1
            widthCounter <<- widthCounter + 1
            if(widthCounter == adaptWidthInterval) {
                for(i in 1:d) {
                    if(widthIndicatorVec[i] == 1) {   # widths that are still adapting
                        if(nExpansions[i] == 0)   nExpansions[i] <<- 1
                        widthAdaptRatio <- nExpansions[i] / (nExpansions[i] + nContracts[i])
                        widthVec[i] <<- widthVec[i] * 2 * widthAdaptRatio
                        adaptWidthInterval <<- 2 * adaptWidthInterval   # double width adapt interval
                        nExpansions[i] <<- 0
                        nContracts[i] <<- 0
                        if(adaptWidthInterval > 16)  # once adapt interval is large enough, determine whether adaptation is finished
                            widthIndicatorVec[i] <<- (abs(widthAdaptRatio - .5) > adaptWidthTolerance)  # equals 1 if adaptation isn't finished
                    }
                }
                allWidthsAdapted <<- 1 - ceiling(mean(widthIndicatorVec))  # equals 1 only if all slice adapt indicators are 0
                widthCounter     <<- 0
            }
            if(adaptWidthMaxIter <= 0)  # alternatively, if max iters have been reached, stop adapting
                allWidthsAdapted <<- 1
        },
        reset = function() {
            gammaMatrix        <<- diag(d)
            empirCov           <<- diag(d)
            widthVec           <<- widthVecOriginal
            nExpansions        <<- rep(0, d)
            nContracts         <<- rep(0, d)
            adaptFactorMaxIter <<- adaptFactorMaxIterOriginal
            factorCounter      <<- 0
            factorTimesAdapted <<- 0
            allWidthsAdapted   <<- 0
            widthCounter       <<- 0
            adaptWidthMaxIter  <<- adaptWidthMaxIterOriginal
            adaptWidthInterval <<- 1
            widthIndicatorVec  <<- rep(1, d)
        }
    )
)

