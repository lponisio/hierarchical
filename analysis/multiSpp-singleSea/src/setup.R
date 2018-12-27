

createEncounterHistory <- function(survey.data){
    ## The detection/non-detection data is reshaped into a three
    ## dimensional array X where the first dimension, j, is the point;
    ## the second dimension, k, is the rep; and the last dimension, i, is
    ## the species.
    survey.data$occupancy <- rep(1, dim(survey.data)[1])

    X = melt(survey.data,
             id.var = c("species", "point", "repetition"),
             measure.var = "occupancy")

    X = cast(X, point ~ repetition ~ species)

    ## Change abundance to presence/absence
    ## Set counts 1 and greater to indicators
    X[which(X > 0)] <- 1

    ## Create all zero encounter histories to add to the detection array X
    ## as part of the data augmentation to account for additional
    ## species (beyond the n observed species).
    X.zero = matrix(0, nrow=dim(X)[1], ncol=dim(X)[2])

    return(list(X=X,X.zero=X.zero))
}

addMissingData <- function(histories, survey.dates){
    ## 'Add' missing data:set X and X.zero for the unsurveyed
    ## repetitions to NA
    X <- histories$X
    X.zero <- histories$X.zero
    for(point in 1:length(unique(survey.data$point))){
        point.index <- which(survey.dates$point == row.names(X)[point])
        missing <- is.na(survey.dates[point.index,][,-(1:1)])
        X[point, missing,] <- NA
        X.zero[point, missing] <- NA
    }
    return(list(X=X,X.zero=X.zero))
}

zeroAugment <- function(histories, n.zeroes){
    ## X.aug is the augmented version of X.  The first n species were
    ## actually observed and the n+1 through n.zeroes species are all
    ## zero encounter histories create an empty 3D array with n + n.zero
    ## species
    X <- histories$X
    X.zero <- histories$X.zero
    X.dim <- dim(X)
    X.dim[3] <- X.dim[3] + n.zeroes
    X.aug <- array( NA, dim = X.dim)

    ## fill in the array with the occurrence data
    X.aug[,,1:dim(X)[3]] <-  X

    ## fill the zero histories
    X.aug[,,-(1:dim(X)[3])] <- rep(X.zero, n.zeroes)
    return(X.aug)
}

siteLevelStandardized <- function(parameter){
    m <- mean(parameter, na.rm = TRUE)
    sd <- sd(parameter, na.rm = TRUE)
    linear <- (parameter - m)/sd
    quadratic <- linear * linear
    return(list(linear=linear,quadratic=quadratic))
}

surveyLevelStandardized <- function(parameter){
    parameter <- as.matrix(parameter)
    m <- mean(parameter, na.rm = TRUE)
    sd <- sd(parameter, na.rm = TRUE)
    linear <- (parameter - m) /  sd
    quadratic <- linear * linear
    linear <- as.matrix(linear)
    quadratic <- as.matrix(quadratic)
    return(list(linear=linear,quadratic=quadratic))
}



reformatData <- function(survey.data,
                         survey.dates,
                         species.groups,
                         habitat,
                         n.zeroes){

    manipulateData <- function(survey.data,
                               survey.dates,
                               species.groups,
                               habitat, n.zeroes){
        histories <- createEncounterHistory(survey.data)
        histories <- addMissingData(histories, survey.dates)
        X.aug <- zeroAugment(histories, n.zeroes)

        num.species <- length(unique(survey.data$species))
        num.points <- length(unique(survey.data$point))
        num.reps <- apply(survey.dates[,-(1:1)],
                          1,
                          function(x) length(which(!is.na(x))))

        ## Create an indicator vector for each assemblage (ground, mid-story)
        ground <- mid <- rep(0, dim(species.groups)[1])
        ground[which(species.groups$group == 1)] <- 1
        mid[which(species.groups$group == 2)] <- 1

        ## Create a vector to indicate which habitat type each point is in
        ## (CATO = 1; FCW  = 0)
        habitat.ind <- rep(0, num.points)
        habitat.ind[grep("CAT", row.names(histories$X))] <- 1

        ## Standardize variables
        ufc <-  siteLevelStandardized(habitat$ufc)
        ufc.linear <- ufc$linear
        ufc.quadratic <- ufc$quadratic

        ba <- siteLevelStandardized(habitat$ba)
        ba.linear <- ba$linear
        ba.quadratic <- ba$quadratic

        date <- surveyLevelStandardized(survey.dates[,-(1:1)])
        date.linear <- date$linear
        date.quadratic <- date$quadratic

        return(list(X.aug = X.aug, num.species = num.species,
                    num.points = num.points,
                    num.reps = num.reps, ground = ground, mid = mid,
                    habitat.ind = habitat.ind, ufc.linear = ufc.linear,
                    ufc.quadratic = ufc.quadratic, ba.linear = ba.linear,
                    ba.quadratic = ba.quadratic, date.linear = date.linear,
                    date.quadratic = date.quadratic))
    }


    return(manipulateData(survey.data,
                          survey.dates,
                          species.groups,
                          habitat,
                          n.zeroes))
}


## prep data for model
prepMutiSpData <- function(survey.data,
                           survey.dates,
                           species.groups,
                           habitat,
                           n.zeroes,
                           remove.zs=TRUE,
                           vectorized=TRUE,
                           hyper.param=TRUE){
    ## reformat data
    data <- reformatData(survey.data,
                         survey.dates,
                         species.groups,
                         habitat,
                         n.zeroes)

    num.species <- data$num.species
    num.points <- data$num.points
    num.reps <- data$num.reps

    ## Z data for whether or not a species was ever observed
    ## zs with 1s as 1s and 0s as NAs
    zs <- apply(data$X, c(1, 3), max, na.rm=TRUE)
    zs[zs == 0] <- NA
    zs[!is.finite(zs)] <- NA

    model.data <- list(Z = zs,
                       X = data$X.aug,
                       ground = data$ground,
                       mid = data$mid,
                       habitat.ind = data$habitat.ind,
                       ufc.linear = data$ufc.linear,
                       ufc.quadratic = data$ufc.quadratic,
                       ba.linear = data$ba.linear,
                       ba.quadratic = data$ba.quadratic,
                       date.linear = data$date.linear,
                       date.quadratic = data$date.quadratic)

    psi.mean.draw <- runif(1, 0.25, 1)

    ## initial values
    omega.draw <- runif(1, num.species/(num.species + n.zeroes), 1)

    ## inital conditions. 1 should be NA, NA should be a 0 or 1
    zinits <- zs
    zinits[zinits == 1] <- 2
    ## zinits[is.na(zinits)] <- 1
    zinits[is.na(zinits)] <- sample(0:1, sum(is.na(zinits)),
                                    replace=TRUE)
    zinits[zinits == 2] <- NA

    if(hyper.param){
        inits <-list(Z=zinits,
                     omega = omega.draw,
                     w = c(rep(1, num.species),
                           rbinom(n.zeroes, size = 1, prob = omega.draw)),
                     u.cato = rnorm(num.species + n.zeroes),
                     v.cato = rnorm(num.species + n.zeroes),
                     u.fcw = rnorm(num.species + n.zeroes) ,
                     v.fcw = rnorm(num.species + n.zeroes),
                     a1 = rnorm(num.species + n.zeroes),
                     a2 = rnorm(num.species + n.zeroes),
                     a3 = rnorm(num.species + n.zeroes),
                     a4 = rnorm(num.species + n.zeroes),
                     b1 = rnorm(num.species + n.zeroes),
                     b2 = rnorm(num.species + n.zeroes))
    } else{
        inits <-list(Z=zinits,
                     omega = omega.draw,
                     w = c(rep(1, num.species),
                           rbinom(n.zeroes, size = 1, prob = omega.draw)),
                     u.cato = rnorm(1),
                     v.cato = rnorm(1),
                     u.fcw = rnorm(1) ,
                     v.fcw = rnorm(1),
                     a1 = rnorm(1),
                     a2 = rnorm(1),
                     a3 = rnorm(1),
                     a4 = rnorm(1),
                     b1 = rnorm(1),
                     b2 = rnorm(1))
    }

    ## constants
    constants <- list(num.species = num.species,
                      num.points = num.points,
                      num.reps = num.reps,
                      n.zeroes = n.zeroes)

    ## for the non-data augmented case
    if(n.zeroes == 0){
        inits[c("w", "omega")] <- NULL
        constants[c("n.zeroes")] <- NULL

    }

    ## these were removed from the model
    model.data[["ground"]] <- NULL
    model.data[["mid"]] <- NULL

    if(vectorized){
        ## Since calculations with date_linear and date_quadratic are now
        ## vectorized, we'll set the NAs to 0
        model.data$date.linear[is.na(model.data$date.linear)] <- 0
        model.data$date.quadratic[is.na(model.data$date.quadratic)] <- 0
    }
    if(remove.zs) {
        ## zs are removed from these models
        model.data[["Z"]] <- NULL
        inits[["Z"]] <- NULL
        ## additional constants and dats for models where z is removed
        constants$max.num.reps <- max(constants$num.reps)
        model.data$onesRow <- matrix(rep(1, constants$max.num.reps),
                                     nrow=1)
    }
    return(list(constants=constants,
                inits=inits,
                data=model.data))
}

## new function for generating initial values for use in nimbleModel(),
## using posterior summary statistics
## -DT
genInits <- function(summary) {
    e <- new.env()
    for(n in names(summary))    eval(parse(text = paste0('e$', n, ' <- summary[\'', n, '\']'))[[1]])
    new_inits <- list()
    for(n in ls(e)) new_inits[[n]] <- e[[n]]
    return(new_inits)
}


