logit <- function(x) {
    log(x/(1 - x))
}

expit <- function(x) {
    exp(x)/(1 + exp(x))
}

genDynamicOccData <- function(nsite = 100,
                              nreps = 10,
                              nyear = 10,
                              psi1 = 0.4,
                              range.p = c(0.2, 0.4),
                              range.phi = c(0.6, 0.8),
                              range.gamma = c(0, 0.1)) {

    ##  Generation and analysis of simulated data for multi season
    ##  occupancy model (Kery and Schaud 2012)

    ## Function to simulate detection/nondetection data for dynamic site-occ model
    ## Annual variation in probabilities of patch survival, colonization and
    ## detection is specified by the bounds of a uniform distribution.

    ## Function arguments:
    ## nsite - Number of sites
    ## nreps - Number of replicate surveys
    ## nyear - Number of years
    ## psi1 - occupancy probability in first year
    ## range.p - bounds of uniform distribution from which annual p drawn
    ## range.psi and range.gamma - same for survival and colonization probability

    ## Set up some required arrays
    site <- 1:nsite					## Sites
    year <- 1:nyear					## Years
    psi <- rep(NA, nyear)				## Occupancy probability
    muZ <- z <- array(dim = c(nsite, nyear))	## Expected and realized occurrence
    y <- array(NA, dim = c(nsite, nreps, nyear))	## Detection histories

    ## Determine initial occupancy and demographic parameters
    psi[1] <- psi1				## Initial occupancy probability
    p <- runif(n = nyear, min = range.p[1], max = range.p[2])
    phi <- runif(n = nyear-1, min = range.phi[1], max = range.phi[2])
    gamma <- runif(n = nyear-1, min = range.gamma[1], max = range.gamma[2])

    ## Generate latent states of occurrence
    ## First year
    z[,1] <- rbinom(nsite, 1, psi[1])		## Initial occupancy state
    ## Later years
    for(site in 1:nsite){				## Loop over sites
        for(year in 2:nyear){				## Loop over years
            muZ[year] <- z[site, year-1]*phi[year-1] +
                (1-z[site, year-1])*gamma[year-1] ## Prob for occ.
            z[site,year] <- rbinom(1, 1, muZ[year])
        }
    }

    ## Generate detection/nondetection data
    for(site in 1:nsite){
        for(year in 1:nyear){
            prob <- z[site, year] * p[year]
            for(rep in 1:nreps){
                y[site, rep ,year] <- rbinom(1, 1, prob)
            }
        }
    }

    ## Compute annual population occupancy
    for (year in 2:nyear){
        psi[year] <- psi[year-1]*phi[year-1] + (1-psi[year-1])*gamma[year-1]
    }

    ## Plot apparent occupancy
    psi.app <- apply(apply(y, c(1,3), max), 2, mean)

    return(list(nsite = nsite,
                nreps = nreps,
                nyear = nyear,
                psi = psi,
                psi.app = psi.app,
                z = z,
                phi = phi,
                gamma = gamma,
                p = p,
                y = y))
}


## prep data for nimble model
prepModDataOcc <- function(data,
                           monitors = c("psi1", "phi",
                                        "gamma",
                                        "p"),
                           include.zs=TRUE){
    ## data zs with 0s set to NAs
    zs <- apply(data$y, c(1, 3), max)
    zs[zs == 0] <- NA

    ## initial condiations, NAs where 1s are in z, and 1s are where NA
    zinits <- zs
    zinits[zinits == 1] <- 2
    zinits[is.na(zinits)] <- 1
    zinits[zinits == 2] <- NA
    inits <- list(z = zinits)

    ## constants
    constants <- list(nsite = dim(data$y)[1],
                      nrep = dim(data$y)[2],
                      nyear = dim(data$y)[3])
    if(include.zs){
        model.data <- list(y = data$y, z = zs)
        inits <- list(z = zinits,
                      phi = data$phi,
                      gamma = data$gamma,
                      p = data$p,
                      psi1= data$psi[1])
    } else{
        model.data <- list(y = data$y)
        inits <- list(gamma = data$gamma,
                      p = data$p,
                      psi1= data$psi[1])
    }
    model.input <- list(data=model.data,
                        monitors=monitors,
                        constants=constants,
                        inits=inits)
    return(model.input)
}
