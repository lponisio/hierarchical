
genDynamicOccData <- function(nsite = 50,
                              nreps = 10,
                              nyear = 15,
                              mu.p,
                              psi1,
                              sigma.p,
                              mu.phi,
                              sigma.phi,
                              mu.gamma,
                              sigma.gamma) {

    ##  Generation and analysis of simulated data for multi season
    ##  occupancy model (adapted from Kery and Schaud 2012)
    ## Set up some required arrays
    site <- 1:nsite					## Sites
    year <- 1:nyear					## Years
    psi <- rep(NA, nyear)				## Occupancy probability
    muZ <- z <- array(dim = c(nsite, nyear))	## Expected and realized occurrence
    y <- array(NA, dim = c(nsite, nreps, nyear))	## Detection histories
    ## Determine initial occupancy and demographic parameters
    psi[1] <- psi1				## Initial occupancy probability
    p <-  rnorm(nyear, mu.p, sigma.p)
    phi <- rnorm(nyear -1, mu.phi, sigma.phi)
    gamma <- rnorm(nyear -1, mu.gamma, sigma.gamma)
    ## Generate latent states of occurrence
    ## First year
    z[,1] <- rbinom(nsite, 1, psi[1])		## Initial occupancy state
    ## Later years
    for(site in 1:nsite){				## Loop over sites
        for(year in 2:nyear){				## Loop over years
            muZ[site, year] <- z[site, year-1]*expit(phi[year-1]) +
                (1-z[site, year-1])*expit(gamma[year-1]) ## Prob for occ.
            z[site,year] <- rbinom(1, 1, muZ[site, year])
        }
    }
    ## Generate detection/nondetection data
    for(site in 1:nsite){
        for(year in 1:nyear){
            prob <- z[site, year]*expit(p[year])
            for(rep in 1:nreps){
                y[site, rep ,year] <- rbinom(1, 1, prob)
            }
        }
    }

    ## Compute annual population occupancy
    for (year in 2:nyear){
        psi[year] <- psi[year-1]*expit(phi[year-1]) +
            (1-psi[year-1])*expit(gamma[year-1])
    }

    return(list(nsite = nsite,
                nreps = nreps,
                nyear = nyear,
                psi = psi,
                z = z,
                phi = phi,
                gamma = gamma,
                p = p,
                y = y,
                mu.p = mu.p,
                sigma.p = sigma.p,
                mu.phi = mu.phi,
                sigma.phi = sigma.phi,
                mu.gamma = mu.gamma,
                sigma.gamma = sigma.gamma))
}


## prep data for nimble model
prepModDataOcc <- function(sim.input,
                           include.zs=TRUE){
    ## data zs with 0s set to NAs
    zs <- apply(sim.input$y, c(1, 3), max, na.rm=TRUE)
    zs[zs == 0] <- NA

    ## initial conditions, NAs where 1s in z, and 1s are where NA
    zinits <- zs
    zinits[zinits == 1] <- 2
    zinits[is.na(zinits)] <- 1
    zinits[zinits == 2] <- NA

    ## constants
    constants <- list(nsite = dim(sim.input$y)[1],
                      nrep = dim(sim.input$y)[2],
                      nyear = dim(sim.input$y)[3])

    model.data <- list(y = sim.input$y, z = zs)
    inits <- list(z = zinits,
                  phi = sim.input$phi,
                  gamma = sim.input$gamma,
                  p = sim.input$p,
                  psi1= sim.input$psi[1],
                  mu.p = sim.input$mu.p,
                  sigma.p = sim.input$sigma.p,
                  mu.phi = sim.input$mu.phi,
                  sigma.phi = sim.input$sigma.phi,
                  mu.gamma = sim.input$mu.gamma,
                  sigma.gamma = sim.input$sigma.gamma,
                  mu.p.mean = expit(sim.input$mu.p),
                  mu.gamma.mean = expit(sim.input$mu.gamma),
                  mu.phi.mean = expit(sim.input$mu.phi))
    if(!include.zs){
        model.data$z <- NULL
        inits$z <- NULL
    }
    model.input <- list(data=model.data,
                        constants=constants,
                        inits=inits)
    return(model.input)
}
