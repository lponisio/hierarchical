
## From section 6.9.2 - for data creation/organization:

if(!exists("DO_PLOT"))
    DO_PLOT <- FALSE

library(AHMbook)
## Code modified to use the SwissTits data set included in the AHMbook package
data(SwissTits)
## str(SwissTits)
## SwissTits$species  ## Available species

## Select Great tit and covariate data from 2013 and
##   drop 4 sites not surveyed in 2013
y0 <- SwissTits$counts[, , '2013', 'Great tit']
( NA.sites <- which(rowSums(is.na(y0)) == 3) ) ## Unsurveyed sites
y <- y0[-NA.sites, ]                 ## Drop them from the count data
tits <- SwissTits$sites[-NA.sites, ] ## Also drop from the site covariates
str(y)  ## Check the matrix of count data
## Get date and duration data for 2013, without the NA.sites rows:
date <- SwissTits$date[-NA.sites, , '2013']
dur <- SwissTits$dur[-NA.sites, , '2013']

## Plot observed data: counts vs survey date (Fig. 6-9)
if(DO_PLOT) matplot(t(date), t(y), type = "l", lwd = 3, lty = 1, frame = F, xlab = "Survey data (1 = April 1)", ylab = "Count of Great Tits")

## Load unmarked, create unmarked data frame and inspect result
library(unmarked)
time <- matrix(rep(as.character(1:3), nrow(y)), ncol = 3, byrow = TRUE)
umf <- unmarkedFramePCount(y = y,
                           siteCovs=data.frame(elev=scale(tits[,"elev"]), forest=scale(tits[,"forest"]), iLength=1/tits[,"rlength"]),
                           obsCovs=list(time = time, date = scale(date), dur = scale(dur)))
summary(umf)                            ## Summarize unmarked data frame
summary(apply(y, 1, max, na.rm = TRUE)) ## Summarize max counts

elev <- umf@siteCovs$elev   ;   elev2 <- elev^2
forest <- umf@siteCovs$forest   ;   forest2 <- forest^2
date <- matrix(umf@obsCovs$date, ncol = 3, byrow = TRUE)
dur <- matrix(umf@obsCovs$dur, ncol = 3, byrow = TRUE)
date[is.na(date)] <- 0   ;   date2 <- date^2
dur[is.na(dur)] <- 0   ;   dur2 <- dur^2
iRoute <- umf@siteCovs$iLength

## Design matrix for abundance model (no intercept)
lamDM <- model.matrix(~ elev + elev2 + forest + forest2 + elev:forest + elev:forest2 + iRoute)[,-1]

## Initial values
Nst <- apply(y, 1, max, na.rm = T) + 1
Nst[is.na(Nst)] <- round(mean(y, na.rm = TRUE))
Nst[Nst == "-Inf"] <- round(mean(y, na.rm = TRUE))
SGT_inits <- function(){ list(N = Nst,
                              beta0 = 0,
                              mean.p = rep(0.5, 3),
                              beta = runif(7, 0, 0),
                              alpha = runif(13, 0, 0)
                              )}


## Bundle data and choose to fit simple ZIP model (model 1)
SGT_data1 <- list(y = y,
                  lamDM = lamDM,
                  elev = elev,
                  date = date,
                  dur = dur,
                  elev2 = elev2,
                  date2 = date2,
                  dur2 = dur2)
## e = 1e-06,
## hlam.on = 0,
## hp.site.on = 0,
## hp.survey.on = 0
## nsite = 263,
## nrep = 3)

constants <- list(nsite = 263,
                  nrep = 3)
## e = 1e-06)



SGT_inits_full <- function(latent, hyper.param){
    ans <- list(N = Nst,
                beta0 = 0,
                mean.p = rep(0.5, 3),
                beta = runif(7, 0, 0),
                alpha = runif(13, 0, 0)
              , phi = 0.5
              , sd.lam = 1
              , sd.p.site = 1
              , sd.p.survey = 1
              , a = rbinom(constants$nsite, prob = 0.5, size = 1)
              , eps.lam = rnorm(constants$nsite)
              , eps.p.site = rnorm(constants$nsite)
              , eps.p.survey = matrix(rnorm(constants$nsite * constants$nrep),
                                      nrow = constants$nsite)
                )
    ans$a[ ans$N > 0 ] <- 1
    if(!latent){
        ans$N <- NULL
    }
    if(!hyper.param){
        ans$sd.lam <- NULL
        ans$sd.p.site <- NULL
        ans$sd.p.survey <- NULL
        ans$eps.lam <- NULL
        ans$eps.p.site <- NULL
        ans$eps.p.survey <- NULL
    }
    return(ans)
}


model.input <- list(data=SGT_data1,
                    constants = constants,
                    inits = SGT_inits_full(latent=latent, hyper.param=hyper.param))
