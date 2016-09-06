## library(devtools)
## install_github("nimble-dev/nimble",
##                ref = "devel",
##                subdir = "packages/nimble")

library(nimble)
library(mcmcplots)
library(igraph)
library(reshape)
source("../all/plotting.R")
source("../all/runNimble.R")
source('../all/sampelrs/sampler_RW_shift.R')
source("src/reformatData.R")
source("src/multispeciesOcc.R")

survey_data <- read.csv("data/occupancy_data.csv")
species_groups <- read.csv("data/species_groups.csv")
survey_dates <- read.csv("data/survey_dates.csv")
habitat <- read.csv("data/habitat.csv")

## reformat data
data <- reformatData(survey_data,
                      survey_dates,
                      species_groups,
                      habitat,
                      n_zeroes)


num_species <- data$num_species
num_points <- data$num_points
num_reps <- data$num_reps

## Z data for whether or not a species was ever observed
## zs with 1s as 1s and 0s as NAs
zs <- apply(data$X, c(1, 3), max)
zs[zs == 0] <- NA

model_data <- list(Z = zs,
                   X = data$X_aug,
                   ground = data$ground,
                   mid = data$mid,
                   habitat_ind = data$habitat_ind,
                   ufc_linear = data$ufc_linear,
                   ufc_quadratic = data$ufc_quadratic,
                   ba_linear = data$ba_linear,
                   ba_quadratic = data$ba_quadratic,
                   date_linear = data$date_linear,
                   date_quadratic = data$date_quadratic)


psi_mean_draw <- runif(1, 0.25, 1)

## initial values
omega_draw <- runif(1, num_species/(num_species + n_zeroes), 1)


## inital conditions with 1s as NAs and Nas as 1s
zinits <- zs
zinits[zinits == 1] <- 2
zinits[is.na(zinits)] <- 1
zinits[zinits == 2] <- NA


inits <-list(Z=zinits,
             omega = omega_draw,
             w = c(rep(1, num_species),
               rbinom(n_zeroes, size = 1, prob = omega_draw)),
             u_cato = rnorm(num_species + n_zeroes),
             v_cato = rnorm(num_species + n_zeroes),
             u_fcw = rnorm(num_species + n_zeroes) ,
             v_fcw = rnorm(num_species + n_zeroes),

             a1 = rnorm(num_species + n_zeroes),
             a2 = rnorm(num_species + n_zeroes),
             a3 = rnorm(num_species + n_zeroes),
             a4 = rnorm(num_species + n_zeroes),
             b1 = rnorm(num_species + n_zeroes),
             b2 = rnorm(num_species + n_zeroes))

## constants
constants <- list(num_species = num_species,
                  num_points = num_points,
                  num_reps = num_reps,
                  n_zeroes = n_zeroes)

## parameters to monitor
monitors <- c('N', 'N_site', 'N_ground', 'N_mid', 'mu_a1','mu_a2',
              'mu_a3','mu_a4', 'sigma_a1','sigma_a2','sigma_a3','sigma_a4',
              'cato_occ_mean', 'fcw_occ_mean', 'cato_det_mean',
              'fcw_det_mean', 'sigma_ucato', 'sigma_vcato', 'sigma_ufcw',
              'sigma_vfcw', 'mu_b1', 'mu_b2', 'sigma_b1', 'sigma_b2')


## for the non-data augmented case
if(n_zeroes == 0){
  inits[c("w", "omega")] <- NULL
  constants[c("n_zeroes")] <- NULL
  monitors <- monitors[!monitors == "N"]
}

## additional constants and dats for models where z is removed
constants$max_num_reps <- max(constants$num_reps)
model_data$onesRow <- matrix(rep(1, constants$max_num_reps), nrow = 1)


## mcmc settings
scale <- 1e2
burnin <- 1e2*scale
niter <- (1e3)*scale
