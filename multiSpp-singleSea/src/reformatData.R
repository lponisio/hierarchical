

createEncounterHistory <- function(survey_data){
  ## The detection/non-detection data is reshaped into a three
  ## dimensional array X where the first dimension, j, is the point;
  ## the second dimension, k, is the rep; and the last dimension, i, is
  ## the species.
  survey_data$occupancy <- rep(1, dim(survey_data)[1])
  
  X = melt(survey_data,
    id.var = c("species", "point", "repetition"),
    measure.var = "occupancy")
  
  X = cast(X, point ~ repetition ~ species)
  
  ## Change abundance to presence/absence
  ## Set counts 1 and greater to indicators
  X[which(X > 0)] <- 1
  
  ## Create all zero encounter histories to add to the detection array X 
  ## as part of the data augmentation to account for additional 
  ## species (beyond the n observed species). 
  X_zero = matrix(0, nrow=dim(X)[1], ncol=dim(X)[2])
  
  return(list(X=X,X_zero=X_zero))
}

addMissingData <- function(histories, survey_dates){
  ## 'Add' missing data:set X and X_zero for the unsurveyed 
  ## repetitions to NA
  X <- histories$X
  X_zero <- histories$X_zero
  for(point in 1:length(unique(survey_data$point))){
    point_index <- which(survey_dates$point == row.names(X)[point])
    missing <- is.na(survey_dates[point_index,][,-(1:1)])
    X[point, missing,] <- NA
    X_zero[point, missing] <- NA
  }
  return(list(X=X,X_zero=X_zero))
}

zeroAugment <- function(histories, n_zeroes){
  ## X_aug is the augmented version of X.  The first n species were
  ## actually observed and the n+1 through n_zeroes species are all
  ## zero encounter histories create an empty 3D array with n + n_zero
  ## species
  X <- histories$X
  X_zero <- histories$X_zero
  X_dim <- dim(X)
  X_dim[3] <- X_dim[3] + n_zeroes
  X_aug <- array( NA, dim = X_dim) 
  
  ## fill in the array with the occurrence data
  X_aug[,,1:dim(X)[3]] <-  X 
  
  ## fill the zero histories
  X_aug[,,-(1:dim(X)[3])] <- rep(X_zero, n_zeroes)
  return(X_aug)
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



reformatData <- function(survey_data,
                         survey_dates,
                         species_groups,
                         habitat,
                         n_zeroes){

  manipulateData <- function(survey_data,
                             survey_dates,
                             species_groups,
                             habitat, n_zeroes){
    histories <- createEncounterHistory(survey_data)
    histories <- addMissingData(histories, survey_dates)
    X_aug <- zeroAugment(histories, n_zeroes)
    
    num_species <- length(unique(survey_data$species))
    num_points <- length(unique(survey_data$point))
    num_reps <- apply(survey_dates[,-(1:1)], 
                      1, 
                      function(x) length(which(!is.na(x))))
    
    ## Create an indicator vector for each assemblage (ground, mid-story)
    ground <- mid <- rep(0, dim(species_groups)[1])
    ground[which(species_groups$group == 1)] <- 1
    mid[which(species_groups$group == 2)] <- 1
    
    ## Create a vector to indicate which habitat type each point is in 
    ## (CATO = 1; FCW  = 0)
    habitat_ind <- rep(0, num_points)
    habitat_ind[grep("CAT", row.names(histories$X))] <- 1
    
    ## Standardize variables
    ufc <-  siteLevelStandardized(habitat$ufc)
    ufc_linear <- ufc$linear
    ufc_quadratic <- ufc$quadratic
    
    ba <- siteLevelStandardized(habitat$ba)
    ba_linear <- ba$linear
    ba_quadratic <- ba$quadratic
    
    date <- surveyLevelStandardized(survey_dates[,-(1:1)])
    date_linear <- date$linear
    date_quadratic <- date$quadratic
    
    return(list(X_aug = X_aug, num_species = num_species,
                num_points = num_points,
                num_reps = num_reps, ground = ground, mid = mid, 
                habitat_ind = habitat_ind, ufc_linear = ufc_linear, 
                ufc_quadratic = ufc_quadratic, ba_linear = ba_linear, 
                ba_quadratic = ba_quadratic, date_linear = date_linear, 
                date_quadratic = date_quadratic)) 
  }


  return(manipulateData(survey_data,
                        survey_dates,
                        species_groups,
                        habitat,
                        n_zeroes))
}
