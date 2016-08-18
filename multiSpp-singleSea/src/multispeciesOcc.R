
## This user-defined distribution iterates over observations at each
## site and each visit and adds up all the probability pieces

dBernDetectionMatrix <- nimbleFunction(
  run = function(x = double(2),
                 occProb = double(1),
                 detectionProb = double(2),
                 num_reps = double(1),
                 log = double(0, default = 0)) {
    returnType(double())
    num_points <- dim(x)[1]
    ans <- 0
    ## new code 5/25/16
    for(j in 1:num_points) {
        probDetectionHistoryGivenOccupied <- 1
        probDetectionHistoryGivenUnoccupied <- 1
        for(k in 1:num_reps[j]) {
            if(x[j,k] == 1) {
                probDetectionHistoryGivenOccupied <- probDetectionHistoryGivenOccupied * detectionProb[j,k]
                probDetectionHistoryGivenUnoccupied <- 0
            } else {
                probDetectionHistoryGivenOccupied <- probDetectionHistoryGivenOccupied * (1-detectionProb[j,k])
            }
        }
        ans <- ans + log(occProb[j] * probDetectionHistoryGivenOccupied + (1-occProb[j]) * probDetectionHistoryGivenUnoccupied)
 ## old code that looks like it was wrong:       
 ##       for(k in 1:num_reps[j]) {
 ##       if(x[j,k] == 1) ans <- ans + log(occProb[j] * detectionProb[j,k])
 ##       else ans <- ans + log(occProb[j] * (1-detectionProb[j,k]) +
 ##                             (1-occProb[j]))
 ##     }
    }
    if(log) return(ans)
    return(exp(ans))
  })

## This is the corresponding "r" function for random simulation, which
## we currently require but won't actually be used in this case
rBernDetectionMatrix <- nimbleFunction(
  run = function(n = integer(), occProb = double(1),
                 detectionProb = double(2), num_reps = double(1)) {
    returnType(double(2))
    declare(ans, double(2))
    num_points <- length(occProb)
    max_num_reps <- dim(detectionProb)[2]
    setSize(ans, num_points, max_num_reps)
    for(j in 1:num_points) {
        ## new code 5/25/16
        if(runif(1,0,1) < occProb[j])
            occupied <- 1
        else
            occupied <- 0

        for(k in 1:num_reps[j]) {
            if(occupied) {
                if(runif(1,0,1) < detectionProb[j, k]) ans[j, k] <- 1
                else ans[j, k] <- 0
            } else {
                ans[j, k] <- 0
            }
        }
        ## old code that looks like it was wrong
##      for(k in 1:num_reps[j]) {
##        if(runif(1,0,1) < occProb[j] * detectionProb[j, k]) ans[j, k] <- 1
##        else ans[j, k] <- 0
##      }
    }
    return(ans)
  })

## This registers the user-provided distribution for use in a model (a
## bit heavy syntax at the moment)
registerDistributions(list(dBernDetectionMatrix = list(
            BUGSdist = "dBernDetectionMatrix(occProb, detectionProb, num_reps)",
            Rdist = "dBernDetectionMatrix(occProb, detectionProb, num_reps)",
            types = c('value = double(2)',
                      'occProb = double(1)',
                      'detectionProb = double(2)',
                      'num_reps = integer(1)'))
            ))





## The next three calls would set up the Bernoulli + detection
## probability distribution if we were doing it one observartion at a
## time.  I wrote this and then decided to replace it with the matrix
## version below:



## dBernDetection <- nimbleFunction(
##     run = function(x = double(), occProb = double(), detectionProb = double(), log = double(0, default = 0)) {
##         returnType(double())
##         if(x == 1) ans <- occProb * detectionProb
##         else ans <- occProb * (1-detectionProb) + (1-occProb)
##         if(log) return(log(ans))
##         return(ans)
##     })
## rBernDetection <- nimbleFunction(
##     run = function(n = integer(), occProb = double(), detectionProb = double()) {
##         returnType(double())
##         if(runif(1,0,1) < occProb * detectionProb) return(1)
##         return(0)
##     })
## registerDistributions(list(dBernDetection = list(
##     BUGSdist = "dBernDetection(occProb, detectionProb)",
##     Rdist = "dBernDetection(occProb, detectionProb)",
##     range = c(0, 1),
##     types = c('value = double()', 'occProb = double(0)', 'detectionProb = double(0)'))
##   ))

