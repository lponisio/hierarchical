
## This user-defined distribution iterates over observations at each
## site and each visit and adds up all the probability pieces

dBernDetectionMatrix <- nimbleFunction(
  run = function(x = double(2),
    occProb = double(1),
    detectionProb = double(2),
    numReps = double(1),
    log = double(0)) {
    returnType(double())
    num_points <- dim(x)[1]
    ans <- 0
    for(j in 1:num_points) {
      probDetectionHistoryGivenOccupied <- 1
      probDetectionHistoryGivenUnoccupied <- 1
      for(k in 1:numReps[j]) {
        if(x[j,k] == 1) {
          probDetectionHistoryGivenOccupied <-
            probDetectionHistoryGivenOccupied * detectionProb[j,k]
          probDetectionHistoryGivenUnoccupied <- 0
        } else {
          probDetectionHistoryGivenOccupied <-
            probDetectionHistoryGivenOccupied * (1-detectionProb[j,k])
        }
      }
      ans <- ans + log(occProb[j] *
                       probDetectionHistoryGivenOccupied +
                       (1-occProb[j]) *
                       probDetectionHistoryGivenUnoccupied)
    }
    if(log) return(ans)
    return(exp(ans))
  })

## corresponding "r" function for random simulation

rBernDetectionMatrix <- nimbleFunction(
  run = function(n = double(), occProb = double(1),
    detectionProb = double(2), numReps = double(1)) {
    returnType(double(2))
    ans <- matrix()
    num_points <- length(occProb)
    max_numReps <- dim(detectionProb)[2]
    setSize(ans, num_points, max_numReps)
    for(j in 1:num_points) {
      if(runif(1,0,1) < occProb[j])
        occupied <- 1
      else
        occupied <- 0

      for(k in 1:numReps[j]) {
        if(occupied) {
          if(runif(1,0,1) < detectionProb[j, k]) ans[j, k] <- 1
          else ans[j, k] <- 0
        } else {
          ans[j, k] <- 0
        }
      }

    }
    return(ans)
  })

## This registers the user-provided distribution for use in a model
registerDistributions(list(dBernDetectionMatrix = list(
                             BUGSdist = "dBernDetectionMatrix(occProb, detectionProb, numReps)",
                             Rdist = "dBernDetectionMatrix(occProb, detectionProb, numReps)",
                             types = c('value = double(2)',
                               'occProb = double(1)',
                               'detectionProb = double(2)',
                               'numReps = double(1)'))
                           ))




