## *********************************************************************
## vectorize Bernoulli draws option 1
## *********************************************************************
## dbern_vec is a user-defined distribution that computes all the
## dbern(psi)'s in one loop

## accepts a scalar p

## we need to take length as an argument because for rbern_vec, there
## is other way to determine the length (it is not the intended use of
## n)

## and the dbern_vec needs to have the save parameter arguments as
## rbern_vec

dbern_vec <- nimbleFunction(
    run = function(x = double(1),
                   p = double(),
                   length = double(),
                   log = double(0, default = 0)) {
        if(length != dim(x)[1]) stop('In dbern_vec, length != length of x')

        ## For the scalar p case, one could also replace the following 6
        ## lines with

        ##sumx <- sum(x)
        ##ans <- log(p) * sumx + log(1-p) * (length-sumx)

        ans <- 0
        logp <- log(p)
        log1mp <- log(1-p)
        for(i in 1:length)
            if(x[i] == 0) ans <- ans + log1mp
            else ans <- ans + logp

        returnType(double())
        if(log) return(ans)
        else return(exp(ans))
    })

## note the n is really a dummy variable for future extensions.  we
## always treat it as 1.  It should not be used as the length needed
rbern_vec <- nimbleFunction(
    run = function(n = double(), p = double(), length = double()) {
        ans <- double(1)
        setSize(ans, length)
        for(i in 1:length) ans[i] <- rbinom(1, size = 1, prob = p)
        returnType(double(1))
        return(ans)
    })

registerDistributions(list(dbern_vec = list(
                               BUGSdist = "dbern_vec(p, length)",
                               Rdist = "dbern_vec(p, length)",
                               range = c(0, 1),
                               types = c('value = double(1)',
                                         'p = double(0)',
                                         'length = double(0)'))
                           ))

## In this case, p provides information on the sizes of x, so we don't
## need separate arguments to provide that information

dbern_matrix <- nimbleFunction(
    run = function(x = double(2),
                   p = double(2),
                   log = double(0, default = 0)) {
        nrow = dim(p)[1]
        ncol = dim(p)[2]
        if(nrow != dim(x)[1]){
            stop('In dbern_matrix, dim[1] != number of rows in x')
        }
        if(ncol != dim(x)[1]){
            stop('In dbern_matrix, dim[1] != number of cols in x')
        }
        ans <- 0
        for(i in 1:nrow)
            for(j in 1:ncol)
                if(x[i,j] == 0) ans <- ans + log(1-p[i,j])
                else ans <- ans + log(p[i,j])

        returnType(double())
        if(log) return(ans)
        else return(exp(ans))
    })

rbern_matrix <- nimbleFunction(
    run = function(n = double(), p = double(2)) {
        ans = double(2)
        nrow = p[1]
        ncol = p[2]
        setSize(ans, nrow, ncol)
        for(i in 1:nrow)
            for(j in 1:ncol)
                ans[i,j] <- rbinom(1, size = 1, prob = p)
        returnType(double(2))
        return(ans)
    })

registerDistributions(list(dbern_matrix = list(
                               BUGSdist = "dbern_matrix(p)",
                               Rdist = "dbern_matrix(p)",
                               range = c(0, 1),
                               types = c('value = double(2)',
                                         'p = double(2)'))
                           ))


## *********************************************************************
## remove Zs, option 3
## *********************************************************************
## DynamicOccupancy removes the z's and muZ's from the model and computes
## the probability of all reps over all years for one site.
dDynamicOccupancy <- nimbleFunction(
    run = function(x = double(2),
                   nrep = double(),
                   psi1 = double(),
                   phi = double(1),
                   gamma = double(1),
                   p = double(1), # should double(2)
                   log = double(0, default = 0)) {

        numObs <- sum(x[,1])
        ## prob of the occupied sites out of the total sites given p
        ProbOccAndCount <- psi1 * dbinom(numObs, size = nrep, p = p[1], log = 0)
        ## prob of the empty sites
        ProbUnoccAndCount <- (1-psi1) * (numObs == 0)
        ## probably of the observed states
        ProbCount <- ProbOccAndCount + ProbUnoccAndCount
        ProbOccGivenCount <- ProbOccAndCount / ProbCount
        ProbOccNextTime <- ProbOccGivenCount * phi[1] +
            (1-ProbOccGivenCount) * gamma[1]
        ll <- log(ProbCount)
        nyears <- dim(x)[2]
        for(t in 2:nyears) {
            numObs <- sum(x[,t])
            ## similar change to this dbinom
            ProbOccAndCount <- ProbOccNextTime *
                dbinom(numObs, size = nrep, p = p[t], log = 0)
            ProbUnoccAndCount <- (1-ProbOccNextTime) * (numObs == 0)
            ProbCount <- ProbOccAndCount + ProbUnoccAndCount
            ProbOccGivenCount <- ProbOccAndCount / ProbCount
            ll <- ll + log(ProbCount)
            if(t < nyears) ProbOccNextTime <- ProbOccGivenCount * phi[t] +
                               (1-ProbOccGivenCount) * gamma[t]
        }
        if(log) return(ll)
        else return(exp(ll))
        returnType(double())
    }
)

## make matching changes here (at least to the input arguments)
## and check UserManual for non-needed "r" functions
rDynamicOccupancy <- nimbleFunction(
    run = function(n = double(),
                   nrep = double(),
                   psi1 = double(),
                   phi = double(1),
                   gamma = double(1),
                   p = double(1),
                   log = double(0, default = 0)) {
        nyear <- length(p)
        ans <- matrix(rbinom(nrep*nyear, 1, prob=p),
                      nrow=nrep, ncol=nyear)
        ## setSize(ans, nrep, nyear)
        returnType(double(2))
        return(ans)
    }
)

registerDistributions(list(
    dDynamicOccupancy = list(
        BUGSdist = "dDynamicOccupancy(nrep, psi1, phi, gamma, p)",
        Rdist = "dDynamicOccupancy(nrep, psi1, phi, gamma, p)",
        types = c('value = double(2)',
                  'nrep = double(0)',
                  'psi1 = double()',
                  'phi = double(1)',
                  'gamma = double(1)',
                  'p = double(1)'))
))
