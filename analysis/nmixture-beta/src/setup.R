


load('~/Dropbox/occupancy/analysis/nmixture-beta/src/data.rdata')


prepData <- function(spp.counts){

    nspp <- dim(spp.counts)[3]
    nsites=dim(spp.counts)[1]

    data.list <- list(counts=spp.counts)
    constants <- list(nspp=nspp,
                      nsites=nsites,
                      nvisits=dim(spp.counts)[2])

    Nst <- apply(spp.counts, c(3,1),  max, na.rm = TRUE) + 1
    Nst[is.na(Nst)] <- round(mean(spp.counts, na.rm = TRUE))
    Nst[Nst == "-Inf"] <- round(mean(spp.counts, na.rm = TRUE))

    inits <- list(mu.a0=runif(1), mu.p0=runif(1),  sigma.a0=runif(1),
                  sigma.p0=runif(1), N=Nst)

    model.input <- list(data=data.list,
                        inits=inits,
                        constants=constants)
    return(model.input)
}
