ms.ms.occ <- nimbleCode({
    ## multi-species priors
    ## detectablility
    mu.p.0     ~ dnorm(0,0.01)
    mu.p.day.1 ~ dnorm(0,0.01)
    mu.p.day.2 ~ dnorm(0,0.01)
    sigma.p.0     ~ dunif(0,10)
    sigma.p.day.1 ~ dunif(0,10)
    sigma.p.day.2 ~ dunif(0,10)

    ## phi/gam random intercepts
    mu.phi.0  ~ dnorm(0,0.01)
    mu.gam.0  ~ dnorm(0,0.01)
    sigma.phi.0 ~ dunif(0,10)
    sigma.gam.0 ~ dunif(0,10)

    ## hedgerow area proximity
    mu.phi.hr.area  ~ dnorm(0,0.01)
    mu.gam.hr.area  ~ dnorm(0,0.01)
    sigma.phi.hr.area ~ dunif(0,10)
    sigma.gam.hr.area ~ dunif(0,10)

    ## semi nat habitat area proximity
    mu.phi.nat.area  ~ dnorm(0,0.01)
    mu.gam.nat.area  ~ dnorm(0,0.01)
    sigma.phi.nat.area ~ dunif(0,10)
    sigma.gam.nat.area ~ dunif(0,10)

    ## floral resource diversity
    mu.phi.fra  ~ dnorm(0,0.01)
    mu.gam.fra  ~ dnorm(0,0.01)
    sigma.phi.fra ~ dunif(0,10)
    sigma.gam.fra ~ dunif(0,10)

    ## diet breadth
    phi.k  ~ dnorm(0,0.01)
    gam.k  ~ dnorm(0,0.01)

    ## body size
    phi.B  ~ dnorm(0,0.01)
    gam.B  ~ dnorm(0,0.01)

    ## interaction between hedgerow proximity and floral resources
    ## (habitat quality)
    phi.hr.area.fra  ~ dnorm(0,0.01)
    gam.hr.area.fra  ~ dnorm(0,0.01)

    ## interaction between semi nat habitat proximity and floral
    ## resources (habitat quality)
    phi.nat.area.fra  ~ dnorm(0,0.01)
    gam.nat.area.fra  ~ dnorm(0,0.01)

    ## interaction between hedgerow proximity and species diet breadth
    phi.hr.area.k  ~ dnorm(0,0.01)
    gam.hr.area.k  ~ dnorm(0,0.01)

    ## interaction between semi nat habitat proximity and species diet
    ## breadth
    phi.nat.area.k  ~ dnorm(0,0.01)
    gam.nat.area.k  ~ dnorm(0,0.01)

    ## interaction between hedgerow proximity and species body size
    phi.hr.area.B  ~ dnorm(0,0.01)
    gam.hr.area.B  ~ dnorm(0,0.01)

    ## interaction between semi nat habitat proximity and body size
    ## breadth
    phi.nat.area.B  ~ dnorm(0,0.01)
    gam.nat.area.B  ~ dnorm(0,0.01)

    ## species-specific  parameters
    for(sp in 1:nsp) {
        ## day
        p.0[sp]     ~ dnorm(mu.p.0,     sd=sigma.p.0)
        p.day.1[sp] ~ dnorm(mu.p.day.1, sd=sigma.p.day.1)
        p.day.2[sp] ~ dnorm(mu.p.day.2, sd=sigma.p.day.2)

        ## species specific intercept
        phi.0[sp] ~ dnorm(mu.phi.0, sd=sigma.phi.0)
        gam.0[sp] ~ dnorm(mu.gam.0, sd=sigma.gam.0)

        ## hedgerow area
        phi.hr.area[sp] ~ dnorm(mu.phi.hr.area,
                                sd=sigma.phi.hr.area)
        gam.hr.area[sp] ~ dnorm(mu.gam.hr.area,
                                sd=sigma.gam.hr.area)

        ## natural habitat
        phi.nat.area[sp] ~ dnorm(mu.phi.nat.area,
                                 sd=sigma.phi.nat.area)
        gam.nat.area[sp] ~ dnorm(mu.gam.nat.area,
                                 sd=sigma.gam.nat.area)

        ## fra
        phi.fra[sp] ~ dnorm(mu.phi.fra,
                            sd=sigma.phi.fra)
        gam.fra[sp] ~ dnorm(mu.gam.fra,
                            sd=sigma.gam.fra)

    }

    for(sp in 1:nsp) {
        for(site in 1:nsite) {

            for(yr in 1:nyear) {
                for(rep in 1:nrep[site,yr,sp]) {
                    logit(p[site,yr,rep,sp]) <-
                        p.0[sp] +
                        p.day.1[sp]*day[site,yr,rep,sp] +
                        p.day.2[sp]*day.2[site,yr,rep,sp]
                }
            }


            ## logit(p[site, 1:nyear, 1:max.nreps, sp]) <- p.0[sp] +
            ##     p.day.1[sp]*day[site, 1:nyear, 1:max.nreps,sp] +
            ##     p.day.2[sp]*day.2[site,1:nyear,1:max.nreps,sp]

            ## start off at the average for each species, site across years
            logit(phi.site.sp.mean[site,sp]) <- mean(phi[site, 1:(nyear-1),sp])
            logit(gam.site.sp.mean[site,sp]) <- mean(gam[site,1:(nyear-1),sp])

            psi.1[site,sp] <- gam.site.sp.mean[site,sp]/
                (1 - phi.site.sp.mean[site,sp] + gam.site.sp.mean[site,sp])

            ## occupancy in year 1
            ## psi is on a logit scale
            psi[site,1,sp] <- psi.1[site,sp]

            ## occupancy in subsequent years
            ## phi and gam are on a linear scale
            for(yr in 1:(nyear-1)) {
                phi[site,yr,sp] <-
                    phi.0[sp] +
                    phi.k*k[sp] +
                    phi.B*B[sp] +
                    phi.hr.area[sp]*HRarea[site] +
                    phi.nat.area[sp]*natural[site] +
                    phi.fra[sp]*fra[site, yr] +
                    phi.hr.area.fra*fra[site, yr]*HRarea[site] +
                    phi.nat.area.fra*fra[site, yr]*natural[site] +
                    phi.hr.area.k*k[sp]*HRarea[site] +
                    phi.nat.area.k*k[sp]*natural[site] +
                    phi.hr.area.B*B[sp]*HRarea[site] +
                    phi.nat.area.B*B[sp]*natural[site]

                gam[site,yr,sp] <-
                    gam.0[sp] +
                    gam.k*k[sp] +
                    gam.B*B[sp] +
                    gam.hr.area[sp]*HRarea[site] +
                    gam.nat.area[sp]*natural[site] +
                    gam.fra[sp]*fra[site, yr] +
                    gam.hr.area.fra*fra[site, yr]*HRarea[site] +
                    gam.nat.area.fra*fra[site, yr]*natural[site] +
                    gam.hr.area.k*k[sp]*HRarea[site] +
                    gam.nat.area.k*k[sp]*natural[site] +
                    gam.hr.area.B*B[sp]*HRarea[site] +
                    gam.nat.area.B*B[sp]*natural[site]
            }

        }
    }



    for(site in 1:nsite) {
        for(sp in 1:nsp) {
            X[site, 1:nyear, 1:max.nreps, sp] ~
                dDynamicOccupancy(nrep=nrep[site, 1:nyear, sp],
                                  psi1=psi[site,1,sp],
                                  phi=expit(phi[site,1:(nyear-1),sp]),
                                  gamma=expit(gam[site,1:(nyear-1),sp]),
                                  p=p[site, 1:nyear, 1:max.nreps, sp])

        }
    }
})

