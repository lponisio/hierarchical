getEffFUN <- function(pattern, save.dir, summary="efficiency", make.plot=TRUE){
    these.files <- list.files(save.dir, pattern=pattern)

    res.list.all <- list()
    for(res in 1:length(these.files)){
        load(file.path(save.dir,
                       these.files[res]))
        ss.ms.samples[[1]] <- rename_MCMC_comparison_method(
            rownames(ss.ms.samples[[1]]$summary),
            gsub(".Rdata", "", these.files[res]),
                                              comparison= ss.ms.samples[[1]])
        res.list.all[[res]] <- ss.ms.samples
    }

    occ.all <-
        do.call(combine_MCMC_comparison_results, unlist(res.list.all,
                                                        recursive=FALSE))

    if(make.plot){
    checkChains(occ.all$MCMCresults$samples,
                f.path = file.path(save.dir,
                                   "../figures/chains/%s.pdf"))
    dir.create(file.path(save.dir, sprintf("../figures/comparisons/%s",
                                         pattern)),
               showWarnings = FALSE)
    make_MCMC_comparison_pages(occ.all,
                           dir=file.path(save.dir,
                                         sprintf("../figures/comparisons/%s",
    pattern)))
    }

    if(summary == "efficiency"){
        out <- occ.all$MCMCresults$efficiency$min
        names(out) <- these.files
    }else if(summary == "mean"){
        out <- apply(occ.all$MCMCresults$summary, 3, function(x)
            x[,"mean"])
        rownames(out) <- these.files
    }
    return(out)
}

getSamplerNames <- function(effs){
    samp.names.prep <- sapply(strsplit(names(effs), "sampler"),
                              function(x) x[2])
    samp.names <-
        sapply(strsplit(samp.names.prep, "_mup"), function(x) x[1])
    mu.names <-
        sapply(strsplit(samp.names.prep, "_mup"), function(x) x[2])
    mu.names <-
        sapply(strsplit(mu.names, ".Rdata"), function(x) x[1])
    return(list(samplers=samp.names, mu=mu.names))
}

plotBar <- function(pattern, effs, adj.names){
    cols <- c("black", "grey")
    these.effs <- effs[grepl(pattern, names(effs))]
    names.effs <- getSamplerNames(these.effs)
    names(cols) <- unique(names.effs$mu)

    bp1 <- barplot(these.effs,
                   names="",
                   col=cols[names.effs$mu],
                   las=1,
                   xlab="", ylab="",
                   ylim=range(c(0,effs)))
    text(bp1, par('usr')[3] - adj.names,
         srt = 45, adj = 1,
         labels = names.effs$samplers,
         xpd = NA,
         cex=1)
}




plotEffSSMS <- function(){
    layout(matrix(1:4, nrow=2))
    par(oma=c(0, 7, 2, 1),
        mar=c(6, 1, 0.5, 3), cex.axis=1.5)
    ## latent, HP
    plotBar("latentTRUE", effsHP, 1)
    legend("topright", legend="a)", bty="n")
    mtext("Min effective sample size per second",
          2, line=6.5, cex=1.5, at=-2)
    mtext("Hyperparameters",
          2, line=4, cex=1.5)
    mtext("Latent states",
          3, line=1, cex=1.5)
    ## latent, no HP
    plotBar("latentTRUE", effsNoHP, 50)
    mtext("No hyperparameters",
          2, line=4, cex=1.5)
    legend("topright", legend="c)", bty="n")
    ## no latent, HP
    plotBar("latentFALSE", effsHP, 1)
    mtext("No latent states",
          3, line=1, cex=1.5)
    legend("topright", legend="b)", bty="n")
    legend("topleft", col=c("black", "grey"), pch=15,
           legend=c("Low detectability", "High detectability"),
           bty="n")
    ## no latent, no HP
    plotBar("latentFALSE", effsNoHP, 50)
    legend("topleft", legend="d)", bty="n")
}

