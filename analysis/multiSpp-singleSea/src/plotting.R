
getEffFUN <- function(pattern, save.dir, summary="efficiency",
                      make.plot=TRUE, adj.xlab=0){

    plotPointsMakeTableWapper <- function() {
        plotPointsMakeTable(occ.all=occ.all, adj.xlab=adj.xlab)
    }

    these.files <- list.files(save.dir, pattern=pattern)

    res.list.all <- list()
    for(res in 1:length(these.files)){
        load(file.path(save.dir,
                       these.files[res]))
        ms.ss.samples[[1]] <- rename_MCMC_comparison_method(
            rownames(ms.ss.samples[[1]]$summary),
            gsub(".Rdata", "", these.files[res]),
            comparison= ms.ss.samples[[1]])
        res.list.all[[res]] <- ms.ss.samples
    }

    occ.all <-
        do.call(combine_MCMC_comparison_results, unlist(res.list.all,
                                                        recursive=FALSE))

    if(make.plot){
        ## checkChains(occ.all$MCMCresults$samples,
        ##             f.path = file.path(save.dir,
        ##                                "../figures/chains/%s.pdf"))
        dir.create(file.path(save.dir,sprintf("../figures/comparisons/%s",
                                              pattern)),
                   showWarnings = FALSE)
        make_MCMC_comparison_pages(occ.all,
                                   dir=file.path(save.dir,
                                                 sprintf("../figures/comparisons/%s",
                                                         pattern)))
    }

    pdf.f(plotPointsMakeTableWapper,
          file= file.path(save.dir,
                          sprintf("../figures/comparisons/msss_%s.pdf",
                                  pattern)),
          height=4, width=8)

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
        sapply(strsplit(samp.names.prep, ".Rdata"), function(x) x[1])
    return(list(samplers=samp.names))
}


plotEffMSSS <- function(){
    layout(matrix(1:4, nrow=2))
    par(oma=c(0, 7, 2, 1),
        mar=c(6, 1, 0.5, 3), cex.axis=1.5)
    ## latent, HP
    plotBar("latentTRUE", effsHP, 0.05)
    legend("topleft", legend="a)", bty="n")
    mtext("Min ESS/second",
          2, line=6.5, cex=1.5, at=-0.25)
    mtext("More hierarchy",
          2, line=4, cex=1.5)
    mtext("Latent states",
          3, line=1, cex=1.5)
    ## latent, no HP
    plotBar("latentTRUE", effsNoHP, 2)
    mtext("Less hierarchy",
          2, line=4, cex=1.5)
    legend("topleft", legend="c)", bty="n")
    ## no latent, HP
    plotBar("latentFALSE", effsHP, 0.05)
    mtext("Latent state integration",
          3, line=1, cex=1.5)
    legend("topleft", legend="b)", bty="n")
    ## no latent, no HP
    plotBar("latentFALSE", effsNoHP, 2)
    legend("topleft", legend="d)", bty="n")
}


plotBar <- function(pattern, effs, adj.names){
    these.effs <- effs[grepl(pattern, names(effs))]
    names.effs <- getSamplerNames(these.effs)

    bp1 <- barplot(these.effs,
                   names="",
                   las=1,
                   xlab="", ylab="", col= "black",
                   ylim=range(c(0,effs)) + c(0, min(effs)))
    text(bp1, par('usr')[3] - adj.names,
         srt = 45, adj = 1,
         labels = names.effs$samplers,
         xpd = NA,
         cex=1)
}
