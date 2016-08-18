
  pdf.f <- function(f, file, ...) {
    cat(sprintf('Writing %s\n', file))
    pdf(file, ...)
    on.exit(dev.off())
    f()
  }

## function to plot values of parameters as a function of
## iterations. Input is a MCMCcompare object

checkChains <- function(all.mods.samps, f.path){
  niter <- dim(samp.all)[3]
  params <- names(samp.all[1,,1])
  mods <- names(samp.all[,1,1])
  if(is.null(mods)) mods <- "nimble"

  lapply(mods, function(z){
    f <- function(){
      layout(matrix(1:4, ncol=2))

      apply(samps[z,,], 1, function(x){
        plot(x, type="l",
             xlab = 'iteration',
             main= params[which(apply(samps[z,,], 1, function(y)
               all(match(x,y))))]
             )
      })
    }

    pdf.f(f,
          file= file.path(sprintf(f.path, z)),
          height=11, width=8.5)
  })
}
