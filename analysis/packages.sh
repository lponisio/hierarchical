#!/usr/bin/env bash

Rscript -e 'install.packages("igraph", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("nimble", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("ggplot2", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("AHMbook", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("unmarked", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("RColorBrewer", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("reshape", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("coda", repos="http://cran.r-project.org")'

## avoids wierd bug in nimble with large models
install_github("nimble-dev/nimble",
               ref = "avoid-protect-stack-overflow",
               subdir = "packages/nimble")
