#!/usr/bin/env bash

Rscript -e 'install.packages("igraph", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("ggplot2", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("AHMbook", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("unmarked", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("RColorBrewer", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("reshape", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("coda", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("devtools", repos="http://cran.r-project.org")'


## within R
library(devtools)
install_github("nimble-dev/nimble", ref = "v0.7.1", subdir = "packages/nimble")
library(nimble)
