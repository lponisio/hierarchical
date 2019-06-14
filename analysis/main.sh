#!/usr/bin/env bash

## All paths are relative to the analysis folder within the github
## repo

## install necessary packages
bash analysis/packages.sh

## first create a folder for saving all of the models and figures

mkdir -p ../hierarchical_saved/multiSpp-multiSea/figures/comparisons
mkdir -p ../hierarchical_saved/multiSpp-multiSea/figures/chains
mkdir -p ../hierarchical_saved/multiSpp-multiSea/saved
mkdir -p ../hierarchical_saved/multiSpp-multiSea/tables

mkdir -p ../hierarchical_saved/singleSpp-multiSea/figures/comparisons
mkdir -p ../hierarchical_saved/singleSpp-multiSea/figures/chains
mkdir -p ../hierarchical_saved/singleSpp-multiSea/saved
mkdir -p ../hierarchical_saved/singleSpp-multiSea/tables

mkdir -p ../hierarchical_saved/multiSpp-singleSea/figures/comparisons
mkdir -p ../hierarchical_saved/multiSpp-singleSea/figures/chains
mkdir -p ../hierarchical_saved/multiSpp-singleSea/saved
mkdir -p ../hierarchical_saved/multiSpp-singleSea/tables

mkdir -p ../hierarchical_saved/nmixture/figures/comparisons
mkdir -p ../hierarchical_saved/nmixture/figures/chains
mkdir -p ../hierarchical_saved/nmixture/saved
mkdir -p ../hierarchical_saved/nmixture/tables

##***************************************************************
## run models
## **************************************************************
## each model run takes quite a bit of memory and time. The first
## argument is whether to run the models, the second is whether to plot
## the chains and ggplot comparison figures
Rscript analysis/singleSpp-multiSea/run.R TRUE TRUE 3e2
Rscript analysis/multiSpp-singleSea/run.R TRUE TRUE 3e2
Rscript analysis/multiSpp-multiSea/run.R TRUE TRUE 3e2
Rscript analysis/nmixture/run.R TRUE TRUE 5e2

## just plotting
Rscript analysis/singleSpp-multiSea/run.R FALSE TRUE 3e2
Rscript analysis/multiSpp-singleSea/run.R FALSE TRUE 3e2
Rscript analysis/multiSpp-multiSea/run.R FALSE TRUE 3e2
Rscript analysis/nmixture/run.R FALSE TRUE 5e2



install_github("nimble-dev/nimble",
               ref = "avoid-protect-stack-overflow",
               subdir = "packages/nimble")
