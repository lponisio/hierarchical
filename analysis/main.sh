#!/usr/bin/env bash

## All paths are relative to the analysis folder within the github
## repo

## install necessary packages
bash analysis/packages.sh

## first create a folder for saving all of the models and figures

##***************************************************************
## run models
## **************************************************************
## each model run takes quite a bit of memory and time. The first
## argument is whether to run the model, the second is whether to plot
## the chains
Rscript analysis/singleSpp-multiSea/run.R TRUE TRUE
Rscript analysis/multiSpp-singleSea/run.R TRUE TRUE
Rscript analysis/multiSpp-multiSea/run.R TRUE TRUE

## just plotting
Rscript analysis/singleSpp-multiSea/run.R FALSE FALSE
Rscript analysis/multiSpp-singleSea/run.R FALSE FALSE
Rscript analysis/multiSpp-multiSea/run.R FALSE FALSE
