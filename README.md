Analysis code for "Factors affecting MCMC performance in complex hierarchical models".  

Analyses were conducted in R 3.6.1 and NIMBLE v0.7.1.  Within the analysis folder, the packages.sh script will install all necessary packages, though the specific version of NIMBLE needs to be installed from within R. The main.sh file will run all models and make the corresponding comparison figures.

Within the analysis folder, each model has a specific folder. Each contains a run.R file, that will run all the flavors of the model and MCMC compared in the study. Within the "src" folder of each model, there is an initialize.R file that loads necaeesary packages and data, and contains functions to srouce for running the models; a prep.R file that preps the data for the model; a model.R file that contains the BUGS models, and a plotting.R file that preps and plots the comparisons figures. Some scripts used in all models (i.e., for plotting chains and making tables) are in the "all" folder.
