setwd("occupancy/analysis/spatial")


rm(list=ls())
source('src/initialize.R')


set.seed(444)
dats <- genSpatialOccData(alpha=0.8, sigma=10, delta=0.1)
model.input <- prepModData(dats$data, dats$y, dats$distance,
                           nsite=500, inits=dats$inits)


rm(list=ls())
source('src/initialize.R')

set.seed(444)
dats <- genSpatialOccData(alpha=0.5, sigma=10, delta=0.1)
model.input <- prepModData(dats$data, dats$y, dats$distance,
                           nsite=500, inits=dats$inits)

rm(list=ls())
source('src/initialize.R')


set.seed(444)
dats <- genSpatialOccData(alpha=0.1, sigma=10, delta=0.1)
model.input <- prepModData(dats$data, dats$y, dats$distance,
                           nsite=500, inits=dats$inits)

rm(list=ls())
source('src/initialize.R')


set.seed(444)
dats <- genSpatialOccData(alpha=0.8, sigma=5, delta=0.1)
model.input <- prepModData(dats$data, dats$y, dats$distance,
                           nsite=500, inits=dats$inits)


rm(list=ls())
source('src/initialize.R')

set.seed(444)
dats <- genSpatialOccData(alpha=0.5, sigma=5, delta=0.1)
model.input <- prepModData(dats$data, dats$y, dats$distance,
                           nsite=500, inits=dats$inits)


rm(list=ls())
source('src/initialize.R')

set.seed(444)
dats <- genSpatialOccData(alpha=0.1, sigma=5, delta=0.1)
model.input <- prepModData(dats$data, dats$y, dats$distance,
                           nsite=500, inits=dats$inits)
