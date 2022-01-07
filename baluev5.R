
# 
# 
## GOAL OF THE WORKSPACE:
##
##                             REDO SYSTEMATICALLY THE FAP MODELLING
##
##
## Change in the selection of data to which the quantile model is fitted:
## 3-wise division:
##  - for training,
##  - for model selection,
##  - for testing the performance.
## Points of the sky will be selected so that the minimal and the maximal number of observations
## are equally well covered, as well as the highest aliasing regions.

## DATAFRAMES HOLDING MOST OF THE RESULTS:
## Estimated location-wise parameters of test statistic distributions, and fractions of significant
## simulated sequences:
##       pars.train.dfr
##       pars.select.dfr
##       pars.test.dfr
## FAPs belonging to all maxima:
##       *fap.train.maxima
##       *fap.test.maxima
## where the * can be combinations of smoothing type and baluev, gev, quantile, and fadm,
## denoting the type of fap.

setwd("/Users/mariasuveges/Documents/SpectralAnalysis/FAP/feasibilityForGaia/OtherFAPs/")


# library(mapproj)
# library(gridBase)
# library(corrplot)
# library(rgl)
# library(numDeriv)
# library(mclust)
# library(evir)
# library(spatial)
library(fields)
library(randomForest)

WORKSPACE.NAME <- "baluev5.RData"
save.image(WORKSPACE.NAME)

colorvec <- c("magenta","powderblue","blueviolet","chocolate","yellow","green","red",
              "orange","violet","blue","deepskyblue","sienna","violetred",
              "aquamarine","royalblue","lightgreen","burlywood","coral","darkgoldenrod",
              "darkgreen","darkorchid","gray","indianred","lightpink","limegreen",
              "olivedrab","maroon","pink","purple","rosybrown","salmon","seagreen",
              "tomato","turquoise","grey")
 


## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------
## Loading previous results and functions
## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------

# # ## Citation from firstlookLine2.R (full-sky dataframe with classes and
# ## the GM models):
# # dfr <- max.power.fullsky.dfr
# # dfr$class7000_12group <- predict(gmxt7000[[12]],
    # # newdata = max.power.fullsky.dfr[, c(3:6,8:19)])$classification
# # dfr$class7000_13group <- predict(gmxt7000[[13]],
    # # newdata = max.power.fullsky.dfr[, c(3:6,8:19)])$classification
# # dfr$class7000_14group <- predict(gmxt7000[[14]],
    # # newdata = max.power.fullsky.dfr[, c(3:6,8:19)])$classification
# # dfr$class7000_15group <- predict(gmxt7000[[15]],
    # # newdata = max.power.fullsky.dfr[, c(3:6,8:19)])$classification
# # dfr$class7000_16group <- predict(gmxt7000[[16]],
    # # newdata = max.power.fullsky.dfr[, c(3:6,8:19)])$classification
# # dfr$class7000_17group <- predict(gmxt7000[[17]],
    # # newdata = max.power.fullsky.dfr[, c(3:6,8:19)])$classification
# # dfr$class7000_18group <- predict(gmxt7000[[18]],
    # # newdata = max.power.fullsky.dfr[, c(3:6,8:19)])$classification
# # dfr$class7000_19group <- predict(gmxt7000[[19]],
    # # newdata = max.power.fullsky.dfr[, c(3:6,8:19)])$classification
# # dfr$class7000_20group <- predict(gmxt7000[[20]],
    # # newdata = max.power.fullsky.dfr[, c(3:6,8:19)])$classification
# # max.power.fullsky7000 <- dfr
# # save(max.power.fullsky7000, file = "max.power.fullsky7000.RObj")
# # save(gmxt7000, file = "gmxt7000.RObj")

# ## The classifiers and the table with the fitted classes
# ## max.power.fullsky7000:
# load(file = "max.power.fullsky7000.RObj")
# ## gmxt7000:
# load(file = "gmxt7000.RObj")


## Necessary functions (coordinate transforms, if any:):
source("cooTransforms.R")


## Noise:

## Leanne's results for the random positions:
traindatadir.name <- "/Users/mariasuveges/Documents/SpectralAnalysis/FAP/feasibilityForGaia/NoiseSimulations011113/NoiseSimulationsRandomPositions/"
testdatadir.name <- "/Users/mariasuveges/Documents/SpectralAnalysis/FAP/feasibilityForGaia/NoiseSimulations011113/NoiseSimulationsRandomMagnitude_RandomPositions_SecondSet/"
## Leanne's results for the line:
datadir.name <- "/Users/mariasuveges/Documents/SpectralAnalysis/FAP/feasibilityForGaia/NoiseSimulations011113/NoiseSimulationsLineGridAll/"
## Leanne's results for the rectangle:
ecldir.name <- "/Users/mariasuveges/Documents/SpectralAnalysis/FAP/feasibilityForGaia/NoiseSimulations011113/NoiseSimulationsEclGridAll/"
newsimdir.name <- "/Users/mariasuveges/Documents/SpectralAnalysis/FAP/feasibilityForGaia/AllSimulations020114and030914/NoiseSimulations_RandomMagnitude_PositionsForNewSimulations20140730/"

## If necessary again:

# ## Time samplings and spectral windows for the whole sky (for backgrounds when necessary):
# roughgriddir.name <- "/Users/mariasuveges/Documents/SpectralAnalysis/AliasGaia/GaiaSpectralWindows/NominalScanningLaw230813/"
# ## Time samplings and spectral windows for the random positions:
# timesamplingdir.randompos.name <- "/Users/mariasuveges/Documents/SpectralAnalysis/AliasGaia/GaiaSpectralWindows/timesamplingsForFAP/"
# ## Time samplings and spectral windows for the line:
# finelinedir.name <- "/Users/mariasuveges/Documents/SpectralAnalysis/AliasGaia/GaiaSpectralWindows/NominalScanningLaw_LineGrid291013/"
# ## Time samplings and spectral windows for the region around the ecliptic:
# fineecldir.name <- "/Users/mariasuveges/Documents/SpectralAnalysis/AliasGaia/GaiaSpectralWindows/NominalScanningLaw_EclGrid291013/"

# ## For spectral windows of the random time samplings: taken from workspace 
# ## /Users/mariasuveges/Documents/SpectralAnalysis/AliasGaia/GaiaSpectralWindows/timesamplingsForFAP/spWin_RandomPositions.RData
# ## the saves were:
# ## save(max.ind.partpg.randompos, max.freq.partpg.randompos, max.power.partpg.randompos, 
# ##    file = "max_aliases_partpg_randompos.RObj")
# ## save(max.ind.fullpg.randompos, max.freq.fullpg.randompos, max.power.fullpg.randompos, 
# ##    file = "max_aliases_fullpg_randompos.RObj")
# load("/Users/mariasuveges/Documents/SpectralAnalysis/AliasGaia/GaiaSpectralWindows/timesamplingsForFAP/max_aliases_fullpg_randompos.RObj")

# ## For spectral windows of the line:
# ## max.ind.line, max.power.line, max.freq.line:
# load(paste(finelinedir.name, "max_aliases.RObj", sep = ""))

# ## For spectral windows of the ecliptic 45° region:
# ## max.ind.ecl, max.power.ecl, max.freq.ecl:
# load(paste(fineecldir.name, "max_aliases.RObj", sep = ""))

## Load the parameters files and the faps themselves from betamodelFitting1.RData
## (contains coordinates, number of obs, alias heights and locations, sums and means of
## spectral window features etc.)
## Objects: pars.ecl.dfr, pars.line.dfr, pars.randompos.dfr,
##    fap.ecl.dfr, fap.line.dfr, fap.randompos.dfr, fap.randompos.check.dfr
load("results_betamodelFitting1/pars.whatever.dfr.RObj")

## Load the similar values from .... for the new simulations, done après coup to check the model
## performances at the extremities of the covariate distribution (high aliasing, low nb of obs)
## (workspace: )
load("results_spWin_NewSimulations20140730/pars.newrand.dfr.RObj")
colnames(pars.ecl.dfr)
colnames(pars.line.dfr)
colnames(pars.randompos.dfr)
colnames(pars.newrand.dfr)
## Give the region name for all regions, and arrange the columns so that I can rbind the frames
## (the new simulations already contain the region indicator):
pars.ecl.dfr$region <- "ecl"
pars.line.dfr$region <- "line"
pars.randompos.dfr$region <- "randompos"
cols.to.use <- c("name","ra","dec","lambda","beta","n","region","at4" ,"at8","at12","at16",
   "at20","at24","at28","at32","at36","at40","at44","at48","at52","at56","at60","at64","at68",              
   "freq.at4","freq.at8","freq.at12","freq.at16","freq.at20","freq.at24","freq.at28","freq.at32",      
   "freq.at36","freq.at40","freq.at44","freq.at48","freq.at52","freq.at56","freq.at60",
   "freq.at64","freq.at68","vart","sum.spwin","highest.spwin",
   "highest2sum.spwin","highest3sum.spwin","highest4sum.spwin","highest5sum.spwin","highest6sum.spwin",
   "highest7sum.spwin","highest8sum.spwin","highest9sum.spwin","sum.of.five.spwin","sum.of.six.spwin",  
   "sum.of.seven.spwin","sum.of.eight.spwin","sum.of.nine.spwin","sum.of.ten.spwin")
pars.all.dfr <- rbind(pars.ecl.dfr[, cols.to.use], pars.line.dfr[, cols.to.use],
   pars.randompos.dfr[, cols.to.use], pars.newrand.dfr[, cols.to.use])
rm(pars.ecl.dfr, pars.line.dfr, pars.randompos.dfr, pars.newrand.dfr)
pars.all.dfr$name <- as.character(pars.all.dfr$name)


## I forgot to change the names of the new simulations to be consistent with the former naming
## (coordinates put together, plus a _gaia extension):
n.tmp <- pars.all.dfr$name[is.element(pars.all.dfr$region, c("lmc","ln.ha"))] 
str <- n.tmp[1]
nameconv.fun <-  function(str)
   {
  	c.tmp <- unlist(strsplit(str, split = "_"))
  	c.tmp[3] <- "gaia"
  	paste(c.tmp[1],c.tmp[2],c.tmp[3], sep = "_")
   }
m.tmp <- sapply(pars.all.dfr$name, nameconv.fun)
pars.all.dfr$name <- m.tmp
rm(n.tmp, m.tmp, c.tmp, str)


## --------------------------------------------------------------------------------------------------------
## Dataframes of all the periodogram maxima of all the simulations
## --------------------------------------------------------------------------------------------------------

## Random positions, first half of the simulations:

randomposdirlist <- list.files(traindatadir.name)
maxima.randompos.dfr <- as.data.frame(matrix(ncol = length(randomposdirlist), nrow = 750))
colnames(maxima.randompos.dfr) <- randomposdirlist
for(ii in 1:length(randomposdirlist))
   {
#    ii <- 2
	dd <- randomposdirlist[ii]
	pp <- list.files(paste(traindatadir.name, dd, sep = ""), full.names = TRUE,
	   pattern = ".dat")
	maxima.randompos.dfr[,dd] <- read.table(pp)[,5]
   }

## Random positions, second half of the simulations:

randomposcheckdirlist <- list.files(testdatadir.name)
identical(randomposcheckdirlist, randomposdirlist)
## All right, TRUE.
## Allow a way to check the number of observations, are they also the same?
n.checkset <- numeric(length(randomposcheckdirlist))
## the number of observations are in the commented-out header in the
##  files:
maxima.randompos.check.dfr <- as.data.frame(matrix(ncol = length(randomposcheckdirlist), nrow = 750))
colnames(maxima.randompos.check.dfr) <- randomposcheckdirlist
for(ii in c(1:67,69:200,202:332,334:465,467:591,593:711,713:720))
   {
#    ii <- 68
	dd <- randomposcheckdirlist[ii]
	pp <- list.files(paste0(testdatadir.name, dd), full.names = TRUE, pattern = ".dat")
    cc <- unlist(strsplit(readLines(pp, n = 2), " "))
    n.checkset[ii] <- as.numeric(cc[length(cc)])
	maxima.randompos.check.dfr[,dd] <- read.table(pp)[,5]
   }
## files unfinished: 68,201,333,466,592,712
identical(n.checkset[-c(68,201,333,466,592,712)], pars.randompos.dfr$n[-c(68,201,333,466,592,712)])
## Ok, true apart from the unfinished files.

## Line:

linedirlist <- list.files(datadir.name)
maxima.line.dfr <- as.data.frame(matrix(ncol = length(linedirlist), nrow = 1500))
colnames(maxima.line.dfr) <- linedirlist
for(ii in 1:length(linedirlist))
   {
#    ii <- 2
	dd <- linedirlist[ii]
	pp <- list.files(paste(datadir.name, dd, sep = ""), full.names = TRUE,
	   pattern = ".dat")
	maxima.line.dfr[,dd] <- read.table(pp)[,5]
   }
# maxima.line.dfr[1:20, 1:20]

## Rectangle:

ecldirlist <- list.files(ecldir.name, pattern = "_gaia")
maxima.ecl.dfr <- as.data.frame(matrix(ncol = length(ecldirlist), nrow = 1500))
colnames(maxima.ecl.dfr) <- ecldirlist
for(ii in 1:length(ecldirlist))
   {
#    ii <- 2
	dd <- ecldirlist[ii]
	pp <- list.files(paste(ecldir.name, dd, sep = ""), full.names = TRUE,
	   pattern = ".dat")
	maxima.ecl.dfr[,dd] <- read.table(pp)[,5]
   }

## Low-nb and high-alias sites:

newsimdirlist <- sapply(list.files(newsimdir.name, pattern = "_gaia"),  function(str) 
   substr(str, start = 18, stop = nchar(str) - 4))
newsimdirlist.long <- list.files(newsimdir.name, pattern = "_gaia", full.names = TRUE)
maxima.newrand.dfr <- as.data.frame(matrix(ncol = length(newsimdirlist), nrow = 1500))
colnames(maxima.newrand.dfr) <- newsimdirlist
for(ii in 1:length(newsimdirlist.long))
   {
#    ii <- 2
	dd <- newsimdirlist.long[ii]
	nn <- newsimdirlist[ii]
	maxima.newrand.dfr[,nn] <- read.table(dd)[,5]
   }
maxima.newrand.dfr[1:10,1:20]
## Looks ok.


## --------------------------------------------------------------------------------------------------------
## Selection of training, selection and test points
## --------------------------------------------------------------------------------------------------------

range(pars.all.dfr$n[pars.all.dfr$region == "randompos"])
range(pars.all.dfr$n[pars.all.dfr$region == "line"])
range(pars.all.dfr$n[pars.all.dfr$region == "ecl"])
range(pars.all.dfr$n[pars.all.dfr$region == "lmc"])
range(pars.all.dfr$n[pars.all.dfr$region == "ln.ha"])

## A rough distribution of the number of time series points:
allskydir <- "/Users/mariasuveges/Documents/SpectralAnalysis/AliasGaia/GaiaSpectralWindows/NominalScanningLaw230813/allsky"
allskydirlist <- list.files(allskydir, pattern = ".0", full.names = TRUE)
nb.obs.allsky <- data.frame(ra = NA, dec = NA, n = NA)
nb.obs.allsky <- nb.obs[FALSE,]
for(ii in 1:length(allskydirlist))
   {
    d.tmp <- list.files(allskydirlist[ii], pattern = "Gaia")
    m.tmp <- sapply(d.tmp, function(chr) 
       {
       	c(ra = as.numeric(unlist(strsplit(chr, split = "_"))[2]),
       	  dec = as.numeric(unlist(strsplit(chr, split = "_"))[3]),
       	  n = as.numeric(unlist(strsplit(unlist(strsplit(chr, split = "_"))[4], split = "\\."))[1]))
       })
    nb.obs.allsky <- rbind(nb.obs.allsky, data.frame(t(m.tmp)))
   }
colnames(nb.obs.allsky) <- c("ra","dec","n")
range(nb.obs.allsky$n)
## [1]  42 243

## Add part of the ecl to the line
## All the random positions, the remaining part of ecl and the new simulations will be in 
## the test set.

## Selection such that a stripe from the rectangle region is included into the training set:
which.ecl <- sort(unique(c(which(abs(pars.all.dfr$ra - 269) <= 0.4 & pars.all.dfr$region == "ecl"))))
## Plus, the whole of the line.
which.line <- sort(which(pars.all.dfr$region == "line"))

pars.all.dfr$set.ind <- NA
pars.all.dfr$set.ind[c(which.line, which.ecl)] <- "train"
 
nrow(pars.all.dfr) - length(c(which.ecl, which.line))
## Divide the non-training part into model selection part and testing part:

which.select <- sort(sample(which(is.na(pars.all.dfr$set.ind)), size = 1283))
pars.all.dfr$set.ind[which.select] <- "select"
pars.all.dfr$set.ind[-c(which.ecl, which.line, which.select)] <- "test"
## All right.

pars.train.dfr <- pars.all.dfr[pars.all.dfr$set.ind == "train",]
pars.select.dfr <- pars.all.dfr[pars.all.dfr$set.ind == "select",]
pars.test.dfr <- pars.all.dfr[pars.all.dfr$set.ind == "test",]

dfr <- cbind(maxima.ecl.dfr, maxima.line.dfr, rbind(maxima.randompos.dfr,
   maxima.randompos.check.dfr), maxima.newrand.dfr)
identical(sort(colnames(dfr)), sort(as.character(pars.all.dfr$name)))
## TRUE, all right.
names.to.omit <- colnames(dfr)[apply(dfr, 2, function(vec) any(is.na(vec)))]
maxima.train.dfr <- dfr[, setdiff(pars.all.dfr$name[pars.all.dfr$set.ind == "train"],
   names.to.omit)]
maxima.select.dfr <- dfr[, setdiff(pars.all.dfr$name[pars.all.dfr$set.ind == "select"],
   names.to.omit)] 
maxima.test.dfr <- dfr[, setdiff(pars.all.dfr$name[pars.all.dfr$set.ind == "test"],
   names.to.omit)]
dim(maxima.train.dfr)
dim(maxima.select.dfr)
dim(maxima.test.dfr)
## Dimensions sound ok.

pars.select.dfr <- pars.select.dfr[!is.element(pars.select.dfr$name, names.to.omit),]
pars.test.dfr <- pars.test.dfr[!is.element(pars.test.dfr$name, names.to.omit),]
pars.all.dfr <- pars.all.dfr[!is.element(pars.all.dfr$name, names.to.omit),]


## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------
## Quantiles
## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------

pars.train.dfr$q95 <- apply(maxima.train.dfr, M = 2, quantile,
    na.rm = TRUE, prob = 0.95)
pars.train.dfr$q99 <- apply(maxima.train.dfr, M = 2, quantile,
    na.rm = TRUE, prob = 0.99)
pars.test.dfr$q95 <- apply(maxima.test.dfr, M = 2, quantile,
    na.rm = TRUE, prob = 0.95)
pars.test.dfr$q99 <- apply(maxima.test.dfr, M = 2, quantile,
    na.rm = TRUE, prob = 0.99)
pars.select.dfr$q95 <- apply(maxima.select.dfr, M = 2, quantile,
    na.rm = TRUE, prob = 0.95)
pars.select.dfr$q99 <- apply(maxima.select.dfr, M = 2, quantile,
    na.rm = TRUE, prob = 0.99)



## --------------------------------------------------------------------------------------------------------
## A plot about the quantile that gives 5% and 1% false positives versus some variables 
## of the locations
## --------------------------------------------------------------------------------------------------------

quartz(height = 7.5, width = 12)
par(mfrow = c(4,5), mar = c(2.5,2.5,2,1), mgp = c(1.6,0.6,0))
#for(ii in c(6,8:26))
for(ii in 42:58)
   {
	plot(pars.train.dfr[,ii], pars.train.dfr$q95,
	   pch = 16, cex = 0.5, main = colnames(pars.train.dfr)[ii], type = "p", ylim = c(0.1,0.6))
	points(pars.select.dfr[,ii],  pars.select.dfr$q95,
	   pch = 17, cex = 0.5, col = "orange")
	points(pars.test.dfr[,ii],  pars.test.dfr$q95,
	   pch = 17, cex = 0.5, col = "deepskyblue")
	points(pars.train.dfr[,ii],  pars.train.dfr$q95,
	   pch = 1, cex = 0.5)
	points(pars.select.dfr[,ii],  pars.select.dfr$q99,
	   pch = 17, cex = 0.5, col = "gold")
	points(pars.test.dfr[,ii],  pars.test.dfr$q99,
	   pch = 17, cex = 0.5, col = "aquamarine")
	points(pars.train.dfr[,ii],  pars.train.dfr$q99,
	   pch = 1, cex = 0.5, col = "grey")
   }
## Fantastic, the addition of the new points does simply not change anything on the aspect
## of the graphs until now... though at aliases not very dominant so far (20,48,60), there is
## apparently a new cloud of points (only yellowish and bluish, not from the greyish training
## set which is composed only of the old ecl and line point sets).
 

## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------
## Baluev FAP 
## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------


## --------------------------------------------------------------------------------------------------------
## Functions
## --------------------------------------------------------------------------------------------------------

baluevfap.epoch.fun <- function(z, t, fmax)
   {
   	N <- length(t)
   	min(1, (1-z)^((N-3)/2) + 
   	(gamma((N-1)/2) / gamma((N-2)/2)) *
   	   sqrt(4*pi*var(t)) * fmax * (1-z)^((N-4)/2) * sqrt(z))
    }

baluevfap.vart.fun <- function(z, vart, fmax, N)
   {
   	min((1-z)^((N-3)/2) + 
   	(gamma((N-1)/2) / gamma((N-2)/2)) *
   	   sqrt(4*pi*vart) * fmax * (1-z)^((N-4)/2) * sqrt(z), 1)
    }


## --------------------------------------------------------------------------------------------------------
## Application
## --------------------------------------------------------------------------------------------------------

baluevfap.train.dfr <- data.frame(matrix(ncol = ncol(maxima.train.dfr),
    nrow =  nrow(maxima.train.dfr)))
colnames(baluevfap.train.dfr) <- colnames(maxima.train.dfr)
for(ii in 1:ncol(maxima.train.dfr))
   {
   	nn <- colnames(maxima.train.dfr)[ii]
   	baluevfap.train.dfr[,nn] <- sapply(maxima.train.dfr[,nn], baluevfap.vart.fun,
   	    fmax = 30,
   	    N = pars.train.dfr$n[pars.train.dfr$name == nn],
   	    vart = pars.train.dfr$vart[pars.train.dfr$name == nn])
   }

baluevfap.select.dfr <- data.frame(matrix(ncol = ncol(maxima.select.dfr),
    nrow =  nrow(maxima.select.dfr)))
colnames(baluevfap.select.dfr) <- colnames(maxima.select.dfr)
for(ii in 1:ncol(maxima.select.dfr))
   {
   	nn <- colnames(maxima.select.dfr)[ii]
   	baluevfap.select.dfr[,nn] <- sapply(maxima.select.dfr[,nn], baluevfap.vart.fun,
   	    fmax = 30,
   	    N = pars.select.dfr$n[pars.select.dfr$name == nn],
   	    vart = pars.select.dfr$vart[pars.select.dfr$name == nn])
   }

baluevfap.test.dfr <- data.frame(matrix(ncol = ncol(maxima.test.dfr),
    nrow =  nrow(maxima.test.dfr)))
colnames(baluevfap.test.dfr) <- colnames(maxima.test.dfr)
for(ii in 1:ncol(maxima.test.dfr))
   {
   	nn <- colnames(maxima.test.dfr)[ii]
   	baluevfap.test.dfr[,nn] <- sapply(maxima.test.dfr[,nn], baluevfap.vart.fun,
   	    fmax = 30,
   	    N = pars.test.dfr$n[pars.test.dfr$name == nn],
   	    vart = pars.test.dfr$vart[pars.test.dfr$name == nn])
   }
save(baluevfap.test.dfr, baluevfap.train.dfr, baluevfap.select.dfr, file = 
   "results_baluev5/baluevfap.all.RObj")
   

## --------------------------------------------------------------------------------------------------------
## Check the number of significant values in the test set
## --------------------------------------------------------------------------------------------------------

pars.train.dfr$fr.baluev.sign05 <- NA
pars.train.dfr$fr.baluev.sign01 <- NA
pars.test.dfr$fr.baluev.sign05 <- NA
pars.test.dfr$fr.baluev.sign01 <- NA
pars.select.dfr$fr.baluev.sign05 <- NA
pars.select.dfr$fr.baluev.sign01 <- NA

## Compute these success rates only on the first 1000 simulations (since F^M and GEV will
## need the other 500 for training)
for(ii in 1:nrow(pars.train.dfr))
   {
   	pars.train.dfr$fr.baluev.sign05[ii] <-
   	   sum(baluevfap.train.dfr[, as.character(pars.train.dfr$name[ii])] < 0.05) / 1500
   	pars.train.dfr$fr.baluev.sign01[ii] <-
   	   sum(baluevfap.train.dfr[, as.character(pars.train.dfr$name[ii])] < 0.01) / 1500
   }
for(ii in 1:nrow(pars.test.dfr))
   {
   	pars.test.dfr$fr.baluev.sign05[ii] <-
   	   sum(baluevfap.test.dfr[, as.character(pars.test.dfr$name[ii])] < 0.05) / 1500
   	pars.test.dfr$fr.baluev.sign01[ii] <-
   	   sum(baluevfap.test.dfr[, as.character(pars.test.dfr$name[ii])] < 0.01) / 1500
   }
for(ii in 1:nrow(pars.select.dfr))
   {
   	pars.select.dfr$fr.baluev.sign05[ii] <-
   	   sum(baluevfap.select.dfr[, as.character(pars.select.dfr$name[ii])] < 0.05) / 1500
   	pars.select.dfr$fr.baluev.sign01[ii] <-
   	   sum(baluevfap.select.dfr[, as.character(pars.select.dfr$name[ii])] < 0.01) / 1500
   }



quartz(height = 7, width = 10)
par(mfrow = c(1,1), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.6,0))

col.x <- "n"
plot(abs(pars.test.dfr[pars.test.dfr$region == "randompos", col.x]), 
   pars.test.dfr$fr.baluev.sign05[pars.test.dfr$region == "randompos"], 
   type = "p", pch = 16, cex = 0.7, ylim = c(0,0.15), 
   xlim = range(abs(pars.test.dfr[,col.x])), main = "Baluev FAP",
   ylab = "Fraction above threshold", xlab = col.x, col = "grey25")
points(abs(pars.test.dfr[pars.test.dfr$region == "randompos", col.x]),
   pars.test.dfr$fr.baluev.sign01[pars.test.dfr$region == "randompos"], 
   type = "p", pch = 16, col = "grey", cex = 0.7)

points(abs(pars.test.dfr[pars.test.dfr$region == "ecl", col.x]), 
   pars.test.dfr$fr.baluev.sign05[pars.test.dfr$region == "ecl"], 
    pch = 16, cex = 0.7, col = "violetred")
points(abs(pars.test.dfr[pars.test.dfr$region == "ecl", col.x]), 
   pars.test.dfr$fr.baluev.sign01[pars.test.dfr$region == "ecl"], 
    pch = 16, cex = 0.7, col = "pink")

points(abs(pars.test.dfr[pars.test.dfr$region == "lmc", col.x]), 
   pars.test.dfr$fr.baluev.sign05[pars.test.dfr$region == "lmc"], 
    pch = 16, cex = 0.7, col = "orange")
points(abs(pars.test.dfr[pars.test.dfr$region == "lmc", col.x]), 
   pars.test.dfr$fr.baluev.sign01[pars.test.dfr$region == "lmc"], 
    pch = 16, cex = 0.7, col = "gold")

points(abs(pars.test.dfr[pars.test.dfr$region == "ln.ha", col.x]), 
   pars.test.dfr$fr.baluev.sign05[pars.test.dfr$region == "ln.ha"], 
    pch = 16, cex = 0.7, col = "deepskyblue")
points(abs(pars.test.dfr[pars.test.dfr$region == "ln.ha", col.x]), 
   pars.test.dfr$fr.baluev.sign01[pars.test.dfr$region == "ln.ha"], 
    pch = 16, cex = 0.7, col = "skyblue")

abline(h = 0.05, lwd = 2)
abline(h = 0.01, lty = 2, lwd = 2)
## This is excellent, no need for better. We will have a bit smaller false alarm
## rate, and so a loss of some variable objects, but so what?

rm(baluevfap.test.dfr, baluevfap.train.dfr, baluevfap.select.dfr)


## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------
## Paltani--Schwarzenberg-Czerny FAP 
## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------


## --------------------------------------------------------------------------------------------------------
## Functions
## --------------------------------------------------------------------------------------------------------

mest.fun <- function(zs, N)
    log(0.5) / log(pbeta(median(zs), 1, (N-3)/2, lower.tail = TRUE)) 

schwczfap.fun <- function(z, mest, N)
    1-pbeta(z, 1, (N-3)/2, lower.tail = TRUE)^mest


## --------------------------------------------------------------------------------------------------------
## Application
## --------------------------------------------------------------------------------------------------------

## Fit the F^M distribution to all maxima at the training set locations
## Workspace baluev1.RData: no way to use estimates based on 500 simulations,
## it is scattered between 0.2 -- 0.13. 

## Estimate M based on 1500 simulations in the training set (and for comparison,
## on the test set too.)  

## For the training set:

for(ii in 1:nrow(pars.train.dfr))
#  for(ii in 5)
   {
   	nn <- as.character(pars.train.dfr$name[ii])
   	nobs <- pars.train.dfr$n[pars.train.dfr$name == nn]
    mest <- mest.fun(zs = maxima.train.dfr[,nn], N = nobs)
    pars.train.dfr[pars.train.dfr$name == nn, "M.full"] <- mest
   }

## For the selection set:

for(ii in 1:nrow(pars.select.dfr))
#  for(ii in 5)
   {
   	nn <- as.character(pars.select.dfr$name[ii])
   	nobs <- pars.select.dfr$n[pars.select.dfr$name == nn]
    mest <- mest.fun(zs = maxima.select.dfr[,nn], N = nobs)
    pars.select.dfr[pars.select.dfr$name == nn, "M.full"] <- mest
   }

## For the test set:

for(ii in 1:nrow(pars.test.dfr))
#  for(ii in 5)
   {
   	nn <- as.character(pars.test.dfr$name[ii])
   	nobs <- pars.test.dfr$n[pars.test.dfr$name == nn]
    mest <- mest.fun(zs = maxima.test.dfr[,nn], N = nobs)
    pars.test.dfr[pars.test.dfr$name == nn, "M.full"] <- mest
   }


## --------------------------------------------------------------------------------------------------------
## A plot about the dependence of M on other parameters
## --------------------------------------------------------------------------------------------------------

quartz(height = 7.5, width = 12)
par(mfrow = c(4,5), mar = c(2.5,2.5,2,1), mgp = c(1.6,0.6,0))
for(ii in c(6,8:24))
#for(ii in 42:58)
   {
	plot(pars.train.dfr[,ii], pars.train.dfr$M.full,
	   pch = 16, cex = 0.5, main = colnames(pars.train.dfr)[ii], type = "p", 
	   ylim = range(c(pars.train.dfr$M.full,pars.test.dfr$M.full,pars.select.dfr$M.full)))
	points(pars.select.dfr[,ii],  pars.select.dfr$M.full,
	   pch = 17, cex = 0.5, col = "orange")
	points(pars.test.dfr[,ii],  pars.test.dfr$M.full,
	   pch = 17, cex = 0.5, col = "deepskyblue")
	points(pars.train.dfr[,ii],  pars.train.dfr$M.full,
	   pch = 1, cex = 0.5)
   }


## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------
## GEV modelling
## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------

## --------------------------------------------------------------------------------------------------------
## Functions
## --------------------------------------------------------------------------------------------------------

## Constraint: the upper limit of the distribution is known to be 1 
boundaryconstr.gev.fun <- function(pars, z)
   {
   	m <- length(z)
   	xi <- pars[1]
   	sig <- pars[2]
   	locsc <- 1 + xi*(z - 1 - sig/xi) / sig
   	locsc[locsc <= 0] <- 1e-120
   	- m*log(sig) - (1 + 1/xi) * sum(log(locsc)) - sum(locsc^(-1/xi))
   }
   
 
## Use all 1500 maxima on the training sample.
## For the training set:

dfr.tmp <- data.frame(matrix(ncol = 7, nrow = nrow(pars.train.dfr)))
colnames(dfr.tmp) <- c("name","xi.full","sig.full","mu.full","xise.full","sigse.full","muse.full")
dfr.tmp$name <- pars.train.dfr$name

for(ii in 1:nrow(pars.train.dfr))
#  for(ii in 5)
   {
   	nn <- as.character(pars.train.dfr$name[ii])
    zest <- maxima.train.dfr[, nn]
    r.tmp <- try(optim(c(-0.01, sqrt(6*var(zest))/pi), boundaryconstr.gev.fun,
	   z = zest, control = list(fnscale = -1, maxit = 10000),
	   hessian = TRUE, method = "Nelder-Mead"))
	s.tmp <- sqrt(diag(-solve(r.tmp$hessian)))
	gradmu.tmp <- c(- r.tmp$par[2]/r.tmp$par[1]^2, 1/r.tmp$par[1])
    dfr.tmp[dfr.tmp$name == nn, 2:7] <-
       if(inherits(r.tmp, "try-error")) {
          rep(NA, 6)
       } else {
          c(r.tmp$par[1], r.tmp$par[2], 1+r.tmp$par[2]/r.tmp$par[1],
#            s.tmp, -s.tmp[1]*r.tmp$par[2]/r.tmp$par[1]^2 + s.tmp[2]/r.tmp$par[1])
            s.tmp, sqrt(-t(gradmu.tmp) %*% solve(r.tmp$hessian) %*% t(t(gradmu.tmp))))
       }
   }

dfr1 <- merge(pars.train.dfr, dfr.tmp)
pars.train.dfr <- dfr1

## How many impossible values (above or below an estimated endoint of distribution)?
sum(apply(dfr1[,68:70], 1, function(vec) any(is.na(vec))))
## all ses could be computed, sounds great


## For the test set (only to check whether the training set is representative):

dfr.tmp <- data.frame(matrix(ncol = 7, nrow = nrow(pars.test.dfr)))
colnames(dfr.tmp) <- c("name","xi.full","sig.full","mu.full","xise.full","sigse.full","muse.full")
dfr.tmp$name <- pars.test.dfr$name

for(ii in 1:nrow(pars.test.dfr))
#  for(ii in 5)
   {
   	nn <- as.character(pars.test.dfr$name[ii])
    zest <- maxima.test.dfr[, nn]
    r.tmp <- try(optim(c(-0.01, sqrt(6*var(zest))/pi), boundaryconstr.gev.fun,
	   z = zest, control = list(fnscale = -1, maxit = 10000),
	   hessian = TRUE, method = "Nelder-Mead"))
	if(inherits(r.tmp, "try-error")) {s.tmp <- NA} else {s.tmp <- sqrt(diag(-solve(r.tmp$hessian)))}
	if(inherits(r.tmp, "try-error")) {gradmu.tmp <- rep(NA, 2)} else {
		gradmu.tmp <- c(- r.tmp$par[2]/r.tmp$par[1]^2, 1/r.tmp$par[1])
		}
    dfr.tmp[dfr.tmp$name == nn, 2:7] <-
       if(inherits(r.tmp, "try-error")) {
          rep(NA, 6)
       } else {
          c(r.tmp$par[1], r.tmp$par[2], 1+r.tmp$par[2]/r.tmp$par[1],
            s.tmp, sqrt(-t(gradmu.tmp) %*% solve(r.tmp$hessian) %*% t(t(gradmu.tmp))))
       }
   }

dfr0 <- merge(pars.test.dfr, dfr.tmp)
## How many impossible values (above or below an estimated endoint of distribution)?
sum(apply(dfr0[,68:70], 1, function(vec) any(is.na(vec))))
## 0, I succeeded to remove them. OK.

pars.test.dfr <- dfr0

## For the model selection set

dfr.tmp <- data.frame(matrix(ncol = 7, nrow = nrow(pars.select.dfr)))
colnames(dfr.tmp) <- c("name","xi.full","sig.full","mu.full","xise.full","sigse.full","muse.full")
dfr.tmp$name <- pars.select.dfr$name

for(ii in 1:nrow(pars.select.dfr))
#  for(ii in 5)
   {
   	nn <- as.character(pars.select.dfr$name[ii])
    zest <- maxima.select.dfr[, nn]
    r.tmp <- try(optim(c(-0.01, sqrt(6*var(zest))/pi), boundaryconstr.gev.fun,
	   z = zest, control = list(fnscale = -1, maxit = 10000),
	   hessian = TRUE, method = "Nelder-Mead"))
	if(inherits(r.tmp, "try-error")) {s.tmp <- NA} else {s.tmp <- sqrt(diag(-solve(r.tmp$hessian)))}
	if(inherits(r.tmp, "try-error")) {gradmu.tmp <- rep(NA, 2)} else {
		gradmu.tmp <- c(- r.tmp$par[2]/r.tmp$par[1]^2, 1/r.tmp$par[1])
		}
    dfr.tmp[dfr.tmp$name == nn, 2:7] <-
       if(inherits(r.tmp, "try-error")) {
          rep(NA, 6)
       } else {
          c(r.tmp$par[1], r.tmp$par[2], 1+r.tmp$par[2]/r.tmp$par[1],
            s.tmp, sqrt(-t(gradmu.tmp) %*% solve(r.tmp$hessian) %*% t(t(gradmu.tmp))))
       }
   }

dfr0 <- merge(pars.select.dfr, dfr.tmp)
## How many impossible values (above or below an estimated endoint of distribution)?
sum(apply(dfr0[,48:50], 1, function(vec) any(is.na(vec))))
## 0, ok.

pars.select.dfr <- dfr0
rm(dfr0,dfr1)


## --------------------------------------------------------------------------------------------------------
## A plot about the dependence of GEV pars on other parameters
## --------------------------------------------------------------------------------------------------------

quartz(height = 7.5, width = 12)
par(mfrow = c(4,5), mar = c(2.5,2.5,2,1), mgp = c(1.6,0.6,0))
#for(ii in c(6,8:24))
for(ii in 42:58)
   {
	plot(pars.train.dfr[,ii], pars.train.dfr$xi.full,
	   pch = 16, cex = 0.5, main = colnames(pars.train.dfr)[ii], type = "p", 
	   ylim = range(c(pars.train.dfr$xi.full,pars.test.dfr$xi.full,pars.select.dfr$xi.full)))
	points(pars.select.dfr[,ii],  pars.select.dfr$xi.full,
	   pch = 17, cex = 0.5, col = "orange")
	points(pars.test.dfr[,ii],  pars.test.dfr$xi.full,
	   pch = 17, cex = 0.5, col = "deepskyblue")
	points(pars.train.dfr[,ii],  pars.train.dfr$xi.full,
	   pch = 1, cex = 0.5)
   }

quartz(height = 7.5, width = 12)
par(mfrow = c(4,5), mar = c(2.5,2.5,2,1), mgp = c(1.6,0.6,0))
#for(ii in c(6,8:24))
for(ii in 42:58)
   {
	plot(pars.train.dfr[,ii], pars.train.dfr$sig.full,
	   pch = 16, cex = 0.5, main = colnames(pars.train.dfr)[ii], type = "p", 
	   ylim = range(c(pars.train.dfr$sig.full,pars.test.dfr$sig.full,pars.select.dfr$sig.full)))
	points(pars.select.dfr[,ii],  pars.select.dfr$sig.full,
	   pch = 17, cex = 0.5, col = "orange")
	points(pars.test.dfr[,ii],  pars.test.dfr$sig.full,
	   pch = 17, cex = 0.5, col = "deepskyblue")
	points(pars.train.dfr[,ii],  pars.train.dfr$sig.full,
	   pch = 1, cex = 0.5)
   }


## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------
## Random Forest regression 
## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------

## --------------------------------------------------------------------------------------------------------
## For M.full
## --------------------------------------------------------------------------------------------------------

## Initial model with all the variables

mfull.start.rf <- randomForest(x = pars.train.dfr[, start.attr], y = pars.train.dfr$M.full,
   xtest = pars.select.dfr[!is.na(pars.select.dfr$M.full), start.attr], 
   ytest = pars.select.dfr[!is.na(pars.select.dfr$M.full), "M.full"],
   importance = TRUE, ntree = 1000, keep.forest = TRUE)

## Look at the results: the mean squared error and pseudo-R^2 on the training set
mfull.start.rf$mse[1000]   ## 75032133 at the end
mfull.start.rf$rsq[1000]   ## 0.8875789

## Look at the results: the mean squared error and pseudo-R^2 on the test set
mfull.start.rf$test$mse    ## 117766220 at the end, bigger
mfull.start.rf$test$rsq    ## ~ 0.8352844 at the end, smaller
## Some deterioration.

## Attribute importance
mfull.start.rf$importance[order(mfull.start.rf$importance[,1]),]

quartz(height = 7.5, width = 7.5)
par(mfcol = c(1,1), mar = c(2.5,10,2,1), mgp = c(1.5,0.6,0))
o.tmp <- order(mfull.start.rf$importance[,1])
barplot(mfull.start.rf$importance[,1][o.tmp], horiz = TRUE, las = 1,
   names.arg = names(mfull.start.rf$importance[,1][o.tmp]), cex.axis = 0.5)

## Based on this, define a restricted attr.set:
mfull.attr1 <- c("sum.of.nine.spwin","highest3sum.spwin","n","vart","sum.of.five.spwin",
  "sum.of.eight.spwin","highest9sum.spwin")

mfull.rf1 <- randomForest(x = pars.train.dfr[, mfull.attr1], y = pars.train.dfr$M.full,
   xtest = pars.select.dfr[!is.na(pars.select.dfr$M.full), mfull.attr1], 
   ytest = pars.select.dfr[!is.na(pars.select.dfr$M.full), "M.full"],
   importance = TRUE, ntree = 1000, keep.forest = TRUE)

## Look at the results: the mean squared error and pseudo-R^2 on the training set
mfull.rf1$mse[1000]   ## 80220732, bigger
mfull.rf1$rsq[1000]   ## 0.8798048, somewhat smaller than full-attribute model

## Look at the results: the mean squared error and pseudo-R^2 on the test set
mfull.rf1$test$mse[1000]    ## 100683448 better than for the full-attribute model
mfull.rf1$test$rsq[1000]    ## 0.8591775, worse than for the full-attribute model

## Attribute importance

quartz(height = 2.5, width = 7.5)
par(mfcol = c(1,1), mar = c(2.5,10,2,1), mgp = c(1.5,0.6,0))
o.tmp <- order(mfull.rf1$importance[,1])
barplot(mfull.rf1$importance[,1][o.tmp], horiz = TRUE, las = 1,
   names.arg = names(mfull.rf1$importance[,1][o.tmp]), cex.axis = 0.5)

## Based on this, check systematically the omission of the variables one by one:

M.rf2var.matr <- matrix(ncol = 4, nrow = length(mfull.attr1)-2)
for(ii in 3:(length(mfull.attr1)))
#for(ii in 1:2)
   {
	r.tmp <- randomForest(x = pars.train.dfr[, mfull.attr1[-ii]], 
	   y = pars.train.dfr$M.full,
	   xtest = pars.select.dfr[, mfull.attr1[-ii]], 
	   ytest = pars.select.dfr[, "M.full"],
	   importance = FALSE, ntree = 1000)
	M.rf2var.matr[ii-2,] <- c(r.tmp$rsq[1000], r.tmp$test$rsq[1000], r.tmp$mse[1000], 
	   r.tmp$test$mse[1000])
   }
## Worst: sum.of.eight.spwin. Omit, and repeat the loop.
mfull.attr2 <- c("sum.of.nine.spwin","highest3sum.spwin","n","vart","sum.of.five.spwin",
  "highest9sum.spwin")
M.rf3var.matr <- matrix(ncol = 4, nrow = length(mfull.attr2)-2)
for(ii in 3:(length(mfull.attr2)))
#for(ii in 1:2)
   {
	r.tmp <- randomForest(x = pars.train.dfr[, mfull.attr2[-ii]], 
	   y = pars.train.dfr$M.full,
	   xtest = pars.select.dfr[, mfull.attr2[-ii]], 
	   ytest = pars.select.dfr[, "M.full"],
	   importance = FALSE, ntree = 1000)
	M.rf3var.matr[ii-2,] <- c(r.tmp$rsq[1000], r.tmp$test$rsq[1000], r.tmp$mse[1000], 
	   r.tmp$test$mse[1000])
   }
## Worst: highest9sum.spwin. Omit, and repeat the loop.
mfull.attr3 <- c("sum.of.nine.spwin","highest3sum.spwin","n","vart","sum.of.five.spwin")
M.rf4var.matr <- matrix(ncol = 4, nrow = length(mfull.attr3)-2)
for(ii in 3:(length(mfull.attr3)))
#for(ii in 1:2)
   {
	r.tmp <- randomForest(x = pars.train.dfr[, mfull.attr3[-ii]], 
	   y = pars.train.dfr$M.full,
	   xtest = pars.select.dfr[, mfull.attr3[-ii]], 
	   ytest = pars.select.dfr[, "M.full"],
	   importance = FALSE, ntree = 1000)
	M.rf4var.matr[ii-2,] <- c(r.tmp$rsq[1000], r.tmp$test$rsq[1000], r.tmp$mse[1000], 
	   r.tmp$test$mse[1000])
   }
## Worst: sum.of.five.spwin. Omit, and repeat the loop.
mfull.attr4 <- c("sum.of.nine.spwin","highest3sum.spwin","n","vart")
M.rf5var.matr <- matrix(ncol = 4, nrow = length(mfull.attr4)-2)
for(ii in 3:(length(mfull.attr4)))
#for(ii in 1:2)
   {
	r.tmp <- randomForest(x = pars.train.dfr[, mfull.attr4[-ii]], 
	   y = pars.train.dfr$M.full,
	   xtest = pars.select.dfr[, mfull.attr4[-ii]], 
	   ytest = pars.select.dfr[, "M.full"],
	   importance = FALSE, ntree = 1000)
	M.rf5var.matr[ii-2,] <- c(r.tmp$rsq[1000], r.tmp$test$rsq[1000], r.tmp$mse[1000], 
	   r.tmp$test$mse[1000])
   }

## Based on this backward elimination, a four-covariate model seems to be a good 
## compromise:
mfull.attr.final <- c("sum.of.nine.spwin","highest3sum.spwin","n","vart")

mfull.rf.final <- randomForest(x = pars.train.dfr[, mfull.attr.final], y = pars.train.dfr$M.full,
   xtest = pars.select.dfr[, mfull.attr.final], 
   ytest = pars.select.dfr[, "M.full"],
   importance = TRUE, ntree = 1000, keep.forest = TRUE)
pred.tmp <- predict(mfull.rf.final, newdata = pars.test.dfr[, mfull.attr.final])

## CONCLUSION:
## The best seems to be the 4-attribute model mfull.rf3 with
## "sum.of.nine.spwin","n","highest3sum.spwin","vart"
## which sounds a bit stupid, but that's what I have.
## Save the model:
save(mfull.rf.final, mfull.attr.final, file = "results_baluev5/fadm.rfmodel.RObj")
rm(mfull.start.rf)

## I have the other predictions in the test set:
mfull.rf.final$test$predicted
mfull.rf.final$predicted

quartz()
plot(pars.train.dfr$M.full, mfull.rf.final$predicted, pch = 16, cex = 0.7, col = "black",
   ylim = range(c(pars.train.dfr$M.full, pars.test.dfr$M.full, pars.select.dfr$M.full, 
   mfull.rf.final$test$predicted, mfull.rf.final$predicted, pred.tmp), na.rm = TRUE))
points(pars.select.dfr$M.full, mfull.rf.final$test$predicted, pch = 16, cex = 0.7, col = "deepskyblue")
points(pars.test.dfr$M.full, pred.tmp, pch = 16, cex = 0.7, col = "orange")
## Some new "branch" seems to come up, but it's not that bad.

pars.train.dfr$M.full.rf <- mfull.rf.final$predicted
pars.select.dfr$M.full.rf <- mfull.rf.final$test$predicted
pars.test.dfr$M.full.rf <- pred.tmp


## --------------------------------------------------------------------------------------------------------
## Compute the FAPs with the estimated Ms
## --------------------------------------------------------------------------------------------------------

fadm.rf.fap.test <- matrix(ncol = ncol(maxima.test.dfr), nrow = nrow(maxima.test.dfr))
dimnames(fadm.rf.fap.test) <- list(NULL, colnames(maxima.test.dfr))
for(ii in 1:ncol(maxima.test.dfr))
   {
   	nn <- as.character(pars.test.dfr$name[ii])
   	nobs <- pars.test.dfr$n[pars.test.dfr$name == nn]
   	mest <- pars.test.dfr$M.full.rf[pars.test.dfr$name == nn]
   	fadm.rf.fap.test[,nn] <- sapply(maxima.test.dfr[,nn], schwczfap.fun, mest = mest, N = nobs)
   }

fadm.rf.fap.select <- matrix(ncol = ncol(maxima.select.dfr), nrow = nrow(maxima.select.dfr))
dimnames(fadm.rf.fap.select) <- list(NULL, colnames(maxima.select.dfr))
for(ii in 1:ncol(maxima.select.dfr))
   {
   	nn <- as.character(pars.select.dfr$name[ii])
   	nobs <- pars.select.dfr$n[pars.select.dfr$name == nn]
   	mest <- pars.select.dfr$M.full.rf[pars.select.dfr$name == nn]
   	fadm.rf.fap.select[,nn] <- sapply(maxima.select.dfr[,nn], schwczfap.fun, mest = mest, N = nobs)
   }

fadm.rf.fap.train <- matrix(ncol = ncol(maxima.train.dfr), nrow = nrow(maxima.train.dfr))
dimnames(fadm.rf.fap.train) <- list(NULL, colnames(maxima.train.dfr))
for(ii in 1:ncol(maxima.train.dfr))
   {
   	nn <- as.character(pars.train.dfr$name[ii])
   	nobs <- pars.train.dfr$n[pars.train.dfr$name == nn]
   	mest <- pars.train.dfr$M.full.rf[pars.train.dfr$name == nn]
	fadm.rf.fap.train[,nn] <- sapply(maxima.train.dfr[,nn], schwczfap.fun, mest = mest, N = nobs)
   }

range(apply(fadm.rf.fap.test, 2, function(vec) sum(vec < 0.05) / length(vec)), na.rm = TRUE)
range(apply(fadm.rf.fap.select, 2, function(vec) sum(vec < 0.05) / length(vec)), na.rm = TRUE)
range(apply(fadm.rf.fap.train, 2, function(vec) sum(vec < 0.05) / length(vec)), na.rm = TRUE)
range(apply(fadm.rf.fap.test, 2, function(vec) sum(vec < 0.01) / length(vec)), na.rm = TRUE)
range(apply(fadm.rf.fap.select, 2, function(vec) sum(vec < 0.01) / length(vec)), na.rm = TRUE)
range(apply(fadm.rf.fap.train, 2, function(vec) sum(vec < 0.01) / length(vec)), na.rm = TRUE)


## --------------------------------------------------------------------------------------------------------
## Check the number of significant values in the test set
## --------------------------------------------------------------------------------------------------------

pars.train.dfr$fr.rffadm.sign05 <- NA
pars.train.dfr$fr.rffadm.sign01 <- NA
pars.test.dfr$fr.rffadm.sign05 <- NA
pars.test.dfr$fr.rffadm.sign01 <- NA
pars.select.dfr$fr.rffadm.sign05 <- NA
pars.select.dfr$fr.rffadm.sign01 <- NA

for(ii in 1:nrow(pars.train.dfr))
   {
   	pars.train.dfr$fr.rffadm.sign05[ii] <-
   	   sum(fadm.rf.fap.train[, as.character(pars.train.dfr$name[ii])] < 0.05) / 1500
   	pars.train.dfr$fr.rffadm.sign01[ii] <-
   	   sum(fadm.rf.fap.train[, as.character(pars.train.dfr$name[ii])] < 0.01) / 1500
   }
for(ii in 1:nrow(pars.test.dfr))
   {
   	pars.test.dfr$fr.rffadm.sign05[ii] <-
   	   sum(fadm.rf.fap.test[, as.character(pars.test.dfr$name[ii])] < 0.05) / 1500
   	pars.test.dfr$fr.rffadm.sign01[ii] <-
   	   sum(fadm.rf.fap.test[, as.character(pars.test.dfr$name[ii])] < 0.01) / 1500
   }
for(ii in 1:nrow(pars.select.dfr))
   {
   	pars.select.dfr$fr.rffadm.sign05[ii] <-
   	   sum(fadm.rf.fap.select[, as.character(pars.select.dfr$name[ii])] < 0.05) / 1500
   	pars.select.dfr$fr.rffadm.sign01[ii] <-
   	   sum(fadm.rf.fap.select[, as.character(pars.select.dfr$name[ii])] < 0.01) / 1500
   }


save(fadm.rf.fap.test, fadm.rf.fap.train, fadm.rf.fap.select, file = "results_baluev5/fadm.rf.fap.all.RObj")
rm(fadm.rf.fap.test, fadm.rf.fap.train, fadm.rf.fap.select)


quartz(height = 7, width = 10)
par(mfrow = c(1,1), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.6,0))

col.x <- "beta"
plot(abs(pars.test.dfr[pars.test.dfr$region == "randompos", col.x]), 
   pars.test.dfr$fr.rffadm.sign05[pars.test.dfr$region == "randompos"], 
   type = "p", pch = 16, cex = 0.7, ylim = c(0,0.15), 
   xlim = range(abs(pars.test.dfr[,col.x])), main = "RF F^M FAP",
   ylab = "Fraction above threshold", xlab = col.x, col = "grey25")
points(abs(pars.test.dfr[pars.test.dfr$region == "randompos", col.x]),
   pars.test.dfr$fr.rffadm.sign01[pars.test.dfr$region == "randompos"], 
   type = "p", pch = 16, col = "grey", cex = 0.7)

points(abs(pars.test.dfr[pars.test.dfr$region == "ecl", col.x]), 
   pars.test.dfr$fr.rffadm.sign05[pars.test.dfr$region == "ecl"], 
    pch = 16, cex = 0.7, col = "violetred")
points(abs(pars.test.dfr[pars.test.dfr$region == "ecl", col.x]), 
   pars.test.dfr$fr.rffadm.sign01[pars.test.dfr$region == "ecl"], 
    pch = 16, cex = 0.7, col = "pink")

points(abs(pars.test.dfr[pars.test.dfr$region == "lmc", col.x]), 
   pars.test.dfr$fr.rffadm.sign05[pars.test.dfr$region == "lmc"], 
    pch = 16, cex = 0.7, col = "orange")
points(abs(pars.test.dfr[pars.test.dfr$region == "lmc", col.x]), 
   pars.test.dfr$fr.rffadm.sign01[pars.test.dfr$region == "lmc"], 
    pch = 16, cex = 0.7, col = "gold")

points(abs(pars.test.dfr[pars.test.dfr$region == "ln.ha", col.x]), 
   pars.test.dfr$fr.rffadm.sign05[pars.test.dfr$region == "ln.ha"], 
    pch = 16, cex = 0.7, col = "deepskyblue")
points(abs(pars.test.dfr[pars.test.dfr$region == "ln.ha", col.x]), 
   pars.test.dfr$fr.rffadm.sign01[pars.test.dfr$region == "ln.ha"], 
    pch = 16, cex = 0.7, col = "skyblue")

abline(h = 0.05, lwd = 2)
abline(h = 0.01, lty = 2, lwd = 2)
## Wow, including all those vart and spectral window variables removed much of the 
## beta-dependent biases! Not all, but a lot....

rm(mfull.rf.final)


## --------------------------------------------------------------------------------------------------------
## Rf for the GEV parameters
## --------------------------------------------------------------------------------------------------------

pars.train.dfr$log.n <- log(pars.train.dfr$n)
pars.train.dfr$log.xi.full <- log(-pars.train.dfr$xi.full)
pars.train.dfr$log.sig.full <- log(pars.train.dfr$sig.full)
pars.train.dfr$log.mu.full <- log(pars.train.dfr$mu.full)
pars.test.dfr$log.n <- log(pars.test.dfr$n)
pars.test.dfr$log.xi.full <- log(-pars.test.dfr$xi.full)
pars.test.dfr$log.sig.full <- log(pars.test.dfr$sig.full)
pars.test.dfr$log.mu.full <- log(pars.test.dfr$mu.full)
pars.select.dfr$log.n <- log(pars.select.dfr$n)
pars.select.dfr$log.xi.full <- log(-pars.select.dfr$xi.full)
pars.select.dfr$log.sig.full <- log(pars.select.dfr$sig.full)
pars.select.dfr$log.mu.full <- log(pars.select.dfr$mu.full)
logstart.attr <- colnames(pars.train.dfr)[c(8:58,74)]

## Initial model with all the variables

xi.start.rf <- randomForest(x = pars.train.dfr[, logstart.attr], y = pars.train.dfr$log.xi.full,
   xtest = pars.select.dfr[, logstart.attr], 
   ytest = pars.select.dfr[, "log.xi.full"],
   importance = TRUE, ntree = 1000, keep.forest = TRUE)

## Look at the results: the mean squared error and pseudo-R^2 on the training set
xi.start.rf$mse[1000]   ## 0.0009540109
xi.start.rf$rsq[1000]   ## 0.9928298

## Look at the results: the mean squared error and pseudo-R^2 on the test set
xi.start.rf$test$mse[1000]    ## 0.00835445 
xi.start.rf$test$rsq[1000]    ## ~ 0.9428968
## Deterioration.

## Attribute importance
quartz(height = 7.5, width = 7.5)
par(mfcol = c(1,1), mar = c(2.5,10,2,1), mgp = c(1.5,0.6,0))
o.tmp <- order(xi.start.rf$importance[,1])
barplot(xi.start.rf$importance[,1][o.tmp], horiz = TRUE, las = 1,
   names.arg = names(xi.start.rf$importance[,1][o.tmp]), cex.axis = 0.5)
## the best are log.n (well, yes), at60, at48, at20, sum.spwin, and at8.

sig.start.rf <- randomForest(x = pars.train.dfr[, logstart.attr], y = pars.train.dfr$log.sig.full,
   xtest = pars.select.dfr[, logstart.attr], 
   ytest = pars.select.dfr[, "log.sig.full"],
   importance = TRUE, ntree = 1000, keep.forest = TRUE)

## Look at the results: the mean squared error and pseudo-R^2 on the training set
sig.start.rf$mse[1000]   ## 0.0007718979
sig.start.rf$rsq[1000]   ## 0.9896197

## Look at the results: the mean squared error and pseudo-R^2 on the test set
sig.start.rf$test$mse[1000]    ## 0.003971609 at the end
sig.start.rf$test$rsq[1000]    ## 0.9517438
## Deterioration.

## Attribute importance
quartz(height = 7.5, width = 7.5)
par(mfcol = c(1,1), mar = c(2.5,10,2,1), mgp = c(1.5,0.6,0))
o.tmp <- order(sig.start.rf$importance[,1])
barplot(sig.start.rf$importance[,1][o.tmp], horiz = TRUE, las = 1,
   names.arg = names(sig.start.rf$importance[,1][o.tmp]), cex.axis = 0.5)
## The same are important for scale as for shape.

## I will compute mu from xi and sigma, not fit it. 

## Let's do some systematic search of reasonable variable.

## xi:

xi.rf2var.matr <- matrix(ncol = 4, nrow = length(logstart.attr)-1)
for(ii in 1:(length(logstart.attr)-1))
#for(ii in 1:2)
   {
	r.tmp <- randomForest(x = pars.train.dfr[, c("log.n", logstart.attr[ii])], 
	   y = pars.train.dfr$log.xi.full,
	   xtest = pars.select.dfr[, c("log.n", logstart.attr[ii])], 
	   ytest = pars.select.dfr[, "log.xi.full"],
	   importance = FALSE, ntree = 1000)
	xi.rf2var.matr[ii,] <- c(r.tmp$rsq[1000], r.tmp$test$rsq[1000], r.tmp$mse[1000], 
	   r.tmp$test$mse[1000])
   }
sort(xi.rf2var.matr[,2], decreasing = TRUE)   ## 6th and 13th: > 0.992
sort(xi.rf2var.matr[,4])   ## 13th < 0.0011
logstart.attr[c(order(xi.rf2var.matr[,2], decreasing = TRUE), length(logstart.attr))]
logstart.attr[c(order(xi.rf2var.matr[,4]), length(logstart.attr))]
## Winners: at52, at24. sum.spwin is somewhere lower (in the first half, still.)

sig.rf2var.matr <- matrix(ncol = 4, nrow = length(logstart.attr)-1)
for(ii in 1:(length(logstart.attr)-1))
#for(ii in 1:2)
   {
	r.tmp <- randomForest(x = pars.train.dfr[, c("log.n", logstart.attr[ii])], 
	   y = pars.train.dfr$log.sig.full,
	   xtest = pars.select.dfr[, c("log.n", logstart.attr[ii])], 
	   ytest = pars.select.dfr[, "log.sig.full"],
	   importance = FALSE, ntree = 1000)
	sig.rf2var.matr[ii,] <- c(r.tmp$rsq[1000], r.tmp$test$rsq[1000], r.tmp$mse[1000], 
	   r.tmp$test$mse[1000])
   }
## test$mse < 0.0009: sum.spwin, at52
## test$rsq > 0.99: sum.spwin, at52
sort(sig.rf2var.matr[,2], decreasing = TRUE)   ## 36th and 13th: > 0.99
sort(sig.rf2var.matr[,4])   ## 13th < 0.0008, 36th < 0.00083
logstart.attr[c(order(sig.rf2var.matr[,2], decreasing = TRUE), length(logstart.attr))]
logstart.attr[c(order(sig.rf2var.matr[,4]), length(logstart.attr))]
## The winner is everywhere at52! I see no reason for this, but it looks like 

## Refit a model with only log.n:

xi.rf2 <- randomForest(x = data.frame(pars.train.dfr[, "log.n"]), y = pars.train.dfr$log.xi.full,
   xtest = data.frame(pars.select.dfr[, "log.n"]), 
   ytest = pars.select.dfr[, "log.xi.full"],
   mtry = 1, importance = TRUE, ntree = 1000)

## Look at the results: the mean squared error and pseudo-R^2 on the training set
xi.rf2$mse[1000]   ## 0.001089166 at the end
xi.rf2$rsq[1000]   ## 0.991814 

## Look at the results: the mean squared error and pseudo-R^2 on the test set
xi.rf2$test$mse[1000]    ## 0.0007993622
xi.rf2$test$rsq[1000]    ## 0.9945363
## Gosh, this is one of the best, by these measures, I could omit the spectral window dependence
## from the model.

sig.rf2 <- randomForest(x = data.frame(pars.train.dfr[, "log.n"]), y = pars.train.dfr$log.sig.full,
   xtest = data.frame(pars.select.dfr[, "log.n"]), 
   ytest = pars.select.dfr[, "log.sig.full"],
   mtry = 1, importance = TRUE, ntree = 1000)

## Look at the results: the mean squared error and pseudo-R^2 on the training set
sig.rf2$mse[1000]   ## 0.001260176
sig.rf2$rsq[1000]   ## 0.9830535

## Look at the results: the mean squared error and pseudo-R^2 on the test set
sig.rf2$test$mse[1000]    ## 0.0009335051
sig.rf2$test$rsq[1000]    ## 0.9886576
sig.rf2var.matr[36,]
## Ok, this is a very little worse without the spectral window dependence. 
## The effect is tiny.

xi.rf.final <- randomForest(x = pars.train.dfr[, c("log.n","at52")],
   y = pars.train.dfr$log.xi.full,
   xtest = pars.select.dfr[, c("log.n","at52")], 
   ytest = pars.select.dfr[, "log.xi.full"],
   importance = TRUE, ntree = 1000, keep.forest = TRUE)
sig.rf.final <- randomForest(x = pars.train.dfr[, c("log.n","at52")],
   y = pars.train.dfr$log.sig.full,
   xtest = pars.select.dfr[, c("log.n","at52")], 
   ytest = pars.select.dfr[, "log.sig.full"],
   importance = TRUE, ntree = 1000, keep.forest = TRUE)

a.tmp <- quantile(pars.train.dfr$at52, pr = c(0.1,0.3,0.5,0.7,0.9))
g.tmp <- gray((pars.train.dfr$at52 - min(pars.train.dfr$at52)) / 
   (max(pars.train.dfr$at52) - min(pars.train.dfr$at52)))
## Plot the fits
quartz(height = 4, width = 8)
par(mfcol = c(1,2), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.6,0))
plot(pars.train.dfr$log.n, pars.train.dfr$xi.full, pch = 16, col = g.tmp, cex = 0.6)
for(ii in 1:5)
   {
	grid.tmp <- data.frame(log.n = log(40:250), at52 = rep(a.tmp[ii], times = 211))
	xi.tmp <- - exp(predict(xi.rf.final, newdata = grid.tmp))
	lines(log(40:250), xi.tmp, col = colorvec[ii])   ## No way!!!!!
   }
plot(pars.train.dfr$log.n, pars.train.dfr$sig.full, pch = 16, col = g.tmp, cex = 0.6)
for(ii in 1:5)
   {
	grid.tmp <- data.frame(log.n = log(40:250), at52 = rep(a.tmp[ii], times = 211))
	sig.tmp <- exp(predict(sig.rf.final, newdata = grid.tmp))
	lines(log(40:250), sig.tmp, col = colorvec[ii])   ## No way!!!!!
   }
## There is a restricted space in which "log.n","sum.spwin" pairs can occur; outside it,
## the rf fits go very wrong, and follow arbitrary paths over empty space. A lot more 
## points, which is a bit further off from the "log.n","sum.spwin" pairs occurring in the
## training set, can be very badly modelled.

save(xi.rf.final, sig.rf.final, file = "results_baluev5/gev.rfmodels.RObj")
rm(xi.start.rf, xi.rf1, xi.rf2, sig.start.rf, sig.rf1, sig.rf2)

pars.train.dfr$xi.full.rf <- -exp(xi.rf.final$predicted)
pars.select.dfr$xi.full.rf <- -exp(xi.rf.final$test$predicted)
pars.test.dfr$xi.full.rf <- -exp(predict(xi.rf.final, newdata = pars.test.dfr[, c("log.n","at52")]))
pars.train.dfr$sig.full.rf <- exp(sig.rf.final$predicted)
pars.select.dfr$sig.full.rf <- exp(sig.rf.final$test$predicted)
pars.test.dfr$sig.full.rf <- exp(predict(sig.rf.final, newdata = pars.test.dfr[, c("log.n","at52")]))
pars.train.dfr$mu.full.rf <- 1+pars.train.dfr$sig.full.rf/pars.train.dfr$xi.full.rf
pars.test.dfr$mu.full.rf <- 1+pars.test.dfr$sig.full.rf/pars.test.dfr$xi.full.rf
pars.select.dfr$mu.full.rf <- 1+pars.select.dfr$sig.full.rf/pars.select.dfr$xi.full.rf

quartz(height = 6, width = 6)
par(mfrow = c(2,2), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.6,0))
plot(pars.train.dfr$xi.full, pars.train.dfr$xi.full.rf, pch = 16, cex = 0.7, col = "black",
   ylim = range(c(pars.train.dfr$xi.full, pars.test.dfr$xi.full, pars.select.dfr$xi.full, 
   pars.train.dfr$xi.full.rf, pars.test.dfr$xi.full.rf, pars.select.dfr$xi.full.rf), na.rm = TRUE),
   xlim = range(c(pars.train.dfr$xi.full, pars.test.dfr$xi.full, pars.select.dfr$xi.full, 
   pars.train.dfr$xi.full.rf, pars.test.dfr$xi.full.rf, pars.select.dfr$xi.full.rf), na.rm = TRUE))
points(pars.test.dfr$xi.full, pars.test.dfr$xi.full.rf, pch = 16, cex = 0.7, col = "deepskyblue")
points(pars.select.dfr$xi.full, pars.select.dfr$xi.full.rf, pch = 16, cex = 0.7, col = "orange")
plot(pars.train.dfr$sig.full, pars.train.dfr$sig.full.rf, pch = 16, cex = 0.7, col = "black",
   ylim = range(c(pars.train.dfr$sig.full, pars.test.dfr$sig.full, pars.select.dfr$sig.full, 
   pars.train.dfr$sig.full.rf, pars.test.dfr$sig.full.rf, pars.select.dfr$sig.full.rf), na.rm = TRUE),
   xlim = range(c(pars.train.dfr$sig.full, pars.test.dfr$sig.full, pars.select.dfr$sig.full, 
   pars.train.dfr$sig.full.rf, pars.test.dfr$sig.full.rf, pars.select.dfr$sig.full.rf), na.rm = TRUE))
points(pars.test.dfr$sig.full, pars.test.dfr$sig.full.rf, pch = 16, cex = 0.7, col = "deepskyblue")
points(pars.select.dfr$sig.full, pars.select.dfr$sig.full.rf, pch = 16, cex = 0.7, col = "orange")
plot(pars.train.dfr$mu.full, pars.train.dfr$mu.full.rf, pch = 16, cex = 0.7, col = "black",
   ylim = range(c(pars.train.dfr$mu.full, pars.test.dfr$mu.full, pars.select.dfr$mu.full, 
   pars.train.dfr$mu.full.rf, pars.test.dfr$mu.full.rf, pars.select.dfr$mu.full.rf), na.rm = TRUE),
   xlim = range(c(pars.train.dfr$mu.full, pars.test.dfr$mu.full, pars.select.dfr$mu.full, 
   pars.train.dfr$mu.full.rf, pars.test.dfr$mu.full.rf, pars.select.dfr$mu.full.rf), na.rm = TRUE))
points(pars.test.dfr$mu.full, pars.test.dfr$mu.full.rf, pch = 16, cex = 0.7, col = "deepskyblue")
points(pars.select.dfr$mu.full, pars.select.dfr$mu.full.rf, pch = 16, cex = 0.7, col = "orange")
## This does not look very good, there are some discrepancies at the extremities.


## --------------------------------------------------------------------------------------------------------
## Compute the individual FAPs at each test site (and possibly train site)
## --------------------------------------------------------------------------------------------------------

rf.gevfap.test.dfr <- data.frame(matrix(ncol = ncol(maxima.test.dfr),
    nrow = 1500))
colnames(rf.gevfap.test.dfr) <- colnames(maxima.test.dfr)

for(ii in 1:nrow(pars.test.dfr))
#  for(ii in 5)
   {
   	nn <- as.character(pars.test.dfr$name[ii])
   	shape <- pars.test.dfr$xi.full.rf[pars.test.dfr$name == nn]
   	scale <- pars.test.dfr$sig.full.rf[pars.test.dfr$name == nn]
   	loc <- pars.test.dfr$mu.full.rf[pars.test.dfr$name == nn]
    zass <- maxima.test.dfr[, nn]
    if(is.na(shape) | is.na(scale) | is.na(loc)) {
    	rf.gevfap.test.dfr[, nn] <- rep(NA, 1500)
       } else {
	    rf.gevfap.test.dfr[, nn] <- 1 - sapply(zass, evd::pgev,
	       shape = shape, loc = loc, scale = scale)
	   }
   }
as.numeric(apply(rf.gevfap.test.dfr, 2, function(vec) sum(is.na(vec))))

rf.gevfap.select.dfr <- data.frame(matrix(ncol = ncol(maxima.select.dfr),
    nrow = 1500))
colnames(rf.gevfap.select.dfr) <- colnames(maxima.select.dfr)

for(ii in 1:nrow(pars.select.dfr))
#  for(ii in 5)
   {
   	nn <- as.character(pars.select.dfr$name[ii])
   	shape <- pars.select.dfr$xi.full.rf[pars.select.dfr$name == nn]
   	scale <- pars.select.dfr$sig.full.rf[pars.select.dfr$name == nn]
   	loc <- pars.select.dfr$mu.full.rf[pars.select.dfr$name == nn]
    zass <- maxima.select.dfr[, nn]
    if(is.na(shape) | is.na(scale) | is.na(loc)) {
    	rf.gevfap.select.dfr[, nn] <- rep(NA, 1500)
       } else {
	    rf.gevfap.select.dfr[, nn] <- 1 - sapply(zass, evd::pgev,
	       shape = shape, loc = loc, scale = scale)
	   }
   }
as.numeric(apply(rf.gevfap.select.dfr, 2, function(vec) sum(is.na(vec))))

rf.gevfap.train.dfr <- data.frame(matrix(ncol = ncol(maxima.train.dfr),
    nrow = 1500))
colnames(rf.gevfap.train.dfr) <- colnames(maxima.train.dfr)

for(ii in 1:nrow(pars.train.dfr))
#  for(ii in 5)
   {
   	nn <- as.character(pars.train.dfr$name[ii])
   	shape <- pars.train.dfr$xi.full.rf[pars.train.dfr$name == nn]
   	scale <- pars.train.dfr$sig.full.rf[pars.train.dfr$name == nn]
   	loc <- pars.train.dfr$mu.full.rf[pars.train.dfr$name == nn]
    zass <- maxima.train.dfr[, nn]
    if(is.na(shape) | is.na(scale) | is.na(loc)) {
    	rf.gevfap.train.dfr[, nn] <- rep(NA, 1500)
       } else {
	    rf.gevfap.train.dfr[, nn] <- 1 - sapply(zass, evd::pgev,
	       shape = shape, loc = loc, scale = scale)
	   }
   }
save(rf.gevfap.test.dfr, file = "results_baluev5/rf.gevfap.test.dfr.RObj")
save(rf.gevfap.select.dfr, file = "results_baluev5/rf.gevfap.select.dfr.RObj")
save(rf.gevfap.train.dfr, file = "results_baluev5/rf.gevfap.train.dfr.RObj")


pars.train.dfr$fr.rfgev95 <- NA
pars.train.dfr$fr.rfgev99 <- NA
pars.test.dfr$fr.rfgev95 <- NA
pars.test.dfr$fr.rfgev99 <- NA
pars.select.dfr$fr.rfgev95 <- NA
pars.select.dfr$fr.rfgev99 <- NA

## Compute the success rates
for(ii in 1:nrow(pars.train.dfr))
   {
   	pars.train.dfr$fr.rfgev95[ii] <-
   	   sum(rf.gevfap.train.dfr[, as.character(pars.train.dfr$name[ii])] < 0.05) / 1500
   	pars.train.dfr$fr.rfgev99[ii] <-
   	   sum(rf.gevfap.train.dfr[, as.character(pars.train.dfr$name[ii])] < 0.01) / 1500
   }
for(ii in 1:nrow(pars.test.dfr))
   {
   	pars.test.dfr$fr.rfgev95[ii] <-
   	   sum(rf.gevfap.test.dfr[, as.character(pars.test.dfr$name[ii])] < 0.05) / 1500
   	pars.test.dfr$fr.rfgev99[ii] <-
   	   sum(rf.gevfap.test.dfr[, as.character(pars.test.dfr$name[ii])] < 0.01) / 1500
   }
for(ii in 1:nrow(pars.select.dfr))
   {
   	pars.select.dfr$fr.rfgev95[ii] <-
   	   sum(rf.gevfap.select.dfr[, as.character(pars.select.dfr$name[ii])] < 0.05) / 1500
   	pars.select.dfr$fr.rfgev99[ii] <-
   	   sum(rf.gevfap.select.dfr[, as.character(pars.select.dfr$name[ii])] < 0.01) / 1500
   }




quartz(height = 7, width = 10)
par(mfrow = c(1,1), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.6,0))

col.x <- "n"
plot(abs(pars.test.dfr[pars.test.dfr$region == "randompos", col.x]), 
   pars.test.dfr$fr.rfgev95[pars.test.dfr$region == "randompos"], 
   type = "p", pch = 16, cex = 0.7, ylim = c(0,0.15), 
   xlim = range(abs(pars.test.dfr[,col.x])), main = "Gev FAP",
   ylab = "Fraction above threshold", xlab = col.x, col = "grey25")
points(abs(pars.test.dfr[pars.test.dfr$region == "randompos", col.x]),
   pars.test.dfr$fr.rfgev99[pars.test.dfr$region == "randompos"], 
   type = "p", pch = 16, col = "grey", cex = 0.7)

points(abs(pars.test.dfr[pars.test.dfr$region == "ecl", col.x]), 
   pars.test.dfr$fr.rfgev95[pars.test.dfr$region == "ecl"], 
    pch = 16, cex = 0.7, col = "violetred")
points(abs(pars.test.dfr[pars.test.dfr$region == "ecl", col.x]), 
   pars.test.dfr$fr.rfgev99[pars.test.dfr$region == "ecl"], 
    pch = 16, cex = 0.7, col = "pink")

points(abs(pars.test.dfr[pars.test.dfr$region == "lmc", col.x]), 
   pars.test.dfr$fr.rfgev95[pars.test.dfr$region == "lmc"], 
    pch = 16, cex = 0.7, col = "orange")
points(abs(pars.test.dfr[pars.test.dfr$region == "lmc", col.x]), 
   pars.test.dfr$fr.rfgev99[pars.test.dfr$region == "lmc"], 
    pch = 16, cex = 0.7, col = "gold")

points(abs(pars.test.dfr[pars.test.dfr$region == "ln.ha", col.x]), 
   pars.test.dfr$fr.rfgev95[pars.test.dfr$region == "ln.ha"], 
    pch = 16, cex = 0.7, col = "deepskyblue")
points(abs(pars.test.dfr[pars.test.dfr$region == "ln.ha", col.x]), 
   pars.test.dfr$fr.rfgev99[pars.test.dfr$region == "ln.ha"], 
    pch = 16, cex = 0.7, col = "skyblue")

abline(h = 0.05, lwd = 2)
abline(h = 0.01, lty = 2, lwd = 2)
## This is very bad. No need to save it, this will not be the fit I use.

## Repeat the thin-plate spline for the gev, now comparing the performance with at12 to
## the best of the attributes, sum.spwin.

rm(rf.gevfap.test.dfr, rf.gevfap.train.dfr, rf.gevfap.select.dfr)

rm(xi.rf.final, sig.rf.final)


## --------------------------------------------------------------------------------------------------------
## For the quantiles
## --------------------------------------------------------------------------------------------------------

## Initial model with all the variables

q95.start.rf <- randomForest(x = pars.train.dfr[, start.attr], y = pars.train.dfr$q95,
   xtest = pars.select.dfr[, start.attr], ytest = pars.select.dfr[, "q95"],
   importance = TRUE, ntree = 1000, keep.forest = TRUE)

q99.start.rf <- randomForest(x = pars.train.dfr[, start.attr], y = pars.train.dfr$q99,
   xtest = pars.select.dfr[, start.attr], ytest = pars.select.dfr[, "q99"],
   importance = TRUE, ntree = 1000, keep.forest = TRUE)

## Look at the results: the mean squared error and pseudo-R^2 on the training set
q95.start.rf$mse[1000]   ## 2.696249e-05
q95.start.rf$rsq[1000]   ## 0.9961472

## Look at the results: the mean squared error and pseudo-R^2 on the test set
q95.start.rf$test$mse[1000]    ##  0.0006496201
q95.start.rf$test$rsq[1000]    ## 0.9144163
## Like for GEV, deterioration: overfitting.

## Look at the results: the mean squared error and pseudo-R^2 on the training set
q99.start.rf$mse[1000]   ## 4.58856e-05
q99.start.rf$rsq[1000]   ## 0.9942436

## Look at the results: the mean squared error and pseudo-R^2 on the test set
q99.start.rf$test$mse[1000]    ## 0.0007156179
q99.start.rf$test$rsq[1000]    ## 0.9175609
## Like for GEV, deterioration: overfitting.

## Attribute importance

quartz(height = 7.5, width = 12.5)
par(mfcol = c(1,2), mar = c(2.5,10,2,1), mgp = c(1.5,0.6,0))
o.tmp <- order(q95.start.rf$importance[,1])
barplot(q95.start.rf$importance[,1][o.tmp], horiz = TRUE, las = 1,
   names.arg = names(q95.start.rf$importance[,1][o.tmp]), cex.axis = 0.5)
o.tmp <- order(q99.start.rf$importance[,1])
barplot(q99.start.rf$importance[,1][o.tmp], horiz = TRUE, las = 1,
   names.arg = names(q99.start.rf$importance[,1][o.tmp]), cex.axis = 0.5)
## the six best are the same as with the gev. No wonder, to my knowledge.

## Systematic forward search:

q95.rf2var.matr <- matrix(ncol = 4, nrow = length(start.attr)-1)
for(ii in 2:(length(start.attr)))
#for(ii in 1:2)
   {
	r.tmp <- randomForest(x = pars.train.dfr[, c("n", start.attr[ii])], 
	   y = pars.train.dfr$q95,
	   xtest = pars.select.dfr[, c("n", start.attr[ii])], 
	   ytest = pars.select.dfr[, "q95"],
	   importance = FALSE, ntree = 1000)
	q95.rf2var.matr[ii-1,] <- c(r.tmp$rsq[1000], r.tmp$test$rsq[1000], r.tmp$mse[1000], 
	   r.tmp$test$mse[1000])
   }
## test$rsq > 0.9955: at48, at60, at20
## test$mse < 3e-5: at48, at60, at16, at36, at20
sort(q95.rf2var.matr[,2], decreasing = TRUE)   ## at24 and at52: > 0.994
sort(q95.rf2var.matr[,4])   ## at24 and at52: < 4.5e-5
start.attr[c(order(q95.rf2var.matr[,2], decreasing = TRUE)+1, 1)]
start.attr[c(order(q95.rf2var.matr[,4])+1, 1)]
## Best two: at24, at52

q99.rf2var.matr <- matrix(ncol = 4, nrow = length(start.attr)-1)
for(ii in 2:(length(start.attr)))
#for(ii in 1:2)
   {
	r.tmp <- randomForest(x = pars.train.dfr[, c("n", start.attr[ii])], 
	   y = pars.train.dfr$q99,
	   xtest = pars.select.dfr[, c("n", start.attr[ii])], 
	   ytest = pars.select.dfr[, "q99"],
	   importance = FALSE, ntree = 1000)
	q99.rf2var.matr[ii-1,] <- c(r.tmp$rsq[1000], r.tmp$test$rsq[1000], r.tmp$mse[1000], 
	   r.tmp$test$mse[1000])
   }
sort(q99.rf2var.matr[,2], decreasing = TRUE)   ## at24 and at52: > 0.994
sort(q99.rf2var.matr[,4])   ## at24 and at52: < 4.5e-5
start.attr[c(order(q99.rf2var.matr[,2], decreasing = TRUE)+1, 1)]
start.attr[c(order(q99.rf2var.matr[,4])+1, 1)]
## Best two: at24, at52. Whoaa.

quan.attr.final <- c("n","at24")

q95.rf.final <- randomForest(x = pars.train.dfr[, quan.attr.final], y = pars.train.dfr$q95,
   xtest = pars.select.dfr[, quan.attr.final], ytest = pars.select.dfr[, "q95"],
   importance = TRUE, ntree = 5000, keep.forest = TRUE)

q99.rf.final <- randomForest(x = pars.train.dfr[, quan.attr.final], y = pars.train.dfr$q99,
   xtest = pars.select.dfr[, quan.attr.final], ytest = pars.select.dfr[, "q99"],
   importance = TRUE, ntree = 5000, keep.forest = TRUE)

pars.train.dfr$q95.rf <- q95.rf.final$predicted
pars.train.dfr$q99.rf <- q99.rf.final$predicted
pars.select.dfr$q95.rf <- q95.rf.final$test$predicted
pars.select.dfr$q99.rf <- q99.rf.final$test$predicted
pars.test.dfr$q95.rf <- predict(q95.rf.final, newdata = pars.test.dfr[,c("n","at24")])
pars.test.dfr$q99.rf <- predict(q99.rf.final, newdata = pars.test.dfr[,c("n","at24")])

save(q95.rf.final,q99.rf.final, file = "results_baluev5/rfmodel.quan.both.RObj")
rm(q95.start.rf,q99.start.rf)


## --------------------------------------------------------------------------------------------------------
## Compute the false alarm fractions at each test site (and possibly train site)
## --------------------------------------------------------------------------------------------------------

## Compute the success rates

for(ii in 1:nrow(pars.train.dfr))
   {
   	pars.train.dfr$fr.above.thr95rf[ii] <-
   	   sum(maxima.train.dfr[, as.character(pars.train.dfr$name[ii])] >
   	   pars.train.dfr$q95.rf[ii]) / 1500
   	pars.train.dfr$fr.above.thr99rf[ii] <-
   	   sum(maxima.train.dfr[, as.character(pars.train.dfr$name[ii])] >
   	   pars.train.dfr$q99.rf[ii]) / 1500
   }
for(ii in 1:nrow(pars.test.dfr))
   {
   	pars.test.dfr$fr.above.thr95rf[ii] <-
   	   sum(maxima.test.dfr[, as.character(pars.test.dfr$name[ii])] >
   	   pars.test.dfr$q95.rf[ii]) / 1500
   	pars.test.dfr$fr.above.thr99rf[ii] <-
   	   sum(maxima.test.dfr[, as.character(pars.test.dfr$name[ii])] >
   	   pars.test.dfr$q99.rf[ii]) / 1500
   }
for(ii in 1:nrow(pars.select.dfr))
   {
   	pars.select.dfr$fr.above.thr95rf[ii] <-
   	   sum(maxima.select.dfr[, as.character(pars.select.dfr$name[ii])] >
   	   pars.select.dfr$q95.rf[ii]) / 1500
   	pars.select.dfr$fr.above.thr99rf[ii] <-
   	   sum(maxima.select.dfr[, as.character(pars.select.dfr$name[ii])] >
   	   pars.select.dfr$q99.rf[ii]) / 1500
   }


quartz(height = 7, width = 10)
par(mfrow = c(1,1), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.6,0))

col.x <- "n"
plot(abs(pars.test.dfr[pars.test.dfr$region == "randompos", col.x]), 
   pars.test.dfr$fr.above.thr95rf[pars.test.dfr$region == "randompos"], 
   type = "p", pch = 16, cex = 0.7, ylim = c(0,0.3), 
   xlim = range(abs(pars.test.dfr[,col.x])), main = "Quantile method",
   ylab = "Fraction above threshold", xlab = col.x, col = "grey25")
points(abs(pars.test.dfr[pars.test.dfr$region == "randompos", col.x]),
   pars.test.dfr$fr.above.thr99rf[pars.test.dfr$region == "randompos"], 
   type = "p", pch = 16, col = "grey", cex = 0.7)

points(abs(pars.test.dfr[pars.test.dfr$region == "ecl", col.x]), 
   pars.test.dfr$fr.above.thr95rf[pars.test.dfr$region == "ecl"], 
    pch = 16, cex = 0.7, col = "violetred")
points(abs(pars.test.dfr[pars.test.dfr$region == "ecl", col.x]), 
   pars.test.dfr$fr.above.thr99rf[pars.test.dfr$region == "ecl"], 
    pch = 16, cex = 0.7, col = "pink")

points(abs(pars.test.dfr[pars.test.dfr$region == "lmc", col.x]), 
   pars.test.dfr$fr.above.thr95rf[pars.test.dfr$region == "lmc"], 
    pch = 16, cex = 0.7, col = "orange")
points(abs(pars.test.dfr[pars.test.dfr$region == "lmc", col.x]), 
   pars.test.dfr$fr.above.thr99rf[pars.test.dfr$region == "lmc"], 
    pch = 16, cex = 0.7, col = "gold")

points(abs(pars.test.dfr[pars.test.dfr$region == "ln.ha", col.x]), 
   pars.test.dfr$fr.above.thr95rf[pars.test.dfr$region == "ln.ha"], 
    pch = 16, cex = 0.7, col = "deepskyblue")
points(abs(pars.test.dfr[pars.test.dfr$region == "ln.ha", col.x]), 
   pars.test.dfr$fr.above.thr99rf[pars.test.dfr$region == "ln.ha"], 
    pch = 16, cex = 0.7, col = "skyblue")

abline(h = 0.05, lwd = 2)
abline(h = 0.01, lty = 2, lwd = 2)

## This is just as bad as the rf fit for GEV. Not necessary to push it further.

rm(a.tmp, col.x, grid.tmp,ii,loc,nn,o.tmp,q95.rf.final,q99.rf.final,r.tmp, rf.gevfap.select.dfr, rf.gevfap.test.dfr,rf.gevfap.train.dfr,scale,shape,sig.rf.final,sig.tmp, xi.rf.final, xi.tmp, zass)
save.image(WORKSPACE.NAME)


## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------
## Models: Smooth spline 
## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------

## --------------------------------------------------------------------------------------------------------
## Smooth spline for the quantiles
## --------------------------------------------------------------------------------------------------------

q95.smoo <- smooth.spline(x = log(pars.train.dfr$n), y = log(pars.train.dfr$q95),
   spar = 0.8)

quartz(height = 7, width = 9)
par(mfcol = c(1,1), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.6,0))
plot(log(pars.train.dfr$n), log(pars.train.dfr$q95), pch = 16, cex = 0.6, xlab = "N", ylab = "0.95-quantile",
   main = "0.95-quantile")
points(log(pars.test.dfr$n), log(pars.test.dfr$q95), pch = 1, cex = 0.6, col = "deepskyblue")
points(log(pars.select.dfr$n), log(pars.select.dfr$q95), pch = 1, cex = 0.6, col = "orangered")
lines(q95.smoo, col = "blueviolet", lwd = 2.5)

q99.smoo <- smooth.spline(x = log(pars.train.dfr$n), y = log(pars.train.dfr$q99),
   spar = 0.8)

quartz(height = 7, width = 9)
par(mfcol = c(1,1), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.6,0))
plot(log(pars.train.dfr$n), log(pars.train.dfr$q99), pch = 16, cex = 0.6, xlab = "N", ylab = "0.99-quantile",
   main = "0.99-quantile")
points(log(pars.test.dfr$n), log(pars.test.dfr$q99), pch = 1, cex = 0.6, col = "deepskyblue")
points(log(pars.select.dfr$n), log(pars.select.dfr$q99), pch = 1, cex = 0.6, col = "orangered")
lines(q99.smoo, col = "blueviolet", lwd = 2.5)
## They look truly nice.
save(q95.smoo, q99.smoo, file = "results_baluev5/smoomodel.quan.both.RObj")


## --------------------------------------------------------------------------------------------------------
## Predict smooth spline 
## --------------------------------------------------------------------------------------------------------

pars.select.dfr$q95.smoo <- exp(predict(q95.smoo, x = log(pars.select.dfr$n))$y)
pars.select.dfr$q99.smoo <- exp(predict(q99.smoo, x = log(pars.select.dfr$n))$y)
quartz()
plot(pars.select.dfr$q95, pars.select.dfr$q95.smoo); abline(c(0,1), col = "orangered", lwd = 2)

pars.test.dfr$q95.smoo <- exp(predict(q95.smoo, x = log(pars.test.dfr$n))$y)
pars.test.dfr$q99.smoo <- exp(predict(q99.smoo, x = log(pars.test.dfr$n))$y)
quartz()
plot(pars.test.dfr$q95, pars.test.dfr$q95.smoo); abline(c(0,1), col = "orangered", lwd = 2)

pars.train.dfr$q95.smoo <- exp(predict(q95.smoo, x = log(pars.train.dfr$n))$y)
pars.train.dfr$q99.smoo <- exp(predict(q99.smoo, x = log(pars.train.dfr$n))$y)
quartz()
plot(pars.train.dfr$q95, pars.train.dfr$q95.smoo); abline(c(0,1), col = "orangered", lwd = 2)

## Looks good. 


## --------------------------------------------------------------------------------------------------------
## Compute the individual FAPs at each test site (and possibly train site)
## --------------------------------------------------------------------------------------------------------


pars.train.dfr$fr.above.thr95smoo <- NA
pars.train.dfr$fr.above.thr99smoo <- NA
pars.test.dfr$fr.above.thr95smoo <- NA
pars.test.dfr$fr.above.thr99smoo <- NA
pars.select.dfr$fr.above.thr95smoo <- NA
pars.select.dfr$fr.above.thr99smoo <- NA

## Compute these success rates only on the first 1000 simulations (since F^M and GEV will
## need the other 500 for training)
for(ii in 1:nrow(pars.train.dfr))
   {
   	pars.train.dfr$fr.above.thr95smoo[ii] <-
   	   sum(maxima.train.dfr[, as.character(pars.train.dfr$name[ii])] >
   	   pars.train.dfr$q95.smoo[ii]) / 1500
   	pars.train.dfr$fr.above.thr99smoo[ii] <-
   	   sum(maxima.train.dfr[, as.character(pars.train.dfr$name[ii])] >
   	   pars.train.dfr$q99.smoo[ii]) / 1500
   }
for(ii in 1:nrow(pars.test.dfr))
   {
   	pars.test.dfr$fr.above.thr95smoo[ii] <-
   	   sum(maxima.test.dfr[, as.character(pars.test.dfr$name[ii])] >
   	   pars.test.dfr$q95.smoo[ii]) / 1500
   	pars.test.dfr$fr.above.thr99smoo[ii] <-
   	   sum(maxima.test.dfr[, as.character(pars.test.dfr$name[ii])] >
   	   pars.test.dfr$q99.smoo[ii]) / 1500
   }
for(ii in 1:nrow(pars.select.dfr))
   {
   	pars.select.dfr$fr.above.thr95smoo[ii] <-
   	   sum(maxima.select.dfr[, as.character(pars.select.dfr$name[ii])] >
   	   pars.select.dfr$q95.smoo[ii]) / 1500
   	pars.select.dfr$fr.above.thr99smoo[ii] <-
   	   sum(maxima.select.dfr[, as.character(pars.select.dfr$name[ii])] >
   	   pars.select.dfr$q99.smoo[ii]) / 1500
   }

sqrt(mean((pars.test.dfr$q95 - pars.test.dfr$q95.smoo)^2))   ## 0.003328131
sqrt(mean((pars.test.dfr$q99 - pars.test.dfr$q99.smoo)^2))   ## 0.004740872

## Show pattern of % of maxima above the fitted quantile
## Workspaces: /Users/mariasuveges/Documents/SpectralAnalysis/FAP/feasibilityForGaia/OtherFAPs/quantileFitting1.R
## 

quartz(height = 7, width = 10)
par(mfrow = c(1,1), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.6,0))

col.x <- "beta"
plot(abs(pars.test.dfr[pars.test.dfr$region == "randompos", col.x]), 
   pars.test.dfr$fr.above.thr95smoo[pars.test.dfr$region == "randompos"], 
   type = "p", pch = 16, cex = 0.7, ylim = c(0,0.3), 
   xlim = range(abs(pars.test.dfr[,col.x])), main = "Quantile method",
   ylab = "Fraction above threshold", xlab = col.x, col = "grey25")
points(abs(pars.test.dfr[pars.test.dfr$region == "randompos", col.x]),
   pars.test.dfr$fr.above.thr99smoo[pars.test.dfr$region == "randompos"], 
   type = "p", pch = 16, col = "grey", cex = 0.7)

points(abs(pars.test.dfr[pars.test.dfr$region == "ecl", col.x]), 
   pars.test.dfr$fr.above.thr95smoo[pars.test.dfr$region == "ecl"], 
    pch = 16, cex = 0.7, col = "violetred")
points(abs(pars.test.dfr[pars.test.dfr$region == "ecl", col.x]), 
   pars.test.dfr$fr.above.thr99smoo[pars.test.dfr$region == "ecl"], 
    pch = 16, cex = 0.7, col = "pink")

points(abs(pars.test.dfr[pars.test.dfr$region == "lmc", col.x]), 
   pars.test.dfr$fr.above.thr95smoo[pars.test.dfr$region == "lmc"], 
    pch = 16, cex = 0.7, col = "orange")
points(abs(pars.test.dfr[pars.test.dfr$region == "lmc", col.x]), 
   pars.test.dfr$fr.above.thr99smoo[pars.test.dfr$region == "lmc"], 
    pch = 16, cex = 0.7, col = "gold")

points(abs(pars.test.dfr[pars.test.dfr$region == "ln.ha", col.x]), 
   pars.test.dfr$fr.above.thr95smoo[pars.test.dfr$region == "ln.ha"], 
    pch = 16, cex = 0.7, col = "deepskyblue")
points(abs(pars.test.dfr[pars.test.dfr$region == "ln.ha", col.x]), 
   pars.test.dfr$fr.above.thr99smoo[pars.test.dfr$region == "ln.ha"], 
    pch = 16, cex = 0.7, col = "skyblue")

abline(h = 0.05, lwd = 2)
abline(h = 0.01, lty = 2, lwd = 2)
## Ok, so indeed, there is the decrease of the performance at the high-aliasing places.
## Certainly needs a variable capturing the high aliases, and so the 2D spline fits.  

rm(q95.smoo, q99.smoo)
save.image(WORKSPACE.NAME)


## --------------------------------------------------------------------------------------------------------
## Smoothing for the GEV parameters
## --------------------------------------------------------------------------------------------------------

## --------------------------------------------------------------------------------------------------------
## Smooth spline for the two parameters
## --------------------------------------------------------------------------------------------------------

xi.smoo <- smooth.spline(x = log(pars.train.dfr$n), y = log(-pars.train.dfr$xi.full),
   spar = 0.8)
#xi.lm <- lm(log(-pars.train.dfr$xi.full) ~ log(pars.train.dfr$n))

quartz(height = 7, width = 9)
par(mfcol = c(1,1), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.6,0))
plot(log(pars.train.dfr$n), log(-pars.train.dfr$xi.full), pch = 16, cex = 0.6, xlab = "N", ylab = "Shape",
   main = "Shape")
points(log(pars.test.dfr$n), log(-pars.test.dfr$xi.full), pch = 1, cex = 0.6, col = "deepskyblue")
points(log(pars.select.dfr$n), log(-pars.select.dfr$xi.full), pch = 1, cex = 0.6, col = "orangered")
lines(xi.smoo, col = "blueviolet", lwd = 2.5)
#lines(log(dfr1$n), xi.lm$fitted, col = "seagreen", lwd = 2.5)
## Ok, this has not become linear.
#rm(xi.lm) 

sig.smoo <- smooth.spline(x = log(pars.train.dfr$n), y = log(pars.train.dfr$sig.full),
   spar = 0.8)

quartz(height = 7, width = 9)
par(mfcol = c(1,1), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.6,0))
plot(log(pars.train.dfr$n), log(pars.train.dfr$sig.full), pch = 16, cex = 0.6, xlab = "N", ylab = "Scale",
   main = "Scale")
points(log(pars.test.dfr$n), log(pars.test.dfr$sig.full), pch = 1, cex = 0.6, col = "deepskyblue")
points(log(pars.select.dfr$n), log(pars.select.dfr$sig.full), pch = 1, cex = 0.6, col = "orangered")
lines(sig.smoo, col = "blueviolet", lwd = 2.5)

save(xi.smoo, sig.smoo, file = "results_baluev5/smoomodel.gev.both.RObj")


## --------------------------------------------------------------------------------------------------------
## Predict smooth spline 
## --------------------------------------------------------------------------------------------------------

pars.test.dfr$xi.smoo <- -exp(predict(xi.smoo, x = log(pars.test.dfr$n))$y)
pars.test.dfr$sig.smoo <- exp(predict(sig.smoo, x = log(pars.test.dfr$n))$y)
pars.test.dfr$mu.smoo <- 1 + pars.test.dfr$sig.smoo / pars.test.dfr$xi.smoo
quartz()
plot(pars.test.dfr$xi.full, pars.test.dfr$xi.smoo) 
abline(c(0,1), col = "orangered", lwd = 2)

pars.select.dfr$xi.smoo <- -exp(predict(xi.smoo, x = log(pars.select.dfr$n))$y)
pars.select.dfr$sig.smoo <- exp(predict(sig.smoo, x = log(pars.select.dfr$n))$y)
pars.select.dfr$mu.smoo <- 1 + pars.select.dfr$sig.smoo / pars.select.dfr$xi.smoo
quartz()
plot(pars.select.dfr$xi.full, pars.select.dfr$xi.smoo) 
abline(c(0,1), col = "orangered", lwd = 2)

pars.train.dfr$xi.smoo <- -exp(predict(xi.smoo, x = log(pars.train.dfr$n))$y)
pars.train.dfr$sig.smoo <- exp(predict(sig.smoo, x = log(pars.train.dfr$n))$y)
pars.train.dfr$mu.smoo <- 1 + pars.train.dfr$sig.smoo / pars.train.dfr$xi.smoo
quartz()
plot(pars.train.dfr$xi.full, pars.train.dfr$xi.smoo) 
abline(c(0,1), col = "orangered", lwd = 2)
## The low end of the shape is very scattered. I doubt about the good modelling.


## --------------------------------------------------------------------------------------------------------
## Compute the individual FAPs at each test site (and possibly train site)
## --------------------------------------------------------------------------------------------------------


smoogevfap.test.dfr <- data.frame(matrix(ncol = ncol(maxima.test.dfr),
    nrow = 1500))
colnames(smoogevfap.test.dfr) <- colnames(maxima.test.dfr)

for(ii in 1:nrow(pars.test.dfr))
#  for(ii in 5)
   {
   	nn <- as.character(pars.test.dfr$name[ii])
   	shape <- pars.test.dfr$xi.smoo[pars.test.dfr$name == nn]
   	scale <- pars.test.dfr$sig.smoo[pars.test.dfr$name == nn]
   	loc <- pars.test.dfr$mu.smoo[pars.test.dfr$name == nn]
    zass <- maxima.test.dfr[, nn]
    smoogevfap.test.dfr[, nn] <- 1 - sapply(zass, evd::pgev,
       shape = shape, loc = loc, scale = scale)
   }

smoogevfap.select.dfr <- data.frame(matrix(ncol = ncol(maxima.select.dfr),
    nrow = 1500))
colnames(smoogevfap.select.dfr) <- colnames(maxima.select.dfr)

for(ii in 1:nrow(pars.select.dfr))
#  for(ii in 5)
   {
   	nn <- as.character(pars.select.dfr$name[ii])
   	shape <- pars.select.dfr$xi.smoo[pars.select.dfr$name == nn]
   	scale <- pars.select.dfr$sig.smoo[pars.select.dfr$name == nn]
   	loc <- pars.select.dfr$mu.smoo[pars.select.dfr$name == nn]
    zass <- maxima.select.dfr[, nn]
    smoogevfap.select.dfr[, nn] <- 1 - sapply(zass, evd::pgev,
       shape = shape, loc = loc, scale = scale)
   }

smoogevfap.train.dfr <- data.frame(matrix(ncol = ncol(maxima.train.dfr),
    nrow = 1500))
colnames(smoogevfap.train.dfr) <- colnames(maxima.train.dfr)

for(ii in 1:nrow(pars.train.dfr))
#  for(ii in 5)
   {
   	nn <- as.character(pars.train.dfr$name[ii])
   	shape <- pars.train.dfr$xi.smoo[pars.train.dfr$name == nn]
   	scale <- pars.train.dfr$sig.smoo[pars.train.dfr$name == nn]
   	loc <- pars.train.dfr$mu.smoo[pars.train.dfr$name == nn]
    zass <- maxima.train.dfr[, nn]
    smoogevfap.train.dfr[, nn] <- 1 - sapply(zass, evd::pgev,
       shape = shape, loc = loc, scale = scale)
   }

save(smoogevfap.train.dfr, smoogevfap.select.dfr, smoogevfap.test.dfr,
   file = "results_baluev5/smoogev.fap.all.RObj")

pars.train.dfr$fr.smoogev95 <- NA
pars.train.dfr$fr.smoogev99 <- NA
pars.test.dfr$fr.smoogev95 <- NA
pars.test.dfr$fr.smoogev99 <- NA
pars.select.dfr$fr.smoogev95 <- NA
pars.select.dfr$fr.smoogev99 <- NA

## Compute these success rates only on the first 1000 simulations (since individual
## F^M and GEV will need the other 500 for training)
for(ii in 1:nrow(pars.test.dfr))
   {
   	pars.test.dfr$fr.smoogev95[ii] <-
   	   sum(smoogevfap.test.dfr[1:1000, as.character(pars.test.dfr$name[ii])] < 0.05) / 1000
   	pars.test.dfr$fr.smoogev99[ii] <-
   	   sum(smoogevfap.test.dfr[1:1000, as.character(pars.test.dfr$name[ii])] < 0.01) / 1000
   }
for(ii in 1:nrow(pars.select.dfr))
   {
   	pars.select.dfr$fr.smoogev95[ii] <-
   	   sum(smoogevfap.select.dfr[1:1000, as.character(pars.select.dfr$name[ii])] < 0.05) / 1000
   	pars.select.dfr$fr.smoogev99[ii] <-
   	   sum(smoogevfap.select.dfr[1:1000, as.character(pars.select.dfr$name[ii])] < 0.01) / 1000
   }
for(ii in 1:nrow(pars.train.dfr))
   {
   	pars.train.dfr$fr.smoogev95[ii] <-
   	   sum(smoogevfap.train.dfr[1:1000, as.character(pars.train.dfr$name[ii])] < 0.05) / 1000
   	pars.train.dfr$fr.smoogev99[ii] <-
   	   sum(smoogevfap.train.dfr[1:1000, as.character(pars.train.dfr$name[ii])] < 0.01) / 1000
   }


quartz(height = 7, width = 10)
par(mfrow = c(1,1), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.6,0))

col.x <- "beta"
plot(abs(pars.test.dfr[pars.test.dfr$region == "randompos", col.x]), 
   pars.test.dfr$fr.smoogev95[pars.test.dfr$region == "randompos"], 
   type = "p", pch = 16, cex = 0.7, ylim = c(0,0.15), 
   xlim = range(abs(pars.test.dfr[,col.x])), main = "Spline GEV",
   ylab = "Fraction above threshold", xlab = col.x, col = "grey25")
points(abs(pars.test.dfr[pars.test.dfr$region == "randompos", col.x]),
   pars.test.dfr$fr.smoogev99[pars.test.dfr$region == "randompos"], 
   type = "p", pch = 16, col = "grey", cex = 0.7)

points(abs(pars.test.dfr[pars.test.dfr$region == "ecl", col.x]), 
   pars.test.dfr$fr.smoogev95[pars.test.dfr$region == "ecl"], 
    pch = 16, cex = 0.7, col = "violetred")
points(abs(pars.test.dfr[pars.test.dfr$region == "ecl", col.x]), 
   pars.test.dfr$fr.smoogev99[pars.test.dfr$region == "ecl"], 
    pch = 16, cex = 0.7, col = "pink")

points(abs(pars.test.dfr[pars.test.dfr$region == "lmc", col.x]), 
   pars.test.dfr$fr.smoogev95[pars.test.dfr$region == "lmc"], 
    pch = 16, cex = 0.7, col = "orange")
points(abs(pars.test.dfr[pars.test.dfr$region == "lmc", col.x]), 
   pars.test.dfr$fr.smoogev99[pars.test.dfr$region == "lmc"], 
    pch = 16, cex = 0.7, col = "gold")

points(abs(pars.test.dfr[pars.test.dfr$region == "ln.ha", col.x]), 
   pars.test.dfr$fr.smoogev95[pars.test.dfr$region == "ln.ha"], 
    pch = 16, cex = 0.7, col = "deepskyblue")
points(abs(pars.test.dfr[pars.test.dfr$region == "ln.ha", col.x]), 
   pars.test.dfr$fr.smoogev99[pars.test.dfr$region == "ln.ha"], 
    pch = 16, cex = 0.7, col = "skyblue")

abline(h = 0.05, lwd = 2)
abline(h = 0.01, lty = 2, lwd = 2)

rm(xi.smoo, sig.smoo, smoogevfap.train.dfr, smoogevfap.select.dfr, smoogevfap.test.dfr)
save.image(WORKSPACE.NAME)


## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------
## Models: thin-plate splines, now try with other variables 
## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------

## --------------------------------------------------------------------------------------------------------
## Thin-plate spline for M
## --------------------------------------------------------------------------------------------------------

pars.train.dfr$log.M.full <- log(pars.train.dfr$M.full)
pars.test.dfr$log.M.full <- log(pars.test.dfr$M.full)
pars.select.dfr$log.M.full <- log(pars.select.dfr$M.full)

## A thin-plate spline from package fields:
surface.M.tps0 <- Tps(pars.train.dfr[,c("log.n", "sum.of.nine.spwin")], Y = pars.train.dfr$log.M.full,
   m = 3, lambda = 0.00001)
surface.M.tps1 <- Tps(pars.train.dfr[,c("log.n", "vart")], Y = pars.train.dfr$log.M.full,
   m = 3, lambda = 0.00001)
surface.M.tps2 <- Tps(pars.train.dfr[,c("log.n", "highest3sum.spwin")], Y = pars.train.dfr$log.M.full,
   m = 3, lambda = 0.00001)
surface.M.tps3 <- Tps(pars.train.dfr[,c("log.n", "sum.spwin")], Y = pars.train.dfr$log.M.full,
   m = 3, lambda = 0.00001)

quartz(width = 7, height = 7)
c.tmp12 <- rgb(green = (pars.train.dfr[,"log.M.full"] - min(pars.train.dfr[,"log.M.full"])) /
   (max(pars.train.dfr[,"log.M.full"]) - min(pars.train.dfr[,"log.M.full"])), 
   blue = 1 - (pars.train.dfr[,"log.M.full"] - min(pars.train.dfr[,"log.M.full"])) /
   (max(pars.train.dfr[,"log.M.full"]) - min(pars.train.dfr[,"log.M.full"])),
   red = 0, alpha = 0.4)
c.tmp34 <- rgb(green = (pars.test.dfr[,"log.M.full"] - min(pars.test.dfr[,"log.M.full"])) /
   (max(pars.test.dfr[,"log.M.full"]) - min(pars.test.dfr[,"log.M.full"])), 
   blue = 1 - (pars.test.dfr[,"log.M.full"] - min(pars.test.dfr[,"log.M.full"])) /
   (max(pars.test.dfr[,"log.M.full"]) - min(pars.test.dfr[,"log.M.full"])),
   red = 0, alpha = 0.4)
c.tmp56 <- rgb(green = (pars.select.dfr[,"log.M.full"] - min(pars.select.dfr[,"log.M.full"])) /
   (max(pars.select.dfr[,"log.M.full"]) - min(pars.select.dfr[,"log.M.full"])), 
   blue = 1 - (pars.select.dfr[,"log.M.full"] - min(pars.select.dfr[,"log.M.full"])) /
   (max(pars.select.dfr[,"log.M.full"]) - min(pars.select.dfr[,"log.M.full"])),
   red = 0, alpha = 0.4)
par(mfrow = c(2,2), mar = c(2.5,2.5,2, 1), mgp = c(1.5,0.6,0))
surface(surface.M.tps0, type = "c")
points(pars.train.dfr$log.n, pars.train.dfr$sum.of.nine.spwin, pch = 15, cex = 0.5, col = c.tmp12)
points(pars.test.dfr$log.n, pars.test.dfr$sum.of.nine.spwin, pch = 15, cex = 0.5, col = c.tmp34)
points(pars.select.dfr$log.n, pars.select.dfr$sum.of.nine.spwin, pch = 15, cex = 0.5, col = c.tmp56)
surface(surface.M.tps1, type = "c")
points(pars.train.dfr$log.n, pars.train.dfr$vart, pch = 15, cex = 0.5, col = c.tmp12)
points(pars.test.dfr$log.n, pars.test.dfr$vart, pch = 15, cex = 0.5, col = c.tmp34)
points(pars.select.dfr$log.n, pars.select.dfr$vart, pch = 15, cex = 0.5, col = c.tmp56)
surface(surface.M.tps2, type = "c")
points(pars.train.dfr$log.n, pars.train.dfr$highest3sum.spwin, pch = 15, cex = 0.5, col = c.tmp12)
points(pars.test.dfr$log.n, pars.test.dfr$highest3sum.spwin, pch = 15, cex = 0.5, col = c.tmp34)
points(pars.select.dfr$log.n, pars.select.dfr$highest3sum.spwin, pch = 15, cex = 0.5, col = c.tmp56)
surface(surface.M.tps3, type = "c")
points(pars.train.dfr$log.n, pars.train.dfr$sum.spwin, pch = 15, cex = 0.5, col = c.tmp12)
points(pars.test.dfr$log.n, pars.test.dfr$sum.spwin, pch = 15, cex = 0.5, col = c.tmp34)
points(pars.select.dfr$log.n, pars.select.dfr$highest3sum.spwin, pch = 15, cex = 0.5, col = c.tmp56)
## vart looks a bit messy, the other two give similar contour lines (they are very correlated.)

## Prediction of M

pred.test.fadm <- data.frame(name = pars.test.dfr$name, log.n = NA,
   vart = NA, highest3sum.spwin = NA, sum.of.nine.spwin = NA, 
   sum.spwin = NA, log.M.full = pars.test.dfr$log.M.full,
   tps0.fitted.log.M.full = NA, tps1.fitted.log.M.full = NA, 
   tps2.fitted.log.M.full = NA, tps3.fitted.log.M.full = NA)
#ii <- 1   
for(ii in seq(1, nrow(pred.test.fadm), by = 2))
   {
	p.tmp <- predict.surface(surface.M.tps0, list(
	   log.n = pars.test.dfr$log.n[ii:(ii+1)],
	   sum.of.nine.spwin = pars.test.dfr$sum.of.nine.spwin[ii:(ii+1)]), extrap = TRUE)
	pred.test.fadm[ii:(ii+1), "log.n"] <- p.tmp$x
	pred.test.fadm[ii:(ii+1), "sum.of.nine.spwin"] <- p.tmp$y
	pred.test.fadm[ii:(ii+1), "tps0.fitted.log.M.full"] <- diag(p.tmp$z)
	p.tmp <- predict.surface(surface.M.tps1, list(
	   log.n = pars.test.dfr$log.n[ii:(ii+1)],
	   vart = pars.test.dfr$vart[ii:(ii+1)]), extrap = TRUE)
	pred.test.fadm[ii:(ii+1), "vart"] <- p.tmp$y
	pred.test.fadm[ii:(ii+1), "tps1.fitted.log.M.full"] <- diag(p.tmp$z)
	p.tmp <- predict.surface(surface.M.tps2, list(
	   log.n = pars.test.dfr$log.n[ii:(ii+1)],
	   highest3sum.spwin = pars.test.dfr$highest3sum.spwin[ii:(ii+1)]), extrap = TRUE)
	pred.test.fadm[ii:(ii+1), "highest3sum.spwin"] <- p.tmp$y
	pred.test.fadm[ii:(ii+1), "tps2.fitted.log.M.full"] <- diag(p.tmp$z)
	p.tmp <- predict.surface(surface.M.tps3, list(
	   log.n = pars.test.dfr$log.n[ii:(ii+1)],
	   sum.spwin = pars.test.dfr$sum.spwin[ii:(ii+1)]), extrap = TRUE)
	pred.test.fadm[ii:(ii+1), "sum.spwin"] <- p.tmp$y
	pred.test.fadm[ii:(ii+1), "tps3.fitted.log.M.full"] <- diag(p.tmp$z)
   }

pred.select.fadm <- data.frame(name = pars.select.dfr$name, log.n = NA,
   vart = NA, sum.spwin = NA, sum.of.nine.spwin = NA, 
   sum.spwin = NA, log.M.full = pars.select.dfr$log.M.full,
   tps0.fitted.log.M.full = NA, tps1.fitted.log.M.full = NA, 
   tps2.fitted.log.M.full = NA, tps3.fitted.log.M.full = NA)
#ii <- 1   
for(ii in seq(1, nrow(pred.select.fadm), by = 2))
   {
	p.tmp <- predict.surface(surface.M.tps0, list(
	   log.n = pars.select.dfr$log.n[ii:(ii+1)],
	   sum.of.nine.spwin = pars.select.dfr$sum.of.nine.spwin[ii:(ii+1)]), extrap = TRUE)
	pred.select.fadm[ii:(ii+1), "log.n"] <- p.tmp$x
	pred.select.fadm[ii:(ii+1), "sum.of.nine.spwin"] <- p.tmp$y
	pred.select.fadm[ii:(ii+1), "tps0.fitted.log.M.full"] <- diag(p.tmp$z)
	p.tmp <- predict.surface(surface.M.tps1, list(
	   log.n = pars.select.dfr$log.n[ii:(ii+1)],
	   vart = pars.select.dfr$vart[ii:(ii+1)]), extrap = TRUE)
	pred.select.fadm[ii:(ii+1), "vart"] <- p.tmp$y
	pred.select.fadm[ii:(ii+1), "tps1.fitted.log.M.full"] <- diag(p.tmp$z)
	p.tmp <- predict.surface(surface.M.tps2, list(
	   log.n = pars.select.dfr$log.n[ii:(ii+1)],
	   highest3sum.spwin = pars.select.dfr$highest3sum.spwin[ii:(ii+1)]), extrap = TRUE)
	pred.select.fadm[ii:(ii+1), "highest3sum.spwin"] <- p.tmp$y
	pred.select.fadm[ii:(ii+1), "tps2.fitted.log.M.full"] <- diag(p.tmp$z)
	p.tmp <- predict.surface(surface.M.tps3, list(
	   log.n = pars.select.dfr$log.n[ii:(ii+1)],
	   sum.spwin = pars.select.dfr$sum.spwin[ii:(ii+1)]), extrap = TRUE)
	pred.select.fadm[ii:(ii+1), "sum.spwin"] <- p.tmp$y
	pred.select.fadm[ii:(ii+1), "tps3.fitted.log.M.full"] <- diag(p.tmp$z)
   }

quartz(height = 4, width = 8)
par(mfrow = c(1,2), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.6,0))
plot(pred.test.fadm$log.M.full, pred.test.fadm$tps0.fitted.log.M.full, pch = 16, cex = 0.5, col = "orangered",
   xlim = range(c(pred.test.fadm$log.M.full, pred.test.fadm$tps0.fitted.log.M.full, 
   pred.test.fadm$tps1.fitted.log.M.full), na.rm = TRUE),
   ylim = range(c(pred.test.fadm$log.M.full, pred.test.fadm$tps0.fitted.log.M.full, 
   pred.test.fadm$tps1.fitted.log.M.full), na.rm = TRUE))
#points(pred.test.fadm$log.M.full, pred.test.fadm$tps1.fitted.log.M.full, pch = 3, cex = 0.5, col = "blue")
#points(pred.test.fadm$log.M.full, pred.test.fadm$tps2.fitted.log.M.full, pch = 3, cex = 0.5, col = "deepskyblue")
points(pred.test.fadm$log.M.full, pred.test.fadm$tps3.fitted.log.M.full, pch = 3, cex = 0.5, col = "green")
#points(pars.test.dfr$log.M.full, log(pars.test.dfr$M.full.rf), pch = 3, cex = 0.5, col = "gold")
abline(c(0,1))
plot(pred.select.fadm$log.M.full, pred.select.fadm$tps0.fitted.log.M.full, pch = 16, cex = 0.5, col = "orangered",
   xlim = range(c(pred.select.fadm$log.M.full, pred.select.fadm$tps0.fitted.log.M.full, 
   pred.select.fadm$tps1.fitted.log.M.full), na.rm = TRUE),
   ylim = range(c(pred.select.fadm$log.M.full, pred.select.fadm$tps0.fitted.log.M.full, 
   pred.select.fadm$tps1.fitted.log.M.full), na.rm = TRUE))
#points(pred.select.fadm$log.M.full, pred.select.fadm$tps1.fitted.log.M.full, pch = 3, cex = 0.5, col = "blue")
#points(pred.select.fadm$log.M.full, pred.select.fadm$tps2.fitted.log.M.full, pch = 3, cex = 0.5, col = "deepskyblue")
points(pred.select.fadm$log.M.full, pred.select.fadm$tps3.fitted.log.M.full, pch = 3, cex = 0.5, col = "green")
#points(pars.select.dfr$log.M.full, log(pars.select.dfr$M.full.rf), pch = 3, cex = 0.5, col = "gold")
abline(c(0,1))
## Neither is very good. The best may be sum.of.nine.spwin, but far less line-like than either for GEV 
## parameters, or for quantiles. Almost as good is highest3sum.spwin. RF misses the lower tail, but 
## otherwise seeems as good as the best of the tps.

## Check the mse:
sqrt(mean((pred.test.fadm$log.M.full - pred.test.fadm$tps0.fitted.log.M.full)^2))   ## 0.1012966
sqrt(mean((pred.test.fadm$log.M.full - pred.test.fadm$tps1.fitted.log.M.full)^2))   ## 0.2158934
sqrt(mean((pred.test.fadm$log.M.full - pred.test.fadm$tps2.fitted.log.M.full)^2))   ## 0.1051636
sqrt(mean((pred.test.fadm$log.M.full - pred.test.fadm$tps3.fitted.log.M.full)^2))   ## 0.1096512
sqrt(mean((pars.test.dfr$log.M.full - log(pars.test.dfr$M.full.rf))^2))             ## 0.07736174
sqrt(mean((pred.select.fadm$log.M.full - pred.select.fadm$tps0.fitted.log.M.full)^2))   ## 0.0971162
sqrt(mean((pred.select.fadm$log.M.full - pred.select.fadm$tps1.fitted.log.M.full)^2))   ## 0.2110031
sqrt(mean((pred.select.fadm$log.M.full - pred.select.fadm$tps2.fitted.log.M.full)^2))   ## 0.09784272
sqrt(mean((pars.select.dfr$log.M.full - log(pars.select.dfr$M.full.rf))^2))             ## 0.07078074
## In general, RandomForest does the best, but it misses specifically the low-M tail, so
## maye despite a larger scatter, tps0 is preferable.

pred.test.fadm$tps0.fitted.M.full <- exp(pred.test.fadm$tps0.fitted.log.M.full)
pred.select.fadm$tps0.fitted.M.full <- exp(pred.select.fadm$tps0.fitted.log.M.full)

pred.train.fadm <- data.frame(name = pars.train.dfr$name, 
   log.n = pars.train.dfr$log.n, vart = pars.train.dfr$at12, 
   sum.of.nine.spwin = pars.train.dfr$sum.of.nine.spwin,
   sum.spwin = pars.train.dfr$sum.spwin,
   log.M.full = pars.train.dfr$log.M.full,
   tps0.fitted.log.M.full = surface.M.tps0$fitted.values,
   tps0.fitted.M.full = exp(surface.M.tps0$fitted.values),
   tps1.fitted.log.M.full = surface.M.tps1$fitted.values,
   tps1.fitted.M.full = exp(surface.M.tps1$fitted.values),
   tps2.fitted.log.M.full = surface.M.tps2$fitted.values,
   tps2.fitted.M.full = exp(surface.M.tps2$fitted.values),
   tps3.fitted.log.M.full = surface.M.tps3$fitted.values,
   tps3.fitted.M.full = exp(surface.M.tps3$fitted.values))


## --------------------------------------------------------------------------------------------------------
## Compute the FAPs with the estimated Ms
## --------------------------------------------------------------------------------------------------------

fadm.tps0.fap.test <- matrix(ncol = ncol(maxima.test.dfr), nrow = nrow(maxima.test.dfr))
dimnames(fadm.tps0.fap.test) <- list(NULL, colnames(maxima.test.dfr))
for(ii in 1:nrow(pred.test.fadm))
   {
   	nn <- as.character(pred.test.fadm$name[ii])
   	nobs <- exp(pred.test.fadm$log.n[as.character(pred.test.fadm$name) == nn])
   	mest <- pred.test.fadm$tps0.fitted.M.full[as.character(pred.test.fadm$name) == nn]
   	fadm.tps0.fap.test[,nn] <- sapply(maxima.test.dfr[,nn], schwczfap.fun, mest = mest, N = nobs)
   }

fadm.tps0.fap.select <- matrix(ncol = ncol(maxima.select.dfr), nrow = nrow(maxima.select.dfr))
dimnames(fadm.tps0.fap.select) <- list(NULL, colnames(maxima.select.dfr))
for(ii in 1:nrow(pred.select.fadm))
   {
   	nn <- as.character(pred.select.fadm$name[ii])
   	nobs <- exp(pred.select.fadm$log.n[as.character(pred.select.fadm$name) == nn])
   	mest <- pred.select.fadm$tps0.fitted.M.full[as.character(pred.select.fadm$name) == nn]
   	fadm.tps0.fap.select[,nn] <- sapply(maxima.select.dfr[,nn], schwczfap.fun, mest = mest, N = nobs)
   }

fadm.tps0.fap.train <- matrix(ncol = ncol(maxima.train.dfr), nrow = nrow(maxima.train.dfr))
dimnames(fadm.tps0.fap.train) <- list(NULL, colnames(maxima.train.dfr))
for(ii in 1:ncol(maxima.train.dfr))
   {
   	nn <- as.character(pred.train.fadm$name[ii])
   	nobs <- exp(pred.train.fadm$log.n[as.character(pred.train.fadm$name) == nn])
   	mest <- pred.train.fadm$tps0.fitted.M.full[as.character(pred.train.fadm$name) == nn]
	fadm.tps0.fap.train[,nn] <- sapply(maxima.train.dfr[,nn], schwczfap.fun, mest = mest, N = nobs)
   }

range(apply(fadm.tps0.fap.test, 2, function(vec) sum(vec < 0.05) / length(vec)), na.rm = TRUE)
range(apply(fadm.tps0.fap.train, 2, function(vec) sum(vec < 0.05) / length(vec)), na.rm = TRUE)
range(apply(fadm.tps0.fap.test, 2, function(vec) sum(vec < 0.01) / length(vec)), na.rm = TRUE)
range(apply(fadm.tps0.fap.train, 2, function(vec) sum(vec < 0.01) / length(vec)), na.rm = TRUE)

## Save all results, all modelling, then keep only the selected model for potential later use.

rm(p.tmp)
dfr.tmp <- merge(pars.train.dfr, pred.train.fadm[,c("name","tps0.fitted.M.full")])
dim(pars.train.dfr)
pars.train.dfr <- dfr.tmp
dfr.tmp <- merge(pars.test.dfr, pred.test.fadm[,c("name","log.M.full","tps0.fitted.M.full")], all.x = TRUE)
dim(pars.test.dfr)
pars.test.dfr <- dfr.tmp
dfr.tmp <- merge(pars.select.dfr, pred.select.fadm[,c("name","log.M.full","tps0.fitted.M.full")], all.x = TRUE)
dim(pars.select.dfr)
pars.select.dfr <- dfr.tmp


## --------------------------------------------------------------------------------------------------------
## Check the number of significant values in the test set
## --------------------------------------------------------------------------------------------------------

pars.train.dfr$fr.tps0fadm.sign05 <- NA
pars.train.dfr$fr.tps0fadm.sign01 <- NA
pars.select.dfr$fr.tps0fadm.sign05 <- NA
pars.select.dfr$fr.tps0fadm.sign01 <- NA
pars.test.dfr$fr.tps0fadm.sign05 <- NA
pars.test.dfr$fr.tps0fadm.sign01 <- NA

for(ii in 1:nrow(pars.train.dfr))
   {
   	pars.train.dfr$fr.tps0fadm.sign05[ii] <-
   	   sum(fadm.tps0.fap.train[, as.character(pars.train.dfr$name[ii])] < 0.05) / 1500
   	pars.train.dfr$fr.tps0fadm.sign01[ii] <-
   	   sum(fadm.tps0.fap.train[, as.character(pars.train.dfr$name[ii])] < 0.01) / 1500
   }
for(ii in 1:nrow(pars.select.dfr))
   {
   	pars.select.dfr$fr.tps0fadm.sign05[ii] <-
   	   sum(fadm.tps0.fap.select[, as.character(pars.select.dfr$name[ii])] < 0.05) / 1500
   	pars.select.dfr$fr.tps0fadm.sign01[ii] <-
   	   sum(fadm.tps0.fap.select[, as.character(pars.select.dfr$name[ii])] < 0.01) / 1500
   }
for(ii in 1:nrow(pars.test.dfr))
   {
   	pars.test.dfr$fr.tps0fadm.sign05[ii] <-
   	   sum(fadm.tps0.fap.test[, as.character(pars.test.dfr$name[ii])] < 0.05) / 1500
   	pars.test.dfr$fr.tps0fadm.sign01[ii] <-
   	   sum(fadm.tps0.fap.test[, as.character(pars.test.dfr$name[ii])] < 0.01) / 1500
   }

## The final model should be the rf, it has less distortions and systematic biases.
## Voilà what I saved, contains mfull.rf.final, mfull.attr.final:
# load(file = "results_baluev5/fadm.rfmodel.RObj")  
## I don't keep it in the workspace, it is big.


quartz(height = 7, width = 10)
par(mfrow = c(1,1), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.6,0))

col.x <- "beta"
plot(abs(pars.test.dfr[pars.test.dfr$region == "randompos", col.x]), 
   pars.test.dfr$fr.tps0fadm.sign05[pars.test.dfr$region == "randompos"], 
   type = "p", pch = 16, cex = 0.7, ylim = c(0,0.15), 
   xlim = range(abs(pars.test.dfr[,col.x])), main = "Thin-plate spline F^M",
   ylab = "Fraction above threshold", xlab = col.x, col = "grey25")
points(abs(pars.test.dfr[pars.test.dfr$region == "randompos", col.x]),
   pars.test.dfr$fr.tps0fadm.sign01[pars.test.dfr$region == "randompos"], 
   type = "p", pch = 16, col = "grey", cex = 0.7)

points(abs(pars.test.dfr[pars.test.dfr$region == "ecl", col.x]), 
   pars.test.dfr$fr.tps0fadm.sign05[pars.test.dfr$region == "ecl"], 
    pch = 16, cex = 0.7, col = "violetred")
points(abs(pars.test.dfr[pars.test.dfr$region == "ecl", col.x]), 
   pars.test.dfr$fr.tps0fadm.sign01[pars.test.dfr$region == "ecl"], 
    pch = 16, cex = 0.7, col = "pink")

points(abs(pars.test.dfr[pars.test.dfr$region == "lmc", col.x]), 
   pars.test.dfr$fr.tps0fadm.sign05[pars.test.dfr$region == "lmc"], 
    pch = 16, cex = 0.7, col = "orange")
points(abs(pars.test.dfr[pars.test.dfr$region == "lmc", col.x]), 
   pars.test.dfr$fr.tps0fadm.sign01[pars.test.dfr$region == "lmc"], 
    pch = 16, cex = 0.7, col = "gold")

points(abs(pars.test.dfr[pars.test.dfr$region == "ln.ha", col.x]), 
   pars.test.dfr$fr.tps0fadm.sign05[pars.test.dfr$region == "ln.ha"], 
    pch = 16, cex = 0.7, col = "deepskyblue")
points(abs(pars.test.dfr[pars.test.dfr$region == "ln.ha", col.x]), 
   pars.test.dfr$fr.tps0fadm.sign01[pars.test.dfr$region == "ln.ha"], 
    pch = 16, cex = 0.7, col = "skyblue")

abline(h = 0.05, lwd = 2)
abline(h = 0.01, lty = 2, lwd = 2)


save(fadm.tps0.fap.test, fadm.tps0.fap.select, fadm.tps0.fap.train, file =
   "results_baluev5/fadm.tps0.fap.both.RObj")
save(surface.M.tps0, surface.M.tps1, surface.M.tps2, surface.M.tps3, 
   file = "results_baluev5/fadm.tpsmodels.RObj")
save(pred.train.fadm, pred.test.fadm, pred.select.fadm, file = "results_baluev5/fadm.tps0.predictions.RObj")
rm(fadm.tps0.fap.test, fadm.tps0.fap.train, fadm.tps0.fap.select, surface.M.tps0, surface.M.tps1, surface.M.tps2,
   pred.test.fadm, pred.select.fadm, pred.train.fadm)

save.image(WORKSPACE.NAME)


## --------------------------------------------------------------------------------------------------------
## For GEV
## --------------------------------------------------------------------------------------------------------

## Initial model with all the variables

## A thin-plate spline from package fields:
surface.xi.tps0 <- Tps(pars.train.dfr[,c("log.n", "sum.spwin")], Y = pars.train.dfr$log.xi.full, m = 3, lambda = 0.00001)
surface.xi.tps1 <- Tps(pars.train.dfr[,c("log.n", "at52")], Y = pars.train.dfr$log.xi.full, m = 3, lambda = 0.00001)
surface.sig.tps0 <- Tps(pars.train.dfr[,c("log.n", "sum.spwin")], Y = pars.train.dfr$log.sig.full, m = 3, lambda = 0.00001)
surface.sig.tps1 <- Tps(pars.train.dfr[,c("log.n", "at52")], Y = pars.train.dfr$log.sig.full, m = 3, lambda = 0.00001)
surface.gev.tps0 <- list(log.xi.full = surface.xi.tps0, log.sig.full = surface.sig.tps0)
surface.gev.tps1 <- list(log.xi.full = surface.xi.tps1, log.sig.full = surface.sig.tps1)
rm(surface.xi.tps0, surface.xi.tps1, surface.sig.tps0,surface.sig.tps1, surface.tps0, surface.tps1)
save(surface.gev.tps0, file = "results_baluev5/surface.gev.tps0.RObj")
save(surface.gev.tps1, file = "results_baluev5/surface.gev.tps1.RObj")

quartz(height = 7, width = 7)
par(mfrow = c(2,2), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.6,0))
for(ii in c("log.xi.full","log.sig.full"))
   {
	c.tmp12 <- rgb(green = (pars.train.dfr[,ii] - min(pars.train.dfr[,ii])) /
	   (max(pars.train.dfr[,ii]) - min(pars.train.dfr[,ii])), 
	   blue = 1 - (pars.train.dfr[,ii] - min(pars.train.dfr[,ii])) /
	   (max(pars.train.dfr[,ii]) - min(pars.train.dfr[,ii])),
	   red = 0, alpha = 0.4)
	c.tmp34 <- rgb(green = (pars.test.dfr[,ii] - min(pars.test.dfr[,ii], na.rm = TRUE)) /
	   (max(pars.test.dfr[,ii], na.rm = TRUE) - min(pars.test.dfr[,ii], na.rm = TRUE)), 
	   blue = 1 - (pars.test.dfr[,ii] - min(pars.test.dfr[,ii], na.rm = TRUE)) /
	   (max(pars.test.dfr[,ii], na.rm = TRUE) - min(pars.test.dfr[,ii], na.rm = TRUE)),
	   red = 0, alpha = 0.4)
	c.tmp56 <- rgb(green = (pars.select.dfr[,ii] - min(pars.select.dfr[,ii], na.rm = TRUE)) /
	   (max(pars.select.dfr[,ii], na.rm = TRUE) - min(pars.select.dfr[,ii], na.rm = TRUE)), 
	   blue = 1 - (pars.select.dfr[,ii] - min(pars.select.dfr[,ii], na.rm = TRUE)) /
	   (max(pars.select.dfr[,ii], na.rm = TRUE) - min(pars.select.dfr[,ii], na.rm = TRUE)),
	   red = 0, alpha = 0.4)
    surface(surface.gev.tps0[[ii]], type = "c")
	points(pars.train.dfr$log.n, pars.train.dfr$sum.spwin, pch = 15, cex = 0.7, col = c.tmp12)
	points(pars.test.dfr$log.n, pars.test.dfr$sum.spwin, pch = 15, cex = 0.7, col = c.tmp34)
	points(pars.select.dfr$log.n, pars.select.dfr$sum.spwin, pch = 15, cex = 0.7, col = c.tmp56)
   }
for(ii in c("log.xi.full","log.sig.full"))
   {
	c.tmp12 <- rgb(green = (pars.train.dfr[,ii] - min(pars.train.dfr[,ii])) /
	   (max(pars.train.dfr[,ii]) - min(pars.train.dfr[,ii])), 
	   blue = 1 - (pars.train.dfr[,ii] - min(pars.train.dfr[,ii])) /
	   (max(pars.train.dfr[,ii]) - min(pars.train.dfr[,ii])),
	   red = 0, alpha = 0.4)
	c.tmp34 <- rgb(green = (pars.test.dfr[,ii] - min(pars.test.dfr[,ii], na.rm = TRUE)) /
	   (max(pars.test.dfr[,ii], na.rm = TRUE) - min(pars.test.dfr[,ii], na.rm = TRUE)), 
	   blue = 1 - (pars.test.dfr[,ii] - min(pars.test.dfr[,ii], na.rm = TRUE)) /
	   (max(pars.test.dfr[,ii], na.rm = TRUE) - min(pars.test.dfr[,ii], na.rm = TRUE)),
	   red = 0, alpha = 0.4)
	c.tmp56 <- rgb(green = (pars.select.dfr[,ii] - min(pars.select.dfr[,ii], na.rm = TRUE)) /
	   (max(pars.select.dfr[,ii], na.rm = TRUE) - min(pars.select.dfr[,ii], na.rm = TRUE)), 
	   blue = 1 - (pars.select.dfr[,ii] - min(pars.select.dfr[,ii], na.rm = TRUE)) /
	   (max(pars.select.dfr[,ii], na.rm = TRUE) - min(pars.select.dfr[,ii], na.rm = TRUE)),
	   red = 0, alpha = 0.4)
    surface(surface.gev.tps1[[ii]], type = "c")
	points(pars.train.dfr$log.n, pars.train.dfr$at52, pch = 15, cex = 0.7, col = c.tmp12)
	points(pars.test.dfr$log.n, pars.test.dfr$at52, pch = 15, cex = 0.7, col = c.tmp34)
	points(pars.select.dfr$log.n, pars.select.dfr$at52, pch = 15, cex = 0.7, col = c.tmp56)
   }


## Prediction of the parameters

pred.test.gev <- data.frame(name = pars.test.dfr$name, log.n = NA, at52 = NA, sum.spwin = NA,
   log.xi.full = pars.test.dfr$log.xi.full,
   tps0.fitted.log.xi.full = NA, tps1.fitted.log.xi.full = NA, 
   log.sig.full = pars.test.dfr$log.sig.full,
   tps0.fitted.log.sig.full = NA, tps1.fitted.log.sig.full = NA,
   tps0.fitted.mu.full = NA, tps1.fitted.mu.full = NA)
#ii <- 1
for(ii in seq(1, nrow(pred.test.gev), by = 2))   ## Fortunately, nrow() is an even number.
   {
   	## for the shape parameter:
	p.tmp <- predict.surface(surface.gev.tps0$log.xi.full, list(
	   log.n = pars.test.dfr$log.n[ii:(ii+1)],
	   sum.spwin = pars.test.dfr$sum.spwin[ii:(ii+1)]), extrap = TRUE)
	pred.test.gev[ii:(ii+1), "log.n"] <- p.tmp$x
	pred.test.gev[ii:(ii+1), "sum.spwin"] <- p.tmp$y
	pred.test.gev[ii:(ii+1), "tps0.fitted.log.xi.full"] <- diag(p.tmp$z)
	p.tmp <- predict.surface(surface.gev.tps1$log.xi.full, list(
	   log.n = pars.test.dfr$log.n[ii:(ii+1)],
	   at52 = pars.test.dfr$at52[ii:(ii+1)]), extrap = TRUE)
	pred.test.gev[ii:(ii+1), "at52"] <- p.tmp$y
	pred.test.gev[ii:(ii+1), "tps1.fitted.log.xi.full"] <- diag(p.tmp$z)
   	## for the scale parameter:
	p.tmp <- predict.surface(surface.gev.tps0$log.sig.full, list(
	   log.n = pars.test.dfr$log.n[ii:(ii+1)],
	   sum.spwin = pars.test.dfr$sum.spwin[ii:(ii+1)]), extrap = TRUE)
	pred.test.gev[ii:(ii+1), "log.n"] <- p.tmp$x
	pred.test.gev[ii:(ii+1), "sum.spwin"] <- p.tmp$y
	pred.test.gev[ii:(ii+1), "tps0.fitted.log.sig.full"] <- diag(p.tmp$z)
	p.tmp <- predict.surface(surface.gev.tps1$log.sig.full, list(
	   log.n = pars.test.dfr$log.n[ii:(ii+1)],
	   at52 = pars.test.dfr$at52[ii:(ii+1)]), extrap = TRUE)
	pred.test.gev[ii:(ii+1), "at52"] <- p.tmp$y
	pred.test.gev[ii:(ii+1), "tps1.fitted.log.sig.full"] <- diag(p.tmp$z)
	pred.test.gev[ii:(ii+1), "tps0.fitted.mu.full"] <- 1 -
	   exp(pred.test.gev[ii:(ii+1), "tps0.fitted.log.sig.full"]) / 
	   exp(pred.test.gev[ii:(ii+1), "tps0.fitted.log.xi.full"])
	pred.test.gev[ii:(ii+1), "tps1.fitted.mu.full"] <- 1 -
	   exp(pred.test.gev[ii:(ii+1), "tps1.fitted.log.sig.full"]) / 
	   exp(pred.test.gev[ii:(ii+1), "tps1.fitted.log.xi.full"])
   }
pred.test.gev$tps0.fitted.xi.full <- - exp(pred.test.gev$tps0.fitted.log.xi.full)
pred.test.gev$tps1.fitted.xi.full <- - exp(pred.test.gev$tps1.fitted.log.xi.full)
pred.test.gev$tps0.fitted.sig.full <- exp(pred.test.gev$tps0.fitted.log.sig.full)
pred.test.gev$tps1.fitted.sig.full <- exp(pred.test.gev$tps1.fitted.log.sig.full)
pred.test.gev$mu.full <- pars.test.dfr$mu.full
pred.test.gev$log.mu.full <- pars.test.dfr$log.mu.full

pred.select.gev <- data.frame(name = pars.select.dfr$name, log.n = NA, at52 = NA, sum.spwin = NA,
   log.xi.full = pars.select.dfr$log.xi.full,
   tps0.fitted.log.xi.full = NA, tps1.fitted.log.xi.full = NA, 
   log.sig.full = pars.select.dfr$log.sig.full,
   tps0.fitted.log.sig.full = NA, tps1.fitted.log.sig.full = NA,
   tps0.fitted.mu.full = NA, tps1.fitted.mu.full = NA)
#ii <- 1
for(ii in seq(1, nrow(pred.select.gev), by = 2))   ## Fortunately, nrow() is an even number.
   {
   	## for the shape parameter:
	p.tmp <- predict.surface(surface.gev.tps0$log.xi.full, list(
	   log.n = pars.select.dfr$log.n[ii:(ii+1)],
	   sum.spwin = pars.select.dfr$sum.spwin[ii:(ii+1)]), extrap = TRUE)
	pred.select.gev[ii:(ii+1), "log.n"] <- p.tmp$x
	pred.select.gev[ii:(ii+1), "sum.spwin"] <- p.tmp$y
	pred.select.gev[ii:(ii+1), "tps0.fitted.log.xi.full"] <- diag(p.tmp$z)
	p.tmp <- predict.surface(surface.gev.tps1$log.xi.full, list(
	   log.n = pars.select.dfr$log.n[ii:(ii+1)],
	   at52 = pars.select.dfr$at52[ii:(ii+1)]), extrap = TRUE)
	pred.select.gev[ii:(ii+1), "at52"] <- p.tmp$y
	pred.select.gev[ii:(ii+1), "tps1.fitted.log.xi.full"] <- diag(p.tmp$z)
   	## for the scale parameter:
	p.tmp <- predict.surface(surface.gev.tps0$log.sig.full, list(
	   log.n = pars.select.dfr$log.n[ii:(ii+1)],
	   sum.spwin = pars.select.dfr$sum.spwin[ii:(ii+1)]), extrap = TRUE)
	pred.select.gev[ii:(ii+1), "log.n"] <- p.tmp$x
	pred.select.gev[ii:(ii+1), "sum.spwin"] <- p.tmp$y
	pred.select.gev[ii:(ii+1), "tps0.fitted.log.sig.full"] <- diag(p.tmp$z)
	p.tmp <- predict.surface(surface.gev.tps1$log.sig.full, list(
	   log.n = pars.select.dfr$log.n[ii:(ii+1)],
	   at52 = pars.select.dfr$at52[ii:(ii+1)]), extrap = TRUE)
	pred.select.gev[ii:(ii+1), "at52"] <- p.tmp$y
	pred.select.gev[ii:(ii+1), "tps1.fitted.log.sig.full"] <- diag(p.tmp$z)
	pred.select.gev[ii:(ii+1), "tps0.fitted.mu.full"] <- 1 -
	   exp(pred.select.gev[ii:(ii+1), "tps0.fitted.log.sig.full"]) / 
	   exp(pred.select.gev[ii:(ii+1), "tps0.fitted.log.xi.full"])
	pred.select.gev[ii:(ii+1), "tps1.fitted.mu.full"] <- 1 -
	   exp(pred.select.gev[ii:(ii+1), "tps1.fitted.log.sig.full"]) / 
	   exp(pred.select.gev[ii:(ii+1), "tps1.fitted.log.xi.full"])
   }
pred.select.gev$tps0.fitted.xi.full <- - exp(pred.select.gev$tps0.fitted.log.xi.full)
pred.select.gev$tps1.fitted.xi.full <- - exp(pred.select.gev$tps1.fitted.log.xi.full)
pred.select.gev$tps0.fitted.sig.full <- exp(pred.select.gev$tps0.fitted.log.sig.full)
pred.select.gev$tps1.fitted.sig.full <- exp(pred.select.gev$tps1.fitted.log.sig.full)
pred.select.gev$mu.full = pars.select.dfr$mu.full
pred.select.gev$log.mu.full <- pars.select.dfr$log.mu.full

pred.train.gev <- data.frame(name = pars.train.dfr$name, 
   log.n = pars.train.dfr$log.n, at52 = pars.train.dfr$at52, 
   sum.spwin = pars.train.dfr$sum.spwin,
   log.xi.full = pars.train.dfr$log.xi.full,
   tps0.fitted.log.xi.full = surface.gev.tps0[["log.xi.full"]]$fitted.values, 
   tps1.fitted.log.xi.full = surface.gev.tps1[["log.xi.full"]]$fitted.values, 
   tps0.fitted.xi.full = -exp(surface.gev.tps0[["log.xi.full"]]$fitted.values), 
   tps1.fitted.xi.full = -exp(surface.gev.tps1[["log.xi.full"]]$fitted.values), 
   log.sig.full = pars.train.dfr$log.sig.full,
   tps0.fitted.log.sig.full = surface.gev.tps0[["log.sig.full"]]$fitted.values, 
   tps1.fitted.log.sig.full = surface.gev.tps1[["log.sig.full"]]$fitted.values,
   tps0.fitted.sig.full = exp(surface.gev.tps0[["log.sig.full"]]$fitted.values), 
   tps1.fitted.sig.full = exp(surface.gev.tps1[["log.sig.full"]]$fitted.values),
   log.mu.full = pars.train.dfr$log.mu.full,
   tps0.fitted.mu.full = NA, 
   tps1.fitted.mu.full = NA)
pred.train.gev$tps0.fitted.mu.full <- 1 -
	   exp(pred.train.gev[, "tps0.fitted.log.sig.full"]) / 
	   exp(pred.train.gev[, "tps0.fitted.log.xi.full"])
pred.train.gev$tps1.fitted.mu.full <- 1 -
	   exp(pred.train.gev[, "tps1.fitted.log.sig.full"]) / 
	   exp(pred.train.gev[, "tps1.fitted.log.xi.full"])

quartz(height = 3.5, width = 10.5)
par(mfrow = c(1,3), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.6,0))
## Shape:
plot(pred.test.gev$log.xi.full, pred.test.gev$tps0.fitted.log.xi.full, pch = 16, cex = 0.5, col = "orangered",
   xlim = range(c(pred.test.gev$log.xi.full,pred.test.gev$tps0.fitted.log.xi.full,
   pred.test.gev$tps1.fitted.log.xi.full), na.rm = TRUE),
   ylim = range(c(pred.test.gev$log.xi.full,pred.test.gev$tps0.fitted.log.xi.full,
   pred.test.gev$tps1.fitted.log.xi.full), na.rm = TRUE))
points(pred.select.gev$log.xi.full, pred.select.gev$tps0.fitted.log.xi.full, pch = 3, cex = 0.5, col = "orangered")
points(pred.test.gev$log.xi.full, pred.test.gev$tps1.fitted.log.xi.full, pch = 3, cex = 0.5, col = "deepskyblue")
points(pred.select.gev$log.xi.full, pred.select.gev$tps1.fitted.log.xi.full, pch = 3, cex = 0.5, col = "deepskyblue")
abline(c(0,1))
## Scale:
plot(pred.test.gev$log.sig.full, pred.test.gev$tps0.fitted.log.sig.full, pch = 16, cex = 0.5, col = "orangered",
   xlim = range(c(pred.test.gev$log.sig.full,pred.test.gev$tps0.fitted.log.sig.full,
   pred.test.gev$tps1.fitted.log.sig.full), na.rm = TRUE),
   ylim = range(c(pred.test.gev$log.sig.full,pred.test.gev$tps0.fitted.log.sig.full,
   pred.test.gev$tps1.fitted.log.sig.full), na.rm = TRUE))
points(pred.select.gev$log.sig.full, pred.select.gev$tps0.fitted.log.sig.full, pch = 3, cex = 0.5, col = "orangered")
points(pred.test.gev$log.sig.full, pred.test.gev$tps1.fitted.log.sig.full, pch = 3, cex = 0.5, col = "deepskyblue")
points(pred.select.gev$log.sig.full, pred.select.gev$tps1.fitted.log.sig.full, pch = 3, cex = 0.5, col = "deepskyblue")
abline(c(0,1))
## Location:
plot(exp(pred.test.gev$log.mu.full), pred.test.gev$tps0.fitted.mu.full, pch = 16, cex = 0.5, col = "orangered",
   xlim = range(c(exp(pred.test.gev$log.mu.full),pred.test.gev$tps0.fitted.mu.full,
   pred.test.gev$tps1.fitted.mu.full), na.rm = TRUE),
   ylim = range(c(exp(pred.test.gev$log.mu.full),pred.test.gev$tps0.fitted.mu.full,
   pred.test.gev$tps1.fitted.mu.full), na.rm = TRUE))
points(exp(pred.select.gev$log.mu.full), pred.select.gev$tps0.fitted.mu.full, pch = 3, cex = 0.5, col = "orangered")
points(exp(pred.test.gev$log.mu.full), pred.test.gev$tps1.fitted.mu.full, pch = 3, cex = 0.5, col = "deepskyblue")
points(exp(pred.select.gev$log.mu.full), pred.select.gev$tps1.fitted.mu.full, pch = 3, cex = 0.5, col = "deepskyblue")
abline(c(0,1))
## The overall aspect is better than with the randomForest fit: there are no systematic mistaken branches
## and wiggles, the line (0,1) is followed, even if it's somewhat thick in its middle.

sqrt(mean((pred.select.gev$log.xi.full - pred.select.gev$tps0.fitted.log.xi.full)^2))
sqrt(mean((pred.select.gev$log.xi.full - pred.select.gev$tps1.fitted.log.xi.full)^2))
sqrt(mean((pred.select.gev$log.sig.full - pred.select.gev$tps0.fitted.log.sig.full)^2))
sqrt(mean((pred.select.gev$log.sig.full - pred.select.gev$tps1.fitted.log.sig.full)^2))
sqrt(mean((pred.select.gev$log.mu.full - log(pred.select.gev$tps0.fitted.mu.full))^2))
sqrt(mean((pred.select.gev$log.mu.full - log(pred.select.gev$tps1.fitted.mu.full))^2))
## In all cases, tps0 wins.

sqrt(mean((pred.train.gev$log.xi.full - pred.train.gev$tps0.fitted.log.xi.full)^2))
sqrt(mean((pred.train.gev$log.xi.full - pred.train.gev$tps1.fitted.log.xi.full)^2))
sqrt(mean((pred.train.gev$log.sig.full - pred.train.gev$tps0.fitted.log.sig.full)^2))
sqrt(mean((pred.train.gev$log.sig.full - pred.train.gev$tps1.fitted.log.sig.full)^2))
sqrt(mean((pred.train.gev$log.mu.full - log(pred.train.gev$tps0.fitted.mu.full))^2))
sqrt(mean((pred.train.gev$log.mu.full - log(pred.train.gev$tps1.fitted.mu.full))^2))
## In all cases, tps0 wins.

## Join the best estimate of all parameters to the pars* dataframes:

dfr.tmp <- merge(pars.train.dfr, pred.train.gev[,c("name","tps0.fitted.xi.full",
   "tps0.fitted.sig.full","tps0.fitted.mu.full")])
dim(pars.train.dfr)
pars.train.dfr <- dfr.tmp
dfr.tmp <- merge(pars.test.dfr, pred.test.gev[,c("name","tps0.fitted.xi.full",
   "tps0.fitted.sig.full","tps0.fitted.mu.full")], all.x = TRUE)
dim(pars.test.dfr)
pars.test.dfr <- dfr.tmp
dfr.tmp <- merge(pars.select.dfr, pred.select.gev[,c("name","tps0.fitted.xi.full",
   "tps0.fitted.sig.full","tps0.fitted.mu.full")], all.x = TRUE)
dim(pars.select.dfr)
pars.select.dfr <- dfr.tmp


## --------------------------------------------------------------------------------------------------------
## Compute the individual FAPs at each site with each method
## --------------------------------------------------------------------------------------------------------

fapcomp.fun <- function(maxima.dfr, pars.dfr, xicol, sigcol, mucol)
   {
	tps.tmp <- data.frame(matrix(ncol = ncol(maxima.dfr),
	    nrow = nrow(maxima.dfr)))
	colnames(tps.tmp) <- colnames(maxima.dfr)
	for(ii in 1:ncol(maxima.dfr))
	#  for(ii in 5)
	   {
	   	nn <- as.character(colnames(maxima.dfr)[ii])
	   	shape <- pars.dfr[pars.dfr$name == nn, xicol]
	   	scale <- pars.dfr[pars.dfr$name == nn, sigcol]
	   	loc <- pars.dfr[pars.dfr$name == nn, mucol]
	    zass <- maxima.dfr[, nn]
	    if(is.na(shape) | is.na(scale) | is.na(loc)) {
	    	tps.tmp[, nn] <- rep(NA, nrow(maxima.dfr))
	       } else {
		    tps.tmp[, nn] <- 1 - sapply(zass, evd::pgev,
		       shape = shape, loc = loc, scale = scale)
		   }
	   }
	tps.tmp
   }

tps0.gevfaps.test.dfr <- fapcomp.fun(
   maxima.dfr = maxima.test.dfr,
   pars.dfr = pred.test.gev,
   xicol = "tps0.fitted.xi.full",
   sigcol = "tps0.fitted.sig.full",
   mucol = "tps0.fitted.mu.full")
save(tps0.gevfaps.test.dfr, file = "results_baluev5/tps0.gevfaps.test.dfr.RObj")

tps1.gevfaps.test.dfr <- fapcomp.fun(
   maxima.dfr = maxima.test.dfr,
   pars.dfr = pred.test.gev,
   xicol = "tps1.fitted.xi.full",
   sigcol = "tps1.fitted.sig.full",
   mucol = "tps1.fitted.mu.full")
save(tps1.gevfaps.test.dfr, file = "results_baluev5/tps1.gevfaps.test.dfr.RObj")

tps0.gevfaps.select.dfr <- fapcomp.fun(
   maxima.dfr = maxima.select.dfr,
   pars.dfr = pred.select.gev,
   xicol = "tps0.fitted.xi.full",
   sigcol = "tps0.fitted.sig.full",
   mucol = "tps0.fitted.mu.full")
save(tps0.gevfaps.select.dfr, file = "results_baluev5/tps0.gevfaps.select.dfr.RObj")

tps1.gevfaps.select.dfr <- fapcomp.fun(
   maxima.dfr = maxima.select.dfr,
   pars.dfr = pred.select.gev,
   xicol = "tps1.fitted.xi.full",
   sigcol = "tps1.fitted.sig.full",
   mucol = "tps1.fitted.mu.full")
save(tps1.gevfaps.select.dfr, file = "results_baluev5/tps1.gevfaps.select.dfr.RObj")

tps0.gevfaps.train.dfr <- fapcomp.fun(
   maxima.dfr = maxima.train.dfr,
   pars.dfr = pred.train.gev,
   xicol = "tps0.fitted.xi.full",
   sigcol = "tps0.fitted.sig.full",
   mucol = "tps0.fitted.mu.full")
save(tps0.gevfaps.train.dfr, file = "results_baluev5/tps0.gevfaps.train.dfr.RObj")

tps1.gevfaps.train.dfr <- fapcomp.fun(
   maxima.dfr = maxima.train.dfr,
   pars.dfr = pred.train.gev,
   xicol = "tps1.fitted.xi.full",
   sigcol = "tps1.fitted.sig.full",
   mucol = "tps1.fitted.mu.full")
save(tps1.gevfaps.train.dfr, file = "results_baluev5/tps1.gevfaps.train.dfr.RObj")

# load("results_baluev4/tps0fitted.gevfaps.test.dfr.RObj")
# load("results_baluev4/tps0computed.gevfaps.test.dfr.RObj")
# load("results_baluev4/tps1fitted.gevfaps.test.dfr.RObj")
# load("results_baluev4/tps1computed.gevfaps.test.dfr.RObj")
# load("results_baluev4/tps0fitted.gevfaps.train.dfr.RObj")
load("results_baluev4/tps0computed.gevfaps.train.dfr.RObj")
# load("results_baluev4/tps1fitted.gevfaps.train.dfr.RObj")
# load("results_baluev4/tps1computed.gevfaps.train.dfr.RObj")

## See the performance on the test set:
pars.test.dfr$fr.tps0fitted.gev95 <- NA
pars.test.dfr$fr.tps0fitted.gev99 <- NA
pars.train.dfr$fr.tps0fitted.gev95 <- NA
pars.train.dfr$fr.tps0fitted.gev99 <- NA
pars.select.dfr$fr.tps0fitted.gev95 <- NA
pars.select.dfr$fr.tps0fitted.gev99 <- NA

## Compute these success rates
for(ii in 1:nrow(pars.select.dfr))
#for(ii in 1:5)
   {
   	pars.select.dfr$fr.tps0fitted.gev95[ii] <-
   	   sum(tps0.gevfaps.select.dfr[, as.character(pars.select.dfr$name[ii])] < 0.05) / 1500
   	pars.select.dfr$fr.tps0fitted.gev99[ii] <-
   	   sum(tps0.gevfaps.select.dfr[, as.character(pars.select.dfr$name[ii])] < 0.01) / 1500
   }
for(ii in 1:nrow(pars.train.dfr))
#for(ii in 1:5)
   {
   	pars.train.dfr$fr.tps0fitted.gev95[ii] <-
   	   sum(tps0.gevfaps.train.dfr[, as.character(pars.train.dfr$name[ii])] < 0.05) / 1500
   	pars.train.dfr$fr.tps0fitted.gev99[ii] <-
   	   sum(tps0.gevfaps.train.dfr[, as.character(pars.train.dfr$name[ii])] < 0.01) / 1500
   }
for(ii in 1:nrow(pars.test.dfr))
#for(ii in 1:5)
   {
   	pars.test.dfr$fr.tps0fitted.gev95[ii] <-
   	   sum(tps0.gevfaps.test.dfr[, as.character(pars.test.dfr$name[ii])] < 0.05) / 1500
   	pars.test.dfr$fr.tps0fitted.gev99[ii] <-
   	   sum(tps0.gevfaps.test.dfr[, as.character(pars.test.dfr$name[ii])] < 0.01) / 1500
   }

## The usual plot:

quartz(height = 7, width = 10)
par(mfrow = c(1,1), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.6,0))

col.x <- "n"
plot(abs(pars.test.dfr[pars.test.dfr$region == "randompos", col.x]), 
   pars.test.dfr$fr.tps0fitted.gev95[pars.test.dfr$region == "randompos"], 
   type = "p", pch = 16, cex = 0.7, ylim = c(0,0.3), 
   xlim = range(abs(pars.test.dfr[,col.x])), main = "Thin-plate spline GEV",
   ylab = "Fraction above threshold", xlab = col.x, col = "grey25")
points(abs(pars.test.dfr[pars.test.dfr$region == "randompos", col.x]),
   pars.test.dfr$fr.tps0fitted.gev99[pars.test.dfr$region == "randompos"], 
   type = "p", pch = 16, col = "grey", cex = 0.7)

points(abs(pars.test.dfr[pars.test.dfr$region == "ecl", col.x]), 
   pars.test.dfr$fr.tps0fitted.gev95[pars.test.dfr$region == "ecl"], 
    pch = 16, cex = 0.7, col = "violetred")
points(abs(pars.test.dfr[pars.test.dfr$region == "ecl", col.x]), 
   pars.test.dfr$fr.tps0fitted.gev99[pars.test.dfr$region == "ecl"], 
    pch = 16, cex = 0.7, col = "pink")

points(abs(pars.test.dfr[pars.test.dfr$region == "lmc", col.x]), 
   pars.test.dfr$fr.tps0fitted.gev95[pars.test.dfr$region == "lmc"], 
    pch = 16, cex = 0.7, col = "orange")
points(abs(pars.test.dfr[pars.test.dfr$region == "lmc", col.x]), 
   pars.test.dfr$fr.tps0fitted.gev99[pars.test.dfr$region == "lmc"], 
    pch = 16, cex = 0.7, col = "gold")

points(abs(pars.test.dfr[pars.test.dfr$region == "ln.ha", col.x]), 
   pars.test.dfr$fr.tps0fitted.gev95[pars.test.dfr$region == "ln.ha"], 
    pch = 16, cex = 0.7, col = "deepskyblue")
points(abs(pars.test.dfr[pars.test.dfr$region == "ln.ha", col.x]), 
   pars.test.dfr$fr.tps0fitted.gev99[pars.test.dfr$region == "ln.ha"], 
    pch = 16, cex = 0.7, col = "skyblue")

abline(h = 0.05, lwd = 2)
abline(h = 0.01, lty = 2, lwd = 2)
## Biased downwards, and most of it at low-n&high-al sites. But on the whole, not bad.  

gev.finaltps.model <- surface.gev.tps0

save(tps0.gevfaps.test.dfr, tps0.gevfaps.train.dfr, tps0.gevfaps.select.dfr,
   tps1.gevfaps.test.dfr, tps1.gevfaps.train.dfr, tps1.gevfaps.select.dfr, file =
   "results_baluev5/gev.tps.fap.all.RObj")
save(surface.gev.tps0, surface.gev.tps1, file = "results_baluev5/fadm.tpsmodels.RObj")
save(pred.train.gev, pred.test.gev, pred.select.gev, file =
   "results_baluev5/gev.tps0.predictions.RObj")
save(gev.finaltps.model, file = "results_baluev5/gev.finaltps.model.RObj")
save(pars.train.dfr, pars.select.dfr, pars.test.dfr, file = 
   "results_baluev5/pars.all.dfr.RObj")

rm(tps0.gevfaps.test.dfr, tps0.gevfaps.train.dfr, tps0.gevfaps.select.dfr,
   tps1.gevfaps.test.dfr, tps1.gevfaps.train.dfr, tps1.gevfaps.select.dfr,
   surface.gev.tps0, surface.gev.tps1,pred.train.gev, pred.test.gev, pred.select.gev)

rm("c.tmp12","c.tmp34","c.tmp56","col.x", dfr.tmp,ii,nn,p.tmp,"pars.ecl.dfr","pars.line.dfr",            
   "pars.newrand.dfr","pars.randompos.dfr")

save.image(WORKSPACE.NAME)


## --------------------------------------------------------------------------------------------------------
## For quantiles
## --------------------------------------------------------------------------------------------------------

## Initial model with all the variables

## A thin-plate spline from package fields:
surface.q95.tps0 <- Tps(pars.train.dfr[,c("log.n", "sum.spwin")], Y = pars.train.dfr$q95, 
   m = 3, lambda = 0.00001)
surface.q95.tps1 <- Tps(pars.train.dfr[,c("log.n", "at48")], Y = pars.train.dfr$q95, 
   m = 3, lambda = 0.00001)
surface.q99.tps0 <- Tps(pars.train.dfr[,c("log.n", "sum.spwin")], Y = pars.train.dfr$q99,
   m = 3, lambda = 0.00001)
surface.q99.tps1 <- Tps(pars.train.dfr[,c("log.n", "at48")], Y = pars.train.dfr$q99,
   m = 3, lambda = 0.00001)
surface.quan.tps0 <- list(q95 = surface.q95.tps0, q99 = surface.q99.tps0)
surface.quan.tps1 <- list(q95 = surface.q95.tps1, q99 = surface.q99.tps1)
rm(surface.q95.tps0,surface.q99.tps0,surface.q95.tps1,surface.q99.tps1)

save(surface.quan.tps0, file = "results_baluev5/surface.quan.tps0.RObj")
save(surface.quan.tps1, file = "results_baluev5/surface.quan.tps1.RObj")

quartz(height = 7, width = 7)
par(mfrow = c(2,2), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.6,0))
for(ii in c("q95","q99"))
   {
	c.tmp12 <- rgb(green = (pars.train.dfr[,ii] - min(pars.train.dfr[,ii])) /
	   (max(pars.train.dfr[,ii]) - min(pars.train.dfr[,ii])), 
	   blue = 1 - (pars.train.dfr[,ii] - min(pars.train.dfr[,ii])) /
	   (max(pars.train.dfr[,ii]) - min(pars.train.dfr[,ii])),
	   red = 0, alpha = 0.4)
	c.tmp34 <- rgb(green = (pars.test.dfr[,ii] - min(pars.test.dfr[,ii], na.rm = TRUE)) /
	   (max(pars.test.dfr[,ii], na.rm = TRUE) - min(pars.test.dfr[,ii], na.rm = TRUE)), 
	   blue = 1 - (pars.test.dfr[,ii] - min(pars.test.dfr[,ii], na.rm = TRUE)) /
	   (max(pars.test.dfr[,ii], na.rm = TRUE) - min(pars.test.dfr[,ii], na.rm = TRUE)),
	   red = 0, alpha = 0.4)
	c.tmp56 <- rgb(green = (pars.select.dfr[,ii] - min(pars.select.dfr[,ii], na.rm = TRUE)) /
	   (max(pars.select.dfr[,ii], na.rm = TRUE) - min(pars.select.dfr[,ii], na.rm = TRUE)), 
	   blue = 1 - (pars.select.dfr[,ii] - min(pars.select.dfr[,ii], na.rm = TRUE)) /
	   (max(pars.select.dfr[,ii], na.rm = TRUE) - min(pars.select.dfr[,ii], na.rm = TRUE)),
	   red = 0, alpha = 0.4)
    surface(surface.quan.tps0[[ii]], type = "c", main = "(log(N), sum.spwin)")
	points(pars.train.dfr$log.n, pars.train.dfr$sum.spwin, pch = 15, cex = 0.7, col = c.tmp12)
	points(pars.test.dfr$log.n, pars.test.dfr$sum.spwin, pch = 15, cex = 0.7, col = c.tmp34)
	points(pars.select.dfr$log.n, pars.select.dfr$sum.spwin, pch = 15, cex = 0.7, col = c.tmp56)
   }
for(ii in c("q95","q99"))
   {
	c.tmp12 <- rgb(green = (pars.train.dfr[,ii] - min(pars.train.dfr[,ii])) /
	   (max(pars.train.dfr[,ii]) - min(pars.train.dfr[,ii])), 
	   blue = 1 - (pars.train.dfr[,ii] - min(pars.train.dfr[,ii])) /
	   (max(pars.train.dfr[,ii]) - min(pars.train.dfr[,ii])),
	   red = 0, alpha = 0.4)
	c.tmp34 <- rgb(green = (pars.test.dfr[,ii] - min(pars.test.dfr[,ii], na.rm = TRUE)) /
	   (max(pars.test.dfr[,ii], na.rm = TRUE) - min(pars.test.dfr[,ii], na.rm = TRUE)), 
	   blue = 1 - (pars.test.dfr[,ii] - min(pars.test.dfr[,ii], na.rm = TRUE)) /
	   (max(pars.test.dfr[,ii], na.rm = TRUE) - min(pars.test.dfr[,ii], na.rm = TRUE)),
	   red = 0, alpha = 0.4)
	c.tmp56 <- rgb(green = (pars.select.dfr[,ii] - min(pars.select.dfr[,ii], na.rm = TRUE)) /
	   (max(pars.select.dfr[,ii], na.rm = TRUE) - min(pars.select.dfr[,ii], na.rm = TRUE)), 
	   blue = 1 - (pars.select.dfr[,ii] - min(pars.select.dfr[,ii], na.rm = TRUE)) /
	   (max(pars.select.dfr[,ii], na.rm = TRUE) - min(pars.select.dfr[,ii], na.rm = TRUE)),
	   red = 0, alpha = 0.4)
    surface(surface.quan.tps1[[ii]], type = "c", main = "(log(N), at48)")
	points(pars.train.dfr$log.n, pars.train.dfr$at48, pch = 15, cex = 0.7, col = c.tmp12)
	points(pars.test.dfr$log.n, pars.test.dfr$at48, pch = 15, cex = 0.7, col = c.tmp34)
	points(pars.select.dfr$log.n, pars.select.dfr$sum.spwin, pch = 15, cex = 0.7, col = c.tmp56)
   }
## The contour lines from the surface fits look almost vertical, as if none of them
## had an effect on the fit.


## Prediction of the parameters

pred.test.quan <- data.frame(name = pars.test.dfr$name, log.n = NA, at48 = NA, sum.spwin = NA,
   q95 = pars.test.dfr$q95, tps0.q95 = NA, tps1.q95 = NA,
   q99 = pars.test.dfr$q99, tps0.q99 = NA, tps1.q99 = NA)
#ii <- 1
for(ii in seq(1, nrow(pred.test.quan), by = 2))   ## Fortunately, nrow() is an even number.
   {
   	## for the shape parameter:
	p.tmp <- predict.surface(surface.quan.tps0$q95, list(
	   log.n = pars.test.dfr$log.n[ii:(ii+1)],
	   sum.spwin = pars.test.dfr$sum.spwin[ii:(ii+1)]), extrap = TRUE)
	pred.test.quan[ii:(ii+1), "log.n"] <- p.tmp$x
	pred.test.quan[ii:(ii+1), "sum.spwin"] <- p.tmp$y
	pred.test.quan[ii:(ii+1), "tps0.q95"] <- diag(p.tmp$z)
	p.tmp <- predict.surface(surface.quan.tps1$q95, list(
	   log.n = pars.test.dfr$log.n[ii:(ii+1)],
	   at48 = pars.test.dfr$at48[ii:(ii+1)]), extrap = TRUE)
	pred.test.quan[ii:(ii+1), "at48"] <- p.tmp$y
	pred.test.quan[ii:(ii+1), "tps1.q95"] <- diag(p.tmp$z)
   	## for the scale parameter:
	p.tmp <- predict.surface(surface.quan.tps0$q99, list(
	   log.n = pars.test.dfr$log.n[ii:(ii+1)],
	   sum.spwin = pars.test.dfr$sum.spwin[ii:(ii+1)]), extrap = TRUE)
	pred.test.quan[ii:(ii+1), "log.n"] <- p.tmp$x
	pred.test.quan[ii:(ii+1), "sum.spwin"] <- p.tmp$y
	pred.test.quan[ii:(ii+1), "tps0.q99"] <- diag(p.tmp$z)
	p.tmp <- predict.surface(surface.quan.tps1$q99, list(
	   log.n = pars.test.dfr$log.n[ii:(ii+1)],
	   at48 = pars.test.dfr$at48[ii:(ii+1)]), extrap = TRUE)
	pred.test.quan[ii:(ii+1), "at48"] <- p.tmp$y
	pred.test.quan[ii:(ii+1), "tps1.q99"] <- diag(p.tmp$z)
   }
pred.test.quan$region <- pars.test.dfr$region
pred.test.quan$beta <- pars.test.dfr$beta

pred.select.quan <- data.frame(name = pars.select.dfr$name, log.n = NA, at48 = NA, sum.spwin = NA,
   q95 = pars.select.dfr$q95, tps0.q95 = NA, tps1.q95 = NA,
   q99 = pars.select.dfr$q99, tps0.q99 = NA, tps1.q99 = NA)
#ii <- 1
for(ii in seq(1, nrow(pred.select.quan), by = 2))   ## Fortunately, nrow() is an even number.
   {
   	## for the shape parameter:
	p.tmp <- predict.surface(surface.quan.tps0$q95, list(
	   log.n = pars.select.dfr$log.n[ii:(ii+1)],
	   sum.spwin = pars.select.dfr$sum.spwin[ii:(ii+1)]), extrap = TRUE)
	pred.select.quan[ii:(ii+1), "log.n"] <- p.tmp$x
	pred.select.quan[ii:(ii+1), "sum.spwin"] <- p.tmp$y
	pred.select.quan[ii:(ii+1), "tps0.q95"] <- diag(p.tmp$z)
	p.tmp <- predict.surface(surface.quan.tps1$q95, list(
	   log.n = pars.select.dfr$log.n[ii:(ii+1)],
	   at48 = pars.select.dfr$at48[ii:(ii+1)]), extrap = TRUE)
	pred.select.quan[ii:(ii+1), "at48"] <- p.tmp$y
	pred.select.quan[ii:(ii+1), "tps1.q95"] <- diag(p.tmp$z)
   	## for the scale parameter:
	p.tmp <- predict.surface(surface.quan.tps0$q99, list(
	   log.n = pars.select.dfr$log.n[ii:(ii+1)],
	   sum.spwin = pars.select.dfr$sum.spwin[ii:(ii+1)]), extrap = TRUE)
	pred.select.quan[ii:(ii+1), "log.n"] <- p.tmp$x
	pred.select.quan[ii:(ii+1), "sum.spwin"] <- p.tmp$y
	pred.select.quan[ii:(ii+1), "tps0.q99"] <- diag(p.tmp$z)
	p.tmp <- predict.surface(surface.quan.tps1$q99, list(
	   log.n = pars.select.dfr$log.n[ii:(ii+1)],
	   at48 = pars.select.dfr$at48[ii:(ii+1)]), extrap = TRUE)
	pred.select.quan[ii:(ii+1), "at48"] <- p.tmp$y
	pred.select.quan[ii:(ii+1), "tps1.q99"] <- diag(p.tmp$z)
   }
pred.select.quan$region <- pars.select.dfr$region
pred.select.quan$beta <- pars.select.dfr$beta

pred.train.quan <- data.frame(name = pars.train.dfr$name, 
   log.n = pars.train.dfr$log.n, at48 = pars.train.dfr$at48, 
   sum.spwin = pars.train.dfr$sum.spwin,
   q95 = pars.train.dfr$q95,
   tps0.q95 = surface.quan.tps0[["q95"]]$fitted.values, 
   tps1.q95 = surface.quan.tps1[["q95"]]$fitted.values, 
   q99 = pars.train.dfr$q99,
   tps0.q99 = surface.quan.tps0[["q99"]]$fitted.values, 
   tps1.q99 = surface.quan.tps1[["q99"]]$fitted.values)
pred.train.quan$region <- pars.train.dfr$region
pred.train.quan$beta <- pars.train.dfr$beta

sqrt(mean((pred.test.quan$q95 - pred.test.quan$tps0.q95)^2))   ## 0.002657046
sqrt(mean((pred.test.quan$q95 - pred.test.quan$tps1.q95)^2))   ## 0.003145592
sqrt(mean((pred.test.quan$q99 - pred.test.quan$tps0.q99)^2))   ## 0.004418725
sqrt(mean((pred.test.quan$q99 - pred.test.quan$tps1.q99)^2))   ## 0.00471334
sqrt(mean((pred.select.quan$q95 - pred.select.quan$tps0.q95)^2))   ## 0.002595238
sqrt(mean((pred.select.quan$q95 - pred.select.quan$tps1.q95)^2))   ## 0.003309012
sqrt(mean((pred.select.quan$q99 - pred.select.quan$tps0.q99)^2))   ## 0.004335668
sqrt(mean((pred.select.quan$q99 - pred.select.quan$tps1.q99)^2))   ## 0.004758469
## Seems that tps0 using sum.spwin finds better the individually estimated quantiles.
## Consistent with the results of GEV.
## Though I should know the df of both fits, I should use those instead of these. 
## Seemingly slight improvement over the spline fits (the mses are similar to the 
## tps1 with at48, so I should check its success rate as a function of beta),
## but possibly this is exactly those few with high aliasing.
save(pred.train.quan, pred.test.quan, pred.select.quan, file = "results_baluev5/pred.all.quan.RObj")
quan.final.model <- surface.quan.tps0

## Spline fits:
# sqrt(mean((pars.test.dfr$q95 - pars.test.dfr$q95.smoo)^2))   ## 0.003328131
# sqrt(mean((pars.test.dfr$q99 - pars.test.dfr$q99.smoo)^2))   ## 0.004740872

## Compute the success rates

for(ii in 1:nrow(pred.train.quan))
   {
   	pred.train.quan$fr.above.thr95tps0[ii] <-
   	   sum(maxima.train.dfr[, as.character(pred.train.quan$name[ii])] >
   	   pred.train.quan$tps0.q95[ii]) / 1500
   	pred.train.quan$fr.above.thr99tps0[ii] <-
   	   sum(maxima.train.dfr[, as.character(pred.train.quan$name[ii])] >
   	   pred.train.quan$tps0.q99[ii]) / 1500
   	pred.train.quan$fr.above.thr95tps1[ii] <-
   	   sum(maxima.train.dfr[, as.character(pred.train.quan$name[ii])] >
   	   pred.train.quan$tps1.q95[ii]) / 1500
   	pred.train.quan$fr.above.thr99tps1[ii] <-
   	   sum(maxima.train.dfr[, as.character(pred.train.quan$name[ii])] >
   	   pred.train.quan$tps1.q99[ii]) / 1500
   }
for(ii in 1:nrow(pred.test.quan))
   {
   	pred.test.quan$fr.above.thr95tps0[ii] <-
   	   sum(maxima.test.dfr[, as.character(pred.test.quan$name[ii])] >
   	   pred.test.quan$tps0.q95[ii]) / 1500
   	pred.test.quan$fr.above.thr99tps0[ii] <-
   	   sum(maxima.test.dfr[, as.character(pred.test.quan$name[ii])] >
   	   pred.test.quan$tps0.q99[ii]) / 1500
   	pred.test.quan$fr.above.thr95tps1[ii] <-
   	   sum(maxima.test.dfr[, as.character(pred.test.quan$name[ii])] >
   	   pred.test.quan$tps1.q95[ii]) / 1500
   	pred.test.quan$fr.above.thr99tps1[ii] <-
   	   sum(maxima.test.dfr[, as.character(pred.test.quan$name[ii])] >
   	   pred.test.quan$tps1.q99[ii]) / 1500
   }
for(ii in 1:nrow(pred.select.quan))
   {
   	pred.select.quan$fr.above.thr95tps0[ii] <-
   	   sum(maxima.select.dfr[, as.character(pred.select.quan$name[ii])] >
   	   pred.select.quan$tps0.q95[ii]) / 1500
   	pred.select.quan$fr.above.thr99tps0[ii] <-
   	   sum(maxima.select.dfr[, as.character(pred.select.quan$name[ii])] >
   	   pred.select.quan$tps0.q99[ii]) / 1500
   	pred.select.quan$fr.above.thr95tps1[ii] <-
   	   sum(maxima.select.dfr[, as.character(pred.select.quan$name[ii])] >
   	   pred.select.quan$tps1.q95[ii]) / 1500
   	pred.select.quan$fr.above.thr99tps1[ii] <-
   	   sum(maxima.select.dfr[, as.character(pred.select.quan$name[ii])] >
   	   pred.select.quan$tps1.q99[ii]) / 1500
   }

## Join the best estimate of all parameters to the pars* dataframes:

dfr.tmp <- merge(pars.train.dfr, pred.train.quan[,c("name","tps0.q95","tps0.q99",
   "fr.above.thr95tps0","fr.above.thr99tps0")])
dim(pars.train.dfr)
dim(dfr.tmp)
pars.train.dfr <- dfr.tmp
dfr.tmp <- merge(pars.test.dfr, pred.test.quan[,c("name","tps0.q95","tps0.q99",
   "fr.above.thr95tps0","fr.above.thr99tps0")], all.x = TRUE)
dim(pars.test.dfr)
dim(dfr.tmp)
dfr.tmp[1:10,]
pars.test.dfr <- dfr.tmp
dfr.tmp <- merge(pars.select.dfr, pred.select.quan[,c("name","tps0.q95","tps0.q99",
   "fr.above.thr95tps0","fr.above.thr99tps0")], all.x = TRUE)
dim(pars.select.dfr)
dim(dfr.tmp)
dfr.tmp[1:10,]
pars.select.dfr <- dfr.tmp


## Plot the false alarm rates against N, like for GEV:

quartz(height = 7, width = 10)
par(mfrow = c(1,1), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.6,0))

col.x <- "beta"
plot(abs(pars.test.dfr[pars.test.dfr$region == "randompos", col.x]), 
   pars.test.dfr$fr.above.thr95tps0[pars.test.dfr$region == "randompos"], 
   type = "p", pch = 16, cex = 0.7, ylim = c(0,0.15), 
   xlim = range(abs(pars.test.dfr[,col.x])), main = "Spline GEV",
   ylab = "Fraction above threshold", xlab = col.x, col = "grey25")
points(abs(pars.test.dfr[pars.test.dfr$region == "randompos", col.x]),
   pars.test.dfr$fr.above.thr99tps0[pars.test.dfr$region == "randompos"], 
   type = "p", pch = 16, col = "grey", cex = 0.7)

points(abs(pars.test.dfr[pars.test.dfr$region == "ecl", col.x]), 
   pars.test.dfr$fr.above.thr95tps0[pars.test.dfr$region == "ecl"], 
    pch = 16, cex = 0.7, col = "violetred")
points(abs(pars.test.dfr[pars.test.dfr$region == "ecl", col.x]), 
   pars.test.dfr$fr.above.thr99tps0[pars.test.dfr$region == "ecl"], 
    pch = 16, cex = 0.7, col = "pink")

points(abs(pars.test.dfr[pars.test.dfr$region == "lmc", col.x]), 
   pars.test.dfr$fr.above.thr95tps0[pars.test.dfr$region == "lmc"], 
    pch = 16, cex = 0.7, col = "orange")
points(abs(pars.test.dfr[pars.test.dfr$region == "lmc", col.x]), 
   pars.test.dfr$fr.above.thr99tps0[pars.test.dfr$region == "lmc"], 
    pch = 16, cex = 0.7, col = "gold")

points(abs(pars.test.dfr[pars.test.dfr$region == "ln.ha", col.x]), 
   pars.test.dfr$fr.above.thr95tps0[pars.test.dfr$region == "ln.ha"], 
    pch = 16, cex = 0.7, col = "deepskyblue")
points(abs(pars.test.dfr[pars.test.dfr$region == "ln.ha", col.x]), 
   pars.test.dfr$fr.above.thr99tps0[pars.test.dfr$region == "ln.ha"], 
    pch = 16, cex = 0.7, col = "skyblue")

abline(h = 0.05, lwd = 2)
abline(h = 0.01, lty = 2, lwd = 2)

## All right. The sum.spwin variable is capable to adjust for the high alias
## effect in the polar regions, when using at48, that is not. Nevertheless, 
## on the line plot, it seems that there are still 2nd-order effects that are
## unexplained.


range(pars.train.dfr$n)
range(pars.test.dfr$n)
range(pars.select.dfr$n)
range(pars.train.dfr$sum.spwin)
range(pars.test.dfr$sum.spwin)
range(pars.select.dfr$sum.spwin)

gev.attr.final <- c("log.n", "sum.spwin")
quan.attr.final <- c("log.n", "sum.spwin")
load(file = "results_baluev5/fadm.rfmodel.RObj")  

save(pars.all.dfr, mfull.rf.final, mfull.attr.final, gev.finaltps.model, gev.attr.final,
   quan.final.model, quan.attr.final, file = "results_baluev5/for_further_work.RObj")

rm(mfull.rf.final, mfull.attr.final)





## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------
## For other uses, the ra and dec of my random positions and my line (EBs for Andrej)
## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------

pt.tmp <- pars.all.dfr[pars.all.dfr$region == "line", c("ra", "dec")]
write.table(pt.tmp, file = "/Users/mariasuveges/Documents/SemiSupervised/ECLsubclass/CompAndrej/line.txt",
    quote = F, sep = " ", row.names = F, col.names = F)
pt.tmp <- pars.all.dfr[pars.all.dfr$region == "randompos", c("ra", "dec")]
write.table(pt.tmp, file = "/Users/mariasuveges/Documents/SemiSupervised/ECLsubclass/CompAndrej/randompositions.txt",
    quote = F, sep = " ", row.names = F, col.names = F)

sum(pars.all.dfr$region == "randompos")














