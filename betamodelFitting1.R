


## GOAL OF THE WORKSPACE:
##
##               MODEL FITS FOR a, b (AND LATER FOR M, shape, scale, location)
##                    USING THE HALF OF THE RANDOM POSITIONS AND THE LINE
##
##                    COMPARE THE MODELS ON THE OTHER HALF OF SIMULATIONS 
##                            AND ON THE DENSELY SAMPLED REGION


## CONCLUSIONS FROM firstlookLine1.RData:
 
## 1. THE BEST (SO FAR) REGRESSION IS OBTAINED BY REGRESSING THE COEFFICIENT
##    b ON THE SUM OF THE 6 HIGHEST SPECTRAL WINDOW PEAKS.
## 2. THE DEPENDENCE IS NOT COMPLETELY SUMMARIZED BY THE REGRESSION ON
##    THIS SUM: IT IS NECESSARY TO INCLUDE THE SKY REGION INTO THE REGRESSION.
## 3. NO DEPENDENCE ON NUMBER OF OBSERVATIONS WHATSOEVER.

## CONCLUSIONS FROM firstlookLine2.RData:
 
## 1. IN THE SPACE SPANNED BY THE HIGHEST SPECTRAL WINDOW PEAKS,
##    THE STRUCTURE IS RATHER LIKE A CURLED-UP TUBE, ITS OUTER END LONG AND
##    THIN (THE POLES), THE INNERMOST END IS A THICK BUNCH (THE ECLIPTIC)
## 2. NO SEPARATE GROUPS, RATHER A CONTINUOUS TRANSITION.
## 3. THESE DATA HAVE NO NOISE, BUT USING THE LLE'S "STRUCTURE-CUTTING"
##    THRESHOLD ALPHA = 0.99 CORRESPONDING TO NO-NOISE CASE RESULTS IN
##    PRACTICALLY NO DIMENSION REDUCTION. STRONG CUTS CAN GIVE NICE PATTERNS
##    WITH LONG BRANCHES COLLECTING THE MOST CHARACTERISTIC SPECTRAL WINDOWS
##    IN LOW-DIMENSIONAL SPACE, BUT NO GROUPS HERE, EITHER. MOST POINTS ARE
##    IN THE MIDDLE, BUNCHED UP.
## ==> SEPARATION INTO BANDS: RATHER ARBITRARY BOUNDARIES, NOT PHYSICAL GROUPS.

## CONCLUSIONS FROM modelfitting1, *1-2 and *2.RData:

## 1. NO POINT TO TRY TO USE THE EM ALGORITHM:
##    BETTER, BUT THE GROUPING IS NOT PREDICTABLE
## 2. USE OF CUTPOINTS OR CLASSIFICATION LEADS TO NEAR-EQUAL QUALITY MODELS
 
##    
## Use half of Leanne's simulations ** AT RANDOM POSITIONS ** to derive the models
##
## Method:
##  - Using the predicted group of the 7-class classification as 'band appartenance',
##    and regarding it as fixed and known, fit glm with different groupings 
##    (no EM algorithm, but use beta distribution), compare BIC/likelihood
##  - using a simple cut point
## Steps:
## 1. Pre-processing (dataframe with all the variables [ra, dec, lambda, beta, n, 
##    variance of time points, spectral window peak heights, fitted single-location
##    beta distribution parameters, etc.])
## 2. Try the groupings (3) - (1,2,4,5,6,7), (2,3) - (1,4,5,6,7) and (1,2,3) -
##    (4,5,6,7).
## 3. Try different cut points.


## DATAFRAME HOLDING MOST OF THE RESULTS:  pars.randompos.dfr
## 
## Columns (a.single, b.single): single-location estimates of the best-fit beta distr.


setwd("/Users/mariasuveges/Documents/SpectralAnalysis/FAP/feasibilityForGaia/OtherFAPs/")


library(mapproj)
library(gridBase)
library(corrplot)
library(rgl)
library(mclust)
library(numDeriv)

WORKSPACE.NAME <- "betamodelFitting1.RData"
save.image(WORKSPACE.NAME)


## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------
## Loading previous results and functions
## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------


## Citation from firstlookLine2.R (full-sky dataframe with classes and
## the GM models):
# dfr <- max.power.fullsky.dfr
# dfr$class7000_12group <- predict(gmxt7000[[12]],
    # newdata = max.power.fullsky.dfr[, c(3:6,8:19)])$classification
# dfr$class7000_13group <- predict(gmxt7000[[13]],
    # newdata = max.power.fullsky.dfr[, c(3:6,8:19)])$classification
# dfr$class7000_14group <- predict(gmxt7000[[14]],
    # newdata = max.power.fullsky.dfr[, c(3:6,8:19)])$classification
# dfr$class7000_15group <- predict(gmxt7000[[15]],
    # newdata = max.power.fullsky.dfr[, c(3:6,8:19)])$classification
# dfr$class7000_16group <- predict(gmxt7000[[16]],
    # newdata = max.power.fullsky.dfr[, c(3:6,8:19)])$classification
# dfr$class7000_17group <- predict(gmxt7000[[17]],
    # newdata = max.power.fullsky.dfr[, c(3:6,8:19)])$classification
# dfr$class7000_18group <- predict(gmxt7000[[18]],
    # newdata = max.power.fullsky.dfr[, c(3:6,8:19)])$classification
# dfr$class7000_19group <- predict(gmxt7000[[19]],
    # newdata = max.power.fullsky.dfr[, c(3:6,8:19)])$classification
# dfr$class7000_20group <- predict(gmxt7000[[20]],
    # newdata = max.power.fullsky.dfr[, c(3:6,8:19)])$classification
# max.power.fullsky7000 <- dfr
# save(max.power.fullsky7000, file = "max.power.fullsky7000.RObj")
# save(gmxt7000, file = "gmxt7000.RObj")

## The classifiers and the table with the fitted classes
## max.power.fullsky7000:
load(file = "max.power.fullsky7000.RObj")
## gmxt7000:
load(file = "gmxt7000.RObj")


## Necessary functions:
## Beta likelihood, EM algorithm for the beta-linear mixture:
source("EMbetaFunctions.R")
source("cooTransform.R")

## Leanne's results for the random positions:
traindatadir.name <- "/Users/mariasuveges/Documents/SpectralAnalysis/FAP/feasibilityForGaia/NoiseSimulations011113/NoiseSimulationsRandomPositions/"
testdatadir.name <- "/Users/mariasuveges/Documents/SpectralAnalysis/FAP/feasibilityForGaia/NoiseSimulations011113/NoiseSimulationsRandomMagnitude_RandomPositions_SecondSet/"
## Leanne's results for the line:
datadir.name <- "/Users/mariasuveges/Documents/SpectralAnalysis/FAP/feasibilityForGaia/NoiseSimulations011113/NoiseSimulationsLineGridAll/"
ecldir.name <- "/Users/mariasuveges/Documents/SpectralAnalysis/FAP/feasibilityForGaia/NoiseSimulations011113/NoiseSimulationsEclGridAll/"

## Time samplings and spectral windows for the whole sky (for backgrounds when necessary):
roughgriddir.name <- "/Users/mariasuveges/Documents/SpectralAnalysis/AliasGaia/GaiaSpectralWindows/NominalScanningLaw230813/"
## Time samplings and spectral windows for the random positions:
timesamplingdir.randompos.name <- "/Users/mariasuveges/Documents/SpectralAnalysis/AliasGaia/GaiaSpectralWindows/timesamplingsForFAP/"
## Time samplings and spectral windows for the line:
finelinedir.name <- "/Users/mariasuveges/Documents/SpectralAnalysis/AliasGaia/GaiaSpectralWindows/NominalScanningLaw_LineGrid291013/"
## Time samplings and spectral windows for the region around the ecliptic:
fineecldir.name <- "/Users/mariasuveges/Documents/SpectralAnalysis/AliasGaia/GaiaSpectralWindows/NominalScanningLaw_EclGrid291013/"

## For spectral windows of the random time samplings: taken from workspace 
## /Users/mariasuveges/Documents/SpectralAnalysis/AliasGaia/GaiaSpectralWindows/timesamplingsForFAP/spWin_RandomPositions.RData
## the saves were:
## save(max.ind.partpg.randompos, max.freq.partpg.randompos, max.power.partpg.randompos, 
##    file = "max_aliases_partpg_randompos.RObj")
## save(max.ind.fullpg.randompos, max.freq.fullpg.randompos, max.power.fullpg.randompos, 
##    file = "max_aliases_fullpg_randompos.RObj")

load("/Users/mariasuveges/Documents/SpectralAnalysis/AliasGaia/GaiaSpectralWindows/timesamplingsForFAP/max_aliases_fullpg_randompos.RObj")



## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------
## Dataframe for all the random locations, holding all estimated or given parameters:
## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------

randomposdirlist <- list.files(traindatadir.name)
pars.randompos.dfr <- data.frame(name = randomposdirlist)
pars.randompos.dfr$ra <- as.numeric(sapply(randomposdirlist, function(chr)
   {
#  	chr <- randomposdirlist[1] 
  	unlist(strsplit(chr, split = "_"))[1]
   }))
pars.randompos.dfr$dec <- as.numeric(sapply(randomposdirlist, function(chr)
   {
#  	chr <- randomposdirlist[1] 
  	unlist(strsplit(chr, split = "_"))[2]
   }))
pars.randompos.dfr$lambda <- NA
pars.randompos.dfr$beta <- NA
pars.randompos.dfr[,4:5] <- t(apply(pars.randompos.dfr[,2:3], 1, function(vec)
   {
   	eq2ecl.fun(vec[1], vec[2])
   }))
pars.randompos.dfr$n <- NA
pars.randompos.dfr[1:10,]


## --------------------------------------------------------------------------------------------------------
## Download the FAPs and fill up the pars.randompos.dfr$n
## --------------------------------------------------------------------------------------------------------

## the number of observations are in the commented-out header in the
## FAP files:

fap.randompos.dfr <- as.data.frame(matrix(ncol = length(randomposdirlist), nrow = 750))
colnames(fap.randompos.dfr) <- randomposdirlist
for(ii in 1:length(randomposdirlist))
   {
#    ii <- 2
	dd <- randomposdirlist[ii]
	pp <- list.files(paste(traindatadir.name, dd, sep = ""), full.names = TRUE,
	   pattern = ".dat")
    cc <- unlist(strsplit(readLines(pp, n = 2), " "))
    pars.randompos.dfr[pars.randompos.dfr$name == dd, "n"] <- as.numeric(cc[length(cc)])
	fap.randompos.dfr[,dd] <- read.table(pp)[,6]
   }
# fap.randompos.dfr[1:20, 1:20]


## --------------------------------------------------------------------------------------------------------
## Get the spectral window values and frequencies at its maximum power 
## --------------------------------------------------------------------------------------------------------

## max.ind.fullpg.randompos, max.power.fullpg.randompos, max.freq.fullpg.randompos:

# dim(max.power.line)
dfr <- as.data.frame(t(max.power.fullpg.randompos))
# dfr$ra <- as.numeric(sapply(row.names(dfr), function(chr)
   # {
  	# unlist(strsplit(chr, split = "_"))[1]
   # }))
# dfr$dec <- as.numeric(sapply(row.names(dfr), function(chr)
   # {
  	# unlist(strsplit(chr, split = "_"))[2]
   # }))
# dfr[1:10,]
colnames(dfr)[3:19] <- paste("at", (1:17)*4, sep = "") 
# is.numeric(dfr$ra)
# is.numeric(dfr$dec)
# intersect(dfr$dec, pars.randompos.dfr$dec)
dfr1 <- merge(pars.randompos.dfr, dfr)
# dfr1[1:10,]
dim(dfr1)
## Don't forget to leave out column spwin.at.20 from the matrix given to classification.
# colnames(max.power.fullsky7000)
## Ok, the names are the same.

dfr <- as.data.frame(t(max.freq.fullpg.randompos))
# dfr[1:5,]
colnames(dfr)[3:19] <- paste("freq.at", (1:17)*4, sep = "") 
# is.numeric(dfr$ra)
# is.numeric(dfr$dec)
# identical(dfr$dec, pars.randompos.dfr$dec)
# identical(dfr$ra, pars.randompos.dfr$ra)
dfr2 <- merge(dfr1, dfr)
# dfr2[1:10,]

pars.randompos.dfr <- dfr2
rm(dfr,dfr1,dfr2,max.power.fullpg.randompos,max.freq.fullpg.randompos,max.ind.fullpg.randompos)
# pars.randompos.dfr[801:810,]


## --------------------------------------------------------------------------------------------------------
## Histograms of FAPs
## --------------------------------------------------------------------------------------------------------

# # fapdensity.randompos <- apply(fap.randompos.dfr, 2, function(vec)
   # hist(vec, breaks = seq(0, 1, by = 0.05), plot = F)$counts/1500)
# dim(fapdensity.randompos)
 

## --------------------------------------------------------------------------------------------------------
## Attach variance of the times
## --------------------------------------------------------------------------------------------------------


pars.randompos.dfr$vart <- NA
d.tmp <- list.files(paste0(timesamplingdir.randompos.name, "/RandomPositions"))
for(nn in pars.randompos.dfr$name) 
  {
   n.tmp <- substr(nn, start = 1, stop = nchar(nn) - 5)
   ii <- grep(d.tmp, pattern = n.tmp, value = TRUE)
   t.tmp <- read.table(paste0(timesamplingdir.randompos.name, "/RandomPositions/", ii))[,1]
   pars.randompos.dfr$vart[as.character(pars.randompos.dfr$name) == nn] <- var(t.tmp)
  }
  

## --------------------------------------------------------------------------------------------------------
## Attach the sums of various alias powers to the pars file 
## --------------------------------------------------------------------------------------------------------


## Put the sums of the highest k peaks into the parameters dataframe:

	dfr <- pars.randompos.dfr
	dfr$sum.spwin <- apply(dfr[,7:23], MAR = 1, FUN = sum)
	dfr$mean.spwinpeaks <- apply(dfr[,c(6:23)], MAR = 1, FUN = function(vec)
       sum(vec[2:18]/vec[1]))
	dfr$total.spwinpeaks <- apply(dfr[,c(6:23)], MAR = 1, FUN = function(vec)
       sum(vec[2:18]*vec[1]))
	dfr$highest.spwin <- apply(dfr[,7:23], MAR = 1, FUN = max)
	dfr$highest2sum.spwin <- apply(dfr[,7:23], MAR = 1, FUN = function(vec)
       {
       	v.tmp <- vec
       	v.tmp <- v.tmp[order(v.tmp, decreasing = TRUE)]
       	sum(v.tmp[1:2])
       	})
	dfr$highest3sum.spwin <- apply(dfr[,7:23], MAR = 1, FUN = function(vec)
       {
       	v.tmp <- vec
       	v.tmp <- v.tmp[order(v.tmp, decreasing = TRUE)]
       	sum(v.tmp[1:3])
       	})
	dfr$highest4sum.spwin <- apply(dfr[,7:23], MAR = 1, FUN = function(vec)
       {
       	v.tmp <- vec
       	v.tmp <- v.tmp[order(v.tmp, decreasing = TRUE)]
       	sum(v.tmp[1:4])
       	})
	dfr$highest5sum.spwin <- apply(dfr[,7:23], MAR = 1, FUN = function(vec)
       {
       	v.tmp <- vec
       	v.tmp <- v.tmp[order(v.tmp, decreasing = TRUE)]
       	sum(v.tmp[1:5])
       	})
	dfr$highest6sum.spwin <- apply(dfr[,7:23], MAR = 1, FUN = function(vec)
       {
       	v.tmp <- vec
       	v.tmp <- v.tmp[order(v.tmp, decreasing = TRUE)]
       	sum(v.tmp[1:6])
       	})
	dfr$highest7sum.spwin <- apply(dfr[,7:23], MAR = 1, FUN = function(vec)
       {
       	v.tmp <- vec
       	v.tmp <- v.tmp[order(v.tmp, decreasing = TRUE)]
       	sum(v.tmp[1:7])
       	})
	dfr$highest8sum.spwin <- apply(dfr[,7:23], MAR = 1, FUN = function(vec)
       {
       	v.tmp <- vec
       	v.tmp <- v.tmp[order(v.tmp, decreasing = TRUE)]
       	sum(v.tmp[1:8])
       	})
	dfr$highest9sum.spwin <- apply(dfr[,7:23], MAR = 1, FUN = function(vec)
       {
       	v.tmp <- vec
       	v.tmp <- v.tmp[order(v.tmp, decreasing = TRUE)]
       	sum(v.tmp[1:9])
       	})
    dfr[1:10,]
    pars.randompos.dfr <- dfr
#	dfr <- pars.randompos.dfr[order(pars.randompos.dfr[,ii]),]

pars.randompos.dfr$sum.of.five.spwin <- 
   pars.randompos.dfr$at12 + pars.randompos.dfr$at16 + pars.randompos.dfr$at28 + 
   pars.randompos.dfr$at40 + pars.randompos.dfr$at68
dim(pars.randompos.dfr)
pars.randompos.dfr$sum.of.six.spwin <- pars.randompos.dfr$at4 + 
   pars.randompos.dfr$at12 + pars.randompos.dfr$at16 + pars.randompos.dfr$at28 + 
   pars.randompos.dfr$at40 + pars.randompos.dfr$at68
pars.randompos.dfr$sum.of.seven.spwin <- pars.randompos.dfr$at4 + 
   pars.randompos.dfr$at12 + pars.randompos.dfr$at16 + pars.randompos.dfr$at28 + 
   pars.randompos.dfr$at40 + pars.randompos.dfr$at52 + pars.randompos.dfr$at68
pars.randompos.dfr$sum.of.eight.spwin <- pars.randompos.dfr$at4 + 
   pars.randompos.dfr$at12 + pars.randompos.dfr$at16 + pars.randompos.dfr$at28 + 
   pars.randompos.dfr$at40 + pars.randompos.dfr$at52 + pars.randompos.dfr$at68 + 
   pars.randompos.dfr$at56
pars.randompos.dfr$sum.of.nine.spwin <- pars.randompos.dfr$at4 + 
   pars.randompos.dfr$at12 + pars.randompos.dfr$at16 + pars.randompos.dfr$at28 + 
   pars.randompos.dfr$at40 + pars.randompos.dfr$at52 + pars.randompos.dfr$at68 + 
   pars.randompos.dfr$at56 + pars.randompos.dfr$at24
pars.randompos.dfr$sum.of.ten.spwin <- pars.randompos.dfr$at4 + 
   pars.randompos.dfr$at12 + pars.randompos.dfr$at16 + pars.randompos.dfr$at28 + 
   pars.randompos.dfr$at40 + pars.randompos.dfr$at52 + pars.randompos.dfr$at68 + 
   pars.randompos.dfr$at56 + pars.randompos.dfr$at24 + pars.randompos.dfr$at44


## Just a check:
## Which frequencies can occur among the highest spectral peaks?

which.highest <- apply(pars.randompos.dfr[, 7:23], MAR = 1, FUN = function(vec)
   {
    o.tmp <- order(vec, decreasing = TRUE)[1:6]
    o.tmp*4
   })

table(which.highest)
  # 4  12  16  24  28  40  44  52  56  68 
# 257 720 652 161 714 698   3 258 249 608 
## So maybe if I want the 6 highest peaks, it will be enough to compute the
## spectral window in these intervals....
## 12, 28, 40, 16, 68

pairs(pars.randompos.dfr[, c(7,9,10,13,16,19,23)], pch = 16, cex = 0.3)
## This will be definitely different from the line: I have almost no representants from
## the highly-aliased pole regions. Line is still better?

## Get the line also


## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------
## Dataframe for all the line locations, holding all estimated or given parameters:
## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------

linedirlist <- list.files(datadir.name)
pars.line.dfr <- data.frame(name = linedirlist)
pars.line.dfr$ra <- as.numeric(sapply(linedirlist, function(chr)
   {
#  	chr <- linedirlist[1] 
  	unlist(strsplit(chr, split = "_"))[1]
   }))
pars.line.dfr$dec <- as.numeric(sapply(linedirlist, function(chr)
   {
#  	chr <- linedirlist[1] 
  	unlist(strsplit(chr, split = "_"))[2]
   }))
pars.line.dfr$lambda <- NA
pars.line.dfr$beta <- NA
pars.line.dfr[,4:5] <- t(apply(pars.line.dfr[,2:3], 1, function(vec)
   {
   	eq2ecl.fun(vec[1], vec[2])
   }))
pars.line.dfr$n <- NA


## --------------------------------------------------------------------------------------------------------
## Download the FAPs and fill up the pars.line.dfr$n
## --------------------------------------------------------------------------------------------------------

## the number of observations are in the commented-oout header in the
## FAP files:

fap.line.dfr <- as.data.frame(matrix(ncol = length(linedirlist), nrow = 1500))
colnames(fap.line.dfr) <- linedirlist
for(ii in 1:length(linedirlist))
   {
#    ii <- 2
	dd <- linedirlist[ii]
	pp <- list.files(paste(datadir.name, dd, sep = ""), full.names = TRUE,
	   pattern = ".dat")
    cc <- unlist(strsplit(readLines(pp, n = 2), " "))
    pars.line.dfr[pars.line.dfr$name == dd, "n"] <- as.numeric(cc[length(cc)])
	fap.line.dfr[,dd] <- read.table(pp)[,6]
   }
# fap.line.dfr[1:20, 1:20]


## --------------------------------------------------------------------------------------------------------
## Get the spectral window values and frequencies of its maximum power on the line
## --------------------------------------------------------------------------------------------------------

## max.ind.line, max.power.line, max.freq.line:

load(paste(finelinedir.name, "max_aliases.RObj", sep = ""))
# dim(max.power.line)
dfr <- as.data.frame(t(max.power.line))
dfr$ra <- as.numeric(sapply(row.names(dfr), function(chr)
   {
#  	chr <- linedirlist[1] 
  	unlist(strsplit(chr, split = "_"))[2]
   }))
dfr$dec <- as.numeric(sapply(row.names(dfr), function(chr)
   {
#  	chr <- linedirlist[1] 
  	unlist(strsplit(chr, split = "_"))[3]
   }))
# dfr[1:10,]
colnames(dfr)[3:19] <- paste("at", (1:17)*4, sep = "") 
# is.numeric(dfr$ra)
# is.numeric(dfr$dec)
# intersect(dfr$dec, pars.line.dfr$dec)
dfr1 <- merge(pars.line.dfr, dfr)
# dfr1[1:10,]
## Don't forget to leave out column spwin.at.20 from the matrix given to classification.
# colnames(max.power.fullsky7000)
## Ok, the names are the same.

dfr <- as.data.frame(t(max.freq.line))
dfr$ra <- as.numeric(sapply(row.names(dfr), function(chr)
   {
#  	chr <- linedirlist[1] 
  	unlist(strsplit(chr, split = "_"))[2]
   }))
dfr$dec <- as.numeric(sapply(row.names(dfr), function(chr)
   {
#  	chr <- linedirlist[1] 
  	unlist(strsplit(chr, split = "_"))[3]
   }))
# dfr[1:5,]
colnames(dfr)[3:19] <- paste("freq.at", (1:17)*4, sep = "") 
# is.numeric(dfr$ra)
# is.numeric(dfr$dec)
# identical(dfr$dec, pars.line.dfr$dec)
# identical(dfr$ra, pars.line.dfr$ra)
dfr2 <- merge(dfr1, dfr)
# dfr2[1:10,]

pars.line.dfr <- dfr2
rm(dfr,dfr1,dfr2,max.power.line,max.freq.line,max.ind.line)
# pars.line.dfr[801:810,]


## --------------------------------------------------------------------------------------------------------
## Histograms of FAPs
## --------------------------------------------------------------------------------------------------------

# # fapdensity.line <- apply(fap.line.dfr, 2, function(vec)
   # hist(vec, breaks = seq(0, 1, by = 0.05), plot = F)$counts/1500)
# dim(fapdensity.line)
 

## --------------------------------------------------------------------------------------------------------
## Attach variance of the times
## --------------------------------------------------------------------------------------------------------

pars.line.dfr$vart <- NA
d.tmp <- list.files(finelinedir.name)[1:900]
for(ii in d.tmp) 
  {
   f.tmp <- list.files(paste(finelinedir.name, ii, sep = ""))
   t.tmp <- read.table(paste(finelinedir.name, ii, "/", f.tmp, sep = ""))[,1]
   pars.line.dfr$vart[as.character(pars.line.dfr$ra) == ii] <- var(t.tmp)
  }
  

## --------------------------------------------------------------------------------------------------------
## Attach the sums of various alias powers to the pars file 
## --------------------------------------------------------------------------------------------------------


## Put the sums of the highest k peaks into the parameters dataframe:

	dfr <- pars.line.dfr
	dfr$sum.spwin <- apply(dfr[,7:23], MAR = 1, FUN = sum)
	dfr$mean.spwin <- apply(dfr[,c(6:23)], MAR = 1, FUN = function(vec)
       sum(vec[2:18]/vec[1]))
	dfr$total.spwin <- apply(dfr[,c(6:23)], MAR = 1, FUN = function(vec)
       sum(vec[2:18]*vec[1]))
	dfr$highest.spwin <- apply(dfr[,7:23], MAR = 1, FUN = max)
	dfr$highest2sum.spwin <- apply(dfr[,7:23], MAR = 1, FUN = function(vec)
       {
       	v.tmp <- vec
       	v.tmp <- v.tmp[order(v.tmp, decreasing = TRUE)]
       	sum(v.tmp[1:2])
       	})
	dfr$highest3sum.spwin <- apply(dfr[,7:23], MAR = 1, FUN = function(vec)
       {
       	v.tmp <- vec
       	v.tmp <- v.tmp[order(v.tmp, decreasing = TRUE)]
       	sum(v.tmp[1:3])
       	})
	dfr$highest4sum.spwin <- apply(dfr[,7:23], MAR = 1, FUN = function(vec)
       {
       	v.tmp <- vec
       	v.tmp <- v.tmp[order(v.tmp, decreasing = TRUE)]
       	sum(v.tmp[1:4])
       	})
	dfr$highest5sum.spwin <- apply(dfr[,7:23], MAR = 1, FUN = function(vec)
       {
       	v.tmp <- vec
       	v.tmp <- v.tmp[order(v.tmp, decreasing = TRUE)]
       	sum(v.tmp[1:5])
       	})
	dfr$highest6sum.spwin <- apply(dfr[,7:23], MAR = 1, FUN = function(vec)
       {
       	v.tmp <- vec
       	v.tmp <- v.tmp[order(v.tmp, decreasing = TRUE)]
       	sum(v.tmp[1:6])
       	})
	dfr$highest7sum.spwin <- apply(dfr[,7:23], MAR = 1, FUN = function(vec)
       {
       	v.tmp <- vec
       	v.tmp <- v.tmp[order(v.tmp, decreasing = TRUE)]
       	sum(v.tmp[1:7])
       	})
	dfr$highest8sum.spwin <- apply(dfr[,7:23], MAR = 1, FUN = function(vec)
       {
       	v.tmp <- vec
       	v.tmp <- v.tmp[order(v.tmp, decreasing = TRUE)]
       	sum(v.tmp[1:8])
       	})
	dfr$highest9sum.spwin <- apply(dfr[,7:23], MAR = 1, FUN = function(vec)
       {
       	v.tmp <- vec
       	v.tmp <- v.tmp[order(v.tmp, decreasing = TRUE)]
       	sum(v.tmp[1:9])
       	})
    dfr[1:10,]
    pars.line.dfr <- dfr
#	dfr <- pars.line.dfr[order(pars.line.dfr[,ii]),]

pars.line.dfr$sum.of.five.spwin <- 
   pars.line.dfr$at12 + pars.line.dfr$at16 + pars.line.dfr$at28 + 
   pars.line.dfr$at40 + pars.line.dfr$at68
pars.line.dfr$sum.of.six.spwin <- pars.line.dfr$at4 + 
   pars.line.dfr$at12 + pars.line.dfr$at16 + pars.line.dfr$at28 + 
   pars.line.dfr$at40 + pars.line.dfr$at68
pars.line.dfr$sum.of.seven.spwin <- pars.line.dfr$at4 + 
   pars.line.dfr$at12 + pars.line.dfr$at16 + pars.line.dfr$at28 + 
   pars.line.dfr$at40 + pars.line.dfr$at52 + pars.line.dfr$at68
pars.line.dfr$sum.of.eight.spwin <- pars.line.dfr$at4 + 
   pars.line.dfr$at12 + pars.line.dfr$at16 + pars.line.dfr$at28 + 
   pars.line.dfr$at40 + pars.line.dfr$at52 + pars.line.dfr$at68 + 
   pars.line.dfr$at56
pars.line.dfr$sum.of.nine.spwin <- pars.line.dfr$at4 + 
   pars.line.dfr$at12 + pars.line.dfr$at16 + pars.line.dfr$at28 + 
   pars.line.dfr$at40 + pars.line.dfr$at52 + pars.line.dfr$at68 + 
   pars.line.dfr$at56 + pars.line.dfr$at24
pars.line.dfr$sum.of.ten.spwin <- pars.line.dfr$at4 + 
   pars.line.dfr$at12 + pars.line.dfr$at16 + pars.line.dfr$at28 + 
   pars.line.dfr$at40 + pars.line.dfr$at52 + pars.line.dfr$at68 + 
   pars.line.dfr$at56 + pars.line.dfr$at24 + pars.line.dfr$at44
dim(pars.line.dfr)


## Just a check:
## Which frequencies can occur among the highest spectral peaks?

which.highest.line <- apply(pars.line.dfr[, 7:23], MAR = 1, FUN = function(vec)
   {
    o.tmp <- order(vec, decreasing = TRUE)[1:6]
    o.tmp*4
   })

table(which.highest.line)
  # 4  12  16  24  28  40  44  52  56  64  68 
# 419 900 803 168 894 882   1 295 339   1 698 
## So maybe if I want the 6 highest peaks, it will be enough to compute the
## spectral window in these intervals....



## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------
## Beta ditribution fits to the first 750 FAPs position-wise
## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------

## --------------------------------------------------------------------------------------------------------
## Random positions
## --------------------------------------------------------------------------------------------------------

dim(fap.randompos.dfr)
## 750 720 

## Use the function fitbeta.fun in EMbetaFunctions.R
betafits.randompos <- apply(fap.randompos.dfr, M = 2, function(vec) try(fitbeta.fun(vec)))
betapar.poswise.randompos <- matrix(ncol = ncol(fap.randompos.dfr), nrow = 4)
dimnames(betapar.poswise.randompos) <- list(c("a.single","a.single.se","b.single","b.single.se"),
   colnames(fap.randompos.dfr))
betapar.poswise.randompos[c("a.single","b.single"),] <- sapply(betafits.randompos, function(lst) 
   if(inherits(lst, "try-error")) rep(NA,2) else lst$par)
betapar.poswise.randompos[c("a.single.se","b.single.se"),] <- sapply(betafits.randompos, function(lst) 
   if(inherits(lst, "try-error")) rep(NA,2) else sqrt(diag(solve(-lst$hessian))))
   
## Attach the estimated beta-distribution coeffs to the pars.randompos.dfr dataframe:

identical(colnames(betapar.poswise.randompos), as.character(pars.randompos.dfr$name))
## True, so it can go by simple merging..
dfr1 <- data.frame(t(betapar.poswise.randompos))
dfr1$name <- colnames(betapar.poswise.randompos)
dfr <- merge(pars.randompos.dfr, dfr1)
# dfr[1:10,]
## looks ok
pars.randompos.dfr <- dfr
# pars.randompos.dfr[1:10,]
rm(dfr, dfr1)

## --------------------------------------------------------------------------------------------------------
## Line
## --------------------------------------------------------------------------------------------------------

## Use the function fitbeta.fun in EMbetaFunctions.R
betafits.line <- apply(fap.line.dfr[1:750,], M = 2, function(vec) try(fitbeta.fun(vec)))
betapar.poswise.line <- matrix(ncol = ncol(fap.line.dfr), nrow = 4)
dimnames(betapar.poswise.line) <- list(c("a.single","a.single.se","b.single","b.single.se"),
   colnames(fap.line.dfr))
betapar.poswise.line[c("a.single","b.single"),] <- sapply(betafits.line, function(lst) 
   if(inherits(lst, "try-error")) rep(NA,2) else lst$par)
betapar.poswise.line[c("a.single.se","b.single.se"),] <- sapply(betafits.line, function(lst) 
   if(inherits(lst, "try-error")) rep(NA,2) else sqrt(diag(solve(-lst$hessian))))

## Attach the estimated beta-distribution coeffs to the pars.line.dfr dataframe:

identical(colnames(betapar.poswise.line), as.character(pars.line.dfr$name))
## True, so it can go by simple merging..
dfr1 <- data.frame(t(betapar.poswise.line))
dfr1$name <- colnames(betapar.poswise.line)
dfr <- merge(pars.line.dfr, dfr1)
# dfr[1:10,]
## looks ok
pars.line.dfr <- dfr
# pars.line.dfr[1:10,]
rm(dfr, dfr1)



## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------
## Classes attached
## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------

## --------------------------------------------------------------------------------------------------------
## Random positions
## --------------------------------------------------------------------------------------------------------

pars.randompos.dfr$class7000_7group <- predict(gmxt7000[[7]],
    newdata = pars.randompos.dfr[, c(7:10,12:23)])$classification
pars.randompos.dfr$class7000_12group <- predict(gmxt7000[[12]],
    newdata = pars.randompos.dfr[, c(7:10,12:23)])$classification
pars.randompos.dfr$class7000_13group <- predict(gmxt7000[[13]],
    newdata = pars.randompos.dfr[, c(7:10,12:23)])$classification
pars.randompos.dfr$class7000_14group <- predict(gmxt7000[[14]],
    newdata = pars.randompos.dfr[, c(7:10,12:23)])$classification
pars.randompos.dfr$class7000_15group <- predict(gmxt7000[[15]],
    newdata = pars.randompos.dfr[, c(7:10,12:23)])$classification
pars.randompos.dfr$class7000_16group <- predict(gmxt7000[[16]],
    newdata = pars.randompos.dfr[, c(7:10,12:23)])$classification
pars.randompos.dfr$class7000_17group <- predict(gmxt7000[[17]],
    newdata = pars.randompos.dfr[, c(7:10,12:23)])$classification
pars.randompos.dfr$class7000_18group <- predict(gmxt7000[[18]],
    newdata = pars.randompos.dfr[, c(7:10,12:23)])$classification
pars.randompos.dfr$class7000_19group <- predict(gmxt7000[[19]],
    newdata = pars.randompos.dfr[, c(7:10,12:23)])$classification
pars.randompos.dfr$class7000_20group <- predict(gmxt7000[[20]],
    newdata = pars.randompos.dfr[, c(7:10,12:23)])$classification

  
table(pars.randompos.dfr$class7000_7group)
## Very, very few in one class (3), only 4. I won't be able to use it as a single class.
## Need to merge it to something else.
  
pars.randompos.dfr$loga.single <- log(pars.randompos.dfr$a.single)
pars.randompos.dfr$logb.single <- log(pars.randompos.dfr$b.single)

## --------------------------------------------------------------------------------------------------------
## Line
## --------------------------------------------------------------------------------------------------------

pars.line.dfr$class7000_7group <- predict(gmxt7000[[7]],
    newdata = pars.line.dfr[, c(7:10,12:23)])$classification
pars.line.dfr$class7000_12group <- predict(gmxt7000[[12]],
    newdata = pars.line.dfr[, c(7:10,12:23)])$classification
pars.line.dfr$class7000_13group <- predict(gmxt7000[[13]],
    newdata = pars.line.dfr[, c(7:10,12:23)])$classification
pars.line.dfr$class7000_14group <- predict(gmxt7000[[14]],
    newdata = pars.line.dfr[, c(7:10,12:23)])$classification
pars.line.dfr$class7000_15group <- predict(gmxt7000[[15]],
    newdata = pars.line.dfr[, c(7:10,12:23)])$classification
pars.line.dfr$class7000_16group <- predict(gmxt7000[[16]],
    newdata = pars.line.dfr[, c(7:10,12:23)])$classification
pars.line.dfr$class7000_17group <- predict(gmxt7000[[17]],
    newdata = pars.line.dfr[, c(7:10,12:23)])$classification
pars.line.dfr$class7000_18group <- predict(gmxt7000[[18]],
    newdata = pars.line.dfr[, c(7:10,12:23)])$classification
pars.line.dfr$class7000_19group <- predict(gmxt7000[[19]],
    newdata = pars.line.dfr[, c(7:10,12:23)])$classification
pars.line.dfr$class7000_20group <- predict(gmxt7000[[20]],
    newdata = pars.line.dfr[, c(7:10,12:23)])$classification

  
table(pars.line.dfr$class7000_7group)
# ## Too few in one class (7), only 17. I will heva to merge it to something else.
  
pars.line.dfr$loga.single <- log(pars.line.dfr$a.single)
pars.line.dfr$logb.single <- log(pars.line.dfr$b.single)


## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------
## Preliminary modelling: variable selection
## Based on grouping (2,3)/(1,4,5,6,7)
## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------

## Optimize wloglik.beta.logpar.fun from the EMbetafunctions.R source file,
## groupwise, and using equal weights (1 for all locations belonging to
## these groups).

## --------------------------------------------------------------------------------------------------------
## Random positions
## --------------------------------------------------------------------------------------------------------

## --------------------------------------------------------------------------------------------------------
## Merging of classes

## Put the new merged classes into the dataframe:

pars.randompos.dfr$fine.class.for.init2[
   is.element(pars.randompos.dfr$class7000_7group, c(2,3))] <- factor(2, levels = c("1","2"))
pars.randompos.dfr$fine.class.for.init2[
   is.element(pars.randompos.dfr$class7000_7group, c(1,4,5,6,7))] <- factor(1, levels = c("1","2"))
  
table(pars.randompos.dfr$fine.class.for.init2)
## 197 in the merged class (2,3)

## --------------------------------------------------------------------------------------------------------
## Model selection and initial parameters for optim by linear regression 

logb.randompos.lm <- vector(17, mode = "list")
names(logb.randompos.lm) <- c("eightvar","eightvar.stepbic","highest",
  "highest2sum","highest3sum","highest4sum","highest5sum","highest6sum","highest7sum",
  "highest8sum","highest9sum","sum.spwin","mean.spwinpeaks","total.spwinpeaks",
  "sum.of.five.spwin","sum.of.six.spwin","sum.of.seven.spwin")

logb.randompos.lm$eightvar <- lm(logb.single ~ (at4 + at12 + at16 + at28 + at40 +
   at52 + at56 + at68)*(fine.class.for.init2), 
   contrasts = contr.treatment, data = pars.randompos.dfr)
logb.randompos.lm$eightvar.stepbic <- stepAIC(logb.randompos.lm$eightvar,
   k = log(nrow(pars.randompos.dfr)), direction = "both")
logb.randompos.lm$highest <- lm(logb.single ~ highest.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
logb.randompos.lm$highest2sum <- lm(logb.single ~ highest2sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
logb.randompos.lm$highest3sum <- lm(logb.single ~ highest3sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
logb.randompos.lm$highest4sum <- lm(logb.single ~ highest4sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
logb.randompos.lm$highest5sum <- lm(logb.single ~ highest5sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
logb.randompos.lm$highest6sum <- lm(logb.single ~ highest6sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
logb.randompos.lm$highest7sum <- lm(logb.single ~ highest7sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
logb.randompos.lm$highest8sum <- lm(logb.single ~ highest8sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
logb.randompos.lm$highest9sum <- lm(logb.single ~ highest9sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
logb.randompos.lm$sum.spwin <- lm(logb.single ~ sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
logb.randompos.lm$mean.spwinpeaks <- lm(logb.single ~ mean.spwinpeaks*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
logb.randompos.lm$total.spwinpeaks <- lm(logb.single ~ total.spwinpeaks*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
logb.randompos.lm$sum.of.five.spwin <- lm(logb.single ~ sum.of.five.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
logb.randompos.lm$sum.of.six.spwin <- lm(logb.single ~ sum.of.six.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
logb.randompos.lm$sum.of.seven.spwin <- lm(logb.single ~ sum.of.seven.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
logb.randompos.lm$sum.of.eight.spwin <- lm(logb.single ~ sum.of.eight.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
logb.randompos.lm$sum.of.nine.spwin <- lm(logb.single ~ sum.of.nine.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
logb.randompos.lm$sum.of.ten.spwin <- lm(logb.single ~ sum.of.ten.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
sapply(logb.randompos.lm, BIC)
summary(logb.randompos.lm$highest6sum)
summary(logb.randompos.lm$highest7sum)
summary(logb.randompos.lm$highest8sum)
summary(logb.randompos.lm$highest9sum)
summary(logb.randompos.lm$sum.of.eight.spwin)
summary(logb.randompos.lm$sum.of.nine.spwin)


loga.randompos.lm <- vector(17, mode = "list")
names(loga.randompos.lm) <- c("eightvar","eightvar.stepbic","highest",
  "highest2sum","highest3sum","highest4sum","highest5sum","highest6sum","highest7sum",
  "highest8sum","highest9sum","sum.spwin","mean.spwinpeaks","total.spwinpeaks",
  "sum.of.five.spwin","sum.of.six.spwin","sum.of.seven.spwin")

loga.randompos.lm$eightvar <- lm(loga.single ~ (at4 + at12 + at16 + at28 + at40 +
   at52 + at56 + at68)*(fine.class.for.init2), 
   contrasts = contr.treatment, data = pars.randompos.dfr)
loga.randompos.lm$eightvar.stepbic <- stepAIC(loga.randompos.lm$eightvar,
   k = log(nrow(pars.randompos.dfr)), direction = "both")
loga.randompos.lm$highest <- lm(loga.single ~ highest.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
loga.randompos.lm$highest2sum <- lm(loga.single ~ highest2sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
loga.randompos.lm$highest3sum <- lm(loga.single ~ highest3sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
loga.randompos.lm$highest4sum <- lm(loga.single ~ highest4sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
loga.randompos.lm$highest5sum <- lm(loga.single ~ highest5sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
loga.randompos.lm$highest6sum <- lm(loga.single ~ highest6sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
loga.randompos.lm$highest7sum <- lm(loga.single ~ highest7sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
loga.randompos.lm$highest8sum <- lm(loga.single ~ highest8sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
loga.randompos.lm$highest9sum <- lm(loga.single ~ highest9sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
loga.randompos.lm$sum.spwin <- lm(loga.single ~ sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
loga.randompos.lm$mean.spwinpeaks <- lm(loga.single ~ mean.spwinpeaks*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
loga.randompos.lm$total.spwinpeaks <- lm(loga.single ~ total.spwinpeaks*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
loga.randompos.lm$sum.of.five.spwin <- lm(loga.single ~ sum.of.five.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
loga.randompos.lm$sum.of.six.spwin <- lm(loga.single ~ sum.of.six.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
loga.randompos.lm$sum.of.seven.spwin <- lm(loga.single ~ sum.of.seven.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
loga.randompos.lm$sum.of.eight.spwin <- lm(loga.single ~ sum.of.eight.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
loga.randompos.lm$sum.of.nine.spwin <- lm(loga.single ~ sum.of.nine.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
loga.randompos.lm$sum.of.ten.spwin <- lm(loga.single ~ sum.of.ten.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
sapply(loga.randompos.lm, BIC)
summary(loga.randompos.lm$highest6sum)
summary(loga.randompos.lm$highest7sum)
summary(loga.randompos.lm$highest8sum)
summary(loga.randompos.lm$highest9sum)
summary(loga.randompos.lm$sum.of.eight.spwin)
summary(loga.randompos.lm$sum.of.nine.spwin)
summary(loga.randompos.lm$sum.of.ten.spwin)

## Looks like here too, the highest6sum, highest7sum, highest8sum, highest9sum are the
## best (BIC-sense),, and that, for both a and b. The class index and its interaction
## with the main variable are significant for b, but not for a.

## Try to include the number of observations: this is determining for the mean of the marginal
## distribution of the periodogram (beta distribution). Is it showing up somewhere in
## maximum distribution?

logb.randompos.lm$highest6sum.n <- lm(logb.single ~ highest6sum.spwin*fine.class.for.init2 + n, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
logb.randompos.lm$highest7sum.n <- lm(logb.single ~ highest7sum.spwin*fine.class.for.init2 + n, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
logb.randompos.lm$highest8sum.n <- lm(logb.single ~ highest8sum.spwin*fine.class.for.init2 + n, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
logb.randompos.lm$highest9sum.n <- lm(logb.single ~ highest9sum.spwin*fine.class.for.init2 + n, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
logb.randompos.lm$sum.of.eight.spwin.n <- lm(logb.single ~ sum.of.eight.spwin*fine.class.for.init2 + n, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
logb.randompos.lm$sum.of.nine.spwin.n <- lm(logb.single ~ sum.of.nine.spwin*fine.class.for.init2 + n, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
logb.randompos.lm$sum.of.ten.spwin.n <- lm(logb.single ~ sum.of.ten.spwin*fine.class.for.init2 + n, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
sapply(logb.randompos.lm, BIC)
summary(logb.randompos.lm$highest6sum.n)
summary(logb.randompos.lm$highest7sum.n)
summary(logb.randompos.lm$highest8sum.n)
summary(logb.randompos.lm$highest9sum.n)
## n is significant, it even almost takes over the role of the class; almost, but not
## completely, it is still significant.

loga.randompos.lm$highest6sum.n <- lm(loga.single ~ highest6sum.spwin*fine.class.for.init2 + n, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
loga.randompos.lm$highest7sum.n <- lm(loga.single ~ highest7sum.spwin*fine.class.for.init2 + n, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
loga.randompos.lm$highest8sum.n <- lm(loga.single ~ highest8sum.spwin*fine.class.for.init2 + n, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
loga.randompos.lm$highest9sum.n <- lm(loga.single ~ highest9sum.spwin*fine.class.for.init2 + n, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
loga.randompos.lm$sum.of.eight.spwin.n <- lm(loga.single ~ sum.of.eight.spwin*fine.class.for.init2 + n, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
loga.randompos.lm$sum.of.nine.spwin.n <- lm(loga.single ~ sum.of.nine.spwin*fine.class.for.init2 + n, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
loga.randompos.lm$sum.of.ten.spwin.n <- lm(loga.single ~ sum.of.ten.spwin*fine.class.for.init2 + n, 
   contrasts = contr.treatment, data = pars.randompos.dfr) 
sapply(loga.randompos.lm, BIC)
summary(loga.randompos.lm$highest6sum.n)
summary(loga.randompos.lm$highest7sum.n)
summary(loga.randompos.lm$highest8sum.n)
summary(loga.randompos.lm$highest9sum.n)
## n is not significant here: apparently, for a, neither n nor the class is not significant

## Total BICs:

total.bic.randompos <- numeric(length(loga.randompos.lm))
for(ii in 1:length(loga.randompos.lm))
   total.bic.randompos[ii] <- BIC(loga.randompos.lm[[ii]]) + BIC(logb.randompos.lm[[ii]])
names(total.bic.randompos) <- names(loga.randompos.lm)


## --------------------------------------------------------------------------------------------------------
## Line
## --------------------------------------------------------------------------------------------------------

## --------------------------------------------------------------------------------------------------------
## Merging of classes

## Put the new merged classes into the dataframe:

pars.line.dfr$fine.class.for.init2[
   is.element(pars.line.dfr$class7000_7group, c(2,3))] <- factor(2, levels = c("1","2"))
pars.line.dfr$fine.class.for.init2[
   is.element(pars.line.dfr$class7000_7group, c(1,4,5,6,7))] <- factor(1, levels = c("1","2"))
  
table(pars.line.dfr$fine.class.for.init2)
## 237 in the merged class (2,3)

## --------------------------------------------------------------------------------------------------------
## Model selection and initial parameters for optim by linear regression 

logb.line.lm <- vector(17, mode = "list")
names(logb.line.lm) <- c("eightvar","eightvar.stepbic","highest",
  "highest2sum","highest3sum","highest4sum","highest5sum","highest6sum","highest7sum",
  "highest8sum","highest9sum","sum.spwin","mean.spwinpeaks","total.spwinpeaks",
  "sum.of.five.spwin","sum.of.six.spwin","sum.of.seven.spwin")

logb.line.lm$eightvar <- lm(logb.single ~ (at4 + at12 + at16 + at28 + at40 +
   at52 + at56 + at68)*(fine.class.for.init2), 
   contrasts = contr.treatment, data = pars.line.dfr)
logb.line.lm$eightvar.stepbic <- stepAIC(logb.line.lm$eightvar,
   k = log(nrow(pars.line.dfr)), direction = "both")
logb.line.lm$highest <- lm(logb.single ~ highest.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
logb.line.lm$highest2sum <- lm(logb.single ~ highest2sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
logb.line.lm$highest3sum <- lm(logb.single ~ highest3sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
logb.line.lm$highest4sum <- lm(logb.single ~ highest4sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
logb.line.lm$highest5sum <- lm(logb.single ~ highest5sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
logb.line.lm$highest6sum <- lm(logb.single ~ highest6sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
logb.line.lm$highest7sum <- lm(logb.single ~ highest7sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
logb.line.lm$highest8sum <- lm(logb.single ~ highest8sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
logb.line.lm$highest9sum <- lm(logb.single ~ highest9sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
logb.line.lm$sum.spwin <- lm(logb.single ~ sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
logb.line.lm$mean.spwinpeaks <- lm(logb.single ~ mean.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
logb.line.lm$total.spwinpeaks <- lm(logb.single ~ total.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
logb.line.lm$sum.of.five.spwin <- lm(logb.single ~ sum.of.five.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
logb.line.lm$sum.of.six.spwin <- lm(logb.single ~ sum.of.six.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
logb.line.lm$sum.of.seven.spwin <- lm(logb.single ~ sum.of.seven.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
logb.line.lm$sum.of.eight.spwin <- lm(logb.single ~ sum.of.eight.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
logb.line.lm$sum.of.nine.spwin <- lm(logb.single ~ sum.of.nine.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
logb.line.lm$sum.of.ten.spwin <- lm(logb.single ~ sum.of.ten.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
sapply(logb.line.lm, BIC)
summary(logb.line.lm$eightvar.stepbic)
summary(logb.line.lm$highest6sum)
summary(logb.line.lm$sum.of.five.spwin)


loga.line.lm <- vector(17, mode = "list")
names(loga.line.lm) <- c("eightvar","eightvar.stepbic","highest",
  "highest2sum","highest3sum","highest4sum","highest5sum","highest6sum","highest7sum",
  "highest8sum","highest9sum","sum.spwin","mean.spwinpeaks","total.spwinpeaks",
  "sum.of.five.spwin","sum.of.six.spwin","sum.of.seven.spwin")

loga.line.lm$eightvar <- lm(loga.single ~ (at4 + at12 + at16 + at28 + at40 +
   at52 + at56 + at68)*(fine.class.for.init2), 
   contrasts = contr.treatment, data = pars.line.dfr)
loga.line.lm$eightvar.stepbic <- stepAIC(loga.line.lm$eightvar,
   k = log(nrow(pars.line.dfr)), direction = "both")
loga.line.lm$highest <- lm(loga.single ~ highest.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
loga.line.lm$highest2sum <- lm(loga.single ~ highest2sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
loga.line.lm$highest3sum <- lm(loga.single ~ highest3sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
loga.line.lm$highest4sum <- lm(loga.single ~ highest4sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
loga.line.lm$highest5sum <- lm(loga.single ~ highest5sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
loga.line.lm$highest6sum <- lm(loga.single ~ highest6sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
loga.line.lm$highest7sum <- lm(loga.single ~ highest7sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
loga.line.lm$highest8sum <- lm(loga.single ~ highest8sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
loga.line.lm$highest9sum <- lm(loga.single ~ highest9sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
loga.line.lm$sum.spwin <- lm(loga.single ~ sum.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
loga.line.lm$mean.spwinpeaks <- lm(loga.single ~ mean.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
loga.line.lm$total.spwinpeaks <- lm(loga.single ~ total.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
loga.line.lm$sum.of.five.spwin <- lm(loga.single ~ sum.of.five.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
loga.line.lm$sum.of.six.spwin <- lm(loga.single ~ sum.of.six.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
loga.line.lm$sum.of.seven.spwin <- lm(loga.single ~ sum.of.seven.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
loga.line.lm$sum.of.eight.spwin <- lm(loga.single ~ sum.of.eight.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
loga.line.lm$sum.of.nine.spwin <- lm(loga.single ~ sum.of.nine.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
loga.line.lm$sum.of.ten.spwin <- lm(loga.single ~ sum.of.ten.spwin*fine.class.for.init2, 
   contrasts = contr.treatment, data = pars.line.dfr) 
sapply(loga.line.lm, BIC)
summary(loga.line.lm$highest6sum)
summary(loga.line.lm$highest7sum)
summary(loga.line.lm$highest8sum)
summary(loga.line.lm$highest9sum)
summary(loga.line.lm$sum.of.nine.spwin)
summary(loga.line.lm$sum.of.ten.spwin)

## Results are similar to the random positions, but the class index and its interaction
## with the main variable are significant now for both, and that, very much. The 50 more
## examples of the pole time samplings might do this.


## Try to include the number of observations: this is determining for the mean of the marginal
## distribution of the periodogram (beta distribution). Is it showing up somewhere in
## maximum distribution?

logb.line.lm$eightvar.stepbic.n <- update(logb.line.lm$eightvar.stepbic, . ~ . + n, 
   contrasts = contr.treatment, data = pars.line.dfr) 
logb.line.lm$highest6sum.n <- lm(logb.single ~ highest6sum.spwin*fine.class.for.init2 + n, 
   contrasts = contr.treatment, data = pars.line.dfr) 
logb.line.lm$highest7sum.n <- lm(logb.single ~ highest7sum.spwin*fine.class.for.init2 + n, 
   contrasts = contr.treatment, data = pars.line.dfr) 
logb.line.lm$highest8sum.n <- lm(logb.single ~ highest8sum.spwin*fine.class.for.init2 + n, 
   contrasts = contr.treatment, data = pars.line.dfr) 
logb.line.lm$highest9sum.n <- lm(logb.single ~ highest9sum.spwin*fine.class.for.init2 + n, 
   contrasts = contr.treatment, data = pars.line.dfr) 
logb.line.lm$sum.of.eight.spwin.n <- lm(logb.single ~ sum.of.eight.spwin*fine.class.for.init2 + n, 
   contrasts = contr.treatment, data = pars.line.dfr) 
logb.line.lm$sum.of.nine.spwin.n <- lm(logb.single ~ sum.of.nine.spwin*fine.class.for.init2 + n, 
   contrasts = contr.treatment, data = pars.line.dfr) 
logb.line.lm$sum.of.ten.spwin.n <- lm(logb.single ~ sum.of.ten.spwin*fine.class.for.init2 + n, 
   contrasts = contr.treatment, data = pars.line.dfr) 
sapply(logb.line.lm, BIC)
summary(logb.line.lm$eightvar.stepbic.n)
summary(logb.line.lm$highest6sum.n)
summary(logb.line.lm$highest7sum.n)
summary(logb.line.lm$highest8sum.n)
summary(logb.line.lm$highest9sum.n)
## n is significant, but here the class also stays significant.

loga.line.lm$eightvar.stepbic.n <- update(loga.line.lm$eightvar.stepbic, . ~ . + n, 
   contrasts = contr.treatment, data = pars.line.dfr) 
loga.line.lm$highest6sum.n <- lm(loga.single ~ highest6sum.spwin*fine.class.for.init2 + n, 
   contrasts = contr.treatment, data = pars.line.dfr) 
loga.line.lm$highest7sum.n <- lm(loga.single ~ highest7sum.spwin*fine.class.for.init2 + n, 
   contrasts = contr.treatment, data = pars.line.dfr) 
loga.line.lm$highest8sum.n <- lm(loga.single ~ highest8sum.spwin*fine.class.for.init2 + n, 
   contrasts = contr.treatment, data = pars.line.dfr) 
loga.line.lm$highest9sum.n <- lm(loga.single ~ highest9sum.spwin*fine.class.for.init2 + n, 
   contrasts = contr.treatment, data = pars.line.dfr) 
loga.line.lm$sum.of.eight.spwin.n <- lm(loga.single ~ sum.of.eight.spwin*fine.class.for.init2 + n, 
   contrasts = contr.treatment, data = pars.line.dfr) 
loga.line.lm$sum.of.nine.spwin.n <- lm(loga.single ~ sum.of.nine.spwin*fine.class.for.init2 + n, 
   contrasts = contr.treatment, data = pars.line.dfr) 
loga.line.lm$sum.of.ten.spwin.n <- lm(loga.single ~ sum.of.ten.spwin*fine.class.for.init2 + n, 
   contrasts = contr.treatment, data = pars.line.dfr) 
sapply(loga.line.lm, BIC)
summary(loga.line.lm$eightvar.stepbic.n)
summary(loga.line.lm$highest6sum.n)
summary(loga.line.lm$highest7sum.n)
summary(loga.line.lm$highest8sum.n)
summary(loga.line.lm$highest9sum.n)
## n improves a little on BIC, and is significant everywhere. The class also stays
## significant.

## Total BICs:

total.bic.line <- numeric(length(loga.line.lm))
for(ii in 1:length(loga.line.lm))
   total.bic.line[ii] <- BIC(loga.line.lm[[ii]]) + BIC(logb.line.lm[[ii]])
names(total.bic.line) <- names(loga.line.lm)

total.bic.line[-21] + total.bic.randompos 
total.bic.line[-21]
total.bic.randompos

## CONCLUSION:
##
## IF I WANT A COMMON MODEL TO PUT ON a AND b, THE BEST CHOICE IS THE highest6sum.n.
## THIS IS BOTH FOR THE LINE AND THE RANDOM POSITIONS.



###########################################################################################################
## --------------------------------------------------------------------------------------------------------
## Fitting the beta log-likelihood
## --------------------------------------------------------------------------------------------------------
###########################################################################################################

## --------------------------------------------------------------------------------------------------------
## Lists containing FAPs and covariates for each location
## --------------------------------------------------------------------------------------------------------

vars.list.randompos <- lapply(pars.randompos.dfr$name, function(chr)
   {
   	f.tmp <- fap.randompos.dfr[, as.character(chr)]
   	f.tmp[f.tmp == 1] <- 1-1e-11
	f.tmp[f.tmp == 0] <- 1e-11
    list(fap = f.tmp, 
        covar = pars.randompos.dfr[pars.randompos.dfr$name == chr, c("highest6sum.spwin","n")])
   })  
sapply(vars.list.randompos, function(lst) length(lst[[1]]))
vars.list.randompos[[1]]

## Take care to use only the first 750 of this!
vars.list.line <- lapply(pars.line.dfr$name, function(chr)
   {
   	f.tmp <- fap.line.dfr[1:750, as.character(chr)]
   	f.tmp[f.tmp == 1] <- 1-1e-11
	f.tmp[f.tmp == 0] <- 1e-11
    list(fap = f.tmp, 
        covar = pars.line.dfr[pars.line.dfr$name == chr, c("highest6sum.spwin","n")])
   })
sapply(vars.list.line, function(lst) length(lst[[1]]))


## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------
## Grouping (2,3)/(1,4,5,6,7)
## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------

## --------------------------------------------------------------------------------------------------------
## Initial values in the format needed for the log-likelihood function
## --------------------------------------------------------------------------------------------------------

lm.loga.randompos.init2 <- vector(2, mode = "list")
for(ii in 1:2)
   {
   	lm.loga.randompos.init2[[ii]] <- lm(loga.single ~ highest6sum.spwin + n,
   	   data = pars.randompos.dfr, 
   	   subset = pars.randompos.dfr$fine.class.for.init2 == factor(ii, levels = c(1,2)))
   }  
lm.logb.randompos.init2 <- vector(2, mode = "list")
for(ii in 1:2)
   {
   	lm.logb.randompos.init2[[ii]] <- lm(logb.single ~ highest6sum.spwin + n,
   	   data = pars.randompos.dfr, 
   	   subset = pars.randompos.dfr$fine.class.for.init2 == factor(ii, levels = c(1,2)))
   }
lapply(lm.loga.randompos.init2, summary) 
lapply(lm.logb.randompos.init2, summary) 

lm.loga.line.init2 <- vector(2, mode = "list")
for(ii in 1:2)
   {
   	lm.loga.line.init2[[ii]] <- lm(loga.single ~ highest6sum.spwin + n,
   	   data = pars.line.dfr, 
   	   subset = pars.line.dfr$fine.class.for.init2 == factor(ii, levels = c(1,2)))
   }  
lm.logb.line.init2 <- vector(2, mode = "list")
for(ii in 1:2)
   {
   	lm.logb.line.init2[[ii]] <- lm(logb.single ~ highest6sum.spwin + n,
   	   data = pars.line.dfr, 
   	   subset = pars.line.dfr$fine.class.for.init2 == factor(ii, levels = c(1,2)))
   }    
lapply(lm.loga.line.init2, summary) 
lapply(lm.logb.line.init2, summary) 
  




c.tmp <- rainbow(2)[as.numeric(pars.randompos.dfr$fine.class.for.init2)]
ctr.tmp <- rainbow(2, alpha = 0.6)[as.numeric(pars.randompos.dfr$fine.class.for.init2)]
c.tmp2 <- rainbow(2)[as.numeric(pars.line.dfr$fine.class.for.init2)]
ctr.tmp2 <- rainbow(2, alpha = 0.6)[as.numeric(pars.line.dfr$fine.class.for.init2)]
new.tmp <- data.frame(highest6sum.spwin = rep(seq(2.5, 5.7, by = 0.005), 2), 
   n = rep(c(60,120), each = length(seq(2.5, 5.7, by = 0.005))))
new.tmp2 <- data.frame(highest6sum.spwin = rep(c(3.2,4.5), each = 181), 
   n = rep(40:220, 2))

quartz(height = 7.5, width = 10)
#dev.set(2)
par(mfrow = c(2,2), mar = c(2,2,2,1), mgp = c(1.8,0.9,0), oma = c(0,0,2,0))
## dependence of a on highest6sum.spwin:
plot(pars.randompos.dfr$highest6sum.spwin, pars.randompos.dfr$loga.single, pch = 16,
   cex = 0.7, col = ctr.tmp, main = "Random pos., log(a) ~ highest6sum.spwin", ylim = c(-0.4,-0.03))
for(ii in 1:2)
  {
   m.tmp <- predict(lm.loga.randompos.init2[[ii]], newdata = new.tmp[1:(nrow(new.tmp)/2),])
   lines(new.tmp[1:(nrow(new.tmp)/2),1], m.tmp, lwd = 2, col = rainbow(2)[ii], lty = 2) 
   m.tmp <- predict(lm.loga.randompos.init2[[ii]], newdata = new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),])
   lines(new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),1], m.tmp, lwd = 2, col = rainbow(2)[ii]) 
  }
legend(4.5, -0.05, bty = "n", legend = c("1, n = 60", "1, n = 120", "2, n = 60", "2, n = 120"), 
   col = c(rainbow(2)[1],rainbow(2)[1],rainbow(2)[2],rainbow(2)[2]), lty = c(2,1,2,1))
## dependence of a on n:
plot(pars.randompos.dfr$n, pars.randompos.dfr$loga.single, pch = 16,
   cex = 0.7, col = ctr.tmp, main = "Random pos., log(a) ~ n")
for(ii in 1:2)
  {
   m.tmp <- predict(lm.loga.randompos.init2[[ii]], newdata = new.tmp2[1:(nrow(new.tmp2)/2),])
   lines(new.tmp2[1:(nrow(new.tmp2)/2),2], m.tmp, lwd = 2, col = rainbow(2)[ii], lty = 2) 
   m.tmp <- predict(lm.loga.randompos.init2[[ii]], newdata = new.tmp2[(nrow(new.tmp2)/2+1):nrow(new.tmp2),])
   lines(new.tmp2[(nrow(new.tmp2)/2+1):nrow(new.tmp2),2], m.tmp, lwd = 2, col = rainbow(2)[ii]) 
  }
legend(160, -0.05, bty = "n", legend = c("1, sum = 3.2", "1, sum = 4.5", "2, sum = 3.2", "2, sum = 4.5"), 
   col = c(rainbow(2)[1],rainbow(2)[1],rainbow(2)[2],rainbow(2)[2]), lty = c(2,1,2,1))

plot(pars.line.dfr$highest6sum.spwin, pars.line.dfr$loga.single, pch = 16,
   cex = 0.7, col = ctr.tmp2, main = "Line, log(a) ~ highest6sum.spwin", ylim = c(-0.4,-0.03))
for(ii in 1:2)
  {
   m.tmp <- predict(lm.loga.line.init2[[ii]], newdata = new.tmp[1:(nrow(new.tmp)/2),])
   lines(new.tmp[1:(nrow(new.tmp)/2),1], m.tmp, lwd = 2, col = rainbow(2)[ii], lty = 2) 
   m.tmp <- predict(lm.loga.line.init2[[ii]], newdata = new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),])
   lines(new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),1], m.tmp, lwd = 2, col = rainbow(2)[ii]) 
  }
legend(4.5, -0.05, bty = "n", legend = c("1, n = 60", "1, n = 120", "2, n = 60", "2, n = 120"), 
   col = c(rainbow(2)[1],rainbow(2)[1],rainbow(2)[2],rainbow(2)[2]), lty = c(2,1,2,1))
## dependence of b on n:
plot(pars.line.dfr$n, pars.line.dfr$loga.single, pch = 16,
   cex = 0.7, col = ctr.tmp2, main = "Line, log(a) ~ n")
for(ii in 1:2)
  {
   m.tmp <- predict(lm.loga.line.init2[[ii]], newdata = new.tmp2[1:(nrow(new.tmp2)/2),])
   lines(new.tmp2[1:(nrow(new.tmp2)/2),2], m.tmp, lwd = 2, col = rainbow(2)[ii], lty = 2) 
   m.tmp <- predict(lm.loga.line.init2[[ii]], newdata = new.tmp2[(nrow(new.tmp2)/2+1):nrow(new.tmp2),])
   lines(new.tmp2[(nrow(new.tmp2)/2+1):nrow(new.tmp2),2], m.tmp, lwd = 2, col = rainbow(2)[ii]) 
  }
legend(160, -0.05, bty = "n", legend = c("1, sum = 3.2", "1, sum = 4.5", "2, sum = 3.2", "2, sum = 4.5"), 
   col = c(rainbow(2)[1],rainbow(2)[1],rainbow(2)[2],rainbow(2)[2]), lty = c(2,1,2,1))
mtext("(2,3)/(4,5,6,7)", outer = TRUE)


quartz(height = 7.5, width = 10)
#dev.set(2)
par(mfrow = c(2,2), mar = c(2,2,2,1), mgp = c(1.8,0.9,0), oma = c(0,0,2,0))
## dependence of a on highest6sum.spwin:
plot(pars.randompos.dfr$highest6sum.spwin, pars.randompos.dfr$logb.single, pch = 16,
   cex = 0.7, col = ctr.tmp, main = "Random pos., log(b) ~ highest6sum.spwin", ylim = c(-0.6,1))
for(ii in 1:2)
  {
   m.tmp <- predict(lm.logb.randompos.init2[[ii]], newdata = new.tmp[1:(nrow(new.tmp)/2),])
   lines(new.tmp[1:(nrow(new.tmp)/2),1], m.tmp, lwd = 2, col = rainbow(2)[ii], lty = 2) 
   m.tmp <- predict(lm.logb.randompos.init2[[ii]], newdata = new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),])
   lines(new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),1], m.tmp, lwd = 2, col = rainbow(2)[ii]) 
  }
legend(4.5, 1, bty = "n", legend = c("1, n = 60", "1, n = 120", "2, n = 60", "2, n = 120"), 
   col = c(rainbow(2)[1],rainbow(2)[1],rainbow(2)[2],rainbow(2)[2]), lty = c(2,1,2,1))
## dependence of a on n:
plot(pars.randompos.dfr$n, pars.randompos.dfr$logb.single, pch = 16, ylim = c(-0.6,1),
   cex = 0.7, col = ctr.tmp, main = "Random pos., log(b) ~ n")
for(ii in 1:2)
  {
   m.tmp <- predict(lm.logb.randompos.init2[[ii]], newdata = new.tmp2[1:(nrow(new.tmp2)/2),])
   lines(new.tmp2[1:(nrow(new.tmp2)/2),2], m.tmp, lwd = 2, col = rainbow(2)[ii], lty = 2) 
   m.tmp <- predict(lm.logb.randompos.init2[[ii]], newdata = new.tmp2[(nrow(new.tmp2)/2+1):nrow(new.tmp2),])
   lines(new.tmp2[(nrow(new.tmp2)/2+1):nrow(new.tmp2),2], m.tmp, lwd = 2, col = rainbow(2)[ii]) 
  }
legend(160, -0.05, bty = "n", legend = c("1, sum = 3.2", "1, sum = 4.5", "2, sum = 3.2", "2, sum = 4.5"), 
   col = c(rainbow(2)[1],rainbow(2)[1],rainbow(2)[2],rainbow(2)[2]), lty = c(2,1,2,1))

plot(pars.line.dfr$highest6sum.spwin, pars.line.dfr$logb.single, pch = 16,
   cex = 0.7, col = ctr.tmp2, main = "Line, log(b) ~ highest6sum.spwin", ylim = c(-0.6,1))
for(ii in 1:2)
  {
   m.tmp <- predict(lm.logb.line.init2[[ii]], newdata = new.tmp[1:(nrow(new.tmp)/2),])
   lines(new.tmp[1:(nrow(new.tmp)/2),1], m.tmp, lwd = 2, col = rainbow(2)[ii], lty = 2) 
   m.tmp <- predict(lm.logb.line.init2[[ii]], newdata = new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),])
   lines(new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),1], m.tmp, lwd = 2, col = rainbow(2)[ii]) 
  }
legend(4.5, -0.05, bty = "n", legend = c("1, n = 60", "1, n = 120", "2, n = 60", "2, n = 120"), 
   col = c(rainbow(2)[1],rainbow(2)[1],rainbow(2)[2],rainbow(2)[2]), lty = c(2,1,2,1))
## dependence of b on n:
plot(pars.line.dfr$n, pars.line.dfr$logb.single, pch = 16,
   cex = 0.7, col = ctr.tmp2, main = "Line, log(b) ~ n")
for(ii in 1:2)
  {
   m.tmp <- predict(lm.logb.line.init2[[ii]], newdata = new.tmp2[1:(nrow(new.tmp2)/2),])
   lines(new.tmp2[1:(nrow(new.tmp2)/2),2], m.tmp, lwd = 2, col = rainbow(2)[ii], lty = 2) 
   m.tmp <- predict(lm.logb.line.init2[[ii]], newdata = new.tmp2[(nrow(new.tmp2)/2+1):nrow(new.tmp2),])
   lines(new.tmp2[(nrow(new.tmp2)/2+1):nrow(new.tmp2),2], m.tmp, lwd = 2, col = rainbow(2)[ii]) 
  }
legend(160, -0.05, bty = "n", legend = c("1, sum = 3.2", "1, sum = 4.5", "2, sum = 3.2", "2, sum = 4.5"), 
   col = c(rainbow(2)[1],rainbow(2)[1],rainbow(2)[2],rainbow(2)[2]), lty = c(2,1,2,1))
mtext("(2,3)/(4,5,6,7)", outer = TRUE)



# colnames(faps.line.dfr)
# identical(as.character(pars.line.dfr$name), colnames(fap.line.dfr))

pars.init.line2 <- matrix(c(lm.loga.line.init2[[1]]$coef,
    lm.logb.line.init2[[1]]$coef,
    lm.loga.line.init2[[2]]$coef,
    lm.logb.line.init2[[2]]$coef), 
    ncol = 6, byrow = TRUE)
pars.init.line2

pars.init.randompos2 <- matrix(c(lm.loga.randompos.init2[[1]]$coef,
    lm.logb.randompos.init2[[1]]$coef,
    lm.loga.randompos.init2[[2]]$coef,
    lm.logb.randompos.init2[[2]]$coef), 
    ncol = 6, byrow = TRUE)
pars.init.randompos2


## --------------------------------------------------------------------------------------------------------
## Optimization of the beta likelihood
## --------------------------------------------------------------------------------------------------------
 
pars.randompos2 <- list()
system.time(
pars.randompos2[[1]] <- optim(pars.init.randompos2[1,], fn = wloglik.beta.logpar.fun, gr = wloglikgrad.logpar.fun,
   vars.list = vars.list.randompos[pars.randompos.dfr$fine.class.for.init2 == 1], 
   weights = rep(1, sum(pars.randompos.dfr$fine.class.for.init2 == 1)),
   control = list(maxit = 10000, fnscale = -1))
  )
system.time(
pars.randompos2[[2]] <- optim(pars.init.randompos2[2,], fn = wloglik.beta.logpar.fun, gr = wloglikgrad.logpar.fun,
   vars.list = vars.list.randompos[pars.randompos.dfr$fine.class.for.init2 == 2], 
   weights = rep(1, sum(pars.randompos.dfr$fine.class.for.init2 == 2)),
   control = list(maxit = 10000, fnscale = -1))
  )
 
pars.line2 <- list()
system.time(
pars.line2[[1]] <- optim(pars.init.line2[1,], fn = wloglik.beta.logpar.fun, gr = wloglikgrad.logpar.fun,
   vars.list = vars.list.line[pars.line.dfr$fine.class.for.init2 == 1], 
   weights = rep(1, sum(pars.line.dfr$fine.class.for.init2 == 1)),
   control = list(maxit = 10000, fnscale = -1))
  )
system.time(
pars.line2[[2]] <- optim(pars.init.line2[2,], fn = wloglik.beta.logpar.fun, gr = wloglikgrad.logpar.fun,
   vars.list = vars.list.line[pars.line.dfr$fine.class.for.init2 == 2], 
   weights = rep(1, sum(pars.line.dfr$fine.class.for.init2 == 2)),
   control = list(maxit = 10000, fnscale = -1))
  )

## --------------------------------------------------------------------------------------------------------
## Plot of the points and the lines  
## --------------------------------------------------------------------------------------------------------

c.tmp <- rgb(blue = (as.numeric(pars.randompos.dfr$fine.class.for.init2)-1), 
   red = (2-as.numeric(pars.randompos.dfr$fine.class.for.init2)),
   green = (1-(pars.randompos.dfr$n-40) / 180)^2, alpha = 1)
c.tmp2 <- rgb(blue = (as.numeric(pars.line.dfr$fine.class.for.init2)-1), 
   red = (2-as.numeric(pars.line.dfr$fine.class.for.init2)),
   green = (1-(pars.line.dfr$n-40) / 180)^2, alpha = 1)
new.tmp <- data.frame(highest6sum.spwin = rep(seq(2.5, 5.7, by = 0.005), 2), 
   n = rep(c(60,120), each = length(seq(2.5, 5.7, by = 0.005))))
csmalln.tmp <- c(rgb(blue = 0, red = 1, green = (1-(60-40) / 180)^2, alpha = 1),
    rgb(blue = 1, red = 0, green = (1-(60-40) / 180)^2, alpha = 1))
clargen.tmp <- c(rgb(blue = 0, red = 1, green = (1-(120-40) / 180)^2, alpha = 1),
    rgb(blue = 1, red = 0, green = (1-(120-40) / 180)^2, alpha = 1))

quartz(height = 7.5, width = 10)
#dev.set(2)
par(mfrow = c(2,2), mar = c(2,2,3,1), mgp = c(1.8,0.9,0), oma = c(0,0,2,0))
## dependence of a on highest6sum.spwin:
plot(pars.randompos.dfr$highest6sum.spwin, pars.randompos.dfr$logb.single, pch = 16,
   cex = 0.7, col = c.tmp, main = "Random positions \n log(b) ~ highest6sum.spwin", 
   ylim = c(-0.6,1), xlab = "", ylab = "")
for(ii in 1:2)
  {
   m.tmp <- predict(lm.logb.randompos.init2[[ii]], newdata = new.tmp[1:(nrow(new.tmp)/2),])
   lines(new.tmp[1:(nrow(new.tmp)/2),1], m.tmp, lwd = 2, col = csmalln.tmp[ii], lty = 2) 
   m.tmp <- predict(lm.logb.randompos.init2[[ii]], newdata = new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),])
   lines(new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),1], m.tmp, lwd = 2, col = clargen.tmp[ii]) 
  }
legend(4.5, 1, bty = "n", legend = c("1, n = 60", "1, n = 120", "2, n = 60", "2, n = 120"), 
   col = c(csmalln.tmp[1],clargen[1], csmalln.tmp[2],clargen[2]), lty = c(2,1,2,1))

plot(pars.line.dfr$highest6sum.spwin, pars.line.dfr$logb.single, pch = 16,
   cex = 0.7, col = c.tmp2, main = "Line \n log(b) ~ highest6sum.spwin", 
   ylim = c(-0.6,1), xlab = "", ylab = "")
for(ii in 1:2)
  {
   m.tmp <- predict(lm.logb.line.init2[[ii]], newdata = new.tmp[1:(nrow(new.tmp)/2),])
   lines(new.tmp[1:(nrow(new.tmp)/2),1], m.tmp, lwd = 2, col = csmalln.tmp[ii], lty = 2) 
   m.tmp <- predict(lm.logb.line.init2[[ii]], newdata = new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),])
   lines(new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),1], m.tmp, lwd = 2, col = clargen.tmp[ii]) 
  }
legend(4.5, 1, bty = "n", legend = c("1, n = 60", "1, n = 120", "2, n = 60", "2, n = 120"), 
   col = c(csmalln.tmp[1],clargen[1], csmalln.tmp[2],clargen[2]), lty = c(2,1,2,1))

plot(pars.randompos.dfr$highest6sum.spwin, pars.randompos.dfr$loga.single, pch = 16,
   cex = 0.7, col = c.tmp, main = "\n log(a) ~ highest6sum.spwin", ylim = c(-0.4,-0.03), 
   xlab = "", ylab = "")
for(ii in 1:2)
  {
   m.tmp <- predict(lm.loga.randompos.init2[[ii]], newdata = new.tmp[1:(nrow(new.tmp)/2),])
   lines(new.tmp[1:(nrow(new.tmp)/2),1], m.tmp, lwd = 2, col = csmalln.tmp[ii], lty = 2) 
   m.tmp <- predict(lm.loga.randompos.init2[[ii]], newdata = new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),])
   lines(new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),1], m.tmp, lwd = 2, col = clargen.tmp[ii]) 
  }
legend(4.5, -0.05, bty = "n", legend = c("1, n = 60", "1, n = 120", "2, n = 60", "2, n = 120"), 
   col = c(csmalln.tmp[1],clargen[1], csmalln.tmp[2],clargen[2]), lty = c(2,1,2,1))

plot(pars.line.dfr$highest6sum.spwin, pars.line.dfr$loga.single, pch = 16,
   cex = 0.7, col = c.tmp2, main = "\n log(a) ~ highest6sum.spwin", ylim = c(-0.4,-0.03), 
   xlab = "", ylab = "")
for(ii in 1:2)
  {
   m.tmp <- predict(lm.loga.line.init2[[ii]], newdata = new.tmp[1:(nrow(new.tmp)/2),])
   lines(new.tmp[1:(nrow(new.tmp)/2),1], m.tmp, lwd = 2, col = csmalln.tmp[ii], lty = 2) 
   m.tmp <- predict(lm.loga.line.init2[[ii]], newdata = new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),])
   lines(new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),1], m.tmp, lwd = 2, col = clargen.tmp[ii]) 
  }
legend(4.5, -0.05, bty = "n", legend = c("1, n = 60", "1, n = 120", "2, n = 60", "2, n = 120"), 
   col = c(csmalln.tmp[1],clargen[1], csmalln.tmp[2],clargen[2]), lty = c(2,1,2,1))
mtext("(2,3)/(4,5,6,7)", outer = TRUE, font = 2)


pars.init.randompos2
pars.randompos2[[1]]$par
pars.randompos2[[2]]$par
pars.init.line2 
pars.line2[[1]]$par
pars.line2[[2]]$par

## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------
## Grouping (1,3)/(2,4,5,6,7)
## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------

## Optimize wloglik.beta.logpar.fun from the EMbetafunctions.R source file,
## groupwise, and using equal weights (1 for all locations belonging to
## these groups).

## --------------------------------------------------------------------------------------------------------
## Merging of classes
## --------------------------------------------------------------------------------------------------------

## Put the new merged classes into the dataframe:

pars.randompos.dfr$fine.class.for.init3[
   is.element(pars.randompos.dfr$class7000_7group, c(1,3))] <- 2
pars.randompos.dfr$fine.class.for.init3[
   is.element(pars.randompos.dfr$class7000_7group, c(2,4,5,6,7))] <- 1  

pars.line.dfr$fine.class.for.init3[
   is.element(pars.line.dfr$class7000_7group, c(1,3))] <- 2
pars.line.dfr$fine.class.for.init3[
   is.element(pars.line.dfr$class7000_7group, c(2,4,5,6,7))] <- 1
  
table(pars.randompos.dfr$fine.class.for.init3)
table(pars.line.dfr$fine.class.for.init3)


## --------------------------------------------------------------------------------------------------------
## Initial values in the format needed for the log-likelihood function
## --------------------------------------------------------------------------------------------------------

lm.loga.randompos.init3 <- vector(2, mode = "list")
for(ii in 1:2)
   {
   	lm.loga.randompos.init3[[ii]] <- lm(loga.single ~ highest6sum.spwin + n,
   	   data = pars.randompos.dfr, 
   	   subset = pars.randompos.dfr$fine.class.for.init3 == factor(ii, levels = c(1,2)))
   }  
lm.logb.randompos.init3 <- vector(2, mode = "list")
for(ii in 1:2)
   {
   	lm.logb.randompos.init3[[ii]] <- lm(logb.single ~ highest6sum.spwin + n,
   	   data = pars.randompos.dfr, 
   	   subset = pars.randompos.dfr$fine.class.for.init3 == factor(ii, levels = c(1,2)))
   }
lapply(lm.loga.randompos.init3, summary) 
lapply(lm.logb.randompos.init3, summary) 

lm.loga.line.init3 <- vector(2, mode = "list")
for(ii in 1:2)
   {
   	lm.loga.line.init3[[ii]] <- lm(loga.single ~ highest6sum.spwin + n,
   	   data = pars.line.dfr, 
   	   subset = pars.line.dfr$fine.class.for.init3 == factor(ii, levels = c(1,2)))
   }  
lm.logb.line.init3 <- vector(2, mode = "list")
for(ii in 1:2)
   {
   	lm.logb.line.init3[[ii]] <- lm(logb.single ~ highest6sum.spwin + n,
   	   data = pars.line.dfr, 
   	   subset = pars.line.dfr$fine.class.for.init3 == factor(ii, levels = c(1,2)))
   }    
lapply(lm.loga.line.init3, summary) 
lapply(lm.logb.line.init3, summary) 
  

## --------------------------------------------------------------------------------------------------------
## Plots of the obtained models: location-wise estimates of the beta parameters
## and the fitted lines
## --------------------------------------------------------------------------------------------------------


## Parameter a:

c.tmp <- rainbow(2)[as.numeric(pars.randompos.dfr$fine.class.for.init3)]
ctr.tmp <- rainbow(2, alpha = 0.6)[as.numeric(pars.randompos.dfr$fine.class.for.init3)]
c.tmp2 <- rainbow(2)[as.numeric(pars.line.dfr$fine.class.for.init3)]
ctr.tmp2 <- rainbow(2, alpha = 0.6)[as.numeric(pars.line.dfr$fine.class.for.init3)]
new.tmp <- data.frame(highest6sum.spwin = rep(seq(2.5, 5.7, by = 0.005), 2), 
   n = rep(c(60,120), each = length(seq(2.5, 5.7, by = 0.005))))
new.tmp2 <- data.frame(highest6sum.spwin = rep(c(3.2,4.5), each = 181), 
   n = rep(40:220, 2))

quartz(height = 7.5, width = 10)
#dev.set(2)
par(mfrow = c(2,2), mar = c(2,2,2,1), mgp = c(1.8,0.9,0), oma = c(0,0,2,0))
## dependence of a on highest6sum.spwin:
plot(pars.randompos.dfr$highest6sum.spwin, pars.randompos.dfr$loga.single, pch = 16,
   cex = 0.7, col = ctr.tmp, main = "Random pos., log(a) ~ highest6sum.spwin", ylim = c(-0.4,-0.03))
for(ii in 1:2)
  {
   m.tmp <- predict(lm.loga.randompos.init3[[ii]], newdata = new.tmp[1:(nrow(new.tmp)/2),])
   lines(new.tmp[1:(nrow(new.tmp)/2),1], m.tmp, lwd = 2, col = rainbow(2)[ii], lty = 2) 
   m.tmp <- predict(lm.loga.randompos.init3[[ii]], newdata = new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),])
   lines(new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),1], m.tmp, lwd = 2, col = rainbow(2)[ii]) 
  }
legend(4.5, -0.05, bty = "n", legend = c("1, n = 60", "1, n = 120", "2, n = 60", "2, n = 120"), 
   col = c(rainbow(2)[1],rainbow(2)[1],rainbow(2)[2],rainbow(2)[2]), lty = c(2,1,2,1))
## dependence of a on n:
plot(pars.randompos.dfr$n, pars.randompos.dfr$loga.single, pch = 16,
   cex = 0.7, col = ctr.tmp, main = "Random pos., log(a) ~ n")
for(ii in 1:2)
  {
   m.tmp <- predict(lm.loga.randompos.init3[[ii]], newdata = new.tmp2[1:(nrow(new.tmp2)/2),])
   lines(new.tmp2[1:(nrow(new.tmp2)/2),2], m.tmp, lwd = 2, col = rainbow(2)[ii], lty = 2) 
   m.tmp <- predict(lm.loga.randompos.init3[[ii]], newdata = new.tmp2[(nrow(new.tmp2)/2+1):nrow(new.tmp2),])
   lines(new.tmp2[(nrow(new.tmp2)/2+1):nrow(new.tmp2),2], m.tmp, lwd = 2, col = rainbow(2)[ii]) 
  }
legend(160, -0.05, bty = "n", legend = c("1, sum = 3.2", "1, sum = 4.5", "2, sum = 3.2", "2, sum = 4.5"), 
   col = c(rainbow(2)[1],rainbow(2)[1],rainbow(2)[2],rainbow(2)[2]), lty = c(2,1,2,1))

plot(pars.line.dfr$highest6sum.spwin, pars.line.dfr$loga.single, pch = 16,
   cex = 0.7, col = ctr.tmp2, main = "Line, log(a) ~ highest6sum.spwin", ylim = c(-0.4,-0.03))
for(ii in 1:2)
  {
   m.tmp <- predict(lm.loga.line.init3[[ii]], newdata = new.tmp[1:(nrow(new.tmp)/2),])
   lines(new.tmp[1:(nrow(new.tmp)/2),1], m.tmp, lwd = 2, col = rainbow(2)[ii], lty = 2) 
   m.tmp <- predict(lm.loga.line.init3[[ii]], newdata = new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),])
   lines(new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),1], m.tmp, lwd = 2, col = rainbow(2)[ii]) 
  }
legend(4.5, -0.05, bty = "n", legend = c("1, n = 60", "1, n = 120", "2, n = 60", "2, n = 120"), 
   col = c(rainbow(2)[1],rainbow(2)[1],rainbow(2)[2],rainbow(2)[2]), lty = c(2,1,2,1))
## dependence of b on n:
plot(pars.line.dfr$n, pars.line.dfr$loga.single, pch = 16,
   cex = 0.7, col = ctr.tmp2, main = "Line, log(a) ~ n")
for(ii in 1:2)
  {
   m.tmp <- predict(lm.loga.line.init3[[ii]], newdata = new.tmp2[1:(nrow(new.tmp2)/2),])
   lines(new.tmp2[1:(nrow(new.tmp2)/2),2], m.tmp, lwd = 2, col = rainbow(2)[ii], lty = 2) 
   m.tmp <- predict(lm.loga.line.init3[[ii]], newdata = new.tmp2[(nrow(new.tmp2)/2+1):nrow(new.tmp2),])
   lines(new.tmp2[(nrow(new.tmp2)/2+1):nrow(new.tmp2),2], m.tmp, lwd = 2, col = rainbow(2)[ii]) 
  }
legend(160, -0.05, bty = "n", legend = c("1, sum = 3.2", "1, sum = 4.5", "2, sum = 3.2", "2, sum = 4.5"), 
   col = c(rainbow(2)[1],rainbow(2)[1],rainbow(2)[2],rainbow(2)[2]), lty = c(2,1,2,1))
mtext("(1,3)/(4,5,6,7)", outer = TRUE)


## Parameter b:

quartz(height = 7.5, width = 10)
#dev.set(2)
par(mfrow = c(2,2), mar = c(2,2,2,1), mgp = c(1.8,0.9,0), oma = c(0,0,2,0))
## dependence of a on highest6sum.spwin:
plot(pars.randompos.dfr$highest6sum.spwin, pars.randompos.dfr$logb.single, pch = 16,
   cex = 0.7, col = ctr.tmp, main = "Random pos., log(b) ~ highest6sum.spwin", ylim = c(-0.6,1))
for(ii in 1:2)
  {
   m.tmp <- predict(lm.logb.randompos.init3[[ii]], newdata = new.tmp[1:(nrow(new.tmp)/2),])
   lines(new.tmp[1:(nrow(new.tmp)/2),1], m.tmp, lwd = 2, col = rainbow(2)[ii], lty = 2) 
   m.tmp <- predict(lm.logb.randompos.init3[[ii]], newdata = new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),])
   lines(new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),1], m.tmp, lwd = 2, col = rainbow(2)[ii]) 
  }
legend(4.5, 1, bty = "n", legend = c("1, n = 60", "1, n = 120", "2, n = 60", "2, n = 120"), 
   col = c(rainbow(2)[1],rainbow(2)[1],rainbow(2)[2],rainbow(2)[2]), lty = c(2,1,2,1))
## dependence of a on n:
plot(pars.randompos.dfr$n, pars.randompos.dfr$logb.single, pch = 16, ylim = c(-0.6,1),
   cex = 0.7, col = ctr.tmp, main = "Random pos., log(b) ~ n")
for(ii in 1:2)
  {
   m.tmp <- predict(lm.logb.randompos.init3[[ii]], newdata = new.tmp2[1:(nrow(new.tmp2)/2),])
   lines(new.tmp2[1:(nrow(new.tmp2)/2),2], m.tmp, lwd = 2, col = rainbow(2)[ii], lty = 2) 
   m.tmp <- predict(lm.logb.randompos.init3[[ii]], newdata = new.tmp2[(nrow(new.tmp2)/2+1):nrow(new.tmp2),])
   lines(new.tmp2[(nrow(new.tmp2)/2+1):nrow(new.tmp2),2], m.tmp, lwd = 2, col = rainbow(2)[ii]) 
  }
legend(160, -0.05, bty = "n", legend = c("1, sum = 3.2", "1, sum = 4.5", "2, sum = 3.2", "2, sum = 4.5"), 
   col = c(rainbow(2)[1],rainbow(2)[1],rainbow(2)[2],rainbow(2)[2]), lty = c(2,1,2,1))
plot(pars.line.dfr$highest6sum.spwin, pars.line.dfr$logb.single, pch = 16,
   cex = 0.7, col = ctr.tmp2, main = "Line, log(b) ~ highest6sum.spwin", ylim = c(-0.6,1))
for(ii in 1:2)
  {
   m.tmp <- predict(lm.logb.line.init3[[ii]], newdata = new.tmp[1:(nrow(new.tmp)/2),])
   lines(new.tmp[1:(nrow(new.tmp)/2),1], m.tmp, lwd = 2, col = rainbow(2)[ii], lty = 2) 
   m.tmp <- predict(lm.logb.line.init3[[ii]], newdata = new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),])
   lines(new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),1], m.tmp, lwd = 2, col = rainbow(2)[ii]) 
  }
legend(4.5, -0.05, bty = "n", legend = c("1, n = 60", "1, n = 120", "2, n = 60", "2, n = 120"), 
   col = c(rainbow(2)[1],rainbow(2)[1],rainbow(2)[2],rainbow(2)[2]), lty = c(2,1,2,1))
## dependence of b on n:
plot(pars.line.dfr$n, pars.line.dfr$logb.single, pch = 16,
   cex = 0.7, col = ctr.tmp2, main = "Line, log(b) ~ n")
for(ii in 1:2)
  {
   m.tmp <- predict(lm.logb.line.init3[[ii]], newdata = new.tmp2[1:(nrow(new.tmp2)/2),])
   lines(new.tmp2[1:(nrow(new.tmp2)/2),2], m.tmp, lwd = 2, col = rainbow(2)[ii], lty = 2) 
   m.tmp <- predict(lm.logb.line.init3[[ii]], newdata = new.tmp2[(nrow(new.tmp2)/2+1):nrow(new.tmp2),])
   lines(new.tmp2[(nrow(new.tmp2)/2+1):nrow(new.tmp2),2], m.tmp, lwd = 2, col = rainbow(2)[ii]) 
  }
legend(160, -0.05, bty = "n", legend = c("1, sum = 3.2", "1, sum = 4.5", "2, sum = 3.2", "2, sum = 4.5"), 
   col = c(rainbow(2)[1],rainbow(2)[1],rainbow(2)[2],rainbow(2)[2]), lty = c(2,1,2,1))
mtext("(1,3)/(4,5,6,7)", outer = TRUE)


pars.init.line3 <- matrix(c(lm.loga.line.init3[[1]]$coef,
    lm.logb.line.init3[[1]]$coef,
    lm.loga.line.init3[[2]]$coef,
    lm.logb.line.init3[[2]]$coef), 
    ncol = 6, byrow = TRUE)
pars.init.line3

pars.init.randompos3 <- matrix(c(lm.loga.randompos.init3[[1]]$coef,
    lm.logb.randompos.init3[[1]]$coef,
    lm.loga.randompos.init3[[2]]$coef,
    lm.logb.randompos.init3[[2]]$coef), 
    ncol = 6, byrow = TRUE)
pars.init.randompos3


## --------------------------------------------------------------------------------------------------------
## Optimization of the beta likelihood
## --------------------------------------------------------------------------------------------------------

pars.randompos3 <- list()
system.time(
pars.randompos3[[1]] <- optim(pars.init.randompos3[1,], fn = wloglik.beta.logpar.fun, gr = wloglikgrad.logpar.fun,
   vars.list = vars.list.randompos[pars.randompos.dfr$fine.class.for.init3 == 1], 
   weights = rep(1, sum(pars.randompos.dfr$fine.class.for.init3 == 1)),
   control = list(maxit = 10000, fnscale = -1))
  )
system.time(
pars.randompos3[[2]] <- optim(pars.init.randompos3[2,], fn = wloglik.beta.logpar.fun, gr = wloglikgrad.logpar.fun,
   vars.list = vars.list.randompos[pars.randompos.dfr$fine.class.for.init3 == 2], 
   weights = rep(1, sum(pars.randompos.dfr$fine.class.for.init3 == 2)),
   control = list(maxit = 10000, fnscale = -1))
  )
 
pars.line3 <- list()
system.time(
pars.line3[[1]] <- optim(pars.init.line3[1,], fn = wloglik.beta.logpar.fun, gr = wloglikgrad.logpar.fun,
   vars.list = vars.list.line[pars.line.dfr$fine.class.for.init3 == 1], 
   weights = rep(1, sum(pars.line.dfr$fine.class.for.init3 == 1)),
   control = list(maxit = 10000, fnscale = -1))
  )
system.time(
pars.line3[[2]] <- optim(pars.init.line3[2,], fn = wloglik.beta.logpar.fun, gr = wloglikgrad.logpar.fun,
   vars.list = vars.list.line[pars.line.dfr$fine.class.for.init3 == 2], 
   weights = rep(1, sum(pars.line.dfr$fine.class.for.init3 == 2)),
   control = list(maxit = 10000, fnscale = -1))
  )
  
## --------------------------------------------------------------------------------------------------------
## Plot of the points and the lines  
## --------------------------------------------------------------------------------------------------------

c.tmp <- rgb(blue = (as.numeric(pars.randompos.dfr$fine.class.for.init3)-1), 
   red = (2-as.numeric(pars.randompos.dfr$fine.class.for.init3)),
   green = (1-(pars.randompos.dfr$n-40) / 180)^2, alpha = 1)
c.tmp2 <- rgb(blue = (as.numeric(pars.line.dfr$fine.class.for.init3)-1), 
   red = (2-as.numeric(pars.line.dfr$fine.class.for.init3)),
   green = (1-(pars.line.dfr$n-40) / 180)^2, alpha = 1)
new.tmp <- data.frame(highest6sum.spwin = rep(seq(2.5, 5.7, by = 0.005), 2), 
   n = rep(c(60,120), each = length(seq(2.5, 5.7, by = 0.005))))
csmalln.tmp <- c(rgb(blue = 0, red = 1, green = (1-(60-40) / 180)^2, alpha = 1),
    rgb(blue = 1, red = 0, green = (1-(60-40) / 180)^2, alpha = 1))
clargen.tmp <- c(rgb(blue = 0, red = 1, green = (1-(120-40) / 180)^2, alpha = 1),
    rgb(blue = 1, red = 0, green = (1-(120-40) / 180)^2, alpha = 1))

quartz(height = 7.5, width = 10)
#dev.set(2)
par(mfrow = c(2,2), mar = c(2,2,3,1), mgp = c(1.8,0.9,0), oma = c(0,0,2,0))
## dependence of a on highest6sum.spwin:
plot(pars.randompos.dfr$highest6sum.spwin, pars.randompos.dfr$logb.single, pch = 16,
   cex = 0.7, col = c.tmp, main = "Random positions \n log(b) ~ highest6sum.spwin", 
   ylim = c(-0.6,1), xlab = "", ylab = "")
for(ii in 1:2)
  {
   m.tmp <- predict(lm.logb.randompos.init3[[ii]], newdata = new.tmp[1:(nrow(new.tmp)/2),])
   lines(new.tmp[1:(nrow(new.tmp)/2),1], m.tmp, lwd = 2, col = csmalln.tmp[ii], lty = 2) 
   m.tmp <- predict(lm.logb.randompos.init3[[ii]], newdata = new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),])
   lines(new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),1], m.tmp, lwd = 2, col = clargen.tmp[ii]) 
  }
legend(4.5, 1, bty = "n", legend = c("1, n = 60", "1, n = 120", "2, n = 60", "2, n = 120"), 
   col = c(csmalln.tmp[1],clargen[1], csmalln.tmp[2],clargen[2]), lty = c(2,1,2,1))

plot(pars.line.dfr$highest6sum.spwin, pars.line.dfr$logb.single, pch = 16,
   cex = 0.7, col = c.tmp2, main = "Line \n log(b) ~ highest6sum.spwin", 
   ylim = c(-0.6,1), xlab = "", ylab = "")
for(ii in 1:2)
  {
   m.tmp <- predict(lm.logb.line.init3[[ii]], newdata = new.tmp[1:(nrow(new.tmp)/2),])
   lines(new.tmp[1:(nrow(new.tmp)/2),1], m.tmp, lwd = 2, col = csmalln.tmp[ii], lty = 2) 
   m.tmp <- predict(lm.logb.line.init3[[ii]], newdata = new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),])
   lines(new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),1], m.tmp, lwd = 2, col = clargen.tmp[ii]) 
  }
legend(4.5, 1, bty = "n", legend = c("1, n = 60", "1, n = 120", "2, n = 60", "2, n = 120"), 
   col = c(csmalln.tmp[1],clargen[1], csmalln.tmp[2],clargen[2]), lty = c(2,1,2,1))

plot(pars.randompos.dfr$highest6sum.spwin, pars.randompos.dfr$loga.single, pch = 16,
   cex = 0.7, col = c.tmp, main = "\n log(a) ~ highest6sum.spwin", ylim = c(-0.4,-0.03), 
   xlab = "", ylab = "")
for(ii in 1:2)
  {
   m.tmp <- predict(lm.loga.randompos.init3[[ii]], newdata = new.tmp[1:(nrow(new.tmp)/2),])
   lines(new.tmp[1:(nrow(new.tmp)/2),1], m.tmp, lwd = 2, col = csmalln.tmp[ii], lty = 2) 
   m.tmp <- predict(lm.loga.randompos.init3[[ii]], newdata = new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),])
   lines(new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),1], m.tmp, lwd = 2, col = clargen.tmp[ii]) 
  }
legend(4.5, -0.05, bty = "n", legend = c("1, n = 60", "1, n = 120", "2, n = 60", "2, n = 120"), 
   col = c(csmalln.tmp[1],clargen[1], csmalln.tmp[2],clargen[2]), lty = c(2,1,2,1))

plot(pars.line.dfr$highest6sum.spwin, pars.line.dfr$loga.single, pch = 16,
   cex = 0.7, col = c.tmp2, main = "\n log(a) ~ highest6sum.spwin", ylim = c(-0.4,-0.03), 
   xlab = "", ylab = "")
for(ii in 1:2)
  {
   m.tmp <- predict(lm.loga.line.init3[[ii]], newdata = new.tmp[1:(nrow(new.tmp)/2),])
   lines(new.tmp[1:(nrow(new.tmp)/2),1], m.tmp, lwd = 2, col = csmalln.tmp[ii], lty = 2) 
   m.tmp <- predict(lm.loga.line.init3[[ii]], newdata = new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),])
   lines(new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),1], m.tmp, lwd = 2, col = clargen.tmp[ii]) 
  }
legend(4.5, -0.05, bty = "n", legend = c("1, n = 60", "1, n = 120", "2, n = 60", "2, n = 120"), 
   col = c(csmalln.tmp[1],clargen[1], csmalln.tmp[2],clargen[2]), lty = c(2,1,2,1))
mtext("(1,3)/(4,5,6,7)", outer = TRUE, font = 2)


## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------
## Grouping (1,2,3)/(4,5,6,7)
## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------

## Optimize wloglik.beta.logpar.fun from the EMbetafunctions.R source file,
## groupwise, and using equal weights (1 for all locations belonging to
## these groups).

## --------------------------------------------------------------------------------------------------------
## Merging of classes

## Put the new merged classes into the dataframe:

pars.randompos.dfr$fine.class.for.init4[
   is.element(pars.randompos.dfr$class7000_7group, c(1,2,3))] <- 2
pars.randompos.dfr$fine.class.for.init4[
   is.element(pars.randompos.dfr$class7000_7group, c(4,5,6,7))] <- 1
 
pars.line.dfr$fine.class.for.init4[
   is.element(pars.line.dfr$class7000_7group, c(1,2,3))] <- 2
pars.line.dfr$fine.class.for.init4[
   is.element(pars.line.dfr$class7000_7group, c(4,5,6,7))] <- 1
  
table(pars.randompos.dfr$fine.class.for.init4)
table(pars.line.dfr$fine.class.for.init4)


## --------------------------------------------------------------------------------------------------------
## Initial values in the format needed for the log-likelihood function
## --------------------------------------------------------------------------------------------------------

lm.loga.randompos.init4 <- vector(2, mode = "list")
for(ii in 1:2)
   {
   	lm.loga.randompos.init4[[ii]] <- lm(loga.single ~ highest6sum.spwin + n,
   	   data = pars.randompos.dfr, 
   	   subset = pars.randompos.dfr$fine.class.for.init4 == factor(ii, levels = c(1,2)))
   }  
lm.logb.randompos.init4 <- vector(2, mode = "list")
for(ii in 1:2)
   {
   	lm.logb.randompos.init4[[ii]] <- lm(logb.single ~ highest6sum.spwin + n,
   	   data = pars.randompos.dfr, 
   	   subset = pars.randompos.dfr$fine.class.for.init4 == factor(ii, levels = c(1,2)))
   }
lapply(lm.loga.randompos.init4, summary) 
lapply(lm.logb.randompos.init4, summary) 

lm.loga.line.init4 <- vector(2, mode = "list")
for(ii in 1:2)
   {
   	lm.loga.line.init4[[ii]] <- lm(loga.single ~ highest6sum.spwin + n,
   	   data = pars.line.dfr, 
   	   subset = pars.line.dfr$fine.class.for.init4 == factor(ii, levels = c(1,2)))
   }  
lm.logb.line.init4 <- vector(2, mode = "list")
for(ii in 1:2)
   {
   	lm.logb.line.init4[[ii]] <- lm(logb.single ~ highest6sum.spwin + n,
   	   data = pars.line.dfr, 
   	   subset = pars.line.dfr$fine.class.for.init4 == factor(ii, levels = c(1,2)))
   }    
lapply(lm.loga.line.init4, summary) 
lapply(lm.logb.line.init4, summary) 
  

c.tmp <- rainbow(2)[as.numeric(pars.randompos.dfr$fine.class.for.init4)]
ctr.tmp <- rainbow(2, alpha = 0.6)[as.numeric(pars.randompos.dfr$fine.class.for.init4)]
c.tmp2 <- rainbow(2)[as.numeric(pars.line.dfr$fine.class.for.init4)]
ctr.tmp2 <- rainbow(2, alpha = 0.6)[as.numeric(pars.line.dfr$fine.class.for.init4)]
new.tmp <- data.frame(highest6sum.spwin = rep(seq(2.5, 5.7, by = 0.005), 2), 
   n = rep(c(60,120), each = length(seq(2.5, 5.7, by = 0.005))))
new.tmp2 <- data.frame(highest6sum.spwin = rep(c(3.2,4.5), each = 181), 
   n = rep(40:220, 2))

quartz(height = 7.5, width = 10)
#dev.set(2)
par(mfrow = c(2,2), mar = c(2,2,2,1), mgp = c(1.8,0.9,0), oma = c(0,0,2,0))
## dependence of a on highest6sum.spwin:
plot(pars.randompos.dfr$highest6sum.spwin, pars.randompos.dfr$loga.single, pch = 16,
   cex = 0.7, col = ctr.tmp, main = "Random pos., log(a) ~ highest6sum.spwin", ylim = c(-0.4,-0.03))
for(ii in 1:2)
  {
   m.tmp <- predict(lm.loga.randompos.init4[[ii]], newdata = new.tmp[1:(nrow(new.tmp)/2),])
   lines(new.tmp[1:(nrow(new.tmp)/2),1], m.tmp, lwd = 2, col = rainbow(2)[ii], lty = 2) 
   m.tmp <- predict(lm.loga.randompos.init4[[ii]], newdata = new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),])
   lines(new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),1], m.tmp, lwd = 2, col = rainbow(2)[ii]) 
  }
legend(4.5, -0.05, bty = "n", legend = c("1, n = 60", "1, n = 120", "2, n = 60", "2, n = 120"), 
   col = c(rainbow(2)[1],rainbow(2)[1],rainbow(2)[2],rainbow(2)[2]), lty = c(2,1,2,1))
## dependence of a on n:
plot(pars.randompos.dfr$n, pars.randompos.dfr$loga.single, pch = 16,
   cex = 0.7, col = ctr.tmp, main = "Random pos., log(a) ~ n")
for(ii in 1:2)
  {
   m.tmp <- predict(lm.loga.randompos.init4[[ii]], newdata = new.tmp2[1:(nrow(new.tmp2)/2),])
   lines(new.tmp2[1:(nrow(new.tmp2)/2),2], m.tmp, lwd = 2, col = rainbow(2)[ii], lty = 2) 
   m.tmp <- predict(lm.loga.randompos.init4[[ii]], newdata = new.tmp2[(nrow(new.tmp2)/2+1):nrow(new.tmp2),])
   lines(new.tmp2[(nrow(new.tmp2)/2+1):nrow(new.tmp2),2], m.tmp, lwd = 2, col = rainbow(2)[ii]) 
  }
legend(160, -0.05, bty = "n", legend = c("1, sum = 3.2", "1, sum = 4.5", "2, sum = 3.2", "2, sum = 4.5"), 
   col = c(rainbow(2)[1],rainbow(2)[1],rainbow(2)[2],rainbow(2)[2]), lty = c(2,1,2,1))

plot(pars.line.dfr$highest6sum.spwin, pars.line.dfr$loga.single, pch = 16,
   cex = 0.7, col = ctr.tmp2, main = "Line, log(a) ~ highest6sum.spwin", ylim = c(-0.4,-0.03))
for(ii in 1:2)
  {
   m.tmp <- predict(lm.loga.line.init4[[ii]], newdata = new.tmp[1:(nrow(new.tmp)/2),])
   lines(new.tmp[1:(nrow(new.tmp)/2),1], m.tmp, lwd = 2, col = rainbow(2)[ii], lty = 2) 
   m.tmp <- predict(lm.loga.line.init4[[ii]], newdata = new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),])
   lines(new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),1], m.tmp, lwd = 2, col = rainbow(2)[ii]) 
  }
legend(4.5, -0.05, bty = "n", legend = c("1, n = 60", "1, n = 120", "2, n = 60", "2, n = 120"), 
   col = c(rainbow(2)[1],rainbow(2)[1],rainbow(2)[2],rainbow(2)[2]), lty = c(2,1,2,1))
## dependence of b on n:
plot(pars.line.dfr$n, pars.line.dfr$loga.single, pch = 16,
   cex = 0.7, col = ctr.tmp2, main = "Line, log(a) ~ n")
for(ii in 1:2)
  {
   m.tmp <- predict(lm.loga.line.init4[[ii]], newdata = new.tmp2[1:(nrow(new.tmp2)/2),])
   lines(new.tmp2[1:(nrow(new.tmp2)/2),2], m.tmp, lwd = 2, col = rainbow(2)[ii], lty = 2) 
   m.tmp <- predict(lm.loga.line.init4[[ii]], newdata = new.tmp2[(nrow(new.tmp2)/2+1):nrow(new.tmp2),])
   lines(new.tmp2[(nrow(new.tmp2)/2+1):nrow(new.tmp2),2], m.tmp, lwd = 2, col = rainbow(2)[ii]) 
  }
legend(160, -0.05, bty = "n", legend = c("1, sum = 3.2", "1, sum = 4.5", "2, sum = 3.2", "2, sum = 4.5"), 
   col = c(rainbow(2)[1],rainbow(2)[1],rainbow(2)[2],rainbow(2)[2]), lty = c(2,1,2,1))
mtext("(1,2,3)/(4,5,6,7)", outer = TRUE)

quartz(height = 7.5, width = 10)
#dev.set(2)
par(mfrow = c(2,2), mar = c(2,2,2,1), mgp = c(1.8,0.9,0), oma = c(0,0,2,0))
## dependence of a on highest6sum.spwin:
plot(pars.randompos.dfr$highest6sum.spwin, pars.randompos.dfr$logb.single, pch = 16,
   cex = 0.7, col = ctr.tmp, main = "Random pos., log(b) ~ highest6sum.spwin", ylim = c(-0.6,1))
for(ii in 1:2)
  {
   m.tmp <- predict(lm.logb.randompos.init4[[ii]], newdata = new.tmp[1:(nrow(new.tmp)/2),])
   lines(new.tmp[1:(nrow(new.tmp)/2),1], m.tmp, lwd = 2, col = rainbow(2)[ii], lty = 2) 
   m.tmp <- predict(lm.logb.randompos.init4[[ii]], newdata = new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),])
   lines(new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),1], m.tmp, lwd = 2, col = rainbow(2)[ii]) 
  }
legend(4.5, 1, bty = "n", legend = c("1, n = 60", "1, n = 120", "2, n = 60", "2, n = 120"), 
   col = c(rainbow(2)[1],rainbow(2)[1],rainbow(2)[2],rainbow(2)[2]), lty = c(2,1,2,1))
## dependence of a on n:
plot(pars.randompos.dfr$n, pars.randompos.dfr$logb.single, pch = 16, ylim = c(-0.6,1),
   cex = 0.7, col = ctr.tmp, main = "Random pos., log(b) ~ n")
for(ii in 1:2)
  {
   m.tmp <- predict(lm.logb.randompos.init4[[ii]], newdata = new.tmp2[1:(nrow(new.tmp2)/2),])
   lines(new.tmp2[1:(nrow(new.tmp2)/2),2], m.tmp, lwd = 2, col = rainbow(2)[ii], lty = 2) 
   m.tmp <- predict(lm.logb.randompos.init4[[ii]], newdata = new.tmp2[(nrow(new.tmp2)/2+1):nrow(new.tmp2),])
   lines(new.tmp2[(nrow(new.tmp2)/2+1):nrow(new.tmp2),2], m.tmp, lwd = 2, col = rainbow(2)[ii]) 
  }
legend(160, -0.05, bty = "n", legend = c("1, sum = 3.2", "1, sum = 4.5", "2, sum = 3.2", "2, sum = 4.5"), 
   col = c(rainbow(2)[1],rainbow(2)[1],rainbow(2)[2],rainbow(2)[2]), lty = c(2,1,2,1))

plot(pars.line.dfr$highest6sum.spwin, pars.line.dfr$logb.single, pch = 16,
   cex = 0.7, col = ctr.tmp2, main = "Line, log(b) ~ highest6sum.spwin", ylim = c(-0.6,1))
for(ii in 1:2)
  {
   m.tmp <- predict(lm.logb.line.init4[[ii]], newdata = new.tmp[1:(nrow(new.tmp)/2),])
   lines(new.tmp[1:(nrow(new.tmp)/2),1], m.tmp, lwd = 2, col = rainbow(2)[ii], lty = 2) 
   m.tmp <- predict(lm.logb.line.init4[[ii]], newdata = new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),])
   lines(new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),1], m.tmp, lwd = 2, col = rainbow(2)[ii]) 
  }
legend(4.5, 1, bty = "n", legend = c("1, n = 60", "1, n = 120", "2, n = 60", "2, n = 120"), 
   col = c(rainbow(2)[1],rainbow(2)[1],rainbow(2)[2],rainbow(2)[2]), lty = c(2,1,2,1))
## dependence of b on n:
plot(pars.line.dfr$n, pars.line.dfr$logb.single, pch = 16,
   cex = 0.7, col = ctr.tmp2, main = "Line, log(b) ~ n")
for(ii in 1:2)
  {
   m.tmp <- predict(lm.logb.line.init4[[ii]], newdata = new.tmp2[1:(nrow(new.tmp2)/2),])
   lines(new.tmp2[1:(nrow(new.tmp2)/2),2], m.tmp, lwd = 2, col = rainbow(2)[ii], lty = 2) 
   m.tmp <- predict(lm.logb.line.init4[[ii]], newdata = new.tmp2[(nrow(new.tmp2)/2+1):nrow(new.tmp2),])
   lines(new.tmp2[(nrow(new.tmp2)/2+1):nrow(new.tmp2),2], m.tmp, lwd = 2, col = rainbow(2)[ii]) 
  }
legend(160, -0.05, bty = "n", legend = c("1, sum = 3.2", "1, sum = 4.5", "2, sum = 3.2", "2, sum = 4.5"), 
   col = c(rainbow(2)[1],rainbow(2)[1],rainbow(2)[2],rainbow(2)[2]), lty = c(2,1,2,1))
mtext("(1,2,3)/(4,5,6,7)", outer = TRUE)


pars.init.line4 <- matrix(c(lm.loga.line.init4[[1]]$coef,
    lm.logb.line.init4[[1]]$coef,
    lm.loga.line.init4[[2]]$coef,
    lm.logb.line.init4[[2]]$coef), 
    ncol = 6, byrow = TRUE)
pars.init.line4

pars.init.randompos4 <- matrix(c(lm.loga.randompos.init4[[1]]$coef,
    lm.logb.randompos.init4[[1]]$coef,
    lm.loga.randompos.init4[[2]]$coef,
    lm.logb.randompos.init4[[2]]$coef), 
    ncol = 6, byrow = TRUE)
pars.init.randompos4


## --------------------------------------------------------------------------------------------------------
## Optimization of the beta likelihood
## --------------------------------------------------------------------------------------------------------

pars.randompos4 <- list()
system.time(
pars.randompos4[[1]] <- optim(pars.init.randompos4[1,], fn = wloglik.beta.logpar.fun, gr = wloglikgrad.logpar.fun,
   vars.list = vars.list.randompos[pars.randompos.dfr$fine.class.for.init4 == 1], 
   weights = rep(1, sum(pars.randompos.dfr$fine.class.for.init4 == 1)),
   control = list(maxit = 10000, fnscale = -1))
  )
system.time(
pars.randompos4[[2]] <- optim(pars.init.randompos4[2,], fn = wloglik.beta.logpar.fun, gr = wloglikgrad.logpar.fun,
   vars.list = vars.list.randompos[pars.randompos.dfr$fine.class.for.init4 == 2], 
   weights = rep(1, sum(pars.randompos.dfr$fine.class.for.init4 == 2)),
   control = list(maxit = 10000, fnscale = -1))
  )
 
pars.line4 <- list()
system.time(
pars.line4[[1]] <- optim(pars.init.line4[1,], fn = wloglik.beta.logpar.fun, gr = wloglikgrad.logpar.fun,
   vars.list = vars.list.line[pars.line.dfr$fine.class.for.init4 == 1], 
   weights = rep(1, sum(pars.line.dfr$fine.class.for.init4 == 1)),
   control = list(maxit = 10000, fnscale = -1))
  )
system.time(
pars.line4[[2]] <- optim(pars.init.line4[2,], fn = wloglik.beta.logpar.fun, gr = wloglikgrad.logpar.fun,
   vars.list = vars.list.line[pars.line.dfr$fine.class.for.init4 == 2], 
   weights = rep(1, sum(pars.line.dfr$fine.class.for.init4 == 2)),
   control = list(maxit = 10000, fnscale = -1))
  )


## --------------------------------------------------------------------------------------------------------
## Comparison of likelihoods
## --------------------------------------------------------------------------------------------------------

## Look at the total likelihoods in all models.
## (no penalization for using here more parameters: this is 2 more, due to one more slope in
## both a and b, than in modelfitting1 and *2)

pars.randompos2[[1]]$value + pars.randompos2[[2]]$value
## 104834.5
pars.randompos3[[1]]$value + pars.randompos3[[2]]$value
## 104761.2
pars.randompos4[[1]]$value + pars.randompos4[[2]]$value
## 104843.1
## Of these, pars4 with grouping (1,2,3)/(4,5,6,7) seems the best.

pars.line2[[1]]$value + pars.line2[[2]]$value
## 142755.5
pars.line3[[1]]$value + pars.line3[[2]]$value
## 142685.7
pars.line4[[1]]$value + pars.line4[[2]]$value
## 142788.3
## Of these, pars4 with grouping (1,2,3)/(4,5,6,7) seems the best.

## All these are still very far from the EM algorithm results in modelfitting1-2.RData, but
## somewhat better than those in modelfitting2.RData.

## Old results of the EM algorithm
# # firsttrial.em$optimization[[1]]$value + firsttrial.em$optimization[[2]]$value
# ## 143234.1
# secondtrial.em$optimization[[1]]$value + secondtrial.em$optimization[[2]]$value
# ## 142707
# loglik.thirdtrial[[10]]
# ## 143141.8
# loglik.fourthtrial[[8]]
# ## 143060.4

## Old results from modelfitting2.RData, using the simple classification...
# pars1[[1]]$value + pars1[[2]]$value
# ## 142244.8
# pars2[[1]]$value + pars2[[2]]$value
# ## 142451.1
# pars3[[1]]$value + pars3[[2]]$value
# ## 142376.4
# pars4[[1]]$value + pars4[[2]]$value
# ## 142546.3
## And a cutpoint at highest6sum = 3.8:
# pars5[[1]]$value + pars5[[2]]$value
# ## 142502.2
## ... and the simple cuts:
# sapply(cutpoint.fits, function(lst) lst[[1]]$value + lst[[2]]$value)
     # # 3.7      3.8      3.9        4      4.1      4.2 
# # 142435.7 142502.2 142497.5 142498.0 142529.4 142446.2 
# sapply(cutfine.fits, function(lst) lst[[1]]$value + lst[[2]]$value)
# ##     4.01     4.02     4.03     4.04     4.05     4.06     4.07     4.08     4.09      4.1     4.11     4.12 
# ## 142501.6 142513.8 142504.9 142516.0 142515.3 142522.1 142528.7 142545.4 142542.7 142529.4 142520.9 142517.2 

  
## --------------------------------------------------------------------------------------------------------
## Plot of the points and the lines  
## --------------------------------------------------------------------------------------------------------

c.tmp <- rgb(blue = (as.numeric(pars.randompos.dfr$fine.class.for.init4)-1), 
   red = (2-as.numeric(pars.randompos.dfr$fine.class.for.init4)),
   green = (1-(pars.randompos.dfr$n-40) / 180)^2, alpha = 1)
c.tmp2 <- rgb(blue = (as.numeric(pars.line.dfr$fine.class.for.init4)-1), 
   red = (2-as.numeric(pars.line.dfr$fine.class.for.init4)),
   green = (1-(pars.line.dfr$n-40) / 180)^2, alpha = 1)
new.tmp <- data.frame(highest6sum.spwin = rep(seq(2.5, 5.7, by = 0.005), 2), 
   n = rep(c(60,120), each = length(seq(2.5, 5.7, by = 0.005))))
csmalln.tmp <- c(rgb(blue = 0, red = 1, green = (1-(60-40) / 180)^2, alpha = 1),
    rgb(blue = 1, red = 0, green = (1-(60-40) / 180)^2, alpha = 1))
clargen.tmp <- c(rgb(blue = 0, red = 1, green = (1-(120-40) / 180)^2, alpha = 1),
    rgb(blue = 1, red = 0, green = (1-(120-40) / 180)^2, alpha = 1))

quartz(height = 7.5, width = 10)
#dev.set(2)
par(mfrow = c(2,2), mar = c(2,2,3,1), mgp = c(1.8,0.9,0), oma = c(0,0,2,0))
## dependence of a on highest6sum.spwin:
plot(pars.randompos.dfr$highest6sum.spwin, pars.randompos.dfr$logb.single, pch = 16,
   cex = 0.7, col = c.tmp, main = "Random positions \n log(b) ~ highest6sum.spwin", 
   ylim = c(-0.6,1), xlab = "", ylab = "")
for(ii in 1:2)
  {
   m.tmp <- predict(lm.logb.randompos.init4[[ii]], newdata = new.tmp[1:(nrow(new.tmp)/2),])
   lines(new.tmp[1:(nrow(new.tmp)/2),1], m.tmp, lwd = 2, col = csmalln.tmp[ii], lty = 2) 
   m.tmp <- predict(lm.logb.randompos.init4[[ii]], newdata = new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),])
   lines(new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),1], m.tmp, lwd = 2, col = clargen.tmp[ii]) 
  }
legend(4.5, 1, bty = "n", legend = c("1, n = 60", "1, n = 120", "2, n = 60", "2, n = 120"), 
   col = c(csmalln.tmp[1],clargen[1], csmalln.tmp[2],clargen[2]), lty = c(2,1,2,1))

plot(pars.line.dfr$highest6sum.spwin, pars.line.dfr$logb.single, pch = 16,
   cex = 0.7, col = c.tmp2, main = "Line \n log(b) ~ highest6sum.spwin", 
   ylim = c(-0.6,1), xlab = "", ylab = "")
for(ii in 1:2)
  {
   m.tmp <- predict(lm.logb.line.init4[[ii]], newdata = new.tmp[1:(nrow(new.tmp)/2),])
   lines(new.tmp[1:(nrow(new.tmp)/2),1], m.tmp, lwd = 2, col = csmalln.tmp[ii], lty = 2) 
   m.tmp <- predict(lm.logb.line.init4[[ii]], newdata = new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),])
   lines(new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),1], m.tmp, lwd = 2, col = clargen.tmp[ii]) 
  }
legend(4.5, 1, bty = "n", legend = c("1, n = 60", "1, n = 120", "2, n = 60", "2, n = 120"), 
   col = c(csmalln.tmp[1],clargen[1], csmalln.tmp[2],clargen[2]), lty = c(2,1,2,1))

plot(pars.randompos.dfr$highest6sum.spwin, pars.randompos.dfr$loga.single, pch = 16,
   cex = 0.7, col = c.tmp, main = "\n log(a) ~ highest6sum.spwin", ylim = c(-0.4,-0.03), 
   xlab = "", ylab = "")
for(ii in 1:2)
  {
   m.tmp <- predict(lm.loga.randompos.init4[[ii]], newdata = new.tmp[1:(nrow(new.tmp)/2),])
   lines(new.tmp[1:(nrow(new.tmp)/2),1], m.tmp, lwd = 2, col = csmalln.tmp[ii], lty = 2) 
   m.tmp <- predict(lm.loga.randompos.init4[[ii]], newdata = new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),])
   lines(new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),1], m.tmp, lwd = 2, col = clargen.tmp[ii]) 
  }
legend(4.5, -0.05, bty = "n", legend = c("1, n = 60", "1, n = 120", "2, n = 60", "2, n = 120"), 
   col = c(csmalln.tmp[1],clargen[1], csmalln.tmp[2],clargen[2]), lty = c(2,1,2,1))

plot(pars.line.dfr$highest6sum.spwin, pars.line.dfr$loga.single, pch = 16,
   cex = 0.7, col = c.tmp2, main = "\n log(a) ~ highest6sum.spwin", ylim = c(-0.4,-0.03), 
   xlab = "", ylab = "")
for(ii in 1:2)
  {
   m.tmp <- predict(lm.loga.line.init4[[ii]], newdata = new.tmp[1:(nrow(new.tmp)/2),])
   lines(new.tmp[1:(nrow(new.tmp)/2),1], m.tmp, lwd = 2, col = csmalln.tmp[ii], lty = 2) 
   m.tmp <- predict(lm.loga.line.init4[[ii]], newdata = new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),])
   lines(new.tmp[(nrow(new.tmp)/2+1):nrow(new.tmp),1], m.tmp, lwd = 2, col = clargen.tmp[ii]) 
  }
legend(4.5, -0.05, bty = "n", legend = c("1, n = 60", "1, n = 120", "2, n = 60", "2, n = 120"), 
   col = c(csmalln.tmp[1],clargen[1], csmalln.tmp[2],clargen[2]), lty = c(2,1,2,1))
mtext("(1,2,3)/(4,5,6,7)", outer = TRUE, font = 2)


## --------------------------------------------------------------------------------------------------------
## Variances of the parameters  
## --------------------------------------------------------------------------------------------------------

## Forgot to put hessian = TRUE in the optimizations.

system.time(
pars.randompos2[[1]]$hessian <- hessian(func = wloglik.beta.logpar.fun, 
   x = pars.randompos2[[1]]$par,
   vars.list = vars.list.randompos[pars.randompos.dfr$fine.class.for.init2 == 1], 
   weights = rep(1, sum(pars.randompos.dfr$fine.class.for.init2 == 1)))
)
system.time(
pars.randompos2[[2]]$hessian <- hessian(func = wloglik.beta.logpar.fun, 
   x = pars.randompos2[[2]]$par,
   vars.list = vars.list.randompos[pars.randompos.dfr$fine.class.for.init2 == 2], 
   weights = rep(1, sum(pars.randompos.dfr$fine.class.for.init2 == 2)))
)

system.time(
pars.randompos3[[1]]$hessian <- hessian(func = wloglik.beta.logpar.fun, 
   x = pars.randompos3[[1]]$par,
   vars.list = vars.list.randompos[pars.randompos.dfr$fine.class.for.init3 == 1], 
   weights = rep(1, sum(pars.randompos.dfr$fine.class.for.init3 == 1)))
)
system.time(
pars.randompos3[[2]]$hessian <- hessian(func = wloglik.beta.logpar.fun, 
   x = pars.randompos3[[2]]$par,
   vars.list = vars.list.randompos[pars.randompos.dfr$fine.class.for.init3 == 2], 
   weights = rep(1, sum(pars.randompos.dfr$fine.class.for.init3 == 2)))
)

system.time(
pars.randompos4[[1]]$hessian <- hessian(func = wloglik.beta.logpar.fun, 
   x = pars.randompos4[[1]]$par,
   vars.list = vars.list.randompos[pars.randompos.dfr$fine.class.for.init4 == 1], 
   weights = rep(1, sum(pars.randompos.dfr$fine.class.for.init4 == 1)))
)
system.time(
pars.randompos4[[2]]$hessian <- hessian(func = wloglik.beta.logpar.fun, 
   x = pars.randompos4[[2]]$par,
   vars.list = vars.list.randompos[pars.randompos.dfr$fine.class.for.init4 == 2], 
   weights = rep(1, sum(pars.randompos.dfr$fine.class.for.init4 == 2)))
)

system.time(
pars.line2[[1]]$hessian <- hessian(func = wloglik.beta.logpar.fun, 
   x = pars.line2[[1]]$par,
   vars.list = vars.list.line[pars.line.dfr$fine.class.for.init2 == 1], 
   weights = rep(1, sum(pars.line.dfr$fine.class.for.init2 == 1)))
)
system.time(
pars.line2[[2]]$hessian <- hessian(func = wloglik.beta.logpar.fun, 
   x = pars.line2[[2]]$par,
   vars.list = vars.list.line[pars.line.dfr$fine.class.for.init2 == 2], 
   weights = rep(1, sum(pars.line.dfr$fine.class.for.init2 == 2)))
)

system.time(
pars.line3[[1]]$hessian <- hessian(func = wloglik.beta.logpar.fun, 
   x = pars.line3[[1]]$par,
   vars.list = vars.list.line[pars.line.dfr$fine.class.for.init3 == 1], 
   weights = rep(1, sum(pars.line.dfr$fine.class.for.init3 == 1)))
)
system.time(
pars.line3[[2]]$hessian <- hessian(func = wloglik.beta.logpar.fun, 
   x = pars.line3[[2]]$par,
   vars.list = vars.list.line[pars.line.dfr$fine.class.for.init3 == 2], 
   weights = rep(1, sum(pars.line.dfr$fine.class.for.init3 == 2)))
)

system.time(
pars.line4[[1]]$hessian <- hessian(func = wloglik.beta.logpar.fun, 
   x = pars.line4[[1]]$par,
   vars.list = vars.list.line[pars.line.dfr$fine.class.for.init4 == 1], 
   weights = rep(1, sum(pars.line.dfr$fine.class.for.init4 == 1)))
)
system.time(
pars.line4[[2]]$hessian <- hessian(func = wloglik.beta.logpar.fun, 
   x = pars.line4[[2]]$par,
   vars.list = vars.list.line[pars.line.dfr$fine.class.for.init4 == 2], 
   weights = rep(1, sum(pars.line.dfr$fine.class.for.init4 == 2)))
)




## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------
## Cut points
## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------

## Optimize wloglik.beta.logpar.fun from the EMbetafunctions.R source file,
## groupwise, and using equal weights (1 for all locations belonging to
## these groups).


## --------------------------------------------------------------------------------------------------------
## Find the best cut on a rough grid: 
## --------------------------------------------------------------------------------------------------------

cuts.rough <- c(3.7,3.8,3.9,4,4.1,4.2,4.3)
cutrough.line.fits <- vector(7, mode = "list")
names(cutrough.line.fits) <- cuts.rough
cutrough.randompos.fits <- vector(7, mode = "list")
names(cutrough.randompos.fits) <- cuts.rough
for(jj in 1:7)
   {
   	## first for the line:
	fine.class.for.init.tmp <- rep(1, 900)
	fine.class.for.init.tmp[pars.line.dfr$highest6sum.spwin >= cuts.rough[jj]] <- 2
	lm.a.tmp <- vector(2, mode = "list")
	for(ii in 1:2)
	   {
	   	lm.a.tmp[[ii]] <- lm(loga.single ~ highest6sum.spwin + n,
	   	   data = pars.line.dfr, subset = fine.class.for.init.tmp == ii)
	   }  
	lm.b.tmp <- vector(2, mode = "list")
	for(ii in 1:2)
	   {
	   	lm.b.tmp[[ii]] <- lm(logb.single ~ highest6sum.spwin + n,
	   	   data = pars.line.dfr, subset = fine.class.for.init.tmp == ii)
	   }    
	pars.init.tmp <- matrix(c(lm.a.tmp[[1]]$coef, lm.b.tmp[[1]]$coef,
	    lm.a.tmp[[2]]$coef, lm.b.tmp[[2]]$coef), ncol = 6, byrow = TRUE)
	sol.tmp <- list()
	s.tmp <- system.time(
	sol.tmp[[1]] <- optim(pars.init.tmp[1,], fn = wloglik.beta.logpar.fun, gr = wloglikgrad.logpar.fun,
	   vars.list = vars.list.line[fine.class.for.init.tmp == 1], 
	   weights = rep(1, sum(fine.class.for.init.tmp == 1)), method = "BFGS",
	   control = list(maxit = 10000, fnscale = -1), hessian = TRUE)
	  )[3]
	cat("cut ", cuts.rough[jj], ", group 1 of line ready; time: ", s.tmp, "\n")
	s.tmp <- system.time(
	sol.tmp[[2]] <- optim(pars.init.tmp[2,], fn = wloglik.beta.logpar.fun, gr = wloglikgrad.logpar.fun,
	   vars.list = vars.list.line[fine.class.for.init.tmp == 2], 
	   weights = rep(1, sum(fine.class.for.init.tmp == 2)), method = "BFGS",
	   control = list(maxit = 10000, fnscale = -1), hessian = TRUE)
	  )[3]
	cat("cut ", cuts.rough[jj], ", group 2 of line ready; time: ", s.tmp, "\n")
	cutrough.line.fits[[jj]] <- sol.tmp

   	## then for the random positions:
	fine.class.for.init.tmp <- rep(1, 720)
	fine.class.for.init.tmp[pars.randompos.dfr$highest6sum.spwin >= cuts.rough[jj]] <- 2
	lm.a.tmp <- vector(2, mode = "list")
	for(ii in 1:2)
	   {
	   	lm.a.tmp[[ii]] <- lm(loga.single ~ highest6sum.spwin + n,
	   	   data = pars.randompos.dfr, subset = fine.class.for.init.tmp == ii)
	   }  
	lm.b.tmp <- vector(2, mode = "list")
	for(ii in 1:2)
	   {
	   	lm.b.tmp[[ii]] <- lm(logb.single ~ highest6sum.spwin + n,
	   	   data = pars.randompos.dfr, subset = fine.class.for.init.tmp == ii)
	   }    
	pars.init.tmp <- matrix(c(lm.a.tmp[[1]]$coef, lm.b.tmp[[1]]$coef,
	    lm.a.tmp[[2]]$coef, lm.b.tmp[[2]]$coef), ncol = 6, byrow = TRUE)
	sol.tmp <- list()
	s.tmp <- system.time(
	sol.tmp[[1]] <- optim(pars.init.tmp[1,], fn = wloglik.beta.logpar.fun, gr = wloglikgrad.logpar.fun,
	   vars.list = vars.list.randompos[fine.class.for.init.tmp == 1], 
	   weights = rep(1, sum(fine.class.for.init.tmp == 1)),  method = "BFGS",
	   control = list(maxit = 10000, fnscale = -1), hessian = TRUE)
	  )[3]
	cat("cut ", cuts.rough[jj], ", group 1 of random positions ready; time: ", s.tmp, "\n")
	s.tmp <- system.time(
	sol.tmp[[2]] <- optim(pars.init.tmp[2,], fn = wloglik.beta.logpar.fun, gr = wloglikgrad.logpar.fun,
	   vars.list = vars.list.randompos[fine.class.for.init.tmp == 2], 
	   weights = rep(1, sum(fine.class.for.init.tmp == 2)), method = "BFGS",
	   control = list(maxit = 10000, fnscale = -1), hessian = TRUE)
	  )[3]
	cat("cut ", cuts.rough[jj], ", group 2 of random positions ready; time: ", s.tmp, "\n")
	cutrough.randompos.fits[[jj]] <- sol.tmp
   }

sapply(cutrough.line.fits, function(lst) lst[[1]]$value + lst[[2]]$value)
     # 3.7      3.8      3.9        4      4.1      4.2      4.3 
# 142757.4 142781.6 142802.8 142836.6 142831.6 142757.8 142748.5 
sapply(cutrough.randompos.fits, function(lst) lst[[1]]$value + lst[[2]]$value)
     # 3.7      3.8      3.9        4      4.1      4.2      4.3 
# 104810.9 104823.2 104898.6 104898.0 104908.2 104880.2 104800.9 
## Still worse than the unpredictable EM algoritm results, but better than anything
## else until now. The difference might be due to changes in the function, that I use
## now lbeta instead of dbeta, or the inclusion of n in the fit? Let's hope that
## it's the inclusion of n.

## Former results from this workspace:
# # pars.randompos2[[1]]$value + pars.randompos2[[2]]$value
# ## 104834.5
# pars.randompos3[[1]]$value + pars.randompos3[[2]]$value
# ## 104761.2
# pars.randompos4[[1]]$value + pars.randompos4[[2]]$value
# ## 104843.1
# ## Of these, pars4 with grouping (1,2,3)/(4,5,6,7) seems the best.

# pars.line2[[1]]$value + pars.line2[[2]]$value
# ## 142755.5
# pars.line3[[1]]$value + pars.line3[[2]]$value
# ## 142685.7
# pars.line4[[1]]$value + pars.line4[[2]]$value
# ## 142788.3
# ## Of these, pars4 with grouping (1,2,3)/(4,5,6,7) seems the best.

## All these are still far from the EM algorithm results in modelfitting1-2.RData, but
## somewhat better than those in modelfitting2.RData.

## Old results of the EM algorithm
# # firsttrial.em$optimization[[1]]$value + firsttrial.em$optimization[[2]]$value
# ## 143234.1
# secondtrial.em$optimization[[1]]$value + secondtrial.em$optimization[[2]]$value
# ## 142707
# loglik.thirdtrial[[10]]
# ## 143141.8
# loglik.fourthtrial[[8]]
# ## 143060.4

## Old results from modelfitting2.RData, using the simple classification...
# pars1[[1]]$value + pars1[[2]]$value
# ## 142244.8
# pars2[[1]]$value + pars2[[2]]$value
# ## 142451.1
# pars3[[1]]$value + pars3[[2]]$value
# ## 142376.4
# pars4[[1]]$value + pars4[[2]]$value
# ## 142546.3
## And a cutpoint at highest6sum = 3.8:
# pars5[[1]]$value + pars5[[2]]$value
# ## 142502.2
## ... and the simple cuts:
# sapply(cutpoint.fits, function(lst) lst[[1]]$value + lst[[2]]$value)
     # # 3.7      3.8      3.9        4      4.1      4.2 
# # 142435.7 142502.2 142497.5 142498.0 142529.4 142446.2 
# sapply(cutfine.fits, function(lst) lst[[1]]$value + lst[[2]]$value)
# ##     4.01     4.02     4.03     4.04     4.05     4.06     4.07     4.08     4.09      4.1     4.11     4.12 
# ## 142501.6 142513.8 142504.9 142516.0 142515.3 142522.1 142528.7 142545.4 142542.7 142529.4 142520.9 142517.2


## OK, SO A CUT AT 4 SEEMS TO BE OK

  
## --------------------------------------------------------------------------------------------------------
## Plot of the points and the lines  
## --------------------------------------------------------------------------------------------------------

cut.tmp <- 4
c.tmp <- rgb(blue = as.numeric(pars.randompos.dfr$highest6sum.spwin >= cut.tmp), 
   red = as.numeric(pars.randompos.dfr$highest6sum.spwin < cut.tmp),
   green = (1-(pars.randompos.dfr$n-40) / 180)^2, alpha = 1)
c.tmp2 <- rgb(blue = as.numeric(pars.line.dfr$highest6sum.spwin >= cut.tmp), 
   red = as.numeric(pars.line.dfr$highest6sum.spwin < cut.tmp),
   green = (1-(pars.line.dfr$n-40) / 180)^2, alpha = 1)
csmalln.tmp <- c(rgb(blue = 0, red = 1, green = (1-(60-40) / 180)^2, alpha = 1),
    rgb(blue = 1, red = 0, green = (1-(60-40) / 180)^2, alpha = 1))
clargen.tmp <- c(rgb(blue = 0, red = 1, green = (1-(120-40) / 180)^2, alpha = 1),
    rgb(blue = 1, red = 0, green = (1-(120-40) / 180)^2, alpha = 1))

quartz(height = 7.5, width = 10)
#dev.set(2)
par(mfrow = c(2,2), mar = c(2,2,3,1), mgp = c(1.8,0.9,0), oma = c(0,0,2,0))
## dependence of a on highest6sum.spwin:
plot(pars.randompos.dfr$highest6sum.spwin, pars.randompos.dfr$logb.single, pch = 16,
   cex = 0.7, col = c.tmp, main = "Random positions \n log(b) ~ highest6sum.spwin", 
   ylim = c(-0.6,1), xlab = "", ylab = "")
for(ii in 1:2)
  {
   abline(c(cutrough.randompos.fits[[4]][[ii]]$par[4] + 60*cutrough.randompos.fits[[4]][[ii]]$par[6],
      cutrough.randompos.fits[[4]][[ii]]$par[5]), col = csmalln.tmp[ii])
   abline(c(cutrough.randompos.fits[[4]][[ii]]$par[4] + 120*cutrough.randompos.fits[[4]][[ii]]$par[6],
      cutrough.randompos.fits[[4]][[ii]]$par[5]),  col = clargen.tmp[ii])
  }
legend(4.5, 1, bty = "n", legend = c("1, n = 60", "1, n = 120", "2, n = 60", "2, n = 120"), 
   col = c(csmalln.tmp[1],clargen.tmp[1], csmalln.tmp[2],clargen.tmp[2]), lty = c(2,1,2,1))

plot(pars.line.dfr$highest6sum.spwin, pars.line.dfr$logb.single, pch = 16,
   cex = 0.7, col = c.tmp2, main = "Line \n log(b) ~ highest6sum.spwin", 
   ylim = c(-0.6,1), xlab = "", ylab = "")
for(ii in 1:2)
  {
   abline(c(cutrough.line.fits[[4]][[ii]]$par[4] + 60*cutrough.line.fits[[4]][[ii]]$par[6],
      cutrough.line.fits[[4]][[ii]]$par[5]), col = csmalln.tmp[ii])
   abline(c(cutrough.line.fits[[4]][[ii]]$par[4] + 120*cutrough.line.fits[[4]][[ii]]$par[6],
      cutrough.line.fits[[4]][[ii]]$par[5]), col = clargen.tmp[ii])
  }
legend(4.5, 1, bty = "n", legend = c("1, n = 60", "1, n = 120", "2, n = 60", "2, n = 120"), 
   col = c(csmalln.tmp[1],clargen.tmp[1], csmalln.tmp[2],clargen.tmp[2]), lty = c(2,1,2,1))

plot(pars.randompos.dfr$highest6sum.spwin, pars.randompos.dfr$loga.single, pch = 16,
   cex = 0.7, col = c.tmp, main = "\n log(a) ~ highest6sum.spwin", ylim = c(-0.4,-0.03), 
   xlab = "", ylab = "")
for(ii in 1:2)
  {
   abline(c(cutrough.randompos.fits[[4]][[ii]]$par[1] + 60*cutrough.randompos.fits[[4]][[ii]]$par[3],
      cutrough.randompos.fits[[4]][[ii]]$par[2]), col = csmalln.tmp[ii])
   abline(c(cutrough.randompos.fits[[4]][[ii]]$par[1] + 120*cutrough.randompos.fits[[4]][[ii]]$par[3],
      cutrough.randompos.fits[[4]][[ii]]$par[2]),  col = clargen.tmp[ii])
  }
legend(4.5, -0.05, bty = "n", legend = c("1, n = 60", "1, n = 120", "2, n = 60", "2, n = 120"), 
   col = c(csmalln.tmp[1],clargen.tmp[1], csmalln.tmp[2],clargen.tmp[2]), lty = c(2,1,2,1))

plot(pars.line.dfr$highest6sum.spwin, pars.line.dfr$loga.single, pch = 16,
   cex = 0.7, col = c.tmp2, main = "\n log(a) ~ highest6sum.spwin", ylim = c(-0.4,-0.03), 
   xlab = "", ylab = "")
for(ii in 1:2)
  {
   abline(c(cutrough.line.fits[[4]][[ii]]$par[1] + 60*cutrough.line.fits[[4]][[ii]]$par[3],
      cutrough.line.fits[[4]][[ii]]$par[2]), col = csmalln.tmp[ii])
   abline(c(cutrough.line.fits[[4]][[ii]]$par[1] + 120*cutrough.line.fits[[4]][[ii]]$par[3],
      cutrough.line.fits[[4]][[ii]]$par[2]), col = clargen.tmp[ii])
  }
legend(4.5, -0.05, bty = "n", legend = c("1, n = 60", "1, n = 120", "2, n = 60", "2, n = 120"), 
   col = c(csmalln.tmp[1],clargen.tmp[1], csmalln.tmp[2],clargen.tmp[2]), lty = c(2,1,2,1))
mtext("Cutpoint at highest6sum.spwin = 4", outer = TRUE, font = 2)

  
## --------------------------------------------------------------------------------------------------------
## Errors on the parameters of the best model, CI of parameters 
## --------------------------------------------------------------------------------------------------------

round(cutrough.line.fits[[4]][[1]]$par - 1.96*sqrt(diag(-solve(cutrough.line.fits[[4]][[1]]$hessian))), 5)
round(cutrough.line.fits[[4]][[1]]$par, 5)
round(cutrough.line.fits[[4]][[1]]$par + 1.96*sqrt(diag(-solve(cutrough.line.fits[[4]][[1]]$hessian))), 5)
round(cutrough.line.fits[[4]][[2]]$par - 1.96*sqrt(diag(-solve(cutrough.line.fits[[4]][[2]]$hessian))), 5)
round(cutrough.line.fits[[4]][[2]]$par, 5)
round(cutrough.line.fits[[4]][[2]]$par + 1.96*sqrt(diag(-solve(cutrough.line.fits[[4]][[2]]$hessian))), 5)

round(cutrough.randompos.fits[[4]][[1]]$par - 1.96*sqrt(diag(-solve(cutrough.randompos.fits[[4]][[1]]$hessian))), 5)
round(cutrough.randompos.fits[[4]][[1]]$par, 5)
round(cutrough.randompos.fits[[4]][[1]]$par + 1.96*sqrt(diag(-solve(cutrough.randompos.fits[[4]][[1]]$hessian))), 5)
round(cutrough.randompos.fits[[4]][[2]]$par - 1.96*sqrt(diag(-solve(cutrough.randompos.fits[[4]][[2]]$hessian))), 5)
round(cutrough.randompos.fits[[4]][[2]]$par, 5)
round(cutrough.randompos.fits[[4]][[2]]$par + 1.96*sqrt(diag(-solve(cutrough.randompos.fits[[4]][[2]]$hessian))), 5)




################################################################################################
## CONCLUSION:
################################################################################################

## THE BEST MODEL FOR WHICH BASIC QUANTITIES ARE RELATIVELY EASY TO COMPUTE FROM THE RAW TIME
## SERIES IS THE CUTPOINT MODEL
##
## CUTPOINT: highest6sum.spwin = 4
##
## FOR MODELFITTING, USE RATHER THE LINE: IT CONTAINS ENOUGH OBSERVATIONS FROM A RARE GROUP
## OF VERY HIGH ALIASES


################################################################################################
## CONTINUATION:
################################################################################################

## 1. TRANSFORM THE PIPELINE FAPS TO HAVE P-VALUES (FIND THE TRANSFORMATION...)
##    1.1 store the parameters (a0,...,b2) of the relations 
##        alpha = exp(a0 + a1 * highest6sum.spwin + a2 * n)
##        beta = exp(b0 + b1 * highest6sum.spwin + b2 * n)
##        Then p_transformed = pbeta(fap; (shape 1, shape2) = (alpha, beta))
## 2. CHECK TWO MODELS ON THE RANDOMPOS AND ON THE OTHER HALF OF THE LINE SIMULATIONS, AND
##    CHECK THEM ALSO ON THE ECLIPTIC GRID:
##    (A) THE BEST CUTPOINT MODEL CONSTRUCTED ON HALF OF THE OBSERVATIONS ON THE LINE
##    (B) THE BEST CLASSIFICATION MODEL CONSTRUCTED ON HALF OF THE OBSERVATIONS ON THE LINE
##    2.1 Compute pars.*.dfr containing the spectral window peaks and n (and vart) for all 
##        locations ==> 
##        parameters of the linear dependence of alpha and beta of the beta distribution ==>
##        alpha and beta parameters of the beta distribution
##    1.2 p_transformed = pbeta(fap; (shape 1, shape2) = (alpha, beta))


###########################################################################################################
## --------------------------------------------------------------------------------------------------------
## Transformation of p-values
## --------------------------------------------------------------------------------------------------------
###########################################################################################################


## --------------------------------------------------------------------------------------------------------
## Storing the models
## --------------------------------------------------------------------------------------------------------

best.cut.line.model <- cutrough.line.fits[[4]]
best.class.line.model <- pars.line4
best.cut.randompos.model <- cutrough.randompos.fits[[4]]
best.class.randompos.model <- pars.randompos4

sqrt(diag(-solve(best.cut.line.model[[1]]$hessian)))
sqrt(diag(-solve(best.class.line.model[[1]]$hessian)))
sqrt(diag(-solve(best.cut.randompos.model[[1]]$hessian)))
sqrt(diag(-solve(best.class.randompos.model[[1]]$hessian)))

sqrt(diag(-solve(best.cut.line.model[[2]]$hessian)))
sqrt(diag(-solve(best.class.line.model[[2]]$hessian)))
sqrt(diag(-solve(best.cut.randompos.model[[2]]$hessian)))
sqrt(diag(-solve(best.class.randompos.model[[2]]$hessian)))

save(best.cut.line.model, best.class.line.model, file = "results_betamodelFitting1/best.line.models.RObj")
save(best.cut.randompos.model, 
   best.class.randompos.model, file = "results_betamodelFitting1/best.randompos.models.RObj")


## --------------------------------------------------------------------------------------------------------
## Transformation function
## --------------------------------------------------------------------------------------------------------

## ptrans.fun: computes the transformed fap given a pipeline fap, the sum of the highest 6 peaks in
##    the spectral window, the number of observations in the time series, and the name of the model
##    to apply
## Arguments:
##  tspars: a vector of length 2 or 3;
##       first: highest6sum, the sum of hte highest 6 peaks in the spectral window
##       second: n, number of observations
##       third: if exists, the class of the location
##  model.to.apply: which model to apply to generate the parameters of the beta distribution at the 
##    location of the pipeline fap. Possible values are best.cut.line.model, best.class.line.model,
##    best.cut.randompos.model, best.class.randompos.model; minimally, it must be a 
##    list of two components, each component having one element named par. For cut-based models,
##    the order of the classes should be: model.to.apply[[1]] for those with highest6sum < 4, 
##    and model.to.apply[[2]] for those with highest6sum >= 4; for classification, this comes from
##    a preliminary classification using the Gaussian mixture gmxt7000, and its classes (1,2,3)
##    are regrouped to 2, and (4,5,6,7), to 1.
model.to.apply <- best.cut.line.model
tspars <- vars.list.randompos[[1]][[2]]
betapars.fun <- function(tspars, model.to.apply)
   {
   	if(length(tspars) < 2 | length(tspars) > 3) {message("Parameters not the correct length \n"); break}
    highest6sum <- tspars[1]
    n <- tspars[2]
   	betapars <- if(length(tspars) == 2) {
   	   if(highest6sum < 4) 
   	   {
   	   	unlist(
   	   	c(exp(model.to.apply[[1]]$par[1] + model.to.apply[[1]]$par[2]*highest6sum + model.to.apply[[1]]$par[3]*n),
   	   	  exp(model.to.apply[[1]]$par[4] + model.to.apply[[1]]$par[5]*highest6sum + model.to.apply[[1]]$par[6]*n)))
   	   	} else {
   	   	unlist(
   	   	c(exp(model.to.apply[[2]]$par[1] + model.to.apply[[2]]$par[2]*highest6sum + model.to.apply[[2]]$par[3]*n),
   	   	  exp(model.to.apply[[2]]$par[4] + model.to.apply[[2]]$par[5]*highest6sum + model.to.apply[[2]]$par[6]*n)))
   	   	}
   	  } else {
   	   if(tspars[3] == 1) 
   	   {
   	   	unlist(
   	   	c(exp(model.to.apply[[1]]$par[1] + model.to.apply[[1]]$par[2]*highest6sum + model.to.apply[[1]]$par[3]*n),
   	   	  exp(model.to.apply[[1]]$par[4] + model.to.apply[[1]]$par[5]*highest6sum + model.to.apply[[1]]$par[6]*n)))
   	   	} else {
   	   	unlist(
   	   	c(exp(model.to.apply[[2]]$par[1] + model.to.apply[[2]]$par[2]*highest6sum + model.to.apply[[2]]$par[3]*n),
   	   	  exp(model.to.apply[[2]]$par[4] + model.to.apply[[2]]$par[5]*highest6sum + model.to.apply[[2]]$par[6]*n)))
   	   	}
   	  }
   	betapars
   }

## ptrans.fun: computes the transformed fap given a pipeline fap, the sum of the highest 6 peaks in
##    the spectral window, the number of observations in the time series, and the name of the model
##    to apply
ptrans.fun <- function(fap, betapars)
   {
   	pbeta(fap, shape1 = betapars[1], shape2 = betapars[2])
   }


## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------
## Lists containing FAPs and covariates for each location to check the models
## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------


## --------------------------------------------------------------------------------------------------------
## Line, the second 750 faps
## --------------------------------------------------------------------------------------------------------

## Take care to use now the second 750 of this!
vars.list.check.line <- lapply(pars.line.dfr$name, function(chr)
   {
   	f.tmp <- fap.line.dfr[751:1500, as.character(chr)]
   	f.tmp[f.tmp == 1] <- 1-1e-11
	f.tmp[f.tmp == 0] <- 1e-11
    list(fap = f.tmp, 
        covar = pars.line.dfr[pars.line.dfr$name == chr, c("highest6sum.spwin","n")])
   })
sapply(vars.list.check.line, function(lst) length(lst[[1]]))
vars.list.check.line[[431]]
names(vars.list.check.line) <- pars.line.dfr$name


## --------------------------------------------------------------------------------------------------------
## RandomPositions, the second 750 faps (in folder 
## /Users/mariasuveges/Documents/SpectralAnalysis/FAP/feasibilityForGaia/NoiseSimulations011113/NoiseSimulationsRandomMagnitude_RandomPositions_SecondSet/)
## --------------------------------------------------------------------------------------------------------

## Download the filenames from the second RandomPosition set, and check if they are the same
## as the first:
randomposcheckdirlist <- list.files(testdatadir.name)
identical(randomposcheckdirlist, randomposdirlist)
## All right, TRUE.
## Allow a way to check the number of observations, are they also the same?
n.checkset <- numeric(length(randomposcheckdirlist))


## the number of observations are in the commented-out header in the
## FAP files:

vars.list.check.line <- as.data.frame(matrix(ncol = length(randomposcheckdirlist), nrow = 750))
colnames(fap.randompos.check.dfr) <- randomposcheckdirlist
for(ii in 713:length(randomposcheckdirlist))
   {
#    ii <- 68
	dd <- randomposcheckdirlist[ii]
	pp <- list.files(paste0(testdatadir.name, dd), full.names = TRUE, pattern = ".dat")
    cc <- unlist(strsplit(readLines(pp, n = 2), " "))
    n.checkset[ii] <- as.numeric(cc[length(cc)])
	fap.randompos.check.dfr[,dd] <- read.table(pp)[,6]
   }
## error for ii = 68, 201, 333, 466, 592, 712:
# Error in `[<-.data.frame`(`*tmp*`, , dd, value = c(0.2161, 0.038837, 0.48784,  : 
  # replacement has **xxx** rows, data has 750
# In addition: Warning message:
# In scan(file, what, nmax, sep, dec, quote, skip, nlines, na.strings,  :
  # number of items read is not a multiple of the number of columns
## Either I have only fewer/0 faps here, or I leave out these locations from testing.
## Do the latter. 

## Check for the number of observations:
which(pars.randompos.dfr$n - n.checkset != 0)
## There are two 99s; these come from ii = 592 and 712, where the files seem to be empty,
## with no fap. It is strange that both should have 99 observations.

vars.list.check.randompos <- lapply(pars.randompos.dfr$name, function(chr)
   {
   	f.tmp <- fap.randompos.check.dfr[, as.character(chr)]
   	f.tmp[f.tmp == 1] <- 1-1e-11
	f.tmp[f.tmp == 0] <- 1e-11
    list(fap = f.tmp, 
        covar = pars.randompos.dfr[pars.randompos.dfr$name == chr, c("highest6sum.spwin","n")])
   })  
names(vars.list.check.randompos) <- pars.randompos.dfr$name
sapply(vars.list.check.randompos, function(lst) length(lst[[1]]))
vars.list.check.randompos[[592]]
## Ok, this consists of 750 NAs.


## --------------------------------------------------------------------------------------------------------
## Ecliptic region (everything for check)
## --------------------------------------------------------------------------------------------------------

## --------------------------------------------------------------------------------------------------------
## Dataframe for all the ecliptic locations, holding all estimated or given parameters:
## --------------------------------------------------------------------------------------------------------

ecldirlist <- list.files(ecldir.name, pattern = "_gaia")
pars.ecl.dfr <- data.frame(name = ecldirlist)
pars.ecl.dfr$ra <- as.numeric(sapply(ecldirlist, function(chr)
   {
#  	chr <- ecldirlist[1] 
  	unlist(strsplit(chr, split = "_"))[1]
   }))
pars.ecl.dfr$dec <- as.numeric(sapply(ecldirlist, function(chr)
   {
#  	chr <- ecldirlist[1] 
  	unlist(strsplit(chr, split = "_"))[2]
   }))
pars.ecl.dfr$lambda <- NA
pars.ecl.dfr$beta <- NA
pars.ecl.dfr[,4:5] <- t(apply(pars.ecl.dfr[,2:3], 1, function(vec)
   {
   	eq2ecl.fun(vec[1], vec[2])
   }))
pars.ecl.dfr$n <- NA
pars.ecl.dfr[1:10,]

## --------------------------------------------------------------------------------------------------------
## Download the FAPs and fill up the pars.ecl.dfr$n
## --------------------------------------------------------------------------------------------------------

## the number of observations are in the commented-oout header in the
## FAP files:

fap.ecl.dfr <- as.data.frame(matrix(ncol = length(ecldirlist), nrow = 1500))
colnames(fap.ecl.dfr) <- ecldirlist
for(ii in 1:length(ecldirlist))
   {
#    ii <- 2
	dd <- ecldirlist[ii]
	pp <- list.files(paste(ecldir.name, dd, sep = ""), full.names = TRUE,
	   pattern = ".dat")
    cc <- unlist(strsplit(readLines(pp, n = 2), " "))
    pars.ecl.dfr[pars.ecl.dfr$name == dd, "n"] <- as.numeric(cc[length(cc)])
	fap.ecl.dfr[,dd] <- read.table(pp)[,6]
   }
# fap.ecl.dfr[1:20, 1:20]


## --------------------------------------------------------------------------------------------------------
## Get the spectral window values and frequencies of its maximum power on the ecl
## --------------------------------------------------------------------------------------------------------

## max.ind.ecl, max.power.ecl, max.freq.ecl:

load(paste(fineecldir.name, "max_aliases.RObj", sep = ""))
# dim(max.power.ecl)
dfr <- as.data.frame(t(max.power.ecl))
# dfr[1:10,]
colnames(dfr)[3:19] <- paste("at", (1:17)*4, sep = "") 
# is.numeric(dfr$ra)
# is.numeric(dfr$dec)
# intersect(dfr$dec, pars.ecl.dfr$dec)
dfr1 <- merge(pars.ecl.dfr, dfr)
# dfr1[1:10,]
# dim(dfr1)
## Don't forget to leave out column spwin.at.20 from the matrix given to classification.
# colnames(max.power.fullsky7000)
## Ok, the names are the same.

dfr <- as.data.frame(t(max.freq.ecl))
# dfr[1:5,]
colnames(dfr)[3:19] <- paste("freq.at", (1:17)*4, sep = "") 
# is.numeric(dfr$ra)
# is.numeric(dfr$dec)
# identical(dfr$dec, pars.ecl.dfr$dec)
# identical(dfr$ra, pars.ecl.dfr$ra)
dfr2 <- merge(dfr1, dfr)
# dfr2[1:5,]

pars.ecl.dfr <- dfr2
rm(dfr,dfr1,dfr2,max.power.ecl,max.freq.ecl,max.ind.ecl)
# pars.ecl.dfr[801:810,]


## --------------------------------------------------------------------------------------------------------
## Attach variance of the times
## --------------------------------------------------------------------------------------------------------

pars.ecl.dfr$vart <- NA
d.tmp <- list.files(paste0(fineecldir.name, "all_samplings/"))
for(ii in d.tmp) 
  {
   tt <- as.numeric(unlist(strsplit(ii, "_"))[2:3])
   t.tmp <- read.table(paste0(fineecldir.name, "all_samplings/", ii))[,1]
   pars.ecl.dfr$vart[pars.ecl.dfr$ra == tt[1] & pars.ecl.dfr$dec == tt[2]] <- var(t.tmp)
  }
  

## --------------------------------------------------------------------------------------------------------
## Attach the sums of various alias powers to the pars file 
## --------------------------------------------------------------------------------------------------------

## Put the sums of the highest k peaks into the parameters dataframe:

	dfr <- pars.ecl.dfr
	dfr$sum.spwin <- apply(dfr[,7:23], MAR = 1, FUN = sum)
	dfr$mean.spwin <- apply(dfr[,c(6:23)], MAR = 1, FUN = function(vec)
       sum(vec[2:18]/vec[1]))
	dfr$total.spwin <- apply(dfr[,c(6:23)], MAR = 1, FUN = function(vec)
       sum(vec[2:18]*vec[1]))
	dfr$highest.spwin <- apply(dfr[,7:23], MAR = 1, FUN = max)
	dfr$highest2sum.spwin <- apply(dfr[,7:23], MAR = 1, FUN = function(vec)
       {
       	v.tmp <- vec
       	v.tmp <- v.tmp[order(v.tmp, decreasing = TRUE)]
       	sum(v.tmp[1:2])
       	})
	dfr$highest3sum.spwin <- apply(dfr[,7:23], MAR = 1, FUN = function(vec)
       {
       	v.tmp <- vec
       	v.tmp <- v.tmp[order(v.tmp, decreasing = TRUE)]
       	sum(v.tmp[1:3])
       	})
	dfr$highest4sum.spwin <- apply(dfr[,7:23], MAR = 1, FUN = function(vec)
       {
       	v.tmp <- vec
       	v.tmp <- v.tmp[order(v.tmp, decreasing = TRUE)]
       	sum(v.tmp[1:4])
       	})
	dfr$highest5sum.spwin <- apply(dfr[,7:23], MAR = 1, FUN = function(vec)
       {
       	v.tmp <- vec
       	v.tmp <- v.tmp[order(v.tmp, decreasing = TRUE)]
       	sum(v.tmp[1:5])
       	})
	dfr$highest6sum.spwin <- apply(dfr[,7:23], MAR = 1, FUN = function(vec)
       {
       	v.tmp <- vec
       	v.tmp <- v.tmp[order(v.tmp, decreasing = TRUE)]
       	sum(v.tmp[1:6])
       	})
	dfr$highest7sum.spwin <- apply(dfr[,7:23], MAR = 1, FUN = function(vec)
       {
       	v.tmp <- vec
       	v.tmp <- v.tmp[order(v.tmp, decreasing = TRUE)]
       	sum(v.tmp[1:7])
       	})
	dfr$highest8sum.spwin <- apply(dfr[,7:23], MAR = 1, FUN = function(vec)
       {
       	v.tmp <- vec
       	v.tmp <- v.tmp[order(v.tmp, decreasing = TRUE)]
       	sum(v.tmp[1:8])
       	})
	dfr$highest9sum.spwin <- apply(dfr[,7:23], MAR = 1, FUN = function(vec)
       {
       	v.tmp <- vec
       	v.tmp <- v.tmp[order(v.tmp, decreasing = TRUE)]
       	sum(v.tmp[1:9])
       	})
    dfr[1:10,]
    pars.ecl.dfr <- dfr
#	dfr <- pars.ecl.dfr[order(pars.ecl.dfr[,ii]),]

pars.ecl.dfr$sum.of.five.spwin <- 
   pars.ecl.dfr$at12 + pars.ecl.dfr$at16 + pars.ecl.dfr$at28 + 
   pars.ecl.dfr$at40 + pars.ecl.dfr$at68
pars.ecl.dfr$sum.of.six.spwin <- pars.ecl.dfr$at4 + 
   pars.ecl.dfr$at12 + pars.ecl.dfr$at16 + pars.ecl.dfr$at28 + 
   pars.ecl.dfr$at40 + pars.ecl.dfr$at68
pars.ecl.dfr$sum.of.seven.spwin <- pars.ecl.dfr$at4 + 
   pars.ecl.dfr$at12 + pars.ecl.dfr$at16 + pars.ecl.dfr$at28 + 
   pars.ecl.dfr$at40 + pars.ecl.dfr$at52 + pars.ecl.dfr$at68
pars.ecl.dfr$sum.of.eight.spwin <- pars.ecl.dfr$at4 + 
   pars.ecl.dfr$at12 + pars.ecl.dfr$at16 + pars.ecl.dfr$at28 + 
   pars.ecl.dfr$at40 + pars.ecl.dfr$at52 + pars.ecl.dfr$at68 + 
   pars.ecl.dfr$at56
pars.ecl.dfr$sum.of.nine.spwin <- pars.ecl.dfr$at4 + 
   pars.ecl.dfr$at12 + pars.ecl.dfr$at16 + pars.ecl.dfr$at28 + 
   pars.ecl.dfr$at40 + pars.ecl.dfr$at52 + pars.ecl.dfr$at68 + 
   pars.ecl.dfr$at56 + pars.ecl.dfr$at24
pars.ecl.dfr$sum.of.ten.spwin <- pars.ecl.dfr$at4 + 
   pars.ecl.dfr$at12 + pars.ecl.dfr$at16 + pars.ecl.dfr$at28 + 
   pars.ecl.dfr$at40 + pars.ecl.dfr$at52 + pars.ecl.dfr$at68 + 
   pars.ecl.dfr$at56 + pars.ecl.dfr$at24 + pars.ecl.dfr$at44
dim(pars.ecl.dfr)
  

## --------------------------------------------------------------------------------------------------------
## Attach the necessary classification results 
## --------------------------------------------------------------------------------------------------------

pars.ecl.dfr$class7000_7group <- predict(gmxt7000[[7]],
    newdata = pars.ecl.dfr[, c(7:10,12:23)])$classification

pars.ecl.dfr$fine.class.for.init4[
   is.element(pars.ecl.dfr$class7000_7group, c(1,2,3))] <- 2
pars.ecl.dfr$fine.class.for.init4[
   is.element(pars.ecl.dfr$class7000_7group, c(4,5,6,7))] <- 1

## Just a check:
## Which frequencies can occur among the highest spectral peaks?

which.highest.ecl <- apply(pars.ecl.dfr[, 7:23], MAR = 1, FUN = function(vec)
   {
    o.tmp <- order(vec, decreasing = TRUE)[1:6]
    o.tmp*4
   })

table(which.highest.ecl)/1589
table(which.highest.line)/900
table(which.highest)/720

vars.list.check.ecl <- lapply(pars.ecl.dfr$name, function(chr)
   {
   	f.tmp <- fap.ecl.dfr[, as.character(chr)]
   	f.tmp[f.tmp == 1] <- 1-1e-11
	f.tmp[f.tmp == 0] <- 1e-11
    list(fap = f.tmp, 
        covar = pars.ecl.dfr[pars.ecl.dfr$name == chr, c("highest6sum.spwin","n")])
   })  
names(vars.list.check.ecl) <- pars.ecl.dfr$name
sapply(vars.list.check.ecl, function(lst) length(lst[[1]]))
vars.list.check.ecl[[1456]]


## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------
## Compute p-values for all faps at all locations
## --------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------

## Just to recall the names:
# best.cut.line.model <- cutrough.line.fits[[4]]
# best.class.line.model <- pars.line4
# best.cut.randompos.model <- cutrough.randompos.fits[[4]]
# best.class.randompos.model <- pars.randompos4

## --------------------------------------------------------------------------------------------------------
## Check all 4 models on the taken-out faps of the line (750/location)
## --------------------------------------------------------------------------------------------------------

## Using the best.cut.line.model (cut-based, fits were done on the first 750 faps
## of the line):
checkline.best.cut.line.model <- matrix(nrow = 750, ncol = 900)
dimnames(checkline.best.cut.line.model) <- list(NULL, pars.line.dfr$name)
for(ii in 1:length(vars.list.check.line))
#for(ii in 1:2)
   {
   	nn <- as.character(pars.line.dfr$name[ii])
   	p.tmp <- unlist(vars.list.check.line[[nn]][[2]])
   	bp.tmp <- betapars.fun(p.tmp, best.cut.line.model)
   	checkline.best.cut.line.model[, nn] <- 
   	    sapply(vars.list.check.line[[nn]][[1]], ptrans.fun, betapars = bp.tmp)
   }
# sum(checkline.best.cut.line.model[,2] > 0.95)

## Using the best.cut.line.model (class-based, fits were done on the first 750 faps
## of the line):
checkline.best.class.line.model <- matrix(nrow = 750, ncol = 900)
dimnames(checkline.best.class.line.model) <- list(NULL, pars.line.dfr$name)
for(ii in 1:length(vars.list.check.line))
#for(ii in 1:2)
   {
   	nn <- as.character(pars.line.dfr$name[ii])
   	p.tmp <- c(unlist(vars.list.check.line[[nn]][[2]]), 
   	   pars.line.dfr$fine.class.for.init4[pars.line.dfr$name == nn])
   	bp.tmp <- betapars.fun(p.tmp, best.class.line.model)
   	checkline.best.class.line.model[, nn] <- 
   	    sapply(vars.list.check.line[[nn]][[1]], ptrans.fun, betapars = bp.tmp)
   }
# sum(checkline.best.class.line.model[,1] > 0.95)

## Using the best.cut.randompos.model (cut-based, fits were done on the first 750 faps
## of random positions all over the sky):
checkline.best.cut.randompos.model <- matrix(nrow = 750, ncol = 900)
dimnames(checkline.best.cut.randompos.model) <- list(NULL, pars.line.dfr$name)
for(ii in 1:length(vars.list.check.line))
#for(ii in 1:2)
   {
   	nn <- as.character(pars.line.dfr$name[ii])
   	p.tmp <- unlist(vars.list.check.line[[nn]][[2]])
   	bp.tmp <- betapars.fun(p.tmp, best.cut.randompos.model)
   	checkline.best.cut.randompos.model[, nn] <- 
   	    sapply(vars.list.check.line[[nn]][[1]], ptrans.fun, betapars = bp.tmp)
   }
# sum(checkline.best.cut.randompos.model[,1] > 0.95)

## Using the best.class.randompos.model (class-based, fits were done on the first 750 faps
## of random positions all over the sky):
checkline.best.class.randompos.model <- matrix(nrow = 750, ncol = 900)
dimnames(checkline.best.class.randompos.model) <- list(NULL, pars.line.dfr$name)
for(ii in 1:length(vars.list.check.line))
#for(ii in 1:2)
   {
   	nn <- as.character(pars.line.dfr$name[ii])
   	p.tmp <- c(unlist(vars.list.check.line[[nn]][[2]]), 
   	   pars.line.dfr$fine.class.for.init4[pars.line.dfr$name == nn])
   	bp.tmp <- betapars.fun(p.tmp, best.class.randompos.model)
   	checkline.best.class.randompos.model[, nn] <- 
   	    sapply(vars.list.check.line[[nn]][[1]], ptrans.fun, betapars = bp.tmp)
   }
# sum(checkline.best.class.randompos.model[,2] > 0.95)
## All right, not that bad.

## --------------------------------------------------------------------------------------------------------
## Check all 4 models on the taken-out faps of the random positions over the sky (750/location)
## --------------------------------------------------------------------------------------------------------

## Using the best.cut.line.model (cut-based, fits were done on the first 750 faps
## of the line):
checkrandompos.best.cut.line.model <- matrix(nrow = 750, ncol = 720)
dimnames(checkrandompos.best.cut.line.model) <- list(NULL, pars.randompos.dfr$name)
for(ii in 1:length(vars.list.check.randompos))
#for(ii in 1:2)
   {
   	nn <- as.character(pars.randompos.dfr$name[ii])
   	p.tmp <- unlist(vars.list.check.randompos[[nn]][[2]])
   	bp.tmp <- betapars.fun(p.tmp, best.cut.line.model)
   	checkrandompos.best.cut.line.model[, nn] <- 
   	    sapply(vars.list.check.randompos[[nn]][[1]], ptrans.fun, betapars = bp.tmp)
   }
# sum(checkrandompos.best.cut.line.model[,1] > 0.95)
# sum(checkrandompos.best.cut.line.model[,2] > 0.95)

## Using the best.cut.line.model (class-based, fits were done on the first 750 faps
## of the line):
checkrandompos.best.class.line.model <- matrix(nrow = 750, ncol = 720)
dimnames(checkrandompos.best.class.line.model) <- list(NULL, pars.randompos.dfr$name)
for(ii in 1:length(vars.list.check.randompos))
#for(ii in 1:2)
   {
   	nn <- as.character(pars.randompos.dfr$name[ii])
   	p.tmp <- c(unlist(vars.list.check.randompos[[nn]][[2]]), 
   	   pars.randompos.dfr$fine.class.for.init4[pars.randompos.dfr$name == nn])
   	bp.tmp <- betapars.fun(p.tmp, best.class.line.model)
   	checkrandompos.best.class.line.model[, nn] <- 
   	    sapply(vars.list.check.randompos[[nn]][[1]], ptrans.fun, betapars = bp.tmp)
   }
# sum(checkrandompos.best.class.line.model[,1] > 0.95)
# sum(checkrandompos.best.class.line.model[,2] > 0.95)

## Using the best.cut.randompos.model (cut-based, fits were done on the first 750 faps
## of random positions all over the sky):
checkrandompos.best.cut.randompos.model <- matrix(nrow = 750, ncol = 720)
dimnames(checkrandompos.best.cut.randompos.model) <- list(NULL, pars.randompos.dfr$name)
for(ii in 1:length(vars.list.check.randompos))
#for(ii in 1:2)
   {
   	nn <- as.character(pars.randompos.dfr$name[ii])
   	p.tmp <- unlist(vars.list.check.randompos[[nn]][[2]])
   	bp.tmp <- betapars.fun(p.tmp, best.cut.randompos.model)
   	checkrandompos.best.cut.randompos.model[, nn] <- 
   	    sapply(vars.list.check.randompos[[nn]][[1]], ptrans.fun, betapars = bp.tmp)
   }
# sum(checkrandompos.best.cut.randompos.model[,1] > 0.95)
# sum(checkrandompos.best.cut.randompos.model[,2] > 0.95)

## Using the best.class.randompos.model (class-based, fits were done on the first 750 faps
## of random positions all over the sky):
checkrandompos.best.class.randompos.model <- matrix(nrow = 750, ncol = 720)
dimnames(checkrandompos.best.class.randompos.model) <- list(NULL, pars.randompos.dfr$name)
for(ii in 1:length(vars.list.check.randompos))
#for(ii in 1:2)
   {
   	nn <- as.character(pars.randompos.dfr$name[ii])
   	p.tmp <- c(unlist(vars.list.check.randompos[[nn]][[2]]), 
   	   pars.randompos.dfr$fine.class.for.init4[pars.randompos.dfr$name == nn])
   	bp.tmp <- betapars.fun(p.tmp, best.class.randompos.model)
   	checkrandompos.best.class.randompos.model[, nn] <- 
   	    sapply(vars.list.check.randompos[[nn]][[1]], ptrans.fun, betapars = bp.tmp)
   }
# sum(checkrandompos.best.class.randompos.model[,1] > 0.95)
# sum(checkrandompos.best.class.randompos.model[,2] > 0.95)

## --------------------------------------------------------------------------------------------------------
## Check all 4 models on all the faps of the ecliptic region (1500/location)
## --------------------------------------------------------------------------------------------------------

## Using the best.cut.line.model (cut-based, fits were done on the first 750 faps
## of the line):
checkecl.best.cut.line.model <- matrix(nrow = 1500, ncol = 1589)
dimnames(checkecl.best.cut.line.model) <- list(NULL, pars.ecl.dfr$name)
for(ii in 1:length(vars.list.check.ecl))
#for(ii in 1:2)
   {
   	nn <- as.character(pars.ecl.dfr$name[ii])
   	p.tmp <- unlist(vars.list.check.ecl[[nn]][[2]])
   	bp.tmp <- betapars.fun(p.tmp, best.cut.line.model)
   	checkecl.best.cut.line.model[, nn] <- 
   	    sapply(vars.list.check.ecl[[nn]][[1]], ptrans.fun, betapars = bp.tmp)
   }
# sum(checkecl.best.cut.line.model[,1] > 0.95)
# sum(checkecl.best.cut.line.model[,2] > 0.95)

## Using the best.cut.line.model (class-based, fits were done on the first 750 faps
## of the line):
checkecl.best.class.line.model <- matrix(nrow = 1500, ncol = 1589)
dimnames(checkecl.best.class.line.model) <- list(NULL, pars.ecl.dfr$name)
for(ii in 1:length(vars.list.check.ecl))
#for(ii in 1:2)
   {
   	nn <- as.character(pars.ecl.dfr$name[ii])
   	p.tmp <- c(unlist(vars.list.check.ecl[[nn]][[2]]), 
   	   pars.ecl.dfr$fine.class.for.init4[pars.ecl.dfr$name == nn])
   	bp.tmp <- betapars.fun(p.tmp, best.class.line.model)
   	checkecl.best.class.line.model[, nn] <- 
   	    sapply(vars.list.check.ecl[[nn]][[1]], ptrans.fun, betapars = bp.tmp)
   }
# sum(checkecl.best.class.line.model[,1] > 0.95)
# sum(checkecl.best.class.line.model[,2] > 0.95)

## Using the best.cut.randompos.model (cut-based, fits were done on the first 750 faps
## of random positions all over the sky):
checkecl.best.cut.randompos.model <- matrix(nrow = 1500, ncol = 1589)
dimnames(checkecl.best.cut.randompos.model) <- list(NULL, pars.ecl.dfr$name)
for(ii in 1:length(vars.list.check.ecl))
#for(ii in 1:2)
   {
   	nn <- as.character(pars.ecl.dfr$name[ii])
   	p.tmp <- unlist(vars.list.check.ecl[[nn]][[2]])
   	bp.tmp <- betapars.fun(p.tmp, best.cut.randompos.model)
   	checkecl.best.cut.randompos.model[, nn] <- 
   	    sapply(vars.list.check.ecl[[nn]][[1]], ptrans.fun, betapars = bp.tmp)
   }
# sum(checkecl.best.cut.randompos.model[,1] > 0.95)
# sum(checkecl.best.cut.randompos.model[,2] > 0.95)

## Using the best.class.randompos.model (class-based, fits were done on the first 750 faps
## of random positions all over the sky):
checkecl.best.class.randompos.model <- matrix(nrow = 1500, ncol = 1589)
dimnames(checkecl.best.class.randompos.model) <- list(NULL, pars.ecl.dfr$name)
for(ii in 1:length(vars.list.check.ecl))
#for(ii in 1:2)
   {
   	nn <- as.character(pars.ecl.dfr$name[ii])
   	p.tmp <- c(unlist(vars.list.check.ecl[[nn]][[2]]), 
   	   pars.ecl.dfr$fine.class.for.init4[pars.ecl.dfr$name == nn])
   	bp.tmp <- betapars.fun(p.tmp, best.class.randompos.model)
   	checkecl.best.class.randompos.model[, nn] <- 
   	    sapply(vars.list.check.ecl[[nn]][[1]], ptrans.fun, betapars = bp.tmp)
   }
# sum(checkecl.best.class.randompos.model[,1] > 0.95)
# sum(checkecl.best.class.randompos.model[,2] > 0.95)



###########################################################################################################
## --------------------------------------------------------------------------------------------------------
## Checking the model
## --------------------------------------------------------------------------------------------------------
###########################################################################################################
   


## --------------------------------------------------------------------------------------------------------
## Histograms of FAPs
## --------------------------------------------------------------------------------------------------------

tr.fapdens.checkline.best.cut.line.model <- 
   apply(checkline.best.cut.line.model, 2, function(vec)
   hist(vec, breaks = seq(0, 1, by = 0.05), plot = F)$counts/750)
tr.fapdens.checkline.best.class.line.model <-
   apply(checkline.best.class.line.model, 2, function(vec)
   hist(vec, breaks = seq(0, 1, by = 0.05), plot = F)$counts/750)
tr.fapdens.checkline.best.cut.randompos.model <-
   apply(checkline.best.cut.randompos.model, 2, function(vec)
   hist(vec, breaks = seq(0, 1, by = 0.05), plot = F)$counts/750)
tr.fapdens.checkline.best.class.randompos.model <-
   apply(checkline.best.class.randompos.model, 2, function(vec)
   hist(vec, breaks = seq(0, 1, by = 0.05), plot = F)$counts/750)

tr.fapdens.checkrandompos.best.cut.line.model <- 
   apply(checkrandompos.best.cut.line.model, 2, function(vec)
   hist(vec, breaks = seq(0, 1, by = 0.05), plot = F)$counts/750)
tr.fapdens.checkrandompos.best.class.line.model <-
   apply(checkrandompos.best.class.line.model, 2, function(vec)
   hist(vec, breaks = seq(0, 1, by = 0.05), plot = F)$counts/750)
tr.fapdens.checkrandompos.best.cut.randompos.model <-
   apply(checkrandompos.best.cut.randompos.model, 2, function(vec)
   hist(vec, breaks = seq(0, 1, by = 0.05), plot = F)$counts/750)
tr.fapdens.checkrandompos.best.class.randompos.model <-
   apply(checkrandompos.best.class.randompos.model, 2, function(vec)
   hist(vec, breaks = seq(0, 1, by = 0.05), plot = F)$counts/750)

tr.fapdens.checkecl.best.cut.line.model <- 
   apply(checkecl.best.cut.line.model, 2, function(vec)
   hist(vec, breaks = seq(0, 1, by = 0.05), plot = F)$counts/1500)
tr.fapdens.checkecl.best.class.line.model <-
   apply(checkecl.best.class.line.model, 2, function(vec)
   hist(vec, breaks = seq(0, 1, by = 0.05), plot = F)$counts/1500)
tr.fapdens.checkecl.best.cut.randompos.model <-
   apply(checkecl.best.cut.randompos.model, 2, function(vec)
   hist(vec, breaks = seq(0, 1, by = 0.05), plot = F)$counts/1500)
tr.fapdens.checkecl.best.class.randompos.model <-
   apply(checkecl.best.class.randompos.model, 2, function(vec)
   hist(vec, breaks = seq(0, 1, by = 0.05), plot = F)$counts/1500)

fapdens.checkline <- 
   apply(fap.line.dfr[751:1500,], 2, function(vec)
   hist(1-vec, breaks = seq(0, 1, by = 0.05), plot = F)$counts/750)
fapdens.checkrandompos <- 
   apply(fap.randompos.check.dfr, 2, function(vec)
   hist(1-vec, breaks = seq(0, 1, by = 0.05), plot = F)$counts/750)
fapdens.checkecl <- 
   apply(fap.ecl.dfr, 2, function(vec)
   hist(1-vec, breaks = seq(0, 1, by = 0.05), plot = F)$counts/1500)
 
## --------------------------------------------------------------------------------------------------------
## A simple visualization: the lines of histograms on top of each other
## --------------------------------------------------------------------------------------------------------

## Independent check on the line:
quartz(height = 7.8, width = 12.5)
par(mfrow = c(2,2), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.6,0), oma = c(0,0,2,0))
matplot(seq(0.025, 0.975, by = 0.05), fapdens.checkline, col = "grey", 
   type = "l", xlab = "1 - FAP", ylab = "Density", ylim = c(0,0.3),
   main = "Model based on cuts at line ") 
matplot(seq(0.025, 0.975, by = 0.05), tr.fapdens.checkline.best.cut.line.model, 
   type = "l", col = 1, add = TRUE)
abline(h = c(0.05, 0.1, 0.15, 0.2), col = "aquamarine", lwd = c(2,1,1,1))
lines(seq(0.025, 0.975, by = 0.05),
   apply(tr.fapdens.checkline.best.cut.line.model, 1, mean), col = "orange",
   lwd = 2, lty = 2) 
matplot(seq(0.025, 0.975, by = 0.05), fapdens.checkline, col = "grey", 
   type = "l", xlab = "1 - FAP", ylab = "Density", ylim = c(0,0.3),
   main = "Model based on class at line ") 
matplot(seq(0.025, 0.975, by = 0.05), tr.fapdens.checkline.best.class.line.model, 
   type = "l", col = 1, add = TRUE) 
abline(h = c(0.05, 0.1, 0.15, 0.2), col = "aquamarine", lwd = c(2,1,1,1))
lines(seq(0.025, 0.975, by = 0.05),
   apply(tr.fapdens.checkline.best.class.line.model, 1, mean), col = "orange",
   lwd = 2, lty = 2) 
matplot(seq(0.025, 0.975, by = 0.05), fapdens.checkline, col = "grey", 
   type = "l", xlab = "1 - FAP", ylab = "Density", ylim = c(0,0.3),
   main = "Model based on cuts at random positions ") 
matplot(seq(0.025, 0.975, by = 0.05), tr.fapdens.checkline.best.cut.randompos.model, 
   type = "l", col = 1, add = TRUE) 
abline(h = c(0.05, 0.1, 0.15, 0.2), col = "aquamarine", lwd = c(2,1,1,1))
lines(seq(0.025, 0.975, by = 0.05),
   apply(tr.fapdens.checkline.best.cut.randompos.model, 1, mean), col = "orange",
   lwd = 2, lty = 2) 
matplot(seq(0.025, 0.975, by = 0.05), fapdens.checkline, col = "grey", 
   type = "l", xlab = "1 - FAP", ylab = "Density", ylim = c(0,0.3),
   main = "Model based on class at random positions") 
matplot(seq(0.025, 0.975, by = 0.05), tr.fapdens.checkline.best.class.randompos.model, 
   type = "l", col = 1, add = TRUE) 
abline(h = c(0.05, 0.1, 0.15, 0.2), col = "aquamarine", lwd = c(2,1,1,1))
lines(seq(0.025, 0.975, by = 0.05),
   apply(tr.fapdens.checkline.best.class.randompos.model, 1, mean), col = "orange",
   lwd = 2, lty = 2) 
mtext("Line", outer = TRUE, font = 2)

## Independent check on the random positions:
quartz(height = 7.8, width = 12.5)
par(mfrow = c(2,2), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.6,0), oma = c(0,0,2,0))
matplot(seq(0.025, 0.975, by = 0.05), fapdens.checkrandompos, col = "grey", 
   type = "l", xlab = "1 - FAP", ylab = "Density", ylim = c(0,0.3),
   main = "Model based on cuts at line ") 
matplot(seq(0.025, 0.975, by = 0.05), tr.fapdens.checkrandompos.best.cut.line.model, 
   type = "l", col = 1, add = TRUE) 
abline(h = c(0.05, 0.1, 0.15, 0.2), col = "aquamarine", lwd = c(2,1,1,1))
lines(seq(0.025, 0.975, by = 0.05),
   apply(tr.fapdens.checkrandompos.best.cut.line.model, 1, mean), col = "orange",
   lwd = 2, lty = 2) 
matplot(seq(0.025, 0.975, by = 0.05), fapdens.checkrandompos, col = "grey", 
   type = "l", xlab = "1 - FAP", ylab = "Density", ylim = c(0,0.3),
   main = "Model based on class at line ") 
matplot(seq(0.025, 0.975, by = 0.05), tr.fapdens.checkrandompos.best.class.line.model, 
   type = "l", col = 1, add = TRUE) 
abline(h = c(0.05, 0.1, 0.15, 0.2), col = "aquamarine", lwd = c(2,1,1,1))
lines(seq(0.025, 0.975, by = 0.05),
   apply(tr.fapdens.checkrandompos.best.class.line.model, 1, mean), col = "orange",
   lwd = 2, lty = 2) 
matplot(seq(0.025, 0.975, by = 0.05), fapdens.checkrandompos, col = "grey", 
   type = "l", xlab = "1 - FAP", ylab = "Density", ylim = c(0,0.3),
   main = "Model based on cuts at random positions ") 
matplot(seq(0.025, 0.975, by = 0.05), tr.fapdens.checkrandompos.best.cut.randompos.model, 
   type = "l", col = 1, add = TRUE) 
abline(h = c(0.05, 0.1, 0.15, 0.2), col = "aquamarine", lwd = c(2,1,1,1))
lines(seq(0.025, 0.975, by = 0.05),
   apply(tr.fapdens.checkrandompos.best.cut.randompos.model, 1, mean), col = "orange",
   lwd = 2, lty = 2) 
matplot(seq(0.025, 0.975, by = 0.05), fapdens.checkrandompos, col = "grey", 
   type = "l", xlab = "1 - FAP", ylab = "Density", ylim = c(0,0.3),
   main = "Model based on class at random positions") 
matplot(seq(0.025, 0.975, by = 0.05), tr.fapdens.checkrandompos.best.class.randompos.model, 
   type = "l", col = 1, add = TRUE) 
abline(h = c(0.05, 0.1, 0.15, 0.2), col = "aquamarine", lwd = c(2,1,1,1))
lines(seq(0.025, 0.975, by = 0.05),
   apply(tr.fapdens.checkrandompos.best.class.randompos.model, 1, mean), col = "orange",
   lwd = 2, lty = 2) 
mtext("Random positions", outer = TRUE, font = 2)

## Independent check on the ecliptic region:
quartz(height = 7.8, width = 12.5)
par(mfrow = c(2,2), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.6,0), oma = c(0,0,2,0))
matplot(seq(0.025, 0.975, by = 0.05), fapdens.checkecl, col = "grey", 
   type = "l", xlab = "1 - FAP", ylab = "Density", ylim = c(0,0.3),
   main = "Model based on cuts at line ") 
matplot(seq(0.025, 0.975, by = 0.05), tr.fapdens.checkecl.best.cut.line.model, 
   type = "l", col = 1, add = TRUE) 
abline(h = c(0.05, 0.1, 0.15, 0.2), col = "aquamarine", lwd = c(2,1,1,1))
lines(seq(0.025, 0.975, by = 0.05),
   apply(tr.fapdens.checkecl.best.cut.line.model, 1, mean), col = "orange",
   lwd = 2, lty = 2) 
matplot(seq(0.025, 0.975, by = 0.05), fapdens.checkecl, col = "grey", 
   type = "l", xlab = "1 - FAP", ylab = "Density", ylim = c(0,0.3),
   main = "Model based on class at line ") 
matplot(seq(0.025, 0.975, by = 0.05), tr.fapdens.checkecl.best.class.line.model, 
   type = "l", col = 1, add = TRUE) 
abline(h = c(0.05, 0.1, 0.15, 0.2), col = "aquamarine", lwd = c(2,1,1,1))
lines(seq(0.025, 0.975, by = 0.05),
   apply(tr.fapdens.checkecl.best.class.line.model, 1, mean), col = "orange",
   lwd = 2, lty = 2) 
matplot(seq(0.025, 0.975, by = 0.05), fapdens.checkecl, col = "grey", 
   type = "l", xlab = "1 - FAP", ylab = "Density", ylim = c(0,0.3),
   main = "Model based on cuts at random positions ") 
matplot(seq(0.025, 0.975, by = 0.05), tr.fapdens.checkecl.best.cut.randompos.model, 
   type = "l", col = 1, add = TRUE) 
abline(h = c(0.05, 0.1, 0.15, 0.2), col = "aquamarine", lwd = c(2,1,1,1))
lines(seq(0.025, 0.975, by = 0.05),
   apply(tr.fapdens.checkecl.best.cut.randompos.model, 1, mean), col = "orange",
   lwd = 2, lty = 2) 
matplot(seq(0.025, 0.975, by = 0.05), fapdens.checkecl, col = "grey", 
   type = "l", xlab = "1 - FAP", ylab = "Density", ylim = c(0,0.3),
   main = "Model based on class at random positions") 
matplot(seq(0.025, 0.975, by = 0.05), tr.fapdens.checkecl.best.class.randompos.model, 
   type = "l", col = 1, add = TRUE) 
abline(h = c(0.05, 0.1, 0.15, 0.2), col = "aquamarine", lwd = c(2,1,1,1))
lines(seq(0.025, 0.975, by = 0.05),
   apply(tr.fapdens.checkecl.best.class.randompos.model, 1, mean), col = "orange",
   lwd = 2, lty = 2) 
mtext("Ecliptic rectangle", outer = TRUE, font = 2)



## --------------------------------------------------------------------------------------------------------
## Plot: a, b coeffs and the spectral window peaks in the ecliptic region
## --------------------------------------------------------------------------------------------------------

## Location-wise beta fits for the ecliptic region
betafits.ecl <- apply(fap.ecl.dfr, M = 2, function(vec) try(fitbeta.fun(vec)))
betapar.poswise.ecl <- matrix(ncol = ncol(fap.ecl.dfr), nrow = 4)
dimnames(betapar.poswise.ecl) <- list(c("a.single","a.single.se","b.single","b.single.se"),
   colnames(fap.ecl.dfr))
betapar.poswise.ecl[c("a.single","b.single"),] <- sapply(betafits.ecl, function(lst) 
   if(inherits(lst, "try-error")) rep(NA,2) else lst$par)
betapar.poswise.ecl[c("a.single.se","b.single.se"),] <- sapply(betafits.ecl, function(lst) 
   if(inherits(lst, "try-error")) rep(NA,2) else sqrt(diag(solve(-lst$hessian))))

## Merge them into the parameter dataframe
identical(colnames(betapar.poswise.ecl), as.character(pars.ecl.dfr$name))
## True, so it can go by simple merging..
dfr1 <- data.frame(t(betapar.poswise.ecl))
dfr1$name <- colnames(betapar.poswise.ecl)
dfr <- merge(pars.ecl.dfr, dfr1)
# dfr[1:10,]
## looks ok
pars.ecl.dfr <- dfr
# pars.ecl.dfr[1:10,]
rm(dfr, dfr1)

## Take the logarithms, after all, these are what I fit: 
pars.ecl.dfr$loga.single <- log(pars.ecl.dfr$a.single)
pars.ecl.dfr$logb.single <- log(pars.ecl.dfr$b.single)

range(pars.ecl.dfr$ra)
range(pars.ecl.dfr$dec)


## Save some interesting results, the pars.*.* and the fap.*.* files:
save(pars.ecl.dfr, pars.line.dfr, pars.randompos.dfr, 
   file = "results_betamodelFitting1/pars.whatever.dfr.RObj")
save(fap.ecl.dfr, fap.line.dfr, fap.randompos.dfr, fap.randompos.check.dfr,
   file = "results_betamodelFitting1/fap.whatever.dfr.RObj")





























