
## This workspace is about:
## 


## TESTING MY IDEA ABOUT USING EXTREME-VALUE STAT + BOOTSTRAP TO CHECK FOR
## SIGNIFICANCE OF PEAKS

## The idea:
##
##   H0:  the observed time series is noise = the peak in the periodogram is 
##        not exceedingly high for a noise sequence
##
##   H1:  the time series has some periodic content 
##

##   1. Find the peak of the periodogram of the observed signal;
##      Find other information about it:
##      - Spectral window: the locations of its smallest (0) values
##      - characteristic width of the peaks 
##   2. Bootstrap the signal (1000 repetitions); for each repetition, 
##      choose randomly one frequency, and by shifting the center of the spectral
##      window there, identify a set of quasi-orthogonal frequencies {F_i} 
##      (the most independent possible!). They should be half as many as there 
##      are observations in the time series. Procedure: iterative selection,
##      at each iteration, add the frequency that has the smallest cumulative
##      correlation with the already selected set.
##   3. Around each quasi-orth frequency, take also the values in a little 
##      window of width equal to the width of the peaks (due to leakage), and take
##      the maximum of the periodogram values at these frequencies. Also, since this
##      might require too much time, take the single value at the selected independent
##      frequency, to compare to the max of the whole peak. Do this for all
##      bootstrap repetition. 
##   4. Do this with the simple near-zeros of the spectral window (not orthogonal
##      system, quasi-orthogonal to only the first selected random freq.)
##   5. Do this also with a set of randomly selected frequencies (Brockwell-Davis:
##      in an evenly-sampled case, m frequencies are asymptotically independent
##      exponential, if m is fixed and n -> Infty).
##      Maybe something like this is true in the uneven sampling, too?
##   6. Do this also with randomly selected  Fourier frequencies; by the sepctral
##      window, at the Fourier frequencies the window values are quite small.
##   7. Fit a GEV for the peaks, and use it for probability levels.


## VARYING SIGNAL-TO-NOISE RATIO

## EVIR PROCEDURE GEV

## Check it on signals with varying relative amplitude:

## Simulations for the detection limit of the extreme-value limit

##    NONPARAMETRIC BOOTSTRAP OF A (SIMULATED) OBSERVED TIME SERIES

##    Simulation:
##      1. sine, times on a grid of around dT ~ 10 min ~ 5e-03 day.
##      2. ECL, same times
##
##    Say observing run is during T = 25 days ==> N = 25/5e-03 = 5000.
##
##    This is on a dT  = 5e-03 day grid ==> upper frequency limit 1/(2*dT) = 1e+02
##    Fourier transform: independent frequencies: k/(N*dT), k = 0,1,....(N-1)/2:
##       (0, 1/25, 2/25, ...., (5000-1)/50)
##    Peak width: 1/T = 0.04/d, if I want to have an oversampling factor of 8,
##    this means a finer grid of dF = 0.005/day instead of dF1 = 0.02 
##    Use the fixed frequency grid dF = 0.005/day 

##    Extreme-value FAP estimates that should be checked:

##    Effect of thinning (corresponding to the number of missing data on the grid):
##       1. compute the spectral windows for the thinned time series on the fixed
##          frequency grid;
##       2. compute the spectra on the orthogonal or near-orthogonal frequencies of a
##          rnadomly chosen frequency and around them in the little windows, on the 
##          fixed grid, using both the (0,25) and (0,100) frequency ranges with 
##          thinned observations: 5000 -> 2000 -> 1000 -> 500 -> 200 -> 100
##          -> 50 -> 25; 
##          take the maximum of these reduced spectra;
##       3. make plots: diagnostic plots are important, as max of blocks of 100,
##          50 and 25 observations are not necessarily follow a good GEV 

## QUESTIONS

## 1. Is it really necessary to compute the (ii-I):(ii+I) periodogram values?
##    Or define heuristically the extremal index as the inverse of the half-width
##    of the peaks, and use it as in EVS.

## 2. More simulations: Do all this also with
##    - same star, quasi-periodic sampling (only nights at one telescope!)
##    - EA eclipsing binary with irregular sampling
##    - EA eclipsing binary with quasi-periodic sampling 


## TO DO:
## 



colorvec <- c("red","powderblue","blue","chocolate","green","pink",
              "orange","violet","deepskyblue","yellow","lightgreen","sienna","violetred",
              "aquamarine","royalblue","blueviolet","burlywood","coral","darkgoldenrod",
              "darkgreen","darkorchid","gray","indianred","lightpink","limegreen",
              "magenta","maroon","olivedrab","purple","rosybrown","salmon","seagreen",
              "tomato","turquoise","grey")
   

WORKSPACE.NAME <- "extrFAP7FirstStar.RData"
save.image(WORKSPACE.NAME)

require(xtable)
#require(evd)
require(evir)



   
## ----------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------
## Sine (FirstStar), random times on a grid of around dT = 10 min = 5e-03 day
## ----------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------


obsTimes <- seq(0, 25, by = 5e-03)

F <- 3.379865
A1 <- 0.3
obsMagTrue <- A1*sin(2*pi*obsTimes*F)
plot(obsTimes, obsMagTrue, type = "l", xlim = c(0,15))


## Add noise: generate random numbers say with gamma distr., as sdev of the noise at each point
## Then generate noise with these as given error bars, and differends 
sdev <- rgamma(5001, shape = 4, scale = 0.0025)

obsMagNoisy <- obsMagTrue + rnorm(5001, mean = 0, sd = 5*sdev)
quartz()
plot(obsTimes, obsMagTrue, type = "n", xlim = c(0,1), ylim = c(-0.4,0.4))
points(obsTimes, obsMagNoisy, pch = 16, col = "deepskyblue", cex = 0.5)
segments(obsTimes, obsMagNoisy + 10*sdev, obsTimes, obsMagNoisy - 10*sdev, 
   col = "deepskyblue")

## Ok, accept obsMagNoisy generated with 5*sdev:
## During 25 days
## 5000 uniformly randomly distributed observations for the first step

FirstStar <- data.frame(time = obsTimes, 
  mag.g = obsMagNoisy,
  mag.error.g = 5*sdev)
FirstStar[1:20,]

FirstStar$sim5000 <- 1
FirstStar$sim2000 <- 0
FirstStar$sim2000[sample(1:5001, size = 2000, replace = FALSE)] <- 1
FirstStar$sim1000 <- 0
FirstStar$sim1000[sample(1:5001, size = 1000, replace = FALSE)] <- 1
FirstStar$sim500 <- 0
FirstStar$sim500[sample(1:5001, size = 500, replace = FALSE)] <- 1
FirstStar$sim200 <- 0
FirstStar$sim200[sample(1:5001, size = 200, replace = FALSE)] <- 1
FirstStar$sim100 <- 0
FirstStar$sim100[sample(1:5001, size = 100, replace = FALSE)] <- 1
FirstStar$sim50 <- 0
FirstStar$sim50[sample(1:5001, size = 50, replace = FALSE)] <- 1
FirstStar$sim25 <- 0
FirstStar$sim25[sample(1:5001, size = 25, replace = FALSE)] <- 1


mean1 <- mean(FirstStar$mag.g)
sd1 <- sqrt(var(FirstStar$mag.g))
qqnorm(FirstStar$mag.g)
abline(c(mean1,sd1), col = 2)

nObs <- c(2000,1000,500,200,100,50,25)

## Change signal-to-noise:

A1 <- 0.2
obsMagTrue <- A1*sin(2*pi*obsTimes*F)
sdev <- rgamma(5001, shape = 4, scale = 0.0025)
obsMagNoisy <- obsMagTrue + rnorm(5001, mean = 0, sd = 5*sdev)
FirstStar$mag.g2 <- obsMagNoisy
FirstStar$mag.error.g2 <- 5*sdev

A1 <- 0.15
obsMagTrue <- A1*sin(2*pi*obsTimes*F)
sdev <- rgamma(5001, shape = 4, scale = 0.0025)
obsMagNoisy <- obsMagTrue + rnorm(5001, mean = 0, sd = 5*sdev)
FirstStar$mag.g3 <- obsMagNoisy
FirstStar$mag.error.g3 <- 5*sdev

A1 <- 0.1
obsMagTrue <- A1*sin(2*pi*obsTimes*F)
sdev <- rgamma(5001, shape = 4, scale = 0.0025)
obsMagNoisy <- obsMagTrue + rnorm(5001, mean = 0, sd = 5*sdev)
FirstStar$mag.g4 <- obsMagNoisy
FirstStar$mag.error.g4 <- 5*sdev

A1 <- 0.05
obsMagTrue <- A1*sin(2*pi*obsTimes*F)
sdev <- rgamma(5001, shape = 4, scale = 0.0025)
obsMagNoisy <- obsMagTrue + rnorm(5001, mean = 0, sd = 5*sdev)
FirstStar$mag.g5 <- obsMagNoisy
FirstStar$mag.error.g5 <- 5*sdev

A1 <- 0.025
obsMagTrue <- A1*sin(2*pi*obsTimes*F)
sdev <- rgamma(5001, shape = 4, scale = 0.0025)
obsMagNoisy <- obsMagTrue + rnorm(5001, mean = 0, sd = 5*sdev)
FirstStar$mag.g6 <- obsMagNoisy
FirstStar$mag.error.g6 <- 5*sdev

quartz(height = 3.5, width = 10)
plot(obsTimes, obsMagTrue, type = "n", xlim = c(0,1), ylim = c(-0.4,0.4))
points(obsTimes, obsMagNoisy, pch = 16, col = "violet", cex = 0.5)
segments(obsTimes, obsMagNoisy + 10*sdev, obsTimes, obsMagNoisy - 10*sdev, 
   col = "violet")

amplitudes <- c(0.3,0.2,0.15,0.1,0.05,0.025)


## Create nighttime observational sequences too:

obsPoints2000 <- c(1:80, 201:280, 401:480, 601:680, 801:880, 
   1001:1080, 1201:1280, 1401:1480, 1601:1680, 1801:1880,
   2001:2080, 2201:2280, 2401:2480, 2601:2680, 2801:2880,
   3001:3080, 3201:3280, 3401:3480, 3601:3680, 3801:3880,
   4001:4080, 4201:4280, 4401:4480, 4601:4680, 4801:4880)
obsPoints1000 <- c(1:40, 201:240, 401:440, 601:640, 801:840, 
   1001:1040, 1201:1240, 1401:1440, 1601:1640, 1801:1840,
   2001:2040, 2201:2240, 2401:2440, 2601:2640, 2801:2840,
   3001:3040, 3201:3240, 3401:3440, 3601:3640, 3801:3840,
   4001:4040, 4201:4240, 4401:4440, 4601:4640, 4801:4840)

## 5000 uniformly randomly distributed observations for the first step,
## then 2000 observations  on the grid covering the nights (80 gridpoints
## from the 200 representing one day), then 1000 covering a shorter night
## (40 gridpoints from the 200), then randomly chosen points (500, 200, 100, 
## 50, 25) from these last shorter nights.

FirstStar[1:20,]

FirstStar$nights2000 <- 0
FirstStar$nights2000[obsPoints2000] <- 1
FirstStar$nights1000 <- 0
FirstStar$nights1000[obsPoints1000] <- 1
FirstStar$nights500 <- 0
FirstStar$nights500[sample(obsPoints1000, size = 500, replace = FALSE)] <- 1
FirstStar$nights200 <- 0
FirstStar$nights200[sample(obsPoints1000, size = 200, replace = FALSE)] <- 1
FirstStar$nights100 <- 0
FirstStar$nights100[sample(obsPoints1000, size = 100, replace = FALSE)] <- 1
FirstStar$nights50 <- 0
FirstStar$nights50[sample(obsPoints1000, size = 50, replace = FALSE)] <- 1
FirstStar$nights25 <- 0
FirstStar$nights25[sample(obsPoints1000, size = 25, replace = FALSE)] <- 1



quartz(height = 3.5, width = 10)
plot(FirstStar$time[FirstStar$nights200 == 1], 
     FirstStar$mag.g4[FirstStar$nights200 == 1], 
     type = "n", xlim = c(0,2.5), ylim = c(-0.6,0.6))
points(FirstStar$time[FirstStar$nights200 == 1], 
     FirstStar$mag.g4[FirstStar$nights200 == 1], 
     pch = 16, col = "violet", cex = 0.5)
segments(FirstStar$time[FirstStar$nights200 == 1], 
     FirstStar$mag.g4[FirstStar$nights200 == 1] + 
     10*FirstStar$mag.error.g4[FirstStar$nights200 == 1], 
     FirstStar$time[FirstStar$nights200 == 1], 
     FirstStar$mag.g4[FirstStar$nights200 == 1] - 
     10*FirstStar$mag.error.g4[FirstStar$nights200 == 1], 
   col = "violet")
## Looks okay.



## ----------------------------------------------------------------------------------------------
## Functions for floating-mean weighted Lomb-Scargle
## ----------------------------------------------------------------------------------------------


demean.fun <- function(vec)  {vec - mean(vec, na.rm = TRUE)}

# "h" is vector of expression values for time points "t".
# SpectralPowerDensity will be evaluated at given TestFrequencies.
# Nindepedent of the TestFrequencies are assumed to be independent.
## I'll compute it without Tau.

## ff: the freq at which to evaluate LS
## tt: the vector of observation times
## hh: the vector of measurements (length equal to tt)
## ww: the weights


FMSPD3.fun <- function(ff, tt, hh, ww = rep(1/length(tt), length(tt)))
	{
      Omega <- 2*pi*ff
      co <- cos(Omega*tt)
      si <- sin(Omega*tt)
#      co^2 + si^2      
      Y <- sum(hh*ww)
      CO <- sum(co*ww)
      SI <- sum(si*ww)
      CChat <- sum(ww*co^2)
      SShat <- sum(ww*si^2)
      CShat <- sum(co*si*ww)
      
      YY <- sum(ww*(hh - Y)^2)
      YC <- sum((hh - Y)*co*ww)
      YS <- sum((hh - Y)*si*ww)
      CC <- CChat - CO^2
      SS <- SShat - SI^2
      CS <- CShat - CO*SI
      
      DD <- CC*SS - CS^2

      (SS*YC^2 + CC*YS^2 - 2*CS*YC*YS) / (YY * DD)
	}


MyFMLS3.fun <- function(t, h, TestFrequencies, Weights = rep(1/length(t), length(t)))
 {
  stopifnot( length(t) == length(h) )
  if (length(t) > 0)
   {
    hResidual    <- h - mean(h)
    SpectralPowerDensity <- sapply(TestFrequencies, FUN = FMSPD3.fun,
        tt = t, hh = hResidual, ww = Weights)
   } else {
    SpectralPowerDensity <- NA
   }
  return(SpectralPowerDensity)
}


## Misi's deeming code, put into sapply form:
deeming1.fun  <- function(ts,y,f) 
    {
     DFT <- sapply(f, FUN = function(freq)
         	sum(y*exp(1i*2.0*pi*freq*ts)))
     1.0/length(ts)*abs(DFT)
    }


   
## ----------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------
## Calculate the periodograms for random thinnings and random nighttimes of FirstStar
## (both with and without weights)
## ----------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------


FreqGrid1 <- seq(0, 100, by = dfreq)

mag.names <- c("mag.g","mag.g2","mag.g3","mag.g4","mag.g5","mag.g6")
mag.error.names <- c("mag.error.g","mag.error.g2","mag.error.g3",
    "mag.error.g4","mag.error.g5","mag.error.g6")

FirstStar.pgram <- vector(6, mode = "list")
names(FirstStar.pgram) <- mag.names
FirstStar.wpgram <- vector(6, mode = "list")
names(FirstStar.pgram) <- mag.names

for(jj in 1:6)
	{
	bootI.tmp <- data.frame(freq = FreqGrid1,
	  sim2000 = NA, sim1000 = NA, sim500 = NA, 
	  sim200 = NA, sim100 = NA, sim50 = NA, sim25 = NA,
	  nights2000 = NA, nights1000 = NA, nights500 = NA, 
	  nights200 = NA, nights100 = NA, nights50 = NA, nights25 = NA)
	for(nnn in c(2000,1000,500,200,100,50,25))
 		{
  		 print(system.time(
   		 {
  		 colnnn <- paste("sim", as.character(nnn), sep = "")
 	     dfr <- FirstStar[FirstStar[,colnnn] == 1,]
 	     ts <- demean.fun(dfr$time)
      	 ps <- demean.fun(dfr[, mag.names[jj]])
      	 ws <- (1/dfr[, mag.error.names[jj]]^2)/sum(1/dfr[, mag.error.names[jj]]^2)
      	 bootI.tmp[,colnnn] <-  MyFMLS3.fun(TestFrequencies = FreqGrid1, t = ts, h = ps,  
         		Weights = ws)
  		 colnnn <- paste("nights", as.character(nnn), sep = "")
 	     dfr <- FirstStar[FirstStar[,colnnn] == 1,]
 	     ts <- demean.fun(dfr$time)
      	 ps <- demean.fun(dfr[, mag.names[jj]])
      	 ws <- (1/dfr[, mag.error.names[jj]]^2)/sum(1/dfr[, mag.error.names[jj]]^2)
      	 bootI.tmp[,colnnn] <-  MyFMLS3.fun(TestFrequencies = FreqGrid1, t = ts, h = ps,  
         		Weights = ws)
 		  }))
		 }
	FirstStar.wpgram[[mag.names[jj]]] <- bootI.tmp
	
	bootI.tmp <- data.frame(freq = FreqGrid1,
	  sim2000 = NA, sim1000 = NA, sim500 = NA, 
	  sim200 = NA, sim100 = NA, sim50 = NA, sim25 = NA,
	  nights2000 = NA, nights1000 = NA, nights500 = NA, 
	  nights200 = NA, nights100 = NA, nights50 = NA, nights25 = NA)
	for(nnn in c(2000,1000,500,200,100,50,25))
 		{
  		 print(system.time(
   		 {
  		 colnnn <- paste("sim", as.character(nnn), sep = "")
 	     dfr <- FirstStar[FirstStar[,colnnn] == 1,]
 	     ts <- demean.fun(dfr$time)
      	 ps <- demean.fun(dfr[, mag.names[jj]])
      	 ws <- rep(1/nnn, nnn)
      	 bootI.tmp[,colnnn] <-  MyFMLS3.fun(TestFrequencies = FreqGrid1, t = ts, h = ps,  
         		Weights = ws)
  		 colnnn <- paste("nights", as.character(nnn), sep = "")
 	     dfr <- FirstStar[FirstStar[,colnnn] == 1,]
 	     ts <- demean.fun(dfr$time)
      	 ps <- demean.fun(dfr[, mag.names[jj]])
     	 ws <- rep(1/nnn, nnn)
      	 bootI.tmp[,colnnn] <-  MyFMLS3.fun(TestFrequencies = FreqGrid1, t = ts, h = ps,  
         		Weights = ws)
 		  }))
		 }
	FirstStar.pgram[[mag.names[jj]]] <- bootI.tmp
 	} 

quartz(height = 7.5, width = 8)
colnnn <- "nights50"
mag.tmp <- "mag.g6"
magerr.tmp <- "mag.error.g6"
#mycol <- "orange"
par(mfrow = c(3,1), mar = c(3,3,2,1))
plot(FirstStar[FirstStar[, colnnn] == 1, "time"], 
  FirstStar[FirstStar[, colnnn] == 1, mag.tmp], xlim = c(0,10),
  ylim = c(-0.5,0.5), pch = 16, cex = 0.5)
segments(FirstStar[FirstStar[, colnnn] == 1, "time"], 
  FirstStar[FirstStar[, colnnn] == 1, mag.tmp] + 2*FirstStar[FirstStar[, colnnn] == 1, magerr.tmp],
  FirstStar[FirstStar[, colnnn] == 1, "time"], 
  FirstStar[FirstStar[, colnnn] == 1, mag.tmp] - 2*FirstStar[FirstStar[, colnnn] == 1, magerr.tmp])
plot(FirstStar.pgram[[mag.tmp]]$freq, FirstStar.pgram[[mag.tmp]][,colnnn], type = "n",
  xlab = "Frequency", ylab = "Chi-squared", main = paste(colnnn, mag.tmp), ylim = c(0,0.5),
  xlim = c(0,100))
abline(v = F, col = mycol, lwd = 2)
abline(h = max(FirstStar.bootI[[mag.tmp]][[colnnn]]["dep.max",]), col = "orange")
points(FirstStar.pgram[[mag.tmp]]$freq, FirstStar.pgram[[mag.tmp]][,colnnn], type = "h")
plot(FirstStar.wpgram[[mag.tmp]]$freq, FirstStar.wpgram[[mag.tmp]][,colnnn], type = "n",
  xlab = "Frequency", ylab = "Chi-squared", main = paste(colnnn, mag.tmp), ylim = c(0,0.5),
  xlim = c(0,100))
abline(v = F, col = mycol, lwd = 2)
points(FirstStar.wpgram[[mag.tmp]]$freq, FirstStar.wpgram[[mag.tmp]][,colnnn], type = "h")
abline(h = max(FirstStar.bootI[[mag.tmp]][[colnnn]]["wdep.max",]), col = "orange")
## Good, looks like I could cover a range where the signal falls below detection level.


length(FirstStar.wpgram)
length(FirstStar.wpgram[mag.names])
FirstStar.wpgram <- FirstStar.wpgram[mag.names]
dim(FirstStar.wpgram[["mag.g5"]])
   
## ----------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------
## Calculate the spectral window for each random thinning and random nighttimes of FirstStar
## (only non-weighted - does the weighted Deeming exist?)
## ----------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------

## Use deeming

dfreq <- 0.0025

FreqGrid1 <- seq(0, 100, by = dfreq)

SpWin <- data.frame(freq = FreqGrid1,
   sim2000 = NA, sim1000 = NA, sim500 = NA, sim200 = NA, 
   sim100 = NA, sim50 = NA, sim25 = NA,
   nights2000 = NA, nights1000 = NA, nights500 = NA, nights200 = NA, 
   nights100 = NA, nights50 = NA, nights25 = NA)
for(cc in 1:7)
   SpWin[,cc+1] <- deeming1.fun(FirstStar$time[FirstStar[,cc+4] == 1],
      rep(1, sum(FirstStar[,cc+4])), FreqGrid1) 
for(cc in 1:7)
   SpWin[,cc+8] <- deeming1.fun(FirstStar$time[FirstStar[,cc+21] == 1],
      rep(1, sum(FirstStar[,cc+21])), FreqGrid1) 
SpWin[1:100,]


quartz(height = 4.5, width = 12)
plot(FreqGrid1, SpWin[,15], type = "h", xlim = c(0,100))


quartz(height = 7.5, width = 12)
par(mfrow = c(4,2), mar = c(2,2,2,1), mgp = c(2,1,0))
for(ii in 1:7)
   plot(SpWin[,1], SpWin[,8+ii], type = "h", ylim = c(0,1),
       xlim = c(0,0.2), main = colnames(SpWin)[8+ii])


##
SpWin5000 <- deeming1.fun(FirstStar$time[-5001],
      rep(1, 5000), FreqGrid1)
plot(FreqGrid1, SpWin5000, type = "h", xlim = c(0,100))
       
## How many near-zero values are there in each spectral window? Are there more
## near-zero locations when I have 2000 observations than when I have only 25?

quantile(SpWin[,2], 2000/40000)
quantile(SpWin[,3], 1000/40000)
quantile(SpWin[,4], 500/40000)
quantile(SpWin[,5], 200/40000)
quantile(SpWin[,6], 100/40000)
quantile(SpWin[,7], 50/40000)
quantile(SpWin[,8], 25/40000)
## Wow,wow,wow! the quantile values are roughly the same...

SpWin[SpWin[,2] < quantile(SpWin[,2], 2000/40000), 1]
SpWin[SpWin[,3] < quantile(SpWin[,3], 1000/40000), 1]
SpWin[SpWin[,4] < quantile(SpWin[,4], 500/40000), 1]
SpWin[SpWin[,5] < quantile(SpWin[,5], 200/40000), 1]
SpWin[SpWin[,6] < quantile(SpWin[,6], 100/40000), 1]
SpWin[SpWin[,7] < quantile(SpWin[,7], 50/40000), 1]
SpWin[SpWin[,8] < quantile(SpWin[,8], 25/40000), 1]

quartz(height = 3.5, width = 12)
par(mar = c(3,3,2,1), mgp = c(2,1,0))
plot(SpWin[,1], SpWin[,8], type = "l", ylim = c(0,0.06),
       xlim = c(0,0.2), xlab = "Frequency [1/day]", ylab = "Sp. window")
for(ii in 1:6)
   lines(SpWin[,1], SpWin[,8-ii], type = "l", col = colorvec[ii])
abline(h = quantile(SpWin[,8], 12/40000), col = 1)
abline(h = quantile(SpWin[,7], 25/40000), col = colorvec[1])
abline(h = quantile(SpWin[,6], 50/40000), col = colorvec[2])
abline(h = quantile(SpWin[,5], 100/40000), col = colorvec[3])
abline(h = quantile(SpWin[,4], 250/40000), col = colorvec[4])
abline(h = quantile(SpWin[,3], 500/40000), col = colorvec[5])
abline(h = quantile(SpWin[,2], 1000/40000), col = colorvec[6])
legend(-0.2, 0.06, bty = "o", bg = "white",
   legend = c(25, 50, 100, 200, 500, 1000, 2000), col = c(1, colorvec[1:6]), lty = 1)

## OK. So, to find approximately uncorrelated frequencies
## for a randomly selected initial frequency:
##  1. I just get one frequency where the shifted spectral window is close enough
##     to zero
##  2. to the spwindow centred at 0 I add a spwindow centered at this second freq
##  3. look for the minimum of the sum of spwindows and select this freq as
##     the third freq
##  4. repeat from 2: add one more shifted spwindow, and select the freq where the
##     sum of 3 is minimal.
## This does not work very quickly.


## Do the same procedure, but with the center at 40001 as the initial
## value, and with respect to that, do the minimum-finding over the
## double 1:80001 interval (to be able to shift it to any random value)
## This is not ideal, part of the quasi-orthogonal frequencies will fall out
## during shifting, so it will not be the ideally orthogonalized.


#samplingCol <- "sim100"
#nFreq <- 220
orthSelect.fun <- function(spWin, samplingCol, nFr)
   {
    dfr0 <- c(spWin[,c(samplingCol)], 
          spWin[40000:1,samplingCol])
    dfr <- c(dfr0[80001:2], dfr0[1:80001])
    ## center: dfr[80001]
#    dfr[39991:40010,]
    ind.tmp <- numeric(nFr) 
#    ind.tmp[1] <- sample(2:40000, size = 1)
    ind.tmp[1] <- 40001
    sp.tmp <- dfr[(40001):(120001)]
    sp.tmp[intersect((ind.tmp[1]-8):(ind.tmp[1]+8), 1:80001)] <- 
          sp.tmp[intersect((ind.tmp[1]-8):(ind.tmp[1]+8), 1:80001)] + 
          rep(100, length(intersect((ind.tmp[1]-8):(ind.tmp[1]+8), 1:80001)))
#    sp.tmp[intersect((ind.tmp[1]-18):(ind.tmp[1]+18), 1:80001)]
    for(ii in 2:(nFr)) 
      {
       ind1 <- which(sp.tmp == min(sp.tmp, na.rm = TRUE))
       ind.tmp[ii] <- if(length(ind1) == 1) ind1 else 
          ind1[sample(1:length(ind1), size = 1)]
#       round(sp.tmp, 4)
       sp.tmp[intersect((ind.tmp[ii]-8):(ind.tmp[ii]+8), 1:80001)] <- 
          sp.tmp[intersect((ind.tmp[ii]-8):(ind.tmp[ii]+8), 1:80001)] + 
          rep(100, length(intersect((ind.tmp[ii]-8):(ind.tmp[ii]+8), 1:80001)))
#       sp.tmp[intersect((ind.tmp[2]-18):(ind.tmp[2]+18), 1:80001)]
#       sp.tmp[intersect((ind.tmp[ii]-18):(ind.tmp[ii]+18), 1:80001)]
       sp.tmp <- sp.tmp +
         dfr[(80001-ind.tmp[ii] + 1):(160001-ind.tmp[ii] + 1)]
      }
    ind.tmp
   }

## Test it:
system.time(huhh <- orthSelect.fun(SpWin, "sim500", 1100))
plot(c(SpWin[40001:2,"sim500"], SpWin[,"sim500"]), type = "h", lwd = 0.2,
   xlim = c(39950,40050))
abline(v = huhh, col = "red", lwd = 0.2)
## Not ideal, one of the found frequencies is very close to the central peak,
## with a very high correlation.

## Compute the systems of quasi-orthogonal values for all thinning:
## there would be half as many independent frequencies in the periodogram 
## than observations in the time series (amplitudes are kept, phases dropped);
## so: for 2000 obs, 1000 freqs; for 1000 obs, 500 freqs;....
## no need 1000 freqs for one max (max-stability would give approx. give the
## same max-distribution, I hope), so I do it only with 500 freqs;
##  

nNames <- colnames(FirstStar)[c(5:11, 22:28)]
nFreq <- c(1100,1100,600,300,150,70,40)  

orthSystFirst.ind <- vector(14, mode = "list") 
names(orthSystFirst.ind) <- nNames
for(jj in 1:7)
 {
  # random thinning:
  ind.tmp <- orthSelect.fun(SpWin, nNames[jj], nFreq[jj])
  ind.tmp <- sort(ind.tmp, decreasing = FALSE)
  orthSystFirst.ind[[jj]] <- ind.tmp - 40001
  # nighttimes:
  ind.tmp <- orthSelect.fun(SpWin, nNames[jj+7], nFreq[jj])
  ind.tmp <- sort(ind.tmp, decreasing = FALSE)
  orthSystFirst.ind[[jj+7]] <- ind.tmp - 40001
 }

orthSystFirst.ind[[7]] + rev(orthSystFirst.ind[[7]])


## Compute the systems of values orthogonal only to the 
## central frequency, but not to each other, for all thinning:

orthFirst.ind <- vector(14, mode = "list") 
names(orthFirst.ind) <- nNames

for(jj in 1:7)
 {
  # random thinning:
  dfr <- SpWin[,c("freq", paste("sim", nObs[jj], sep = ""))]
  ind.tmp <- numeric(nFreq[jj]/2) 
  for(ii in 1:(nFreq[jj]/2))
    {
     ind.tmp[ii] <- which(dfr[,2] == min(dfr[,2], na.rm = TRUE))
     dfr[(ind.tmp[ii]-8):(ind.tmp[ii]+8),] <- 100
#     ii <- ii+1
    }
  orthFirst.ind[[jj]] <- sort(ind.tmp, decreasing = FALSE)
  # nighttimes:
  dfr <- SpWin[,c("freq", paste("nights", nObs[jj], sep = ""))]
  ind.tmp <- numeric(nFreq[jj]/2) 
  for(ii in 1:(nFreq[jj]/2))
    {
     ind.tmp[ii] <- which(dfr[,2] == min(dfr[,2], na.rm = TRUE))
     dfr[(ind.tmp[ii]-8):(ind.tmp[ii]+8),] <- 100
#     ii <- ii+1
    }
  orthFirst.ind[[jj+7]] <- sort(ind.tmp, decreasing = FALSE)
 }
#orthFirst.ind <- orth.ind 
 
nrow(SpWin[orth.ind[[6]],c(1,7)])
SpWin[orth.ind[["nights100"]],c("freq","nights100")]
SpWin[1:10,]

nSamplings <- c("syst.max", "dep.max", "fou.max", "rand.max",
  "wsyst.max", "wdep.max", "wfou.max", "wrand.max")




## ----------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------
## Nonparametric bootstrap of FirstStar, random thinning (might fail to describe the tail)
## ----------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------


## Peak width: 1/T = 0.04/d, if I want to have an oversampling factor of 8,
## this means the finer grid of dF = 0.005/day instead of dF1 = 0.02. I use 
## even more, the finest grid 

mag.names <- c("mag.g","mag.g2","mag.g3","mag.g4","mag.g5","mag.g6")
mag.error.names <- c("mag.error.g","mag.error.g2","mag.error.g3",
    "mag.error.g4","mag.error.g5","mag.error.g6")

nTestFreq <- c(500,500,250,100,50,25,13,500,500,250,100,50,25,13)  
#nObs <- c(nObs,nObs)

FirstStar.bootI <- vector(6, mode = "list")
names(FirstStar.bootI) <- mag.names

for(jj in 1:6)
	{
	bootI.tmp <- vector(14, mode = "list")
	names(bootI.tmp) <- nNames
	for(ii in 1:14)
 		{
  		 print(system.time(
   		 {
 		 nnn <- nObs[ii]
  		 colnnn <- nNames[ii]
 	     dfr <- FirstStar[FirstStar[,colnnn] == 1,]
 	     ts <- demean.fun(dfr$time)

 	     bootI.tmp[[nNames[ii]]] <- replicate(1000,
		   {
  	  	    CentralInd <- sample(2:40000, size = 1)
  		 
#  		    ## Select the frequencies orthogonal to this central freq:
#  		    ind.tmp <- c(-rev(orthFirst.ind[[colnnn]]), 0, orthFirst.ind[[colnnn]])
#  		    IndGrid0 <- CentralInd + ind.tmp
# 		    IndGrid0 <- IndGrid0[(IndGrid0 > 0 & IndGrid0 < 40001)]
# 		    IndGrid0 <- if(length(IndGrid0) <= nTestFreq[ii])
# 	           IndGrid0 else sort(sample(IndGrid0, size = nTestFreq[ii]))
# 	        FreqGrid0 <- FreqGrid1[IndGrid0]
  		 
  		    ## Select now the quasi-orthogonal frequencies:
  		    IndGridSyst <- CentralInd + orthSystFirst.ind[[colnnn]]
  			IndGridSyst <- IndGridSyst[(IndGridSyst > 0 & IndGridSyst < 40001)]
 	     	IndGridSyst <- if(length(IndGridSyst) <= nTestFreq[ii]) 
 	        	 IndGridSyst else sort(sample(IndGridSyst, size = nTestFreq[ii]))
  		 	FreqGridSyst <- FreqGrid1[IndGridSyst]
  		 
  		    ## Select Fourier frequencies:
  			FreqGridFou <- (1:2499)/25
 	     	FreqGridFou <- sort(sample(FreqGridFou, size = nTestFreq[ii]))
   		 
#  		 	## Select peak-width intervals around these quasi-orth frequencies:
#  		 	IndGrid <- sapply(IndGridSyst, function(kk)
# 	        	(kk-8):(kk+8), simplify = FALSE)
# 	     	IndGrid <- IndGrid[sapply(IndGrid, function(vec)
# 	        	!any(vec < 2 | vec >= 40001))]
#  			FreqGridDep <- lapply(IndGrid, function(vec)
#  		    	FreqGrid1[vec])
  		    
  		 	## Select peak-width intervals around these quasi-orth frequencies:
  		 	FreqGridDep <- sapply(FreqGridFou, function(ff)
 	        	ff + (-8:8)/25, simplify = FALSE)
 	     	FreqGridDep <- FreqGridDep[sapply(FreqGridDep, function(vec)
 	        	!any(vec < 0.04 | vec >= 100))]
            FreqGridDep <- sort(unlist(FreqGridDep))
            
  			## Select now totally randomly nTestFreq frequencies, regardless of orthogonality
  		 	## or independence:
  		 	IndGridRand <- sort(sample(2:40000, size = nTestFreq[ii]))
  		  	FreqGridRand <- FreqGrid1[IndGridRand]
  		  	
            ## The bootstrap:  		 
      	    ind <- sample(1:nrow(dfr), replace = TRUE)
      	    ps <- demean.fun(dfr[ind, mag.names[jj]])
      	    ## Unweighted:
            ws <- rep(1/nnn, nnn)
            frspDep.tmp <- sapply(FreqGridDep, MyFMLS3.fun, t = ts, h = ps,  
         		Weights = ws)
      	    frspFou.tmp <- sapply(FreqGridFou, MyFMLS3.fun, t = ts, h = ps,  
         		Weights = ws)
     	    frspSyst.tmp <- sapply(FreqGridSyst, MyFMLS3.fun, t = ts, h = ps,  
         		Weights = ws)
            frspRand.tmp <- sapply(FreqGridRand, MyFMLS3.fun, t = ts, h = ps,  
         		Weights = ws)
            ## Weighted:
            ws <- (1/dfr[ind, mag.error.names[jj]]^2)/sum(1/dfr[ind, mag.error.names[jj]]^2)
            wfrspDep.tmp <- sapply(FreqGridDep, MyFMLS3.fun, t = ts, h = ps,  
         		Weights = ws)
      	    wfrspFou.tmp <- sapply(FreqGridFou, MyFMLS3.fun, t = ts, h = ps,  
         		Weights = ws)
     	    wfrspSyst.tmp <- sapply(FreqGridSyst, MyFMLS3.fun, t = ts, h = ps,  
         		Weights = ws)
            wfrspRand.tmp <- sapply(FreqGridRand, MyFMLS3.fun, t = ts, h = ps,  
         		Weights = ws)
		    c(
	          syst.max = max(frspSyst.tmp[!is.na(frspSyst.tmp) & abs(frspSyst.tmp) != Inf]),
	          dep.max = max(frspDep.tmp[!is.na(frspDep.tmp) & abs(frspDep.tmp) != Inf]), 
              fou.max = max(frspFou.tmp[!is.na(frspFou.tmp) & abs(frspFou.tmp) != Inf]),
  		      rand.max = max(frspRand.tmp[!is.na(frspRand.tmp) & abs(frspRand.tmp) != Inf]),
	          wsyst.max = max(wfrspSyst.tmp[!is.na(wfrspSyst.tmp) & abs(wfrspSyst.tmp) != Inf]),
	          wdep.max = max(wfrspDep.tmp[!is.na(wfrspDep.tmp) & abs(wfrspDep.tmp) != Inf]), 
              wfou.max = max(wfrspFou.tmp[!is.na(wfrspFou.tmp) & abs(wfrspFou.tmp) != Inf]),
  		      wrand.max = max(wfrspRand.tmp[!is.na(wfrspRand.tmp) & abs(wfrspRand.tmp) != Inf])
  		      )
		   }) 
 		 }))
	  }
	FirstStar.bootI[[mag.names[jj]]] <- bootI.tmp
	cat(mag.names[jj], " ready \n")
 	} 
save.image(WORKSPACE.NAME)

lapply(FirstStar.bootI, function(lst)
  sapply(lst, function(mat) apply(mat, M = 1, F = mean)))


## primpeaks.fgev: first argument: which signal-to-noise ratio;
##   second argument: which sampling pattern / nb of observations;
##   third argument: which peak selection pattern in the periodogram.
primpeaks.fgev <- lapply(FirstStar.bootI, function(lst)
  lapply(lst,  function(mat)
    {
     res1.tmp <- vector(8, mode = "list")
     names(res1.tmp) <- c("syst.max","dep.max","fou.max","rand.max", 
        "wsyst.max","wdep.max","wfou.max","wrand.max") 
     for(ii in 1:8)
       {
       	res1.tmp[[ii]] <- try(gev(mat[ii,]))
       	}
      res1.tmp
      }))
primpeaks.fgev[[6]]

## Extract the parameters, the standard errors, and the endpoint when it is 
## an upper endpoint:

primpeaks.gevpar <- lapply(primpeaks.fgev, function(lst)
   lapply(lst, function(lst2)
        sapply(lst2, function(lst1) if(!inherits(lst1, "try-error")) lst1$par.ests else rep(NA, 3))))

primpeaks.gevse <- lapply(primpeaks.fgev, function(lst)
   lapply(lst, function(lst2)
        sapply(lst2, function(lst1) if(!inherits(lst1, "try-error")) lst1$par.ses else rep(NA, 3))))

primpeaks.gevend <- lapply(primpeaks.gevpar, function(lst)
   sapply(lst, function(mat)
       apply(mat, MAR = 2, FUN = function(vec) 
          if(vec[1] >= 0 | is.na(vec[1])) NA else {vec[3] - vec[2]/vec[1]})))
## This looks incredibly ugly:
##  - the endpoints are very small for high number of observations (= high number 
##    of test frequencies), and increasing as the number of observations decrease;
##  - the endpoints become more scattered as the number of observations decrease
##    (but at least this is something expected.)

endse.fun <- function(pars, varcov)
    {
     fn <- c(pars[2]/pars[1]^2, -1/pars[1], 1)
     sqrt(fn %*% varcov %*% fn)
    }

primpeaks.gevendse <- lapply(primpeaks.fgev, function(lst)
   sapply(lst, function(lst2)
        sapply(lst2, function(lst1) 
            if(!inherits(lst1, "try-error")) endse.fun(lst1$par.ests, lst1$varcov) else NA)))

primpeaks.gevpar$mag.g6$

## -----------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------
## Lots of plots
## -----------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------


## -----------------------------------------------------------------------------------
## Differences of the weighted and non-weighted  shape estimates:

ii <- 5
range(sapply(primpeaks.gevpar[[ii]], function(mat) mat[1,6] - mat[1,2]))
range(sapply(primpeaks.gevpar[[ii]], function(mat) mat[1,5] - mat[1,1]))
range(sapply(primpeaks.gevpar[[ii]], function(mat) mat[1,7] - mat[1,3]))
range(sapply(primpeaks.gevpar[[ii]], function(mat) mat[1,8] - mat[1,4]))
#cand2$ndata[is.element(cand2$name, names(primpeaks.gevpar))]
ii <- 6
par(mfrow = c(1,2), mar = c(3,3,2,1), mgp = c(2,1,0))
plot(log(nObs[1:7]), sapply(primpeaks.gevpar[[ii]][1:7], function(mat) mat[1,5] - mat[1,1]),
   ylim = c(-0.31, 0), main = paste(mag.names[ii], 
   ", Uniform random thinning", sep = ""), ylab = "Difference")
points(log(nObs[1:7]), sapply(primpeaks.gevpar[[ii]][1:7], function(mat) mat[1,6] - mat[1,2]),
   pch = 2, col = "deepskyblue")
points(log(nObs[1:7]), sapply(primpeaks.gevpar[[ii]][1:7], function(mat) mat[1,7] - mat[1,3]),
   pch = 3, col = "red")
points(log(nObs[1:7]), sapply(primpeaks.gevpar[[ii]][1:7], function(mat) mat[1,8] - mat[1,4]),
   pch = 4, col = "orange")
abline(h = 0)

plot(log(nObs[8:14]), sapply(primpeaks.gevpar[[ii]][8:14], function(mat) mat[1,5] - mat[1,1]),
   ylim = c(-0.31, 0), main = paste(mag.names[ii], 
   ", Nighttime random thinning", sep = ""), ylab = "Difference")
points(log(nObs[8:14]), sapply(primpeaks.gevpar[[ii]][8:14], function(mat) mat[1,6] - mat[1,2]),
   pch = 2, col = "deepskyblue")
points(log(nObs[8:14]), sapply(primpeaks.gevpar[[ii]][8:14], function(mat) mat[1,7] - mat[1,3]),
   pch = 3, col = "red")
points(log(nObs[8:14]), sapply(primpeaks.gevpar[[ii]][8:14], function(mat) mat[1,8] - mat[1,4]),
   pch = 4, col = "orange")
abline(h = 0)

## Ugly big difference. The weights themselves look okay this time, but the difference
## is similar to what I found when I made an error in the weight definition. 
## however, the difference between weighted and non-weighted appears to decrease
## when SNR decreases - the observed ts is close to a noise sequence.

## -----------------------------------------------------------------------------------
## Weighted and non-weighted shapes, scales and locations against nObs:

ii <- 6
cols.tmp <- 1:7
#cols.tmp <- 8:14
#dev.set(2)
quartz(height = 6, width = 7.2)
par(mfrow = c(3,2), mar = c(3,3,2,1), mgp = c(2,1,0))
## Shapes:
plot(log(nObs[cols.tmp])-0.03,  sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[1,1]),
   ylim = c(-0.3,0.1), main = paste(if(sum(cols.tmp == 1:7) == 7) 
   "Uniform thinning, " else "Nighttime thinning, ", mag.names[ii], 
   ", Non-weighted", sep = ""), ylab = "Non-weighted shape", 
   xlab = "Log(N)", pch = 21, col = 1, bg = 0)
segments(log(nObs[cols.tmp])-0.03, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[1,1]) + 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[1,1]),
   log(nObs[cols.tmp])-0.03, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[1,1]) - 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[1,1]),
   col = 1, lty = 1)
points(log(nObs[cols.tmp])-0.01,  sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[1,2]), 
   pch = 22, col = "deepskyblue", bg = 0)
segments(log(nObs[cols.tmp])-0.01, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[1,2]) + 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[1,2]),
   log(nObs[cols.tmp])-0.01, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[1,2]) - 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[1,2]),
   col = "deepskyblue", lty = 1)
points(log(nObs[cols.tmp])+0.01,  sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[1,3]), 
   pch = 23, col = "red", bg = 0)
segments(log(nObs[cols.tmp])+0.01, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[1,3]) + 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[1,3]),
   log(nObs[cols.tmp])+0.01, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[1,3]) - 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[1,3]),
   col = "red", lty = 1)
points(log(nObs[cols.tmp])+0.03,  sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[1,4]), 
   pch = 24, col = "orange", bg = 0)
segments(log(nObs[cols.tmp])+0.03, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[1,4]) + 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[1,4]),
   log(nObs[cols.tmp])+0.03, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[1,4]) - 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[1,4]),
   col = "orange", lty = 1)

plot(log(nObs[cols.tmp])-0.03,  sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[1,5]),
   ylim = c(-0.3,0.1), main = paste(if(cols.tmp == 1:7) 
   "Uniform thinning, " else "Nighttime thinning, ", mag.names[ii], 
   ", Weighted", sep = ""), ylab = "Weighted shape", 
   xlab = "Log(N)", pch = 21, col = 1, bg = 0)
segments(log(nObs[cols.tmp])-0.03, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[1,5]) + 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[1,5]),
   log(nObs[cols.tmp])-0.03, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[1,5]) - 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[1,5]),
   col = 1, lty = 1)
points(log(nObs[cols.tmp])-0.01,  sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[1,6]), 
   pch = 22, col = "deepskyblue", bg = 0)
segments(log(nObs[cols.tmp])-0.01, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[1,6]) + 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[1,6]),
   log(nObs[cols.tmp])-0.01, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[1,6]) - 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[1,6]),
   col = "deepskyblue", lty = 1)
points(log(nObs[cols.tmp])+0.01,  sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[1,7]), 
   pch = 23, col = "red", bg = 0)
segments(log(nObs[cols.tmp])+0.01, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[1,7]) + 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[1,7]),
   log(nObs[cols.tmp])+0.01, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[1,7]) - 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[1,7]),
   col = "red", lty = 1)
points(log(nObs[cols.tmp])+0.03,  sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[1,8]), 
   pch = 24, col = "orange", bg = 0)
segments(log(nObs[cols.tmp])+0.03, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[1,8]) + 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[1,8]),
   log(nObs[cols.tmp])+0.03, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[1,8]) - 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[1,8]),
   col = "orange", lty = 1)

## Scales
plot(log(nObs[cols.tmp])-0.03,  sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[2,1]),
   ylim = c(0,0.125), main = "", ylab = "Non-weighted scale", 
   xlab = "Log(N)", pch = 21, col = 1, bg = 0)
segments(log(nObs[cols.tmp])-0.03, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[2,1]) + 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[2,1]),
   log(nObs[cols.tmp])-0.03, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[2,1]) - 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[2,1]),
   col = 1, lty = 1)
points(log(nObs[cols.tmp])-0.01,  sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[2,2]), 
   pch = 22, col = "deepskyblue", bg = 0)
segments(log(nObs[cols.tmp])-0.01, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[2,2]) + 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[2,2]),
   log(nObs[cols.tmp])-0.01, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[2,2]) - 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[2,2]),
   col = "deepskyblue", lty = 1)
points(log(nObs[cols.tmp])+0.01,  sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[2,3]), 
   pch = 23, col = "red", bg = 0)
segments(log(nObs[cols.tmp])+0.01, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[2,3]) + 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[2,3]),
   log(nObs[cols.tmp])+0.01, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[2,3]) - 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[2,3]),
   col = "red", lty = 1)
points(log(nObs[cols.tmp])+0.03,  sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[2,4]), 
   pch = 24, col = "orange", bg = 0)
segments(log(nObs[cols.tmp])+0.03, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[2,4]) + 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[2,4]),
   log(nObs[cols.tmp])+0.03, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[2,4]) - 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[2,4]),
   col = "orange", lty = 1)

plot(log(nObs[cols.tmp])-0.03,  sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[2,5]),
   ylim = c(0,0.125), main = "", ylab = "Weighted scale", 
   xlab = "Log(N)", pch = 21, col = 1, bg = 0)
segments(log(nObs[cols.tmp])-0.03, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[2,5]) + 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[2,5]),
   log(nObs[cols.tmp])-0.03, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[2,5]) - 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[2,5]),
   col = 1, lty = 1)
points(log(nObs[cols.tmp])-0.01,  sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[2,6]), 
   pch = 22, col = "deepskyblue", bg = 0)
segments(log(nObs[cols.tmp])-0.01, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[2,6]) + 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[2,6]),
   log(nObs[cols.tmp])-0.01, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[2,6]) - 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[2,6]),
   col = "deepskyblue", lty = 1)
points(log(nObs[cols.tmp])+0.01,  sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[2,7]), 
   pch = 23, col = "red", bg = 0)
segments(log(nObs[cols.tmp])+0.01, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[2,7]) + 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[2,7]),
   log(nObs[cols.tmp])+0.01, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[2,7]) - 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[2,7]),
   col = "red", lty = 1)
points(log(nObs[cols.tmp])+0.03,  sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[2,8]), 
   pch = 24, col = "orange", bg = 0)
segments(log(nObs[cols.tmp])+0.03, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[2,8]) + 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[2,8]),
   log(nObs[cols.tmp])+0.03, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[2,8]) - 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[2,8]),
   col = "orange", lty = 1)

## Location:
plot(log(nObs[cols.tmp])-0.03,  sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[3,1]),
   ylim = c(0,0.5), main = "", ylab = "Non-weighted location", 
   xlab = "Log(N)", pch = 21, col = 1, bg = 0)
segments(log(nObs[cols.tmp])-0.03, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[3,1]) + 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[3,1]),
   log(nObs[cols.tmp])-0.03, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[3,1]) - 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[3,1]),
   col = 1, lty = 1)
points(log(nObs[cols.tmp])-0.01,  sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[3,2]), 
   pch = 22, col = "deepskyblue", bg = 0)
segments(log(nObs[cols.tmp])-0.01, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[3,2]) + 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[3,2]),
   log(nObs[cols.tmp])-0.01, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[3,2]) - 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[3,2]),
   col = "deepskyblue", lty = 1)
points(log(nObs[cols.tmp])+0.01,  sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[3,3]), 
   pch = 23, col = "red", bg = 0)
segments(log(nObs[cols.tmp])+0.01, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[3,3]) + 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[3,3]),
   log(nObs[cols.tmp])+0.01, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[3,3]) - 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[3,3]),
   col = "red", lty = 1)
points(log(nObs[cols.tmp])+0.03,  sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[3,4]), 
   pch = 24, col = "orange", bg = 0)
segments(log(nObs[cols.tmp])+0.03, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[3,4]) + 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[3,4]),
   log(nObs[cols.tmp])+0.03, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[3,4]) - 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[3,4]),
   col = "orange", lty = 1)

plot(log(nObs[cols.tmp])-0.03,  sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[3,5]),
   ylim = c(0,0.5), main = "", ylab = "Weighted location", 
   xlab = "Log(N)", pch = 21, col = 1, bg = 0)
segments(log(nObs[cols.tmp])-0.03, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[3,5]) + 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[3,5]),
   log(nObs[cols.tmp])-0.03, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[3,5]) - 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[3,5]),
   col = 1, lty = 1)
points(log(nObs[cols.tmp])-0.01,  sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[3,6]), 
   pch = 22, col = "deepskyblue", bg = 0)
segments(log(nObs[cols.tmp])-0.01, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[3,6]) + 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[3,6]),
   log(nObs[cols.tmp])-0.01, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[3,6]) - 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[3,6]),
   col = "deepskyblue", lty = 1)
points(log(nObs[cols.tmp])+0.01,  sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[3,7]), 
   pch = 23, col = "red", bg = 0)
segments(log(nObs[cols.tmp])+0.01, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[3,7]) + 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[3,7]),
   log(nObs[cols.tmp])+0.01, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[3,7]) - 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[3,7]),
   col = "red", lty = 1)
points(log(nObs[cols.tmp])+0.03,  sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[3,8]), 
   pch = 24, col = "orange", bg = 0)
segments(log(nObs[cols.tmp])+0.03, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[3,8]) + 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[3,8]),
   log(nObs[cols.tmp])+0.03, 
   sapply(primpeaks.gevpar[[ii]][cols.tmp], function(mat) mat[3,8]) - 
   2*sapply(primpeaks.gevse[[ii]][cols.tmp], function(mat) mat[3,8]),
   col = "orange", lty = 1)

## Nighttime/uniform thinning (no difference between the two):
## very strong signal (mag.g) + very few observations: very negative shape 
##     compensated with large location and large scale. The shape goes towards 
##     0 when nObs increases. Dependent sampling gives more negative
##     shapes and larger locations than the others for nObs <= 200, but
##     approximately the same scale. For nObs > 200, shape estimates from 
##     dependent sampling is approximately same as the others.
## less strong (mag.g2, .g3) signal: the overall appearance is the same, but
##     the shapes become less negative. However, the difference between
##     dependent sampling and the others is still present, the same
##     way and elsewhere too: more negative shapes and larger locations
##     for  nObs <= 200.
## weak signal (mag.g4-): the nObs-dependence of the parameters is much less
##     clear, they vary a lot (shape for nObs = 50 is almost the highest of
##     all nObs). Overall, shape goes still less negative with increasing nObs, 
##     scale and location decrease. The difference between the dependent
##     sampling and the others seems to decrease, but only slightly, and
##     not sure if that's just some fluctuation due to the particular
##     realized time thinning, or indeed something that should happen.
## even weakening signal (mag.g5): most of the nObs-dependence at very
##     low nObs seems to be due to simple natural estimation fluctuations,
##     because now the shape at nObs = 50 is the most negative. The 
##     dependent sampling yields still the most negative shapes and 
##     the highest locations; no difference in the scale for weighted
##     estimation, some (smaller scale) in non-weighted.
## almost no signal (mag.g6): from the weighted shapes, only the dependent
##     sampling shows the nObs-dependence seen until now. The non-weighted
##     ones all do. 
## All this is present both for weighted and non-weighted estimation,
## but the non-weighted is much closer to 0 for all type of frequency
## sampling, dependent or others.
## For non-weighted, the difference between dependent or other selection
## is less pronounced. The non-weighted seems to show as well more consistent, 
## less varying dependence on nObs. The weighted estimates, especially for
## low nObs, show much more fluctuations.

## Summary:
##  - Patterns of nObs-dependence: general pattern:
##    ° shape: more negative for small nObs, increasing towards 0 when nObs increase
##    ° scale: steadily decreasing when nObs increase
##    ° location: steadily decreasing when nObs increase
##    * Some difference in this overall pattern for weighted and non-weighted estimation:
##      ° shape is much more negative for weighted spectrum estimation than for non-
##        weighted, for non-weighted it is quite close to 0 (no endpoint)
##      ° scale is smaller for non-weighted than for weighted
##      ° location is smaller for non-weighted than for weighted
##      The differences between weighted and non-weighted estimates of scale and location
##      tend to decrease a little when SNR decreases, but not much. Shape is too
##      unstable to tell.
##    * Difference between dependent and other frequency sampling methods:
##      ° shape: the estimate from dependent sampling is much more negative than the
##        others when using weighted periodogram estimation for nObs <= 200; this 
##        difference disappears for higher nObs using weighted estimation, and is
##        not present in non-weighted estimation
##      ° scale: not much difference, the dependent sampling estimate is a bit smaller
##        only for nObs = 25
##      ° location: real persistent difference, for low nObs it is much higher than
##        the other samplings, both for weighted and nonweighted estimates too; the 
##        difference disappears for nObs > 1000

## What does all this tell me, if anything at all?

##  - Weighted/nonweighted: locations are higher for weighted than for non-weighted,
##    shapes are more negative ==> the GEV is shifted upwards, and the probability
##    to find very high periodogram values is cut more sharply, no probability to
##    find larger values than a limit. Non-weighted seems therefore much more prone
##    to producing false high peaks!
##  - Dependent sampling/others:
##    ° Shift of location when using dependent sampling: logical, I go even upwards in a 
##    window, I can obtain only larger values, the location will be shifted upwards.
##    ° Scale: for nObs = 25, it is slightly smaller with dependent sampling, but not
##    for any other sampling or larger nObs. Okay, the extension ("scatter") of the 
##    GEV is the same, but pushed upwards according to the location.
##    ° Shape: no difference between dependent sampling and the others when using
##    nonweighted spectral estimation, all are close to 0, so there is still relatively
##    high probability of getting high peaks in a non-weighted estimation, whatever 
##    the frequency sampling method is, even if we climbed upwards within a peak.
##    The GEV distribution has no limit, it stays close to exponential tail even
##    with this shift! 

## Sg I don't understand: why the location and the scale decreases with increasing 
## nObs?

## mag.g2 and mag.g6 are good illustrations for the talk, I can present



## -----------------------------------------------------------------------------------
## Estimates against signal-to-noise ratio:

smp.tmp <- "sim2000"
sapply(primpeaks.gevpar, function(lst1, smp.tmp) lst1[[smp.tmp]][1,2], smp.tmp = smp.tmp)

## For uniform thinning:
quartz()
par(mfrow = c(3,2), mar = c(3,3,2,1), mgp = c(2,1,0))

plot(amplitudes, 
   sapply(primpeaks.gevpar, function(lst1) lst1[["sim200"]][1,2]), type = "n",
   ylim = c(-0.35,0.1), main = paste("Non-weighted, dependent sampling", 
   sep = ""), ylab = "Non-weighted shape", 
   xlab = "Signal amplitude", pch = 15, col = colorvec[1])
for(jj in 1:7)
    lines(amplitudes, 
       sapply(primpeaks.gevpar, function(lst1, smp.tmp) lst1[[smp.tmp]][1,2], 
       smp.tmp = colnames(FirstStar)[jj+4]), type = "o", pch = 16, 
       col = colorvec[jj], lwd = 2)
legend(0.025, -0.2, legend = c("nObs = 2000","nObs = 1000","nObs = 500",
   "nObs = 200"), lty = 1, col = colorvec[1:4],
   pch = 16, bty = "n")
legend(0.14, -0.2, legend = c("nObs = 100","nObs = 50","nObs = 25"), 
   lty = 1, col = colorvec[5:7], pch = 16, bty = "n")

plot(amplitudes, 
   sapply(primpeaks.gevpar, function(lst1, smp.tmp) lst1[[smp.tmp]][1,6], 
   smp.tmp = "sim200"), type = "n",
   ylim = c(-0.35,0.1), main = paste("Weighted, dependent sampling", 
   sep = ""), ylab = "Weighted shape", 
   xlab = "Signal amplitude", pch = 15, col = colorvec[1])
for(jj in 1:7)
lines(amplitudes, 
   sapply(primpeaks.gevpar, function(lst1, smp.tmp) lst1[[smp.tmp]][1,6], 
   smp.tmp = colnames(FirstStar)[jj+4]), type = "o", pch = 16, 
   col = colorvec[jj], lwd = 2)

plot(amplitudes, 
   sapply(primpeaks.gevpar, function(lst1, smp.tmp) lst1[[smp.tmp]][2,2], 
   smp.tmp = "sim200"), type = "n",
   ylim = c(0,0.11), main = "", ylab = "Non-weighted scale", 
   xlab = "Signal amplitude", pch = 15, col = colorvec[1])
for(jj in 1:7)
lines(amplitudes, 
   sapply(primpeaks.gevpar, function(lst1, smp.tmp) lst1[[smp.tmp]][2,2], 
   smp.tmp = colnames(FirstStar)[jj+4]), type = "o", pch = 16, 
   col = colorvec[jj], lwd = 2)

plot(amplitudes, 
   sapply(primpeaks.gevpar, function(lst1, smp.tmp) lst1[[smp.tmp]][2,6], 
   smp.tmp = "sim200"), type = "n",
   ylim = c(0,0.11), main = "", ylab = "Non-weighted scale", 
   xlab = "Signal amplitude", pch = 15, col = colorvec[1])
for(jj in 1:7)
lines(amplitudes, 
   sapply(primpeaks.gevpar, function(lst1, smp.tmp) lst1[[smp.tmp]][2,6], 
   smp.tmp = colnames(FirstStar)[jj+4]), type = "o", pch = 16, 
   col = colorvec[jj], lwd = 2)

plot(amplitudes, 
   sapply(primpeaks.gevpar, function(lst1, smp.tmp) lst1[[smp.tmp]][3,2], 
   smp.tmp = "sim200"), type = "n",
   ylim = c(0,0.6), main = "", ylab = "Non-weighted location", 
   xlab = "Signal amplitude", pch = 15, col = colorvec[1])
for(jj in 1:7)
lines(amplitudes, 
   sapply(primpeaks.gevpar, function(lst1, smp.tmp) lst1[[smp.tmp]][3,2], 
   smp.tmp = colnames(FirstStar)[jj+4]), type = "o", pch = 16, 
   col = colorvec[jj], lwd = 2)

plot(amplitudes, 
   sapply(primpeaks.gevpar, function(lst1, smp.tmp) lst1[[smp.tmp]][3,6], 
   smp.tmp = "sim200"), type = "n",
   ylim = c(0,0.6), main = "", ylab = "Non-weighted location", 
   xlab = "Signal amplitude", pch = 15, col = colorvec[1])
for(jj in 1:7)
lines(amplitudes, 
   sapply(primpeaks.gevpar, function(lst1, smp.tmp) lst1[[smp.tmp]][3,6], 
   smp.tmp = colnames(FirstStar)[jj+4]), type = "o", pch = 16, 
   col = colorvec[jj], lwd = 2)


## -----------------------------------------------------------------------------------
## Weighted and non-weighted endpoints against nObs:

ii <- 3
#dev.set(2)
quartz(height = 3, width = 7.2)
par(mfrow = c(1,2), mar = c(3,3,2,1), mgp = c(2,1,0))
plot(log(nObs[1:7])-0.035, primpeaks.gevend[[ii]][1,1:7],
   ylim = c(0, 6), main = paste(mag.names[ii], 
   ", Non-weighted", sep = ""), ylab = "Non-weighted endpoint", 
   xlab = "Log(N)", pch = 21, col = 1, bg = 0)
segments(log(nObs[1:7])-0.035, 
   primpeaks.gevend[[ii]][1,1:7] + 2*primpeaks.gevendse[[ii]][1,1:7],
   log(nObs[1:7])-0.035, 
   primpeaks.gevend[[ii]][1,1:7] - 2*primpeaks.gevendse[[ii]][1,1:7],
   col = 1, lty = 2)
points(log(nObs[1:7])-0.025, primpeaks.gevend[[ii]][2,1:7], 
   pch = 22, col = "deepskyblue", bg = 0)
segments(log(nObs[1:7])-0.025, 
   primpeaks.gevend[[ii]][2,1:7] + 2*primpeaks.gevendse[[ii]][2,1:7],
   log(nObs[1:7])-0.025, 
   primpeaks.gevend[[ii]][2,1:7] - 2*primpeaks.gevendse[[ii]][2,1:7],
   col = "deepskyblue", lty = 2)
points(log(nObs[1:7])-0.015, primpeaks.gevend[[ii]][3,1:7], 
   pch = 23, col = "red", bg = 0)
segments(log(nObs[1:7])-0.015, 
   primpeaks.gevend[[ii]][3,1:7] + 2*primpeaks.gevendse[[ii]][3,1:7],
   log(nObs[1:7])-0.015, 
   primpeaks.gevend[[ii]][3,1:7] - 2*primpeaks.gevendse[[ii]][3,1:7],
   col = "red", lty = 2)
points(log(nObs[1:7])-0.005, primpeaks.gevend[[ii]][4,1:7],
   pch = 24, col = "orange", bg = 0)
segments(log(nObs[1:7])-0.005, 
   primpeaks.gevend[[ii]][4,1:7] + 2*primpeaks.gevendse[[ii]][4,1:7],
   log(nObs[1:7])-0.005, 
   primpeaks.gevend[[ii]][4,1:7] - 2*primpeaks.gevendse[[ii]][4,1:7],
   col = "orange", lty = 2)
points(log(nObs[8:14])+0.005, primpeaks.gevend[[ii]][1,8:14],
   pch = 21, col = "black", bg = "black")
segments(log(nObs[8:14])+0.005, 
   primpeaks.gevend[[ii]][1,8:14] + 2*primpeaks.gevendse[[ii]][1,8:14],
   log(nObs[8:14])+0.005, 
   primpeaks.gevend[[ii]][1,8:14] - 2*primpeaks.gevendse[[ii]][1,8:14],
   col = 1, lty = 1)
points(log(nObs[8:14])+0.015, primpeaks.gevend[[ii]][2,8:14],
   pch = 22, col = "deepskyblue", bg = "deepskyblue")
segments(log(nObs[8:14])+0.015, 
   primpeaks.gevend[[ii]][2,8:14] + 2*primpeaks.gevendse[[ii]][2,8:14],
   log(nObs[8:14])+0.015, 
   primpeaks.gevend[[ii]][2,8:14] - 2*primpeaks.gevendse[[ii]][2,8:14],
   col = "deepskyblue", lty = 1)
points(log(nObs[8:14])+0.025, primpeaks.gevend[[ii]][3,8:14],
   pch = 23, col = "red", bg = "red")
segments(log(nObs[8:14])+0.025, 
   primpeaks.gevend[[ii]][3,8:14] + 2*primpeaks.gevendse[[ii]][3,8:14],
   log(nObs[8:14])+0.025, 
   primpeaks.gevend[[ii]][3,8:14] - 2*primpeaks.gevendse[[ii]][3,8:14],
   col = "red", lty = 1)
points(log(nObs[8:14])+0.035, primpeaks.gevend[[ii]][4,8:14],
   pch = 24, col = "orange", bg = "orange")
segments(log(nObs[8:14])+0.035, 
   primpeaks.gevend[[ii]][4,8:14] + 2*primpeaks.gevendse[[ii]][4,8:14],
   log(nObs[8:14])+0.035, 
   primpeaks.gevend[[ii]][4,8:14] - 2*primpeaks.gevendse[[ii]][4,8:14],
   col = "orange", lty = 1)
#abline(v = log(nObs[8:14]), col = "grey70", lty = 3)

plot(log(nObs[1:7])-0.035, primpeaks.gevend[[ii]][5,1:7],
   ylim = c(0, 6), main = paste(mag.names[ii], 
   ", Weighted", sep = ""), ylab = "Weighted endpoint", 
   xlab = "Log(N)", pch = 21, col = 1, bg = 0)
segments(log(nObs[1:7])-0.035, 
   primpeaks.gevend[[ii]][5,1:7] + 2*primpeaks.gevendse[[ii]][5,1:7],
   log(nObs[1:7])-0.035, 
   primpeaks.gevend[[ii]][5,1:7] - 2*primpeaks.gevendse[[ii]][5,1:7],
   col = 1, lty = 2)
points(log(nObs[1:7])-0.025, primpeaks.gevend[[ii]][6,1:7], 
   pch = 22, col = "deepskyblue", bg = 0)
segments(log(nObs[1:7])-0.025, 
   primpeaks.gevend[[ii]][6,1:7] + 2*primpeaks.gevendse[[ii]][6,1:7],
   log(nObs[1:7])-0.025, 
   primpeaks.gevend[[ii]][6,1:7] - 2*primpeaks.gevendse[[ii]][6,1:7],
   col = "deepskyblue", lty = 2)
points(log(nObs[1:7])-0.015, primpeaks.gevend[[ii]][7,1:7], 
   pch = 23, col = "red", bg = 0)
segments(log(nObs[1:7])-0.015, 
   primpeaks.gevend[[ii]][7,1:7] + 2*primpeaks.gevendse[[ii]][7,1:7],
   log(nObs[1:7])-0.015, 
   primpeaks.gevend[[ii]][7,1:7] - 2*primpeaks.gevendse[[ii]][7,1:7],
   col = "red", lty = 2)
points(log(nObs[1:7])-0.005, primpeaks.gevend[[ii]][8,1:7],
   pch = 24, col = "orange", bg = 0)
segments(log(nObs[1:7])-0.005, 
   primpeaks.gevend[[ii]][8,1:7] + 2*primpeaks.gevendse[[ii]][8,1:7],
   log(nObs[1:7])-0.005, 
   primpeaks.gevend[[ii]][8,1:7] - 2*primpeaks.gevendse[[ii]][8,1:7],
   col = "orange", lty = 2)
points(log(nObs[8:14])+0.005, primpeaks.gevend[[ii]][5,8:14],
   pch = 21, col = "black", bg = "black")
segments(log(nObs[8:14])+0.005, 
   primpeaks.gevend[[ii]][5,8:14] + 2*primpeaks.gevendse[[ii]][5,8:14],
   log(nObs[8:14])+0.005, 
   primpeaks.gevend[[ii]][5,8:14] - 2*primpeaks.gevendse[[ii]][5,8:14],
   col = 1, lty = 1)
points(log(nObs[8:14])+0.015, primpeaks.gevend[[ii]][6,8:14],
   pch = 22, col = "deepskyblue", bg = "deepskyblue")
segments(log(nObs[8:14])+0.015, 
   primpeaks.gevend[[ii]][6,8:14] + 2*primpeaks.gevendse[[ii]][6,8:14],
   log(nObs[8:14])+0.015, 
   primpeaks.gevend[[ii]][6,8:14] - 2*primpeaks.gevendse[[ii]][6,8:14],
   col = "deepskyblue", lty = 1)
points(log(nObs[8:14])+0.025, primpeaks.gevend[[ii]][7,8:14],
   pch = 23, col = "red", bg = "red")
segments(log(nObs[8:14])+0.025, 
   primpeaks.gevend[[ii]][7,8:14] + 2*primpeaks.gevendse[[ii]][7,8:14],
   log(nObs[8:14])+0.025, 
   primpeaks.gevend[[ii]][7,8:14] - 2*primpeaks.gevendse[[ii]][7,8:14],
   col = "red", lty = 1)
points(log(nObs[8:14])+0.035, primpeaks.gevend[[ii]][8,8:14],
   pch = 24, col = "orange", bg = "orange")
segments(log(nObs[8:14])+0.035, 
   primpeaks.gevend[[ii]][8,8:14] + 2*primpeaks.gevendse[[ii]][8,8:14],
   log(nObs[8:14])+0.035, 
   primpeaks.gevend[[ii]][8,8:14] - 2*primpeaks.gevendse[[ii]][8,8:14],
   col = "orange", lty = 1)
#abline(v = log(nObs[8:14]), col = "grey70", lty = 3)
## Something contrary to the impression from the GEV parameter plotting: 
##  - the endpoints from the weighted estimation have much smaller errors,
##    they show a much steadier tendency against the number of observations,
##    and for a given nObs, they are almost equal in all types of esti-
##    mation (the GEV parameters showed more scatter and less steady
##    behaviour in the weighted version than in the non-weighted one, so
##    it was practically the contrary)
##  - the endpoint of the distribution decreases when the number of
##    the data points increase, both with weighted or non-weighted least
##    squares.

## -----------------------------------------------------------------------------------
## The time series, the periodograms and the maximum of the maxima:

quartz(height = 7.5, width = 8)
colnnn <- "sim25"
mag.tmp <- "mag.g5"
magerr.tmp <- "mag.error.g5"
#mycol <- "orange"
par(mfrow = c(3,1), mar = c(3,3,2,1))
plot(FirstStar[FirstStar[, colnnn] == 1, "time"], 
  FirstStar[FirstStar[, colnnn] == 1, mag.tmp], xlim = c(0,25),
  ylim = c(-0.5,0.5), pch = 16, cex = 0.5)
segments(FirstStar[FirstStar[, colnnn] == 1, "time"], 
  FirstStar[FirstStar[, colnnn] == 1, mag.tmp] + 2*FirstStar[FirstStar[, colnnn] == 1, magerr.tmp],
  FirstStar[FirstStar[, colnnn] == 1, "time"], 
  FirstStar[FirstStar[, colnnn] == 1, mag.tmp] - 2*FirstStar[FirstStar[, colnnn] == 1, magerr.tmp])
plot(FirstStar.pgram[[mag.tmp]]$freq, FirstStar.pgram[[mag.tmp]][,colnnn], type = "n",
  xlab = "Frequency", ylab = "Chi-squared", main = paste(colnnn, mag.tmp), 
#  ylim = c(0,max(max(FirstStar.bootI[[mag.tmp]][[colnnn]]["dep.max",], na.rm = TRUE), 
#  max(FirstStar.pgram[[mag.tmp]][,colnnn], na.rm = TRUE))),
  ylim = c(0,1),
  xlim = c(0,100))
abline(v = F, col = "violetred", lwd = 2)
abline(h = max(FirstStar.bootI[[mag.tmp]][[colnnn]]["dep.max",]), col = "violetred")
points(FirstStar.pgram[[mag.tmp]]$freq, FirstStar.pgram[[mag.tmp]][,colnnn], type = "h")
plot(FirstStar.wpgram[[mag.tmp]]$freq, FirstStar.wpgram[[mag.tmp]][,colnnn], type = "n",
  xlab = "Frequency", ylab = "Chi-squared", main = paste(colnnn, mag.tmp), 
#  ylim = c(0,max(max(FirstStar.bootI[[mag.tmp]][[colnnn]]["wdep.max",], na.rm = TRUE), 
#  max(FirstStar.wpgram[[mag.tmp]][,colnnn], na.rm = TRUE))),
  ylim = c(0,1),
  xlim = c(0,100))
abline(v = F, col = "violetred", lwd = 2)
points(FirstStar.wpgram[[mag.tmp]]$freq, FirstStar.wpgram[[mag.tmp]][,colnnn], type = "h")
abline(h = max(FirstStar.bootI[[mag.tmp]][[colnnn]]["wdep.max",]), col = "violetred")





## -----------------------------------------------------------------------------------
## Diagnostic plots

length( primpeaks.fgev[[1]])

primpeaks.quan <- lapply(primpeaks.gevpar, function(lst)
  lapply(lst, function(mat)
   apply(mat, MAR = 2, function(parvec) 
    c(q.95 = qgev(0.95, xi = parvec[1], mu = parvec[3], sigma = parvec[2]),
      q.99 = qgev(0.99, xi = parvec[1], mu = parvec[3], sigma = parvec[2]),
      q.999 = qgev(0.999, xi = parvec[1], mu = parvec[3], sigma = parvec[2]),
      q.9995 = qgev(0.9995, xi = parvec[1], mu = parvec[3], sigma = parvec[2]),
      q.9999 = qgev(0.9999, xi = parvec[1], mu = parvec[3], sigma = parvec[2])))))
## First look only at the numbers:
##  - weighted quantiles are in general higher than the corresponding non-weighted ones
##  - the dependent quantiles are in general higher than the others with the
##    same weighting strategy.


## Quantile plots for the weighted dependent case:

ii <- 6
quartz(height = 5, width = 10)
par(mfrow = c(2,4),  mar = c(3,3,2,1), mgp = c(2,1,0))
for(jj in 1:7)
   {
#   	devAskNewPage(ask = TRUE)
    plot(qgev(ppoints(1000), xi = primpeaks.gevpar[[ii]][[jj]][1,6],
         sigma = primpeaks.gevpar[[ii]][[jj]][2,6], 
         mu = primpeaks.gevpar[[ii]][[jj]][3,6]),
         sort(FirstStar.bootI[[ii]][[jj]][6,]), pch = 16, cex = 0.4,
         main = paste(names(primpeaks.gevpar)[ii], ", ",
         names(primpeaks.gevpar[[ii]])[jj], ", wdep", sep = ""),
         xlab = "Theoretical quantiles", ylab = "Sample")
    abline(c(0,1), col = "red")
  	}

## The majority seems to be good model. More or less strong systematic 
## deviation from the fitted GEV is present in:
## nights200 mag.g6
## nights200 mag.g5
## nights50 mag.g5
## nights500 mag.g4
## nights200 mag.g4
## nights25 mag.g4
## nights25 mag.g3
## nights25 mag.g
## sim1000 mag.g
## sim200 mag.g
## sim25 mag.g
## sim1000 mag.g3
## sim25 mag.g3
## sim1000 mag.g4
## sim50 mag.g4
## sim200 mag.g5
## around 20% of all the 7*6*2 cases...

## Quantile plots for the non-weighted dependent case:

ii <- 6
quartz(height = 5, width = 10)
par(mfrow = c(2,4),  mar = c(3,3,2,1), mgp = c(2,1,0))
for(jj in 1:7)
   {
#   	devAskNewPage(ask = TRUE)
    plot(qgev(ppoints(1000), xi = primpeaks.gevpar[[ii]][[jj]][1,2],
         sigma = primpeaks.gevpar[[ii]][[jj]][2,2], 
         mu = primpeaks.gevpar[[ii]][[jj]][3,2]),
         sort(FirstStar.bootI[[ii]][[jj]][2,]), pch = 16, cex = 0.4,
         main = paste(names(primpeaks.gevpar)[ii], ", ",
         names(primpeaks.gevpar[[ii]])[jj], ", dep", sep = ""),
         xlab = "Theoretical quantiles", ylab = "Sample")
    abline(c(0,1), col = "pink")
  	}

## The majority seems to be good model. More or less strong systematic 
## deviation from the fitted GEV is present in:
## nights1000 mag.g4
## nights50 mag.g4
## nights200 mag.g2
## nights2000 mag.g
## nights100 mag.g
## sim100 mag.g2
## sim200 mag.g2
## sim2000 mag.g5
## sim100 mag.g5
## sim25 mag.g5
## ... and less strong deviations than in the weighted case.
## The impression is that these are in general better models than the
## weighted dependent models.


## Quantile plots for the non-weighted random-fr. case:

ii <- 1
quartz(height = 5, width = 10)
par(mfrow = c(2,4),  mar = c(3,3,2,1), mgp = c(2,1,0))
for(jj in 1:7)
   {
#   	devAskNewPage(ask = TRUE)
    plot(qgev(ppoints(1000), xi = primpeaks.gevpar[[ii]][[jj]][1,4],
         sigma = primpeaks.gevpar[[ii]][[jj]][2,4], 
         mu = primpeaks.gevpar[[ii]][[jj]][3,4]),
         sort(FirstStar.bootI[[ii]][[jj]][4,]), pch = 16, cex = 0.4,
         main = paste(names(primpeaks.gevpar)[ii], ", ",
         names(primpeaks.gevpar[[ii]])[jj], ", rand", sep = ""),
         xlab = "Theoretical quantiles", ylab = "Sample")
    abline(c(0,1), col = "blue")
  	}

## The majority seems to be good model. More or less strong systematic 
## deviation from the fitted GEV is present in:
## nights1000 mag.g6
## nights2000 mag.g5
## nights1000 mag.g5
## nights2000 mag.g4
## sim200 mag.g3
## These are almost as nice as the non-weighted dependent sampling. In 
## their tail, there is a bit more flapping-around going on. I cannot
## compare directly to the weighted dependent one, my memories are too 
## vague.


## Quantile plots for the weighted random-fr. case:

ii <- 1
quartz(height = 5, width = 10)
par(mfrow = c(2,4),  mar = c(3,3,2,1), mgp = c(2,1,0))
for(jj in 1:7)
   {
#   	devAskNewPage(ask = TRUE)
    plot(qgev(ppoints(1000), xi = primpeaks.gevpar[[ii]][[jj]][1,8],
         sigma = primpeaks.gevpar[[ii]][[jj]][2,8], 
         mu = primpeaks.gevpar[[ii]][[jj]][3,8]),
         sort(FirstStar.bootI[[ii]][[jj]][8,]), pch = 16, cex = 0.4,
         main = paste(names(primpeaks.gevpar)[ii], ", ",
         names(primpeaks.gevpar[[ii]])[jj], ", wrand", sep = ""),
         xlab = "Theoretical quantiles", ylab = "Sample")
    abline(c(0,1), col = "deepskyblue")
  	}

## The majority seems to be good model. More or less strong systematic 
## deviation from the fitted GEV is present in:
## nights100 mag.g3
## nights500 mag.g4
## nights200 mag.g4
## nights200 mag.g5
## nights50 mag.g5
## nights200 mag.g6
## sim200 mag.g6
## sim25 mag.g6
## sim200 mag.g5
## sim1000 mag.g4
## sim200 mag.g
## Approx. the same goodness as the wdep case, maybe a bit better fits?



col.tmp <- c("brown", "violetred", "deeppink", "gold")

## Huge quantile plots to compare the weighted and non-weighted estimates:

## A common range for all sample maxima can be c(0,0.9)
## The weighted estimates give always much larger spans. Why? Is this natural,
## or again it's me who screwed it up?
jj <- 2
for(jj in 1:14)
{
quartz(height = 4.2, width = 12.5)
par(mfcol = c(2,6),  mar = c(3,3,2,1), mgp = c(2,1,0))
#axl.tmp <- range(c(unlist(FirstStar.bootI[[1]][[jj]][1:4,]),
#   unlist(FirstStar.bootI[[2]][[jj]][1:4,]),
#   unlist(FirstStar.bootI[[3]][[jj]][1:4,]),
#   unlist(FirstStar.bootI[[4]][[jj]][1:4,]),
#   unlist(FirstStar.bootI[[5]][[jj]][1:4,]),
#   unlist(FirstStar.bootI[[6]][[jj]][1:4,])))
#waxl.tmp <- range(c(unlist(FirstStar.bootI[[1]][[jj]][5:8,]),
#   unlist(FirstStar.bootI[[2]][[jj]][5:8,]),
#   unlist(FirstStar.bootI[[3]][[jj]][5:8,]),
#   unlist(FirstStar.bootI[[4]][[jj]][5:8,]),
#   unlist(FirstStar.bootI[[5]][[jj]][5:8,]),
#   unlist(FirstStar.bootI[[6]][[jj]][5:8,])))
for(ii in 1:6)
   {
   	axl.tmp <- range(unlist(FirstStar.bootI[[ii]][[jj]][1:4,]))
    plot(sort(FirstStar.bootI[[1]][[1]][8,]),
        sort(FirstStar.bootI[[1]][[1]][8,]), ylim = axl.tmp,
        xlim = axl.tmp, main = paste(names(primpeaks.gevpar)[ii], ", ",
    	names(primpeaks.gevpar[[ii]])[jj], ", non-weighted", sep = ""),
    	xlab = "Theoretical quantiles", ylab = "Sample", type = "n")
	for(kk in 1:4)
   		{
#   	devAskNewPage(ask = TRUE)
   		 abline(c(0,1), col = "black")
    	 points(qgev(ppoints(1000), xi = primpeaks.gevpar[[ii]][[jj]][1,kk],
         	sigma = primpeaks.gevpar[[ii]][[jj]][2,kk], 
         	mu = primpeaks.gevpar[[ii]][[jj]][3,kk]),
         	sort(FirstStar.bootI[[ii]][[jj]][kk,]), pch = 16, cex = 1-0.12*kk,
         	col = col.tmp[kk])
 	 	}
	legend(axl.tmp[1], axl.tmp[2], legend = c("Syst", "Dep", "Fou", "Rand"), col = col.tmp,
    	pt.cex = 1-0.12*(1:4), pch = 16, bty = "n")

   	waxl.tmp <- range(unlist(FirstStar.bootI[[ii]][[jj]][5:8,]))
	plot(sort(FirstStar.bootI[[1]][[1]][8,]),
    	sort(FirstStar.bootI[[1]][[1]][8,]), ylim = waxl.tmp,
    	xlim = waxl.tmp, main = paste(names(primpeaks.gevpar)[ii], ", ",
    	names(primpeaks.gevpar[[ii]])[jj], ", weighted", sep = ""),
    	xlab = "Theoretical quantiles", ylab = "Sample", type = "n")
	for(kk in 1:4)
	   {
#   	devAskNewPage(ask = TRUE)
	    abline(c(0,1), col = "black")
	    points(qgev(ppoints(1000), xi = primpeaks.gevpar[[ii]][[jj]][1,kk+4],
	         sigma = primpeaks.gevpar[[ii]][[jj]][2,kk+4], 
	         mu = primpeaks.gevpar[[ii]][[jj]][3,kk+4]),
	         sort(FirstStar.bootI[[ii]][[jj]][kk+4,]), pch = 16, cex = 1-0.12*kk,
	         col = col.tmp[kk])
	  	}
	legend(waxl.tmp[1], waxl.tmp[2], legend = c("Syst", "Dep", "Fou", "Rand"), col = col.tmp,
    	pt.cex = 1-0.12*(1:4), pch = 16, bty = "n")
   }
}

## Impression: 
## ugly fits: g4,g5 nights50 weighted   (lower tail)
##       g2,g3 nights100 weighted  (lower tail)
##       g4 nights100 weighted  (upper tail)
##       g4 nights200 weighted  (lower tail)
##       g5 nights200 weighted  (upper tail)
##       g4 nights500 weighted  (lower tail)
##       g4 nights1000 weighted  (lower tail) ??? so many mag.g4, with all night thinnings?
##       g5 sim200 weighted  (lower tail)
## All these more serious distortions occur in the weighted estimates.
## Apart from these cases, the modelling of the upper tail do not show 
## great quality differences. 




## -----------------------------------------------------------------------------------
## The periodograms and the quantiles of the maxima:

quartz(height = 7.5, width = 8)
colnnn <- "nights25"
mag.tmp <- "mag.g6"
magerr.tmp <- "mag.error.g6"
#mycol <- "orange"
par(mfrow = c(4,2), mar = c(3,3,2,1))
for(ii in 1:4)
   {
    plot(FirstStar.pgram[[mag.tmp]]$freq, 
       FirstStar.pgram[[mag.tmp]][,colnnn], type = "n",
       xlab = "Frequency", ylab = "Chi-squared", 
       main = paste(colnnn, mag.tmp, nSamplings[ii]), 
       ylim = c(0,max(max(FirstStar.bootI[[mag.tmp]][[colnnn]][nSamplings[ii],], na.rm = TRUE), 
       max(FirstStar.pgram[[mag.tmp]][,colnnn], na.rm = TRUE),
       primpeaks.quan[[mag.tmp]][[colnnn]][5,ii+4])),
#  ylim = c(0,0.2),
       xlim = c(0,100))
    abline(v = F, col = "violetred", lty = 2)
    points(FirstStar.pgram[[mag.tmp]]$freq, FirstStar.pgram[[mag.tmp]][,colnnn], type = "h")
    abline(h = primpeaks.quan[[mag.tmp]][[colnnn]][,ii], col = c("blue", "deepskyblue",
       "green", "orange", "red"))
    plot(FirstStar.wpgram[[mag.tmp]]$freq, 
       FirstStar.wpgram[[mag.tmp]][,colnnn], type = "n",
       xlab = "Frequency", ylab = "Chi-squared", 
       main = paste(colnnn, mag.tmp, nSamplings[ii+4]), 
       ylim = c(0,max(max(FirstStar.bootI[[mag.tmp]][[colnnn]][nSamplings[ii+4],], na.rm = TRUE), 
       max(FirstStar.wpgram[[mag.tmp]][,colnnn], na.rm = TRUE),
       primpeaks.quan[[mag.tmp]][[colnnn]][5,ii+4])),
#  ylim = c(0,0.2),
       xlim = c(0,100))
    abline(v = F, col = "violetred", lty = 2)
    points(FirstStar.wpgram[[mag.tmp]]$freq, FirstStar.wpgram[[mag.tmp]][,colnnn], type = "h")
    abline(h = primpeaks.quan[[mag.tmp]][[colnnn]][,ii+4], col = c("blue", "deepskyblue",
       "green", "orange", "red"))
   }





## -----------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------
## Successively higher GEV fits
## -----------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------



## asymp.fgev: first argument: how many maxima were used for the model;
##   second: which signal-to-noise ratio;
##   third argument: which sampling pattern / nb of observations;
##   fourth argument: which peak selection pattern in the periodogram.


asympt.fgev <- vector(5, mode = "list")
names(asympt.fgev) <- c("max.all","max.500","max.250","max.120","max.60")

asympt.fgev[[1]] <- primpeaks.fgev

for(jj in 2:5)
  asympt.fgev[[jj]] <- lapply(FirstStar.bootI, function(lst)
    lapply(lst,  function(mat)
      {
       res1.tmp <- vector(8, mode = "list")
       names(res1.tmp) <- c("syst.max","dep.max","fou.max","rand.max", 
          "wsyst.max","wdep.max","wfou.max","wrand.max") 
       for(ii in 1:8)
         {
          nMax <- if(jj < 5) 1000 else 992
          mat1 <- matrix(mat[ii,1:nMax], ncol = 2^(jj-1))
          vec <- apply(mat1, MAR = 1, FUN = max)
       	  res1.tmp[[ii]] <- try(gev(vec))
       	  }
       res1.tmp
       }))
       
length(asympt.fgev[[2]])
asympt.fgev[[5]][[6]][[7]][["wdep.max"]]
asympt.fgev[[5]][[6]][[7]][["wfou.max"]]
asympt.fgev[[1]][[6]][[7]][["wdep.max"]]
asympt.fgev[[1]][[6]][[7]][["wfou.max"]]
       
primpeaks.fgev[[6]]


## Extract the parameters, the standard errors, and the endpoint when it is 
## an upper endpoint:

asympt.gevpar <- lapply(asympt.fgev, function(lst3)
  lapply(lst3, function(lst)
     lapply(lst, function(lst2)
        sapply(lst2, function(lst1) 
           if(!inherits(lst1, "try-error")) lst1$par.ests else rep(NA, 3)))))

asympt.gevse <- lapply(asympt.fgev, function(lst3)
  lapply(lst3, function(lst)
     lapply(lst, function(lst2)
        sapply(lst2, function(lst1) 
          if(!inherits(lst1, "try-error")) lst1$par.ses else rep(NA, 3)))))

asympt.gevend <- lapply(asympt.gevpar, function(lst3)
  lapply(lst3, function(lst)
     sapply(lst, function(mat)
       apply(mat, MAR = 2, FUN = function(vec) 
          if(vec[1] >= 0 | is.na(vec[1])) NA else {vec[3] - vec[2]/vec[1]}))))

asympt.gevendse <- lapply(asympt.fgev, function(lst3)
  lapply(lst3, function(lst)
     sapply(lst, function(lst2)
        sapply(lst2, function(lst1) 
            if(!inherits(lst1, "try-error")) endse.fun(lst1$par.ests, lst1$varcov) else NA))))


## -----------------------------------------------------------------------------------
## Weighted and non-weighted fourier and dependent shapes, scales and
## locations against nObs:

#ii <- 6
color.tmp <- c("black","blue","deepskyblue","green","orange")
rng.tmp <- matrix(c(-0.6, 0.3,
                    0, 0.14,
                    0, 0.6), ncol = 2, byrow = TRUE)
partype.tmp <- c("Shape","Scale","Location")
#ss <- 3   ## frequency sampling type                   

cols.tmp <- 8:14
cols.tmp <- 1:7

#dev.set(2)
for(ii in 1:6)
{
 quartz(height = 7, width = 12)
 par(mfcol = c(3,4), mar = c(3,3,2,1), mgp = c(2,1,0))
# for(ss in c(3,7,2,6))   ## ss: index through frequency sampling type
 for(ss in c(1,5,4,8))   ## ss: index through frequency sampling type
  for(ll in 1:3)         ## ll: index through the gev parameters, xi = 1, sig = 2, mu = 3
   {
    par1.tmp <- sapply(asympt.gevpar, function(lst) 
     {res.tmp <- numeric(7)
      for(jj in 1:7)
#        res.tmp[jj] <- lst[[ii]][[jj+7]][ll,ss]  ## if nighttime
        res.tmp[jj] <- lst[[ii]][[jj]][ll,ss]  ## if uniform
      res.tmp
     })
    se1.tmp <- sapply(asympt.gevse, function(lst) 
     {res.tmp <- numeric(7)
      for(jj in 1:7)
#        res.tmp[jj] <- lst[[ii]][[jj+7]][ll,ss]  ## if nighttime
        res.tmp[jj] <- lst[[ii]][[jj]][ll,ss]  ## if uniform
      res.tmp
     })
    plot(log(nObs[cols.tmp]), par1.tmp[,1],
      ylim = rng.tmp[ll, ], main = paste(mag.names[ii], 
      ", ", nSamplings[ss], if(sum(cols.tmp == 1:7) == 7) 
      ", Uniform thinning" else ", Nighttime thinning",  sep = ""), 
      ylab = partype.tmp[ll], xlab = "Log(N)", pch = 21, col = 1, bg = 0, type = "n")
    for(kk in 1:5)
     {
      points(log(nObs[cols.tmp])+(kk-3)*0.04, par1.tmp[,kk], 
         pch = 16, col = color.tmp[kk])
      segments(log(nObs[cols.tmp])+(kk-3)*0.04, par1.tmp[,kk] + 2*se1.tmp[,kk],
         log(nObs[cols.tmp])+(kk-3)*0.04, par1.tmp[,kk] - 2*se1.tmp[,kk],
         col = color.tmp[kk], lty = 1)
     }
   }
}

##  - Weighted estimates are more variable than
##    the corresponding non-weighted estimates for all numbers of maxima used.
##  - weighted shape estimates are more negative than non-weighted ones.
##  Behaviour as a function of used number of maxima:
##  Shape:
##    estimates are stable using 1000, 500, 250 maxima (and often for 125 maxima,
##    too), for all types of frequency sampling and for weighted and non-weighted.
##    It is different when using 125 or 62 maxima, deviates into both positive and 
##    negative direction, but has a large uncertainty, the CIs usually contain all
##    the other estimates.  
##  Scale:
##    The same overall pattern, decreasing scale with increasing nObs can be
##    observed. The values with dependent sampling decrease if fewer maxima were
##    used, but stay the same with fourier sampling.
##  Location:
##    The same overall pattern, decreasing location with increasing nObs can be
##    observed, but the values increase if fewer maxima are used for the gev fit.
##    This is both for weighted and un-weighted, both for dependent or fourier sampling. 
##  Nighttime and uniform thinning behaves the same way.



## -----------------------------------------------------------------------------------
## Huge quantile plots to compare the weighted and non-weighted estimates:

ss <- 6
ii <- 1
for(ii in 1:6)
 {
  quartz(height = 7.5, width = 12.5)
  par(mfcol = c(3,5),  mar = c(3,3,2,1), mgp = c(2,1,0))
  for(jj in 1:14)
   {
   	axl.tmp <- range(unlist(FirstStar.bootI[[ii]][[jj]][ss,]))
    plot(sort(FirstStar.bootI[[1]][[1]][8,]),
        sort(FirstStar.bootI[[1]][[1]][8,]), ylim = axl.tmp,
        xlim = axl.tmp, main = paste(names(primpeaks.gevpar)[ii], ", ",
    	names(primpeaks.gevpar[[ii]])[jj], ", ", nSamplings[ss], sep = ""),
    	xlab = "Theoretical quantiles", ylab = "Sample", type = "n")
	for(kk in 1:5)
   		{
#   	devAskNewPage(ask = TRUE)
         nMax <- if(kk < 5) 1000 else 992
         mat1 <- matrix(FirstStar.bootI[[ii]][[jj]][ss,1:nMax], ncol = 2^(kk-1))
         vec <- apply(mat1, MAR = 1, FUN = max)
   		 abline(c(0,1), col = "black")
    	 points(qgev(ppoints(length(vec)), xi = asympt.gevpar[[kk]][[ii]][[jj]][1,ss],
         	sigma = asympt.gevpar[[kk]][[ii]][[jj]][2,ss], 
         	mu = asympt.gevpar[[kk]][[ii]][[jj]][3,ss]),
         	sort(vec), pch = 16, cex = 1.3-0.15*kk,
         	col = color.tmp[kk])
 	 	}
	legend(axl.tmp[1], axl.tmp[2], legend = c("nMax = 1000", "nMax = 500", 
	    "nMax = 250", "nMax = 125", "nMax = 62"), col = color.tmp,
    	pt.cex = 1.3-0.15*(1:5), pch = 16, bty = "n")
   }
 }


## -----------------------------------------------------------------------------------
## High quantiles and their stability across nMax

asympt.quan <- lapply(asympt.gevpar, function(lst)
 lapply(lst, function(lst1)
  lapply(lst1, function(mat)
   apply(mat, MAR = 2, function(parvec) 
    c(q.95 = qgev(0.95, xi = parvec[1], mu = parvec[3], sigma = parvec[2]),
      q.99 = qgev(0.99, xi = parvec[1], mu = parvec[3], sigma = parvec[2]),
      q.999 = qgev(0.999, xi = parvec[1], mu = parvec[3], sigma = parvec[2]),
      q.9995 = qgev(0.9995, xi = parvec[1], mu = parvec[3], sigma = parvec[2]),
      q.9999 = qgev(0.9999, xi = parvec[1], mu = parvec[3], sigma = parvec[2])))))



quartz(height = 7.5, width = 8)
colnnn <- "sim50"
mag.tmp <- "mag.g6"
magerr.tmp <- "mag.error.g6"  
#smp.tmp <- 2
#mycol <- "orange"
par(mfrow = c(4,2), mar = c(3,3,2,1))
for(smp.tmp in 1:4)
   {
    plot(FirstStar.pgram[[mag.tmp]]$freq, 
       FirstStar.pgram[[mag.tmp]][,colnnn], type = "n",
       xlab = "Frequency", ylab = "Chi-squared", 
       main = paste(colnnn, mag.tmp, nSamplings[smp.tmp]), 
       ylim = c(0,max(max(FirstStar.bootI[[mag.tmp]][[colnnn]][nSamplings[smp.tmp],], 
       na.rm = TRUE),# max(FirstStar.pgram[[mag.tmp]][,colnnn], na.rm = TRUE),
       unlist(asympt.quan[[3]][[mag.tmp]][[colnnn]]))),
#  ylim = c(0,0.2),
       xlim = c(0,100))
    abline(v = F, col = "violetred", lty = 2)
    points(FirstStar.pgram[[mag.tmp]]$freq, FirstStar.pgram[[mag.tmp]][,colnnn],
       type = "h", col = "grey45")
    abline(h = sapply(asympt.quan, function(lst)
       lst[[mag.tmp]][[colnnn]][2,smp.tmp]), col = c("black","blue", "deepskyblue",
       "green", "orange"))
    plot(FirstStar.wpgram[[mag.tmp]]$freq, 
       FirstStar.wpgram[[mag.tmp]][,colnnn], type = "n",
       xlab = "Frequency", ylab = "Chi-squared", 
       main = paste(colnnn, mag.tmp, nSamplings[smp.tmp+4]), 
       ylim = c(0,max(max(FirstStar.bootI[[mag.tmp]][[colnnn]][nSamplings[smp.tmp+4],],
       na.rm = TRUE),# max(FirstStar.wpgram[[mag.tmp]][,colnnn], na.rm = TRUE),
       unlist(asympt.quan[[3]][[mag.tmp]][[colnnn]]))),
#  ylim = c(0,0.2),
       xlim = c(0,100))
    abline(v = F, col = "violetred", lty = 2)
    points(FirstStar.wpgram[[mag.tmp]]$freq, FirstStar.wpgram[[mag.tmp]][,colnnn],
       type = "h", col = "grey45")
    abline(h = sapply(asympt.quan, function(lst)
       lst[[mag.tmp]][[colnnn]][2,smp.tmp+4]), col = c("black","blue", "deepskyblue",
       "green", "orange"))
   }


## One more check that looks necessary (unfortunately, it will be able 
## to kill all that I have done):
##  - repeat this bootstrap for pure noise (a mag.g7....)
##  - in some 200 case, calculate the whole 40001 spectrum value, 
##    and see how many times the peak exceeds the estimated 1% quantile!
200*40001
## that's 8 million... ok, ask Hubert.

## Another necessary check: check the stability of these GEVs based on all
## 1000 repetitions against the block size = number of test frequencies



## -----------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------
## Stability against block size:
## at how many frequencies to compute the periodogram?
## -----------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------
 


mag.names2 <- c("mag.g3","mag.g5","mag.g6")
mag.error.names2 <- c("mag.error.g3","mag.error.g5","mag.error.g6")

nTestFreq2 <- matrix(c(500,400,300,200,100,50,125,100,75,50,25,13), ncol = 2)  
nObs2 <- c(100,25)
nNames2 <- c("nights100", "nights25")

FirstStar.bootII <- vector(3, mode = "list")
names(FirstStar.bootII) <- mag.names2

for(jj in 1:3)
  {
   bootI.tmp <- vector(2, mode = "list")
   names(bootI.tmp) <- nNames2
   for(ii in 1:2)
 	{
  	 print(system.time(
   	 {
      nnn <- nObs2[ii]
  	  colnnn <- nNames2[ii]
      dfr <- FirstStar[FirstStar[,colnnn] == 1,]
 	  ts <- demean.fun(dfr$time)

      bootIII.tmp <- replicate(1000,
		   {
  	  	    CentralInd <- sample(2:40000, size = 1)
  		   		 
  		    ## Select Fourier frequencies:
  			FreqGridFou <- (1:2499)/25
 	     	FreqGridFou <- sample(FreqGridFou, size = nTestFreq2[1,ii])
   		 
  		 	## Select peak-width intervals around the Fourier frequencies:
  		 	FreqGridDep <- sapply(FreqGridFou, function(ff)
 	        	ff + (-8:8)*dfreq, simplify = FALSE)
 	     	FreqGridDep <- FreqGridDep[sapply(FreqGridDep, function(vec)
 	        	!any(vec < 0.04 | vec >= 100))]
#            FreqGridDep <- unlist(FreqGridDep)
              		  	
            ## The bootstrap:  		 
      	    ind <- sample(1:nrow(dfr), replace = TRUE)
      	    ps <- demean.fun(dfr[ind, mag.names2[jj]])
      	    ## Unweighted:
            ws <- rep(1/nnn, nnn)
            frspDep.tmp <- lapply(FreqGridDep, function(lst) 
                sapply(lst, MyFMLS3.fun, t = ts, h = ps, Weights = ws))
      	    frspFou.tmp <- sapply(FreqGridFou, MyFMLS3.fun, t = ts, h = ps,  
         		Weights = ws)
            ## Weighted:
            ws <- (1/dfr[ind, mag.error.names2[jj]]^2)/sum(1/dfr[ind, mag.error.names2[jj]]^2)
            wfrspDep.tmp <- lapply(FreqGridDep, function(lst) 
                sapply(lst, MyFMLS3.fun, t = ts, h = ps, Weights = ws))
      	    wfrspFou.tmp <- sapply(FreqGridFou, MyFMLS3.fun, t = ts, h = ps,  
         		Weights = ws)

            bootII.tmp <- matrix(ncol = 4, nrow = 6)
            dimnames(bootII.tmp) <- list(nTestFreq2[,ii], 
                c("dep.max","fou.max","wdep.max","wfou.max"))
            for(kk in 1:6)
              {
               frspFou2.tmp <- frspFou.tmp[1:nTestFreq2[kk,ii]] 
               frspDep2.tmp <- unlist(frspDep.tmp[1:min(nTestFreq2[kk,ii],
                 length(frspDep.tmp))])
               wfrspFou2.tmp <- wfrspFou.tmp[1:nTestFreq2[kk,ii]] 
               wfrspDep2.tmp <- unlist(wfrspDep.tmp[1:min(nTestFreq2[kk,ii],
                 length(wfrspDep.tmp))]) 
		       bootII.tmp[as.character(nTestFreq2[kk,ii]), ] <- c(
	             max(frspDep2.tmp[!is.na(frspDep2.tmp) & abs(frspDep2.tmp) != Inf]), 
                 max(frspFou2.tmp[!is.na(frspFou2.tmp) & abs(frspFou2.tmp) != Inf]),
	             max(wfrspDep2.tmp[!is.na(wfrspDep2.tmp) & abs(wfrspDep2.tmp) != Inf]), 
                 max(wfrspFou2.tmp[!is.na(wfrspFou2.tmp) & abs(wfrspFou2.tmp) != Inf])
  		         )
  		       }
  		   bootII.tmp
		   }, simplify = FALSE) 
 		 }))
 	   prb <- vector(6, mode = "list")
       names(prb) <- nTestFreq2[,ii]
 	   for(kk in 1:6)
 	       prb[[as.character(nTestFreq2[kk,ii])]] <- sapply(bootIII.tmp, function(mat) mat[kk, ])
 	   bootI.tmp[[nNames2[ii]]] <- prb
	  }
	FirstStar.bootII[[mag.names2[jj]]] <- bootI.tmp
	cat(mag.names2[jj], " ready \n")
    save.image(WORKSPACE.NAME)
 	} 
save.image(WORKSPACE.NAME)

FirstStar.bootIV <- vector(3, mode = "list")
names(FirstStar.bootIV) <- mag.names2
nTestFreq3 <- c(200, 300, 400, 500, 800)

for(jj in 1:3)
  {
   ii <- 2
 	{
  	 print(system.time(
   	 {
      nnn <- nObs2[ii]
  	  colnnn <- nNames2[ii]
      dfr <- FirstStar[FirstStar[,colnnn] == 1,]
 	  ts <- demean.fun(dfr$time)

      bootIII.tmp <- replicate(1000,
		   {
  	  	    CentralInd <- sample(2:40000, size = 1)
  		   		 
  		    ## Select Fourier frequencies:
  			FreqFou.ind <- 16*(1:2499) + 1
  			FreqFou.rnd <- sample(FreqFou.ind, size = nTestFreq3[kk])
  			FreqGridFou <- FreqGrid1[FreqFou.rnd]
   		 
  		 	## Select peak-width intervals around the Fourier frequencies:
  		 	FreqDep.ind <- sapply(FreqFou.rnd, function(ff)
 	        	ff + (-8:8), simplify = FALSE)
 	     	FreqDep.ind2 <- FreqDep.ind[sapply(FreqDep.ind, function(vec)
 	        	!any(vec < 2 | vec > 40000))]
 	     	FreqGridDep <- lapply(FreqDep.ind2, function(vec)
 	     	     FreqGrid1[vec])
#            FreqGridDep <- sort(unlist(FreqGridDep))
              		  	
            ## The bootstrap:  		 
      	    ind <- sample(1:nrow(dfr), replace = TRUE)
      	    ps <- demean.fun(dfr[ind, mag.names2[jj]])
      	    ## Unweighted:
            ws <- rep(1/nnn, nnn)
            frspDep.tmp <- lapply(FreqGridDep, function(lst) 
                sapply(lst, MyFMLS3.fun, t = ts, h = ps, Weights = ws))
      	    frspFou.tmp <- sapply(FreqGridFou, MyFMLS3.fun, t = ts, h = ps,  
         		Weights = ws)
            ## Weighted:
            ws <- (1/dfr[ind, mag.error.names2[jj]]^2)/sum(1/dfr[ind, mag.error.names2[jj]]^2)
            wfrspDep.tmp <- lapply(FreqGridDep, function(lst) 
                sapply(lst, MyFMLS3.fun, t = ts, h = ps, Weights = ws))
      	    wfrspFou.tmp <- sapply(FreqGridFou, MyFMLS3.fun, t = ts, h = ps,  
         		Weights = ws)

            bootII.tmp <- matrix(ncol = 4, nrow = length(nTestFreq3))
            dimnames(bootII.tmp) <- list(nTestFreq3, 
                c("dep.max","fou.max","wdep.max","wfou.max"))
            for(kk in 1:length(nTestFreq3))
              {
               frspFou2.tmp <- frspFou.tmp[1:nTestFreq3[kk]] 
               frspDep2.tmp <- unlist(frspDep.tmp[1:min(nTestFreq3[kk],
                 length(frspDep.tmp))])
               wfrspFou2.tmp <- wfrspFou.tmp[1:nTestFreq3[kk]] 
               wfrspDep2.tmp <- unlist(wfrspDep.tmp[1:min(nTestFreq3[kk],
                 length(wfrspDep.tmp))]) 
		       bootII.tmp[as.character(nTestFreq3[kk]), ] <- c(
	             max(frspDep2.tmp[!is.na(frspDep2.tmp) & abs(frspDep2.tmp) != Inf]), 
                 max(frspFou2.tmp[!is.na(frspFou2.tmp) & abs(frspFou2.tmp) != Inf]),
	             max(wfrspDep2.tmp[!is.na(wfrspDep2.tmp) & abs(wfrspDep2.tmp) != Inf]), 
                 max(wfrspFou2.tmp[!is.na(wfrspFou2.tmp) & abs(wfrspFou2.tmp) != Inf])
  		         )
  		       }
  		   bootII.tmp
		   }, simplify = FALSE) 
 		 }))
 	   prb <- vector(6, mode = "list")
       names(prb) <- nTestFreq3
 	   for(kk in 1:length(nTestFreq3))
 	       prb[[as.character(nTestFreq3[kk])]] <- sapply(bootIII.tmp, function(mat) mat[kk, ])
 	   bootI.tmp[[nNames2[ii]]] <- prb
	  }
	FirstStar.bootIV[[mag.names2[jj]]] <- bootI.tmp
	cat(mag.names2[jj], " ready \n")
    save.image(WORKSPACE.NAME)
 	} 
save.image(WORKSPACE.NAME)


names(FirstStar.bootIV)
names(FirstStar.bootIV[[1]])
names(FirstStar.bootIV[[2]][[1]])

FirstStar.bootIV[[1]][[2]] <- FirstStar.bootIV[[1]][[2]][c("200","300","400","500","800")]
FirstStar.bootIV[[2]][[2]] <- FirstStar.bootIV[[2]][[2]][c("200","300","400","500","800")]
FirstStar.bootIV[[3]][[2]] <- FirstStar.bootIV[[3]][[2]][c("200","300","400","500","800")]

FirstStar.bootII.backup <- FirstStar.bootII

names(FirstStar.bootII)
names(FirstStar.bootII[[1]])
names(FirstStar.bootII[[2]][["nights100"]])

FirstStar.bootII[[1]][["nights25"]] <- c(FirstStar.bootII[[1]][["nights25"]], FirstStar.bootIV[[1]][["nights25"]])
FirstStar.bootII[[2]][["nights25"]] <- c(FirstStar.bootII[[2]][["nights25"]], FirstStar.bootIV[[2]][["nights25"]])
FirstStar.bootII[[3]][["nights25"]] <- c(FirstStar.bootII[[3]][["nights25"]], FirstStar.bootIV[[3]][["nights25"]])


growblock.fgev <- lapply(FirstStar.bootII, function(lst)
  lapply(lst,  function(lst1)
    lapply(lst1, function(mat)  
      apply(mat, MAR = 1, FUN = function(vec)
       	try(gev(vec))))))
       
#length(asympt.fgev[[2]])
#asympt.fgev[[5]][[6]][[7]][["wdep.max"]]
#asympt.fgev[[5]][[6]][[7]][["wfou.max"]]
#asympt.fgev[[1]][[6]][[7]][["wdep.max"]]
#asympt.fgev[[1]][[6]][[7]][["wfou.max"]]
       

## Extract the parameters, the standard errors, and the endpoint when it is 
## an upper endpoint:

growblock.gevpar <- lapply(growblock.fgev, function(lst3)
  lapply(lst3, function(lst)
     lapply(lst, function(lst2)
        sapply(lst2, function(lst1) 
           if(!inherits(lst1, "try-error")) lst1$par.ests else rep(NA, 3)))))

growblock.gevse <- lapply(growblock.fgev, function(lst3)
  lapply(lst3, function(lst)
     lapply(lst, function(lst2)
        sapply(lst2, function(lst1) 
          if(!inherits(lst1, "try-error")) lst1$par.ses else rep(NA, 3)))))

growblock.gevend <- lapply(growblock.gevpar, function(lst3)
  lapply(lst3, function(lst)
     sapply(lst, function(mat)
       apply(mat, MAR = 2, FUN = function(vec) 
          if(vec[1] >= 0 | is.na(vec[1])) NA else {vec[3] - vec[2]/vec[1]}))))

growblock.gevendse <- lapply(growblock.fgev, function(lst3)
  lapply(lst3, function(lst)
     sapply(lst, function(lst2)
        sapply(lst2, function(lst1) 
            if(!inherits(lst1, "try-error")) endse.fun(lst1$par.ests, lst1$varcov) else NA))))

growblock.quan <- lapply(growblock.gevpar, function(lst)
 lapply(lst, function(lst1)
  lapply(lst1, function(mat)
   apply(mat, MAR = 2, function(parvec) 
    c(q.95 = qgev(0.95, xi = parvec[1], mu = parvec[3], sigma = parvec[2]),
      q.99 = qgev(0.99, xi = parvec[1], mu = parvec[3], sigma = parvec[2]),
      q.995 = qgev(0.995, xi = parvec[1], mu = parvec[3], sigma = parvec[2]),
      q.999 = qgev(0.999, xi = parvec[1], mu = parvec[3], sigma = parvec[2]),
      q.9999 = qgev(0.9999, xi = parvec[1], mu = parvec[3], sigma = parvec[2]))))))

rlevel.gevnoplot <- function (out, k.blocks = 20, add = FALSE, ...) 
{
    par.ests <- out$par.ests
    mu <- par.ests["mu"]
    sigma <- par.ests["sigma"]
    if (!("xi" %in% names(par.ests))) 
        stop("Use this function after a GEV rather than a Gumbel fit")
    else xi <- par.ests["xi"]
    pp <- 1/k.blocks
    v <- qgev((1 - pp), xi, mu, sigma)
    as.numeric(v)
}
 
rlevel.mine <- function(out, k.blocks = 20, add = FALSE, ...) 
   {
   	r1.tmp <- try(rlevel.gevnoplot(out = out, k.blocks = k.blocks, add = add, ...))
   	r2.tmp <- if(inherits(r1.tmp, "try-error")) rep(NA, 3) else r1.tmp
   	r2.tmp
   }
rlevel.mine1 <- function(out, k.blocks = 20, add = FALSE, ...) 
   {
   	r1.tmp <- try(rlevel.gev(out = out, k.blocks = k.blocks, add = add, ...))
   	r2.tmp <- if(inherits(r1.tmp, "try-error")) rep(NA, 3) else r1.tmp
   	r2.tmp
   }

gev.mine <- function(data, block = NA, init.par = "built-in", ...) 
{
    n.all <- NA
    if (!is.na(block)) {
        n.all <- length(data)
        if (is.character(block)) {
            times <- as.POSIXlt(attributes(data)$times)
            if (block %in% c("semester", "quarter")) {
                sem <- quart <- times$mon
                sem[sem %in% 0:5] <- quart[quart %in% 0:2] <- 0
                sem[sem %in% 6:11] <- quart[quart %in% 3:5] <- 1
                quart[quart %in% 6:8] <- 2
                quart[quart %in% 9:11] <- 3
            }
            grouping <- switch(block, semester = paste(times$year, 
                sem), quarter = paste(times$year, quart), month = paste(times$year, 
                times$mon), year = times$year, stop("unknown time period"))
            data <- tapply(data, grouping, max)
        }
        else {
            data <- as.numeric(data)
            nblocks <- (length(data)%/%block) + 1
            grouping <- rep(1:nblocks, rep(block, nblocks))[1:length(data)]
            data <- tapply(data, grouping, max)
        }
    }
    data <- as.numeric(data)
    n <- length(data)
    sigma0 <- if(init.par[1] == "built-in") sqrt(6 * var(data))/pi else init.par[2]
    mu0 <- if(init.par[1] == "built-in") mean(data) - 0.57722 * sigma0 else init.par[3]
    xi0 <- if(init.par[1] == "built-in") 0.1 else init.par[1]
    theta <- c(xi0, sigma0, mu0)
    negloglik <- function(theta, tmp) {
        y <- 1 + (theta[1] * (tmp - theta[3]))/theta[2]
        if ((theta[2] < 0) || (min(y) < 0)) 
            out <- 1e+06
        else {
            term1 <- length(tmp) * logb(theta[2])
            term2 <- sum((1 + 1/theta[1]) * logb(y))
            term3 <- sum(y^(-1/theta[1]))
            out <- term1 + term2 + term3
        }
        out
    }
    fit <- optim(theta, negloglik, hessian = TRUE, ..., tmp = data)
    if (fit$convergence) 
        warning("optimization may not have succeeded")
    par.ests <- fit$par
    varcov <- solve(fit$hessian)
    par.ses <- sqrt(diag(varcov))
    out <- list(n.all = n.all, n = n, data = data, block = block, 
        par.ests = par.ests, par.ses = par.ses, varcov = varcov, 
        converged = fit$convergence, nllh.final = fit$value)
    names(out$par.ests) <- c("xi", "sigma", "mu")
    names(out$par.ses) <- c("xi", "sigma", "mu")
    class(out) <- "gev"
    out
}

#blocklen.tmp <- 500 * 17
#gevfit.tmp <- growblock.fgev[[3]][[1]][["500"]][["dep.max"]]
#rlevel.mine(gevfit.tmp, k.blocks = 100*40000/blocklen.tmp)


## -----------------------------------------------------------------------------------
## CIs by bootstrap for R = 1000
## -----------------------------------------------------------------------------------

  
growblockboot.fgev <- lapply(FirstStar.bootII, function(lst)
  lapply(lst,  function(lst1)
    lapply(lst1, function(mat)  
      lapply(data.frame(t(mat)), FUN = function(vec)
        replicate(1000, { 
          vec1 <- sample(vec, replace = TRUE)
       	  try(gev(vec1))
       	  }
       	  )))))
growblockboot.gevpar <- lapply(growblockboot.fgev, function(lst)
  lapply(lst,  function(lst1)
    lapply(lst1, function(mat)  
      lapply(data.frame(t(mat)), FUN = function(vec)
        replicate(1000, { 
          vec1 <- sample(vec, replace = TRUE)
       	  res.tmp <- try(gev(vec1))
          if(!inherits(res.tmp, "try-error")) res.tmp$par.ests else rep(NA, 3)
       	  }
       	  )))))

names(FirstStar.bootII[[ii]][[jj]][[kk]])     	  
growblockboot.gevpar <- vector(length(FirstStar.bootII), mode = "list")
names(growblockboot.gevpar) <- names(FirstStar.bootII)
growblockfullboot.quan <- vector(length(FirstStar.bootII), mode = "list")
names(growblockfullboot.quan) <- names(FirstStar.bootII)
system.time(
for(ii in 1:length(FirstStar.bootII))
   {
    v1.tmp <- vector(length(FirstStar.bootII[[ii]]), mode = "list")
    names(v1.tmp) <- names(FirstStar.bootII[[ii]])
    w1.tmp <- vector(length(FirstStar.bootII[[ii]]), mode = "list")
    names(w1.tmp) <- names(FirstStar.bootII[[ii]])
    for(jj in 1:length(FirstStar.bootII[[ii]]))
       {
        v2.tmp <- vector(length(FirstStar.bootII[[ii]][[jj]]), mode = "list")
        names(v2.tmp) <- names(FirstStar.bootII[[ii]][[jj]])
        w2.tmp <- vector(length(FirstStar.bootII[[ii]][[jj]]), mode = "list")
        names(w2.tmp) <- names(FirstStar.bootII[[ii]][[jj]])
        nameset <- which(is.element(names(FirstStar.bootII[[ii]][[jj]]), c("100","500")))
#        for(kk in nameset)
        for(kk in 1:length(FirstStar.bootII[[ii]][[jj]]))
           {
            v3.tmp <- vector(nrow(FirstStar.bootII[[ii]][[jj]][[kk]]), mode = "list")
            names(v3.tmp) <- dimnames(FirstStar.bootII[[ii]][[jj]][[kk]])[[1]]
            w3.tmp <- vector(nrow(FirstStar.bootII[[ii]][[jj]][[kk]]), mode = "list")
            names(w3.tmp) <- dimnames(FirstStar.bootII[[ii]][[jj]][[kk]])[[1]]
            nn.tmp <- names(FirstStar.bootII[[ii]][[jj]])[kk]
     	    blocklen.tmp <- as.numeric(nn.tmp)
            # for dependent maxima:
            res.tmp <- replicate(1000, {
            	   vec1 <- sample(FirstStar.bootII[[ii]][[jj]][[kk]]["dep.max",], replace = TRUE);
            	   try(gev.mine(vec1, init.par = growblock.gevpar[[ii]][[jj]][[kk]][,"dep.max"]))
                   }, simplify = FALSE)
            v3.tmp[["dep.max"]] <- sapply(res.tmp, FUN = function(lst)
                   if(!inherits(lst, "try-error")) lst$par.ests else rep(NA, 3))
            w3.tmp[["dep.max"]] <- sapply(res.tmp, FUN = function(lst)
                   if(inherits(lst, "try-error")) rep(NA, 3) else 
                   c(qboot.95 = rlevel.mine(lst, k.blocks = 20*40000/(17*blocklen.tmp)),
                   qboot.99 = rlevel.mine(lst, k.blocks = 100*40000/(17*blocklen.tmp)),
                   qboot.995 = rlevel.mine(lst, k.blocks = 200*40000/(17*blocklen.tmp))))
            # for weighted dependent maxima:
            res.tmp <- replicate(1000, {
            	   vec1 <- sample(FirstStar.bootII[[ii]][[jj]][[kk]]["wdep.max",], replace = TRUE);       	               try(gev.mine(vec1, init.par = growblock.gevpar[[ii]][[jj]][[kk]][,"wdep.max"]))
                   }, simplify = FALSE)
            v3.tmp[["wdep.max"]] <- sapply(res.tmp, FUN = function(lst)
                   if(!inherits(lst, "try-error")) lst$par.ests else rep(NA, 3))
            w3.tmp[["wdep.max"]] <- sapply(res.tmp, FUN = function(lst)
                   if(inherits(lst, "try-error")) rep(NA, 3) else 
                   c(qboot.95 = rlevel.mine(lst, k.blocks = 20*40000/(17*blocklen.tmp)),
                   qboot.99 = rlevel.mine(lst, k.blocks = 100*40000/(17*blocklen.tmp)),
                   qboot.995 = rlevel.mine(lst, k.blocks = 200*40000/(17*blocklen.tmp))))
            # for fourier maxima:
            ei <- 1/16
            res.tmp <- replicate(1000, { 
            	   vec1 <- sample(FirstStar.bootII[[ii]][[jj]][[kk]]["fou.max",], replace = TRUE);       	               try(gev.mine(vec1, init.par = growblock.gevpar[[ii]][[jj]][[kk]][,"fou.max"]))
                   }, simplify = FALSE)
            v3.tmp[["fou.max"]] <- sapply(res.tmp, FUN = function(lst)
                   if(!inherits(lst, "try-error")) lst$par.ests else rep(NA, 3))
            w3.tmp[["fou.max"]] <- sapply(res.tmp, FUN = function(lst)
                   if(inherits(lst, "try-error")) rep(NA, 3) else 
                   c(qboot.95 = rlevel.mine(lst, k.blocks = 20*ei*40000/blocklen.tmp),
                   qboot.99 = rlevel.mine(lst, k.blocks = 100*ei*40000/blocklen.tmp),
                   qboot.995 = rlevel.mine(lst, k.blocks = 200*ei*40000/blocklen.tmp)))
            # for weighted fourier maxima:
            res.tmp <- replicate(1000, { 
            	   vec1 <- sample(FirstStar.bootII[[ii]][[jj]][[kk]]["wfou.max",], replace = TRUE);      	               try(gev.mine(vec1, init.par = growblock.gevpar[[ii]][[jj]][[kk]][,"wfou.max"]))
                   }, simplify = FALSE)
            v3.tmp[["wfou.max"]] <- sapply(res.tmp, FUN = function(lst)
                   if(!inherits(lst, "try-error")) lst$par.ests else rep(NA, 3))
            w3.tmp[["wfou.max"]] <- sapply(res.tmp, FUN = function(lst)
                   if(inherits(lst, "try-error")) rep(NA, 3) else 
                   c(qboot.95 = rlevel.mine(lst, k.blocks = 20*ei*40000/blocklen.tmp),
                   qboot.99 = rlevel.mine(lst, k.blocks = 100*ei*40000/blocklen.tmp),
                   qboot.995 = rlevel.mine(lst, k.blocks = 200*ei*40000/blocklen.tmp)))
            v2.tmp[[kk]] <- v3.tmp
            w2.tmp[[kk]] <- w3.tmp
#            save.image("extrFAP7FirstStar-tmp.RData")
            cat("ii=",ii, "jj=", jj, "kk=",kk, "ready \n")
           }
        v1.tmp[[jj]] <- v2.tmp
        w1.tmp[[jj]] <- w2.tmp
        cat("ii=",ii, "jj=", jj, "ready \n")
       }
    growblockboot.gevpar[[ii]] <- v1.tmp
    growblockfullboot.quan[[ii]] <- w1.tmp
    cat("ii=",ii, "ready \n")
   }
)

   
names(growblockfullboot.quan[["mag.g6"]][["nights25"]])       
names(growblockfullboot.quan[["mag.g6"]][["nights25"]])       

names(growblockfullboot.quan[["mag.g6"]][["nights25"]][["dep.max"]])       
growblockfullboot.quan[["mag.g6"]][["nights100"]][["100"]][["wfou.max"]]

growblockboot.quan <- lapply(growblockboot.gevpar, function(lst)  ## appl to magnitudes
 lapply(lst, function(lst1)   ## appl to ts length
  lapply(lst1, function(lst2)   ## appl to ts length
   lapply(lst2, function(mat)   ## appl to psearch method
    apply(mat, MAR = 2, function(parvec) 
     c(q.95 = qgev(0.95, xi = parvec[1], mu = parvec[3], sigma = parvec[2]),
       q.99 = qgev(0.99, xi = parvec[1], mu = parvec[3], sigma = parvec[2]),
       q.995 = qgev(0.995, xi = parvec[1], mu = parvec[3], sigma = parvec[2]),
       q.999 = qgev(0.999, xi = parvec[1], mu = parvec[3], sigma = parvec[2]),
       q.9999 = qgev(0.9999, xi = parvec[1], mu = parvec[3], sigma = parvec[2])))))))
growblockboot.quan[["mag.g6"]][["nights25"]][["100"]][["dep.max"]]



## ---------------------------------------------------------------------------------------------------
## High quantiles in the reduced periodogram and their stability with the number of test frequencies
## ---------------------------------------------------------------------------------------------------


quartz(height = 7.5, width = 4)
colnnn <- "nights25"
mag.tmp <- "mag.g6"
#magerr.tmp <- "mag.error.g5"  
#smp.tmp <- 2
#mycol <- "orange"
par(mfrow = c(4,1), mar = c(3,3,2,1))
for(smp.tmp in c(2,3))
   {
   	# unweighted:
    plot(FirstStar.pgram[[mag.tmp]]$freq, 
       FirstStar.pgram[[mag.tmp]][,colnnn], type = "n",
       xlab = "Frequency", ylab = "Chi-squared", 
       main = paste(colnnn, mag.tmp, nSamplings[smp.tmp]), 
#       ylim = c(0,max(max(FirstStar.bootII[[mag.tmp]][[colnnn]][[]][nSamplings[smp.tmp],], 
#       na.rm = TRUE), max(FirstStar.pgram[[mag.tmp]][,colnnn], na.rm = TRUE),
#       unlist(growblock.quan[[mag.tmp]][[colnnn]]))),
       ylim = c(0,1),
       xlim = c(0,100))
    abline(v = F, col = "red", lty = 2)
    points(FirstStar.pgram[[mag.tmp]]$freq, FirstStar.pgram[[mag.tmp]][,colnnn],
       type = "h", col = "grey45")
    abline(h = sapply(growblock.quan[[mag.tmp]][[colnnn]][c("50","100","200","300","400","500")], function(mat)
       mat[2,nSamplings[smp.tmp]]), col = c("black","blue", "deepskyblue",
       "green", "orange","red"))

    ## weighted: 
    plot(FirstStar.wpgram[[mag.tmp]]$freq, 
       FirstStar.wpgram[[mag.tmp]][,colnnn], type = "n",
       xlab = "Frequency", ylab = "Chi-squared", 
       main = paste(colnnn, mag.tmp, nSamplings[smp.tmp+4]), 
#       ylim = c(0,max(max(FirstStar.bootII[[mag.tmp]][[colnnn]][[1]][nSamplings[smp.tmp+4],], 
#       na.rm = TRUE), max(FirstStar.wpgram[[mag.tmp]][,colnnn], na.rm = TRUE),
#       unlist(growblock.quan[[mag.tmp]][[colnnn]]))),
       ylim = c(0,1),
       xlim = c(0,100))
    abline(v = F, col = "red", lty = 2)
    points(FirstStar.wpgram[[mag.tmp]]$freq, FirstStar.wpgram[[mag.tmp]][,colnnn],
       type = "h", col = "grey45")
    abline(h = sapply(growblock.quan[[mag.tmp]][[colnnn]][c("50","100","200","300","400","500")], function(mat)
       mat[2,nSamplings[smp.tmp+4]]), col = c("black","blue", "deepskyblue",
       "green", "orange","red"))
   }
## Looks good: the three or four uppermost 0.99 quantiles are relatively stable.
## I should put some CI on them, I hope they would strongly overlap.

rlevel.gev(growblock.fgev[["mag.g6"]][["nights25"]][["500"]][["dep.max"]], k.blocks = 100)

## The CI does not always work.


## ----------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------
## Compute high quantiles for the complete spectrum
## ----------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------

##  - for the dependent sampling, it is just the obtained 
##    F(q)^(number of blocks), that is, 
##    F(q)^(nb of all freqs in the full periodogram / nb of checked freqs for the maxima)
##  - for the simple random sampling, it is 
##    F(q)^(nb of all freqs in the full periodogram / nb of checked freqs for the maxima) 

growblockfull.quan2 <- lapply(growblock.fgev, function(lst)
# lst <- growblock.fgev[[1]] 
  lapply(lst, function(lst1)
    {
#     lst1 <- lst[[2]]
     res.dep <- matrix(ncol = length(lst1), nrow = 18)
     dimnames(res.dep) <- list(c("qlb.95","q.95","qub.95","qlbboot.95","qboot.95","qubboot.95",
          "qlb.99","q.99","qub.99","qlbboot.99","qboot.99","qubboot.99",
          "qlb.995","q.995","qub.995","qlbboot.995","qboot.995","qubboot.995"), names(lst1))
     res.wdep <- matrix(ncol = length(lst1), nrow = 18)
     dimnames(res.wdep) <- list(c("qlb.95","q.95","qub.95","qlbboot.95","qboot.95","qubboot.95",
          "qlb.99","q.99","qub.99","qlbboot.99","qboot.99","qubboot.99",
          "qlb.995","q.995","qub.995","qlbboot.995","qboot.995","qubboot.995"), names(lst1))
     res.fou <- matrix(ncol = length(lst1), nrow = 18)
     dimnames(res.fou) <- list(c("qlb.95","q.95","qub.95","qlbboot.95","qboot.95","qubboot.95",
          "qlb.99","q.99","qub.99","qlbboot.99","qboot.99","qubboot.99",
          "qlb.995","q.995","qub.995","qlbboot.995","qboot.995","qubboot.995"), names(lst1))
     res.wfou <- matrix(ncol = length(lst1), nrow = 18)
     dimnames(res.wfou) <- list(c("qlb.95","q.95","qub.95","qlbboot.95","qboot.95","qubboot.95",
          "qlb.99","q.99","qub.99","qlbboot.99","qboot.99","qubboot.99",
          "qlb.995","q.995","qub.995","qlbboot.995","qboot.995","qubboot.995"), names(lst1))
     for(ii in 1:length(lst1))
       {
        ## for the dependent testfreq selections:
   	    blocklen.tmp <- as.numeric(names(lst1))[ii] * 17
   	    gevfit.tmp <- lst1[[ii]][["dep.max"]]
     	res.dep[c("qlb.95","q.95","qub.95"), ii] <- rlevel.mine1(gevfit.tmp, k.blocks = 20*40000/blocklen.tmp)
    	res.dep[c("qlb.99","q.99","qub.99"), ii] <- rlevel.mine1(gevfit.tmp, k.blocks = 100*40000/blocklen.tmp)
    	res.dep[c("qlb.995","q.995","qub.995"), ii] <- rlevel.mine1(gevfit.tmp, k.blocks = 200*40000/blocklen.tmp)

   	    gevfit.tmp <- lst1[[ii]][["wdep.max"]]
     	res.wdep[c("qlb.95","q.95","qub.95"), ii] <- rlevel.mine1(gevfit.tmp, k.blocks = 20*40000/blocklen.tmp)
    	res.wdep[c("qlb.99","q.99","qub.99"), ii] <- rlevel.mine1(gevfit.tmp, k.blocks = 100*40000/blocklen.tmp)
    	res.wdep[c("qlb.995","q.995","qub.995"), ii] <- rlevel.mine1(gevfit.tmp, k.blocks = 200*40000/blocklen.tmp)

        ## for the simple random testfreq selections:
        ei <- 1/16
#        ei <- 1/8
   	    blocklen.tmp <- as.numeric(names(lst1))[ii]
   	    gevfit.tmp <- lst1[[ii]][["fou.max"]]
     	res.fou[c("qlb.95","q.95","qub.95"), ii] <- rlevel.mine1(gevfit.tmp, k.blocks = 20*ei*40000/blocklen.tmp)
    	res.fou[c("qlb.99","q.99","qub.99"), ii] <- rlevel.mine1(gevfit.tmp, k.blocks = 100*ei*40000/blocklen.tmp)
    	res.fou[c("qlb.995","q.995","qub.995"), ii] <- rlevel.mine1(gevfit.tmp, k.blocks = 200*ei*40000/blocklen.tmp)

   	    gevfit.tmp <- lst1[[ii]][["wfou.max"]]
     	res.wfou[c("qlb.95","q.95","qub.95"), ii] <- rlevel.mine1(gevfit.tmp, k.blocks = 20*ei*40000/blocklen.tmp)
    	res.wfou[c("qlb.99","q.99","qub.99"), ii] <- rlevel.mine1(gevfit.tmp, k.blocks = 100*ei*40000/blocklen.tmp)
    	res.wfou[c("qlb.995","q.995","qub.995"), ii] <- rlevel.mine1(gevfit.tmp, k.blocks = 200*ei*40000/blocklen.tmp)
    	}
     list(dep.max = res.dep, wdep.max = res.wdep, fou.max = res.fou, wfou.max = res.wfou)
    }))

names(growblockfull.quan2[[3]][[2]][[1]])

growblockfull.quan2[["mag.g6"]][["nights25"]][["dep.max"]]
growblockfull.quan[["mag.g6"]][["nights25"]][["dep.max"]]

names(growblockboot.gevpar[["mag.g6"]][["nights25"]][["dep.max"]])       
growblockboot.gevpar[["mag.g6"]][["nights25"]][["100"]][["dep.max"]]

for(ii in 1:length(growblockfullboot.quan))
   for(jj in 1:length(growblockfullboot.quan[[ii]]))
      for(kk in 1:length(growblockfullboot.quan[[ii]][[jj]]))
           {
           	nn.tmp <- names(growblockfullboot.quan[[ii]][[jj]])[kk]
           	growblockfull.quan2[[ii]][[jj]][["dep.max"]][c("qlbboot.95","qboot.95","qubboot.95"), nn.tmp] <- 
           	   quantile(growblockfullboot.quan[[ii]][[jj]][[nn.tmp]][["dep.max"]]["qboot.95",], prob = 
           	   c(0.025,0.5,0.975), na.rm = TRUE)
           	growblockfull.quan2[[ii]][[jj]][["dep.max"]][c("qlbboot.99","qboot.99","qubboot.99"), nn.tmp] <- 
           	   quantile(growblockfullboot.quan[[ii]][[jj]][[nn.tmp]][["dep.max"]]["qboot.99",], prob = 
           	   c(0.025,0.5,0.975), na.rm = TRUE)
           	growblockfull.quan2[[ii]][[jj]][["dep.max"]][c("qlbboot.995","qboot.995","qubboot.995"), nn.tmp] <- 
           	   quantile(growblockfullboot.quan[[ii]][[jj]][[nn.tmp]][["dep.max"]]["qboot.995",], prob = 
           	   c(0.025,0.5,0.975), na.rm = TRUE)
           	growblockfull.quan2[[ii]][[jj]][["wdep.max"]][c("qlbboot.95","qboot.95","qubboot.95"), nn.tmp] <- 
           	   quantile(growblockfullboot.quan[[ii]][[jj]][[nn.tmp]][["wdep.max"]]["qboot.95",], prob = 
           	   c(0.025,0.5,0.975), na.rm = TRUE)
           	growblockfull.quan2[[ii]][[jj]][["wdep.max"]][c("qlbboot.99","qboot.99","qubboot.99"), nn.tmp] <- 
           	   quantile(growblockfullboot.quan[[ii]][[jj]][[nn.tmp]][["wdep.max"]]["qboot.99",], prob = 
           	   c(0.025,0.5,0.975), na.rm = TRUE)
           	growblockfull.quan2[[ii]][[jj]][["wdep.max"]][c("qlbboot.995","qboot.995","qubboot.995"), nn.tmp] <- 
           	   quantile(growblockfullboot.quan[[ii]][[jj]][[nn.tmp]][["wdep.max"]]["qboot.995",], prob = 
           	   c(0.025,0.5,0.975), na.rm = TRUE)
           	growblockfull.quan2[[ii]][[jj]][["fou.max"]][c("qlbboot.95","qboot.95","qubboot.95"), nn.tmp] <- 
           	   quantile(growblockfullboot.quan[[ii]][[jj]][[nn.tmp]][["fou.max"]]["qboot.95",], prob = 
           	   c(0.025,0.5,0.975), na.rm = TRUE)
           	growblockfull.quan2[[ii]][[jj]][["fou.max"]][c("qlbboot.99","qboot.99","qubboot.99"), nn.tmp] <- 
           	   quantile(growblockfullboot.quan[[ii]][[jj]][[nn.tmp]][["fou.max"]]["qboot.99",], prob = 
           	   c(0.025,0.5,0.975), na.rm = TRUE)
           	growblockfull.quan2[[ii]][[jj]][["fou.max"]][c("qlbboot.995","qboot.995","qubboot.995"), nn.tmp] <- 
           	   quantile(growblockfullboot.quan[[ii]][[jj]][[nn.tmp]][["fou.max"]]["qboot.995",], prob = 
           	   c(0.025,0.5,0.975), na.rm = TRUE)
           	growblockfull.quan2[[ii]][[jj]][["wfou.max"]][c("qlbboot.95","qboot.95","qubboot.95"), nn.tmp] <- 
           	   quantile(growblockfullboot.quan[[ii]][[jj]][[nn.tmp]][["wfou.max"]]["qboot.95",], prob = 
           	   c(0.025,0.5,0.975), na.rm = TRUE)
           	growblockfull.quan2[[ii]][[jj]][["wfou.max"]][c("qlbboot.99","qboot.99","qubboot.99"), nn.tmp] <- 
           	   quantile(growblockfullboot.quan[[ii]][[jj]][[nn.tmp]][["wfou.max"]]["qboot.99",], prob = 
           	   c(0.025,0.5,0.975), na.rm = TRUE)
           	growblockfull.quan2[[ii]][[jj]][["wfou.max"]][c("qlbboot.995","qboot.995","qubboot.995"), nn.tmp] <- 
           	   quantile(growblockfullboot.quan[[ii]][[jj]][[nn.tmp]][["wfou.max"]]["qboot.995",], prob = 
           	   c(0.025,0.5,0.975), na.rm = TRUE)
          	}




## -----------------------------------------------------------------------------------
## High quantiles and their stability with the number of test frequencies,
## L = 100, 300, 500 only, with the profile likelihood CIs
## -----------------------------------------------------------------------------------


quartz(width = 12, height = 6)
colnnn <- "nights25"
mag.tmp <- "mag.g5"
nTestFr.tmp <- c("100", "300", "500")
par(mfrow = c(2,2), mar = c(3,3,2,1), mgp = c(2,0.75,0))

## unweighted:
plot(FirstStar.pgram[[mag.tmp]]$freq, 
   FirstStar.pgram[[mag.tmp]][,colnnn], type = "n",
   xlab = "Frequency", ylab = "Chi-squared", 
   main = paste(colnnn, mag.tmp, "Unweighted random"), 
   ylim = c(0.55,0.85),
   xlim = c(0,100))
#abline(v = F, col = "yellow")
points(FirstStar.pgram[[mag.tmp]]$freq, FirstStar.pgram[[mag.tmp]][,colnnn],
   type = "h", col = "grey45")
for(ii in 1:3)
    abline(h = growblockfull.quan[[mag.tmp]][[colnnn]][["fou.max"]][c("q.99","qlb.99","qub.99"),
       nTestFr.tmp[ii]], col = colorvec[ii], lty = c(1,2,2))
legend(80,0.5, bty = "n", legend = c("M = 100","M = 300","M = 500"), col = colorvec[1:3],
   lty = 1)

plot(FirstStar.pgram[[mag.tmp]]$freq, 
   FirstStar.pgram[[mag.tmp]][,colnnn], type = "n",
   xlab = "Frequency", ylab = "Chi-squared", 
   main = paste(colnnn, mag.tmp, "Unweighted dependent"), 
   ylim = c(0.55,0.85),
   xlim = c(0,100))
#abline(v = F, col = "yellow")
points(FirstStar.pgram[[mag.tmp]]$freq, FirstStar.pgram[[mag.tmp]][,colnnn],
   type = "h", col = "grey45")
for(ii in 1:3)
    abline(h = growblockfull.quan[[mag.tmp]][[colnnn]][["dep.max"]][c("q.99","qlb.99","qub.99"),
       nTestFr.tmp[ii]], col = colorvec[ii], lty = c(1,2,2))

## weighted:
plot(FirstStar.wpgram[[mag.tmp]]$freq, 
   FirstStar.wpgram[[mag.tmp]][,colnnn], type = "n",
   xlab = "Frequency", ylab = "Chi-squared", 
   main = paste(colnnn, mag.tmp, "Weighted random"), 
   ylim = c(0.55,0.85),
   xlim = c(0,100))
#abline(v = F, col = "yellow")
points(FirstStar.wpgram[[mag.tmp]]$freq, FirstStar.wpgram[[mag.tmp]][,colnnn],
   type = "h", col = "grey45")
for(ii in 1:3)
    abline(h = growblockfull.quan[[mag.tmp]][[colnnn]][["wfou.max"]][c("q.99","qlb.99","qub.99"),
       nTestFr.tmp[ii]], col = colorvec[ii], lty = c(1,2,2))
#legend(80,0.5, bty = "n", legend = c("M = 100","M = 300","M = 500"), col = colorvec[1:3],
#   lty = 1)

plot(FirstStar.wpgram[[mag.tmp]]$freq, 
   FirstStar.wpgram[[mag.tmp]][,colnnn], type = "n",
   xlab = "Frequency", ylab = "Chi-squared", 
   main = paste(colnnn, mag.tmp, "Weighted dependent"), 
#       ylim = c(0,max(max(FirstStar.bootII[[mag.tmp]][[colnnn]][[]][nSamplings[smp.tmp],], 
#       na.rm = TRUE), max(FirstStar.pgram[[mag.tmp]][,colnnn], na.rm = TRUE),
#       unlist(growblock.quan[[mag.tmp]][[colnnn]]))),
#   ylim = c(0.2,0.5),
   ylim = c(0.55,0.85),
   xlim = c(0,100))
#abline(v = F, col = "yellow")
points(FirstStar.wpgram[[mag.tmp]]$freq, FirstStar.wpgram[[mag.tmp]][,colnnn],
   type = "h", col = "grey45")
for(ii in 1:3)
    abline(h = growblockfull.quan[[mag.tmp]][[colnnn]][["wdep.max"]][c("q.99","qlb.99","qub.99"),
       nTestFr.tmp[ii]], col = colorvec[ii], lty = c(1,2,2))
    
## Some inspection of the 0.99 quantile lines: 
##  - For nights100, the quantiles estimated by the same frequency sampling but
##    different M are fairly together; only the M = 100 is sometimes a little off,
##    but the confidence intervals of different Ms are largely overlapping.
##  - For nights25, they are more scattered, for g3, there could be no CI estimated.
##    They are the least scattered for the lowest amplitude.
##  - For all, the quantiles obtained with dependent sampling pattern are slightly
##    higher than those obtained with extremal index theory + random pattern.

## Unfortunately, the profile likelihood is not always able to give a CI. Bootstrap 
## may be used?


## -----------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------
## Check of the high quantiles: 
## -----------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------



save.image("extrFAPFirstStar-backup.RData")

#mag.names2 <- c("mag.g3","mag.g5","mag.g6")
#mag.error.names2 <- c("mag.error.g3","mag.error.g5","mag.error.g6")

FirstStar.bootIII <- vector(3, mode = "list")
names(FirstStar.bootIII) <- mag.names2

for(jj in 2:3)
	{
	bootI.tmp <- vector(2, mode = "list")
	names(bootI.tmp) <- nNames2
	for(ii in 1:2)
 		{
  		 print(system.time(
   		 {
 		 nnn <- nObs2[ii]
  		 colnnn <- nNames2[ii]
 	     dfr <- FirstStar[FirstStar[,colnnn] == 1,]
 	     ts <- demean.fun(dfr$time)

 	     bootI.tmp[[nNames2[ii]]] <- replicate(400,
		   {  		  	
            ## The bootstrap:  		 
      	    ind <- sample(1:nrow(dfr), replace = TRUE, size = nrow(dfr))
      	    ps <- demean.fun(dfr[ind, mag.names2[jj]])
      	    ## Unweighted:
            ws <- rep(1/nnn, nnn)
            frsp.tmp <- sapply(FreqGrid1, MyFMLS3.fun, t = ts, h = ps,  
         		Weights = ws)
            ## Weighted:
            ws <- (1/dfr[ind, mag.error.names2[jj]]^2)/sum(1/dfr[ind, mag.error.names2[jj]]^2)
            wfrsp.tmp <- sapply(FreqGrid1, MyFMLS3.fun, t = ts, h = ps,  
         		Weights = ws)
		    c(
	          max = max(frsp.tmp[!is.na(frsp.tmp) & abs(frsp.tmp) != Inf]),
	          wmax = max(wfrsp.tmp[!is.na(wfrsp.tmp) & abs(wfrsp.tmp) != Inf])
  		      )
		   }) 
 		 }))
	  }
	FirstStar.bootIII[[mag.names2[jj]]] <- bootI.tmp
    save.image(WORKSPACE.NAME)
	cat(mag.names2[jj], " ready \n")
 	} 

save.image(WORKSPACE.NAME)


save.image("extrFAPFirstStar-backup.RData")

#mag.names2 <- c("mag.g3","mag.g5","mag.g6")
#mag.error.names2 <- c("mag.error.g3","mag.error.g5","mag.error.g6")

FirstStar.bootV <- vector(3, mode = "list")
names(FirstStar.bootV) <- mag.names2

for(jj in 1:3)
	{
	bootI.tmp <- vector(2, mode = "list")
	names(bootI.tmp) <- nNames2
	for(ii in 1:2)
 		{
  		 print(system.time(
   		 {
 		 nnn <- nObs2[ii]
  		 colnnn <- nNames2[ii]
 	     dfr <- FirstStar[FirstStar[,colnnn] == 1,]
 	     ts <- demean.fun(dfr$time)

 	     bootI.tmp[[nNames2[ii]]] <- replicate(1600,
		   {  		  	
            ## The bootstrap:  		 
      	    ind <- sample(1:nrow(dfr), replace = TRUE, size = nrow(dfr))
      	    ps <- demean.fun(dfr[ind, mag.names2[jj]])
      	    ## Unweighted:
            ws <- rep(1/nnn, nnn)
            frsp.tmp <- sapply(FreqGrid1, MyFMLS3.fun, t = ts, h = ps,  
         		Weights = ws)
            ## Weighted:
            ws <- (1/dfr[ind, mag.error.names2[jj]]^2)/sum(1/dfr[ind, mag.error.names2[jj]]^2)
            wfrsp.tmp <- sapply(FreqGrid1, MyFMLS3.fun, t = ts, h = ps,  
         		Weights = ws)
		    c(
	          max = max(frsp.tmp[!is.na(frsp.tmp) & abs(frsp.tmp) != Inf]),
	          wmax = max(wfrsp.tmp[!is.na(wfrsp.tmp) & abs(wfrsp.tmp) != Inf])
  		      )
		   }) 
 		 }))
	  }
	FirstStar.bootV[[mag.names2[jj]]] <- bootI.tmp
    save.image(WORKSPACE.NAME)
	cat(mag.names2[jj], " ready \n")
 	} 

save.image(WORKSPACE.NAME)

FirstStar.bootIII.backup <- FirstStar.bootIII
FirstStar.bootIII[["mag.g3"]][["nights100"]] <- 
  cbind(FirstStar.bootIII[["mag.g3"]][["nights100"]], 
        FirstStar.bootV[["mag.g3"]][["nights100"]]) 
FirstStar.bootIII[["mag.g5"]][["nights100"]] <- 
  cbind(FirstStar.bootIII[["mag.g5"]][["nights100"]], 
        FirstStar.bootV[["mag.g5"]][["nights100"]]) 
FirstStar.bootIII[["mag.g6"]][["nights100"]] <- 
  cbind(FirstStar.bootIII[["mag.g6"]][["nights100"]], 
        FirstStar.bootV[["mag.g6"]][["nights100"]]) 
FirstStar.bootIII[["mag.g3"]][["nights25"]] <- 
  cbind(FirstStar.bootIII[["mag.g3"]][["nights25"]], 
        FirstStar.bootV[["mag.g3"]][["nights25"]]) 
FirstStar.bootIII[["mag.g5"]][["nights25"]] <- 
  cbind(FirstStar.bootIII[["mag.g5"]][["nights25"]], 
        FirstStar.bootV[["mag.g5"]][["nights25"]]) 
FirstStar.bootIII[["mag.g6"]][["nights25"]] <- 
  cbind(FirstStar.bootIII[["mag.g6"]][["nights25"]], 
        FirstStar.bootV[["mag.g6"]][["nights25"]]) 


## Longer time series:

## 95%, non-weighted LS, dependent:
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g3"]][["nights100"]][["dep.max"]]["q.95",nn] > 
      FirstStar.bootIII[["mag.g3"]][["nights100"]][1,])/2000)
## M = 500: 94.95%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g5"]][["nights100"]][["dep.max"]]["q.95",nn] > 
      FirstStar.bootIII[["mag.g5"]][["nights100"]][1,])/2000)
## M = 500: 97.65%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g6"]][["nights100"]][["dep.max"]]["q.95",nn] > 
      FirstStar.bootIII[["mag.g6"]][["nights100"]][1,])/2000)
## M = 500: 96.15%

## 99%, non-weighted LS, dependent:
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g3"]][["nights100"]][["dep.max"]]["q.99",nn] > 
      FirstStar.bootIII[["mag.g3"]][["nights100"]][1,])/2000)
## M = 500: 98.85%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g5"]][["nights100"]][["dep.max"]]["q.99",nn] > 
      FirstStar.bootIII[["mag.g5"]][["nights100"]][1,])/2000)
## M = 500: 99.6%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g6"]][["nights100"]][["dep.max"]]["q.99",nn] > 
      FirstStar.bootIII[["mag.g6"]][["nights100"]][1,])/2000)
## M = 500: 99%

## 99.5%, non-weighted LS, dependent:
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g3"]][["nights100"]][["dep.max"]]["q.995",nn] > 
      FirstStar.bootIII[["mag.g3"]][["nights100"]][1,])/2000)
## M = 500: 99.6%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g5"]][["nights100"]][["dep.max"]]["q.995",nn] > 
      FirstStar.bootIII[["mag.g5"]][["nights100"]][1,])/2000)
## M = 500: 99.7%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g6"]][["nights100"]][["dep.max"]]["q.995",nn] > 
      FirstStar.bootIII[["mag.g6"]][["nights100"]][1,])/2000)
## M = 500: 99.5%
## In these simulation and sampling combinations, lower Ms tend to estimate too
## high quantiles.

## 95%, non-weighted LS, random:
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g3"]][["nights100"]][["fou.max"]]["q.95",nn] > 
      FirstStar.bootIII[["mag.g3"]][["nights100"]][1,])/2000)
## M = 500: 88%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g5"]][["nights100"]][["fou.max"]]["q.95",nn] > 
      FirstStar.bootIII[["mag.g5"]][["nights100"]][1,])/2000)
## M = 500: 94.45%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g6"]][["nights100"]][["fou.max"]]["q.95",nn] > 
      FirstStar.bootIII[["mag.g6"]][["nights100"]][1,])/2000)
## M = 500: 92.05%

## 99%, non-weighted LS, random:
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g3"]][["nights100"]][["fou.max"]]["q.99",nn] > 
      FirstStar.bootIII[["mag.g3"]][["nights100"]][1,])/2000)
## M = 500: 96.5%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g5"]][["nights100"]][["fou.max"]]["q.99",nn] > 
      FirstStar.bootIII[["mag.g5"]][["nights100"]][1,])/2000)
## M = 500: 98.6%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g6"]][["nights100"]][["fou.max"]]["q.99",nn] > 
      FirstStar.bootIII[["mag.g6"]][["nights100"]][1,])/2000)
## M = 500: 97.65%

## 99.5%, non-weighted LS, random:
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g3"]][["nights100"]][["fou.max"]]["q.995",nn] > 
      FirstStar.bootIII[["mag.g3"]][["nights100"]][1,])/2000)
## M = 500: 97.65%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g5"]][["nights100"]][["fou.max"]]["q.995",nn] > 
      FirstStar.bootIII[["mag.g5"]][["nights100"]][1,])/2000)
## M = 500: 99.4%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g6"]][["nights100"]][["fou.max"]]["q.995",nn] > 
      FirstStar.bootIII[["mag.g6"]][["nights100"]][1,])/2000)
## M = 500: 98.75%
## In these simulation and sampling combinations, lower Ms still tend to estimate 
## higher quantiles, but in this case, it is better: this method strongly under-
## estimates the quantiles... Pity, this would have been feasible for Gaia.


## 95%, weighted LS, dependent:
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g3"]][["nights100"]][["wdep.max"]]["q.95",nn] > 
      FirstStar.bootIII[["mag.g3"]][["nights100"]][2,])/2000)
## M = 500: 97.6%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g5"]][["nights100"]][["wdep.max"]]["q.95",nn] > 
      FirstStar.bootIII[["mag.g5"]][["nights100"]][2,])/2000)
## M = 500: 96.75%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g6"]][["nights100"]][["wdep.max"]]["q.95",nn] > 
      FirstStar.bootIII[["mag.g6"]][["nights100"]][2,])/2000)
## M = 500: 97.5%

## 99%, weighted LS, dependent:
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g3"]][["nights100"]][["wdep.max"]]["q.99",nn] > 
      FirstStar.bootIII[["mag.g3"]][["nights100"]][2,])/2000)
## M = 500: 99.5%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g5"]][["nights100"]][["wdep.max"]]["q.99",nn] > 
      FirstStar.bootIII[["mag.g5"]][["nights100"]][2,])/2000)
## M = 500: 99.1%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g6"]][["nights100"]][["wdep.max"]]["q.99",nn] > 
      FirstStar.bootIII[["mag.g6"]][["nights100"]][2,])/2000)
## M = 500: 99.05%

## 99.5%, weighted LS, dependent:
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g3"]][["nights100"]][["wdep.max"]]["q.995",nn] > 
      FirstStar.bootIII[["mag.g3"]][["nights100"]][2,])/2000)
## M = 500: 99.75%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g5"]][["nights100"]][["wdep.max"]]["q.995",nn] > 
      FirstStar.bootIII[["mag.g5"]][["nights100"]][2,])/2000)
## M = 500: 99.3%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g6"]][["nights100"]][["wdep.max"]]["q.995",nn] > 
      FirstStar.bootIII[["mag.g6"]][["nights100"]][2,])/2000)
## M = 500: 99.6%
## In these simulation and sampling combinations, lower Ms tend to estimate 
## higher quantiles. Especially for g3, it is over-estimating the 
## the quantiles; the over-estimation is smaller for the more noise-like 
## time series, which is good, the performance is better and the estimates
## more realistic for the cases when we indeed want to distinguish...
## Why it is g3 for which it overestimates? Strong signal ==> heavier tail
## ==> stronger long-memory character in the periodogram (?) ==> worse 
## approximation when I use independence in the final step? 

## 95%, weighted LS, random:
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g3"]][["nights100"]][["wfou.max"]]["q.95",nn] > 
      FirstStar.bootIII[["mag.g3"]][["nights100"]][2,])/2000)
## M = 500: 95.85%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g5"]][["nights100"]][["wfou.max"]]["q.95",nn] > 
      FirstStar.bootIII[["mag.g5"]][["nights100"]][2,])/2000)
## M = 500: 90.65%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g6"]][["nights100"]][["wfou.max"]]["q.95",nn] > 
      FirstStar.bootIII[["mag.g6"]][["nights100"]][2,])/2000)
## M = 500: 94.55%

## 99%, weighted LS, random:
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g3"]][["nights100"]][["wfou.max"]]["q.99",nn] > 
      FirstStar.bootIII[["mag.g3"]][["nights100"]][2,])/2000)
## M = 500: 98.55%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g5"]][["nights100"]][["wfou.max"]]["q.99",nn] > 
      FirstStar.bootIII[["mag.g5"]][["nights100"]][2,])/2000)
## M = 500: 95.25%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g6"]][["nights100"]][["wfou.max"]]["q.99",nn] > 
      FirstStar.bootIII[["mag.g6"]][["nights100"]][2,])/2000)
## M = 500: 98.1%

## 99.5%, weighted LS, random:
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g3"]][["nights100"]][["wfou.max"]]["q.995",nn] > 
      FirstStar.bootIII[["mag.g3"]][["nights100"]][2,])/2000)
## M = 500: 99.35%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g5"]][["nights100"]][["wfou.max"]]["q.995",nn] > 
      FirstStar.bootIII[["mag.g5"]][["nights100"]][2,])/2000)
## M = 500: 96.45%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g6"]][["nights100"]][["wfou.max"]]["q.995",nn] > 
      FirstStar.bootIII[["mag.g6"]][["nights100"]][2,])/2000)
## M = 500: 98.65%
## In these simulation and sampling combinations, lower Ms still tend to estimate 
## higher quantiles. This might be feasible for Gaia, it has lower computational
## requirements than the dependent sampling, and still not too bad.



## Shorter time series:


## 95%, non-weighted LS, dependent:
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g3"]][["nights25"]][["dep.max"]]["q.95",nn] > 
      FirstStar.bootIII[["mag.g3"]][["nights25"]][1,])/2000)
## M = 500: 95.85%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g5"]][["nights25"]][["dep.max"]]["q.95",nn] > 
      FirstStar.bootIII[["mag.g5"]][["nights25"]][1,])/2000)
## M = 500: 95.3%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g6"]][["nights25"]][["dep.max"]]["q.95",nn] > 
      FirstStar.bootIII[["mag.g6"]][["nights25"]][1,])/2000)
## M = 500: 96.35%

## 99%, non-weighted LS, dependent:
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g3"]][["nights25"]][["dep.max"]]["q.99",nn] > 
      FirstStar.bootIII[["mag.g3"]][["nights25"]][1,])/2000)
## M = 500: 98.55%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g5"]][["nights25"]][["dep.max"]]["q.99",nn] > 
      FirstStar.bootIII[["mag.g5"]][["nights25"]][1,])/2000)
## M = 500: 98.9%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g6"]][["nights25"]][["dep.max"]]["q.99",nn] > 
      FirstStar.bootIII[["mag.g6"]][["nights25"]][1,])/2000)
## M = 500: 99.2%

## 99.5%, non-weighted LS, dependent:
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g3"]][["nights25"]][["dep.max"]]["q.995",nn] > 
      FirstStar.bootIII[["mag.g3"]][["nights25"]][1,])/2000)
## M = 500: 99.1%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g5"]][["nights25"]][["dep.max"]]["q.995",nn] > 
      FirstStar.bootIII[["mag.g5"]][["nights25"]][1,])/2000)
## M = 500: 99.5%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g6"]][["nights25"]][["dep.max"]]["q.995",nn] > 
      FirstStar.bootIII[["mag.g6"]][["nights25"]][1,])/2000)
## M = 500: 99.45%
## Smaller M tend to perform better, here again they give higher quantiles, which
## is better for this case.

## 95%, non-weighted LS, random:
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g3"]][["nights25"]][["fou.max"]]["q.95",nn] > 
      FirstStar.bootIII[["mag.g3"]][["nights25"]][1,])/2000)
## M = 500: 92.9%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g5"]][["nights25"]][["fou.max"]]["q.95",nn] > 
      FirstStar.bootIII[["mag.g5"]][["nights25"]][1,])/2000)
## M = 500: 90.8%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g6"]][["nights25"]][["fou.max"]]["q.95",nn] > 
      FirstStar.bootIII[["mag.g6"]][["nights25"]][1,])/2000)
## M = 500: 94.45%

## 99%, non-weighted LS, random:
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g3"]][["nights25"]][["fou.max"]]["q.99",nn] > 
      FirstStar.bootIII[["mag.g3"]][["nights25"]][1,])/2000)
## M = 500: 97.85%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g5"]][["nights25"]][["fou.max"]]["q.99",nn] > 
      FirstStar.bootIII[["mag.g5"]][["nights25"]][1,])/2000)
## M = 500: 97.55%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g6"]][["nights25"]][["fou.max"]]["q.99",nn] > 
      FirstStar.bootIII[["mag.g6"]][["nights25"]][1,])/2000)
## M = 500: 98.9%

## 99.5%, non-weighted LS, random:
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g3"]][["nights25"]][["fou.max"]]["q.995",nn] > 
      FirstStar.bootIII[["mag.g3"]][["nights25"]][1,])/2000)
## M = 500: 98.7%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g5"]][["nights25"]][["fou.max"]]["q.995",nn] > 
      FirstStar.bootIII[["mag.g5"]][["nights25"]][1,])/2000)
## M = 500: 98.6%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g6"]][["nights25"]][["fou.max"]]["q.995",nn] > 
      FirstStar.bootIII[["mag.g6"]][["nights25"]][1,])/2000)
## M = 500: 99.2%
## Here too, this method  under-estimates the quantiles, though seems to me that
## less than for longer time series. 


## 95%, weighted LS, dependent:
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g3"]][["nights25"]][["wdep.max"]]["q.95",nn] > 
      FirstStar.bootIII[["mag.g3"]][["nights25"]][2,])/2000)
## M = 500: 95.3%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g5"]][["nights25"]][["wdep.max"]]["q.95",nn] > 
      FirstStar.bootIII[["mag.g5"]][["nights25"]][2,])/2000)
## M = 500: 94.95%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g6"]][["nights25"]][["wdep.max"]]["q.95",nn] > 
      FirstStar.bootIII[["mag.g6"]][["nights25"]][2,])/2000)
## M = 500: 97.15%

## 99%, weighted LS, dependent:
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g3"]][["nights25"]][["wdep.max"]]["q.99",nn] > 
      FirstStar.bootIII[["mag.g3"]][["nights25"]][2,])/2000)
## M = 500: 97.1%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g5"]][["nights25"]][["wdep.max"]]["q.99",nn] > 
      FirstStar.bootIII[["mag.g5"]][["nights25"]][2,])/2000)
## M = 500: 98.3%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g6"]][["nights25"]][["wdep.max"]]["q.99",nn] > 
      FirstStar.bootIII[["mag.g6"]][["nights25"]][2,])/2000)
## M = 500: 99.25%

## 99.5%, weighted LS, dependent:
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g3"]][["nights25"]][["wdep.max"]]["q.995",nn] > 
      FirstStar.bootIII[["mag.g3"]][["nights25"]][2,])/2000)
## M = 500: 97.5%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g5"]][["nights25"]][["wdep.max"]]["q.995",nn] > 
      FirstStar.bootIII[["mag.g5"]][["nights25"]][2,])/2000)
## M = 500: 98.95%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g6"]][["nights25"]][["wdep.max"]]["q.995",nn] > 
      FirstStar.bootIII[["mag.g6"]][["nights25"]][2,])/2000)
## M = 500: 99.7%


## 95%, weighted LS, random:
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g3"]][["nights25"]][["wfou.max"]]["q.95",nn] > 
      FirstStar.bootIII[["mag.g3"]][["nights25"]][2,])/2000)
## M = 500: 90.55%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g5"]][["nights25"]][["wfou.max"]]["q.95",nn] > 
      FirstStar.bootIII[["mag.g5"]][["nights25"]][2,])/2000)
## M = 500: 90.8%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g6"]][["nights25"]][["wfou.max"]]["q.95",nn] > 
      FirstStar.bootIII[["mag.g6"]][["nights25"]][2,])/2000)
## M = 500: 94%

## 99%, weighted LS, random:
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g3"]][["nights25"]][["wfou.max"]]["q.99",nn] > 
      FirstStar.bootIII[["mag.g3"]][["nights25"]][2,])/2000)
## M = 500: 95.95%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g5"]][["nights25"]][["wfou.max"]]["q.99",nn] > 
      FirstStar.bootIII[["mag.g5"]][["nights25"]][2,])/2000)
## M = 500: 95.65%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g6"]][["nights25"]][["wfou.max"]]["q.99",nn] > 
      FirstStar.bootIII[["mag.g6"]][["nights25"]][2,])/2000)
## M = 500: 98%

## 99.5%, weighted LS, random:
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g3"]][["nights25"]][["wfou.max"]]["q.995",nn] > 
      FirstStar.bootIII[["mag.g3"]][["nights25"]][2,])/2000)
## M = 500: 96.3%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g5"]][["nights25"]][["wfou.max"]]["q.995",nn] > 
      FirstStar.bootIII[["mag.g5"]][["nights25"]][2,])/2000)
## M = 500: 96.65%
for(nn in c("50","100","200","300","400","500"))
  print(sum(growblockfull.quan[["mag.g6"]][["nights25"]][["wfou.max"]]["q.995",nn] > 
      FirstStar.bootIII[["mag.g6"]][["nights25"]][2,])/2000)
## M = 500: 98.85%
## Underestimation. Pity.


## Plots 

nQuan <- c("q.95","q.99","q.995")
vQuan <- c(0.95, 0.99, 0.995)
felirat.tmp <- 
quartz(width = 9, height = 4)
#par(mfcol = c(3,3), mar = c(3,3,2,1), mgp = c(2,0.7,0))
par(mfcol = c(2,3), mar = c(3,3,2,1), mgp = c(2,0.7,0))
ii <- 1
for(jj in 1:3)
#  for(kk in 1:3)
  for(kk in 2:3)
    {
     plot(c(50,100,200,300,400,500), 
         sapply(c("50","100","200","300","400","500"), function(ss)
            sum(growblockfull.quan[[mag.names2[jj]]][[nNames2[ii]]][["wfou.max"]][nQuan[kk], ss] > 
         FirstStar.bootIII[[mag.names2[jj]]][[nNames2[ii]]][2,])/2000),
         type = "n", ylim = c(0.9,1),# ylog = TRUE, log = "y",
         ylab = "Log[F(q)]", xlab = "M", main = paste(mag.names2[jj], ",   ",
         vQuan[kk], " nominal,   ", nObs2[ii], " obs.", sep = ""))
     abline(h = vQuan[kk], lwd = 2)
     points(c(50,100,200,300,400,500), 
         sapply(c("50","100","200","300","400","500"), function(ss)
            sum(growblockfull.quan[[mag.names2[jj]]][[nNames2[ii]]][["wdep.max"]][nQuan[kk], ss] > 
         FirstStar.bootIII[[mag.names2[jj]]][[nNames2[ii]]][2,])/2000),
         type = "o", col = "darkred", pch = 15)
     points(c(50,100,200,300,400,500), 
         sapply(c("50","100","200","300","400","500"), function(ss)
            sum(growblockfull.quan[[mag.names2[jj]]][[nNames2[ii]]][["wfou.max"]][nQuan[kk], ss] > 
         FirstStar.bootIII[[mag.names2[jj]]][[nNames2[ii]]][2,])/2000),
         type = "o", col = "blue", pch = 16)
     points(c(50,100,200,300,400,500), 
         sapply(c("50","100","200","300","400","500"), function(ss)
            sum(growblockfull.quan[[mag.names2[jj]]][[nNames2[ii]]][["dep.max"]][nQuan[kk], ss] > 
         FirstStar.bootIII[[mag.names2[jj]]][[nNames2[ii]]][1,])/2000),
         type = "o", col = "red", pch = 17)
     points(c(50,100,200,300,400,500), 
         sapply(c("50","100","200","300","400","500"), function(ss)
            sum(growblockfull.quan[[mag.names2[jj]]][[nNames2[ii]]][["fou.max"]][nQuan[kk], ss] > 
         FirstStar.bootIII[[mag.names2[jj]]][[nNames2[ii]]][1,])/2000),
         type = "o", col = "deepskyblue", pch = 18)
    }
legend(200,0.97, legend = c("Weighted dependent","Weighted random",
  "Non-weighted dependent","Non-weighted random"), col = c("darkred","blue","red",
  "deepskyblue"), lty = 1, pch = c(15,16,17,18), bty = "n")

## Conclusions:
##  - No point in ruining the nice effect with showing the 0.95 level. For the others:
##  - Unfortunately, the dependent sampling is better, either with or without weights, 
##    the simple random sampling + use of theta seems to be unstable, now this, then
##    that is strongly biased.
##  - Sugar in the bitter: dep.sampling good is even at low number of tested frequencies
##    (50 or 100 is already okay), the lines seem to change almost nothing.
##  - Non-weighted dependent seems to be overall better, because the weighted 
##    dependent is bad for higher amplitudes combined with few observations; otherwise,
##    weighted dependent seems to be slightly better. 
##  - The non-weighted random is generally better among the random ones, but it still
##    tends to underestimate the quantile, and the underestimation is the worst for
##    high-amplitude, longer time series (which is not critical at all, as the signal is 
##    very strong with frequent observations, and so the peak is very big.)
growblock.quan[[1]][["nights25"]][["125"]]
bootI.tmp[["nights25"]]



## -----------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------
## Bootstrap CI for the high quantiles (M = 500)
## -----------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------



names(FirstStar.bootII)
names(FirstStar.bootII[[1]])
names(FirstStar.bootII[[1]][[1]])
names(FirstStar.bootII[[1]][[2]])

print(system.time(
gblboot.quan <- lapply(FirstStar.bootII[mag.names2], function(lst)
  lapply(lst, function(lst1)
    lapply(lst1, function(lst2)
     {
      res2.tmp <- vector(4, mode = "list")
      names(res2.tmp) <- c("dep.max","fou.max","wdep.max","wfou.max")
      for(ii in 1:4)
       	{
         vec <- lst2[ii,]
       	 res1.tmp <- replicate(1000, 
       	   {
       	   	vec1 <- sample(vec, replace = TRUE)
       	    res.tmp <- try(gev(vec1))
       	    parvec <- if(inherits(res.tmp, "try-error")) rep(NA, 3) else res.tmp$par.ests
       	    res1.tmp <- c(q95 = qgev(0.95, xi = parvec[1], mu = parvec[3], sigma = parvec[2]),
       	      q99 = qgev(0.99, xi = parvec[1], mu = parvec[3], sigma = parvec[2]),
       	      q995 = qgev(0.995, xi = parvec[1], mu = parvec[3], sigma = parvec[2]),
       	      q999 = qgev(0.999, xi = parvec[1], mu = parvec[3], sigma = parvec[2]))
       	   })
         res2.tmp[[ii]] <- apply(res1.tmp, MAR = 1, FUN = quantile, 
       	    probs = c(0.025, 0.5, 0.975), na.rm = TRUE) 
      	}
       res2.tmp
       })))

))


rlevel.gev(growblock.fgev[[1]][[1]][[4]][[2]], k.blocks = 5)

## -----------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------
## Stability of the fits against # of repetitions
## 0.99 quantile of the periodogram subset fits, CI from profile likelihood
## -----------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------
 

names(FirstStar.bootII[[1]][[2]])

R.ind <- c(200, 400, 600, 800, 1000)
randR.ind <- vector(4, mode = "list")
randR.ind[[5]] <- 1:1000
for(ii in 1:4)
   randR.ind[[ii]] <- sample(1:1000, size = R.ind[ii])

numRep.fgev <- vector(5, mode = "list")
names(numRep.fgev) <- c("R.200","R.400","R.600","R.800","R.1000")
for(ii in 1:4)
  {
   numRep.fgev[[ii]] <- lapply(FirstStar.bootII, function(lst)
     lapply(lst,  function(lst1)
       lapply(lst1, function(mat)  
         apply(mat, MAR = 1, FUN = function(vec)
       	   try(gev(vec[randR.ind[[ii]]]))))))
#       	   try(gev(vec[1:R.ind[ii]]))))))
  }   
numRep.fgev[["R.1000"]] <- growblock.fgev

numRep.gevpar <- vector(5, mode = "list")
names(numRep.gevpar) <- c("R.200","R.400","R.600","R.800","R.1000")
for(ii in 1:5)
  {
   numRep.gevpar[[ii]] <- lapply(numRep.fgev[[ii]], function(lst3)
     lapply(lst3, function(lst)
        lapply(lst, function(lst2)
           sapply(lst2, function(lst1) 
              if(!inherits(lst1, "try-error")) lst1$par.ests else rep(NA, 3)))))
   }
numRep.gevpar[["R.1000"]] <- growblock.gevpar

numRep.gevse <- vector(5, mode = "list")
names(numRep.gevse) <- c("R.200","R.400","R.600","R.800","R.1000")
for(ii in 1:5)
  {
   numRep.gevse[[ii]] <- lapply(numRep.fgev[[ii]], function(lst3)
    lapply(lst3, function(lst)
      lapply(lst, function(lst2)
        sapply(lst2, function(lst1) 
          if(!inherits(lst1, "try-error")) lst1$par.ses else rep(NA, 3)))))
   }
numRep.gevse[["R.1000"]] <- growblock.gevse

numRep.quan <- vector(5, mode = "list")
names(numRep.quan) <- c("R.200","R.400","R.600","R.800","R.1000")
for(ii in 1:5)
  {
   numRep.quan[[ii]] <- lapply(numRep.gevpar[[ii]], function(lst)
    lapply(lst, function(lst1)
     lapply(lst1, function(mat)
      apply(mat, MAR = 2, function(parvec) 
       c(q.95 = qgev(0.95, xi = parvec[1], mu = parvec[3], sigma = parvec[2]),
         q.99 = qgev(0.99, xi = parvec[1], mu = parvec[3], sigma = parvec[2]),
         q.995 = qgev(0.995, xi = parvec[1], mu = parvec[3], sigma = parvec[2]),
         q.999 = qgev(0.999, xi = parvec[1], mu = parvec[3], sigma = parvec[2]),
         q.9999 = qgev(0.9999, xi = parvec[1], mu = parvec[3], sigma = parvec[2]))))))
   }
numRep.quan[["R.1000"]] <- growblock.quan

numRep.rlev <- vector(5, mode = "list")
names(numRep.rlev) <- c("R.200","R.400","R.600","R.800","R.1000")
for(ii in 1:5)
  {
   numRep.rlev[[ii]] <- lapply(numRep.fgev[[ii]], function(lst3)
     lapply(lst3, function(lst)
        lapply(lst, function(lst2)
           sapply(lst2, function(lst1) 
              if(!inherits(lst1, "try-error")) rlevel.gev(lst1, k.blocks = 100) else rep(NA, 3)))))
   }



## Shapes, scales, locations and quantiles against R.ind:

names(numRep.gevpar)
#[1] "R.200"  "R.400"  "R.600"  "R.800"  "R.1000"
names(numRep.gevpar[[1]])
#[1] "mag.g3" "mag.g5" "mag.g6"
names(numRep.gevpar[[1]][[1]])
#[1] "nights100" "nights25" 
names(numRep.gevpar[[1]][[1]][[1]])
#[1] "500" "400" "300" "200" "100" "50" 
dimnames(numRep.gevpar[[1]][[1]][[1]][[1]])
#[[1]]
#[1] "xi"    "sigma" "mu"   
#[[2]]
#[1] "dep.max"  "fou.max"  "wdep.max" "wfou.max"


rng.tmp <- matrix(c(-0.15, 0.05,
                    0, 0.1,
                    0.2, 0.5), ncol = 2, byrow = TRUE)
rng.tmp <- matrix(c(-0.15, 0.05,
                    0.01, 0.03,
                    0.12, 0.15), ncol = 2, byrow = TRUE)
partype.tmp <- c("Shape","Scale","Location")

cols.tmp <- "nights100"
nMax.tmp <- 1   ## For nights100, this means 500 maxima, for nights25, 125 

#dev.set(2)
quartz(height = 7, width = 9.3)
par(mfrow = c(3,4), mar = c(3,3,2,1), mgp = c(2,1,0))
for(ii in 1:3)            ## ii: index through magnitude
  {
   par2.tmp <- vector(3, mode = "list")
   se2.tmp <- vector(3, mode = "list")
   quan2.tmp <- vector(3, mode = "list")
   for(kk in 1:3) 
    {
     par1.tmp <- matrix(ncol = 5, nrow = 3)
     se1.tmp <- matrix(ncol = 5, nrow = 3)
     quan1.tmp <- numeric(5)
     for(jj in 1:5)         ## jj: index through the number of repetitions
        {
        par1.tmp[,jj] <- numRep.gevpar[[jj]][[ii]][[cols.tmp]][[kk]][,1]  
        se1.tmp[,jj] <- numRep.gevse[[jj]][[ii]][[cols.tmp]][[kk]][,1]  
        quan1.tmp[jj] <- numRep.quan[[jj]][[ii]][[cols.tmp]][[kk]][2,1]  
        }
      par2.tmp[[kk]] <- par1.tmp
      se2.tmp[[kk]] <- se1.tmp
      quan2.tmp[[kk]] <- quan1.tmp
     }
    for(ll in 1:3)        ## ll: index through the gev parameters
      {
       plot(R.ind, par1.tmp[ll,], xlim = c(100,1100),
         ylim = rng.tmp[ll, ], main = paste(mag.names2[ii], sep = ""), 
         ylab = partype.tmp[ll], xlab = "R", type = "n")
       points(R.ind-30, par2.tmp[[1]][ll,], pch = 16, col = "black")       
       segments(R.ind-30, par2.tmp[[1]][ll,] + 2*se2.tmp[[1]][ll,],
         R.ind-30, par2.tmp[[1]][ll,] - 2*se2.tmp[[1]][ll,], col = "black")
       points(R.ind, par2.tmp[[2]][ll,], pch = 16, col = "red")       
       segments(R.ind, par2.tmp[[2]][ll,] + 2*se2.tmp[[2]][ll,],
         R.ind, par2.tmp[[2]][ll,] - 2*se2.tmp[[2]][ll,], col = "red")
       points(R.ind+30, par2.tmp[[3]][ll,], pch = 16, col = "deepskyblue")       
       segments(R.ind+30, par2.tmp[[3]][ll,] + 2*se2.tmp[[3]][ll,],
         R.ind+30, par2.tmp[[3]][ll,] - 2*se2.tmp[[3]][ll,], col = "deepskyblue")
      }
    plot(R.ind-30, quan2.tmp[[1]], xlim = c(100,1100),
         ylim = c(0,1), main = paste(mag.names2[ii], sep = ""), 
         ylab = "0.99 quantile", xlab = "R", pch = "-", cex = 2)
    points(R.ind, quan2.tmp[[2]], pch = "-", col = "red", cex = 2)
    points(R.ind+30, quan2.tmp[[3]], pch = "-", col = "deepskyblue", cex = 2)
   }
  
 

## -----------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------
## Stability of the fits against # of repetitions
## 0.95,0.99,0.995 quantile of the whole periodogram, CI from bootstrap
## -----------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------
 


names(FirstStar.bootII)     	  
names(FirstStar.bootII[[ii]])     	  
names(FirstStar.bootII[[ii]][[jj]])     	  

growrepboot.gevpar <- vector(5, mode = "list")
names(growrepboot.gevpar) <- c("R.200","R.400","R.600","R.800","R.1000")
growrepfullboot.quan <- vector(5, mode = "list")
names(growrepfullboot.quan) <- c("R.200","R.400","R.600","R.800","R.1000")
system.time(
#for(ll in 1:4)   ## ll through the nb of repetitions
for(ll in 5)   ## ll through the nb of repetitions
	{
	growrepboot.gevpar.tmp <- vector(length(FirstStar.bootII), mode = "list")
	names(growrepboot.gevpar.tmp) <- names(FirstStar.bootII)
	growrepfullboot.quan.tmp <- vector(length(FirstStar.bootII), mode = "list")
	names(growrepfullboot.quan.tmp) <- names(FirstStar.bootII)
	for(ii in 1:length(FirstStar.bootII))    # ii over the magnitudes
#	for(ii in 1)
	   {
	    v1.tmp <- vector(length(FirstStar.bootII[[ii]]), mode = "list")
	    names(v1.tmp) <- names(FirstStar.bootII[[ii]])
	    w1.tmp <- vector(length(FirstStar.bootII[[ii]]), mode = "list")
	    names(w1.tmp) <- names(FirstStar.bootII[[ii]])
	    for(jj in 1:length(FirstStar.bootII[[ii]]))    # jj over the ts lengths
#	    for(jj in 1)
	       {
	        v2.tmp <- vector(length(FirstStar.bootII[[ii]][[jj]]), mode = "list")
	        names(v2.tmp) <- names(FirstStar.bootII[[ii]][[jj]])
	        w2.tmp <- vector(length(FirstStar.bootII[[ii]][[jj]]), mode = "list")
	        names(w2.tmp) <- names(FirstStar.bootII[[ii]][[jj]])
	        nameset <- which(is.element(names(FirstStar.bootII[[ii]][[jj]]), c("100","500")))
#	        for(kk in nameset)
	        for(kk in 1:length(FirstStar.bootII[[ii]][[jj]]))    # loop over the block size
	           {
	            v3.tmp <- vector(nrow(FirstStar.bootII[[ii]][[jj]][[kk]]), mode = "list")
	            names(v3.tmp) <- dimnames(FirstStar.bootII[[ii]][[jj]][[kk]])[[1]]
	            w3.tmp <- vector(nrow(FirstStar.bootII[[ii]][[jj]][[kk]]), mode = "list")
	            names(w3.tmp) <- dimnames(FirstStar.bootII[[ii]][[jj]][[kk]])[[1]]
	            nn.tmp <- names(FirstStar.bootII[[ii]][[jj]])[kk]
	     	    blocklen.tmp <- as.numeric(nn.tmp)
	            # for dependent maxima:
	            RR <- 1000
	            res.tmp <- replicate(RR, {
	            	   vec1 <- sample(FirstStar.bootII[[ii]][[jj]][[kk]]["dep.max", randR.ind[[ll]]], replace = TRUE);
	            	   try(gev.mine(vec1, init.par = growblock.gevpar[[ii]][[jj]][[kk]][,"dep.max"]))
	                   }, simplify = FALSE)
	            v3.tmp[["dep.max"]] <- sapply(res.tmp, FUN = function(lst)
	                   if(!inherits(lst, "try-error")) lst$par.ests else rep(NA, 3))
	            w3.tmp[["dep.max"]] <- sapply(res.tmp, FUN = function(lst)
	                   if(inherits(lst, "try-error")) rep(NA, 3) else 
	                   c(qboot.95 = rlevel.mine(lst, k.blocks = 20*40000/(17*blocklen.tmp)),
	                   qboot.99 = rlevel.mine(lst, k.blocks = 100*40000/(17*blocklen.tmp)),
	                   qboot.995 = rlevel.mine(lst, k.blocks = 200*40000/(17*blocklen.tmp))))
	            # for weighted dependent maxima:
	            res.tmp <- replicate(RR, {
	            	   vec1 <- sample(FirstStar.bootII[[ii]][[jj]][[kk]]["wdep.max", randR.ind[[ll]]], replace = TRUE);       		               try(gev.mine(vec1, init.par = growblock.gevpar[[ii]][[jj]][[kk]][,"wdep.max"]))
	                   }, simplify = FALSE)
	            v3.tmp[["wdep.max"]] <- sapply(res.tmp, FUN = function(lst)
	                   if(!inherits(lst, "try-error")) lst$par.ests else rep(NA, 3))
	            w3.tmp[["wdep.max"]] <- sapply(res.tmp, FUN = function(lst)
	                   if(inherits(lst, "try-error")) rep(NA, 3) else 
	                   c(qboot.95 = rlevel.mine(lst, k.blocks = 20*40000/(17*blocklen.tmp)),
	                   qboot.99 = rlevel.mine(lst, k.blocks = 100*40000/(17*blocklen.tmp)),
	                   qboot.995 = rlevel.mine(lst, k.blocks = 200*40000/(17*blocklen.tmp))))
	            # for fourier maxima:
	            ei <- 1/16
	            res.tmp <- replicate(RR, { 
	            	   vec1 <- sample(FirstStar.bootII[[ii]][[jj]][[kk]]["fou.max", randR.ind[[ll]]], replace = TRUE);       		               try(gev.mine(vec1, init.par = growblock.gevpar[[ii]][[jj]][[kk]][,"fou.max"]))
	                   }, simplify = FALSE)
	            v3.tmp[["fou.max"]] <- sapply(res.tmp, FUN = function(lst)
	                   if(!inherits(lst, "try-error")) lst$par.ests else rep(NA, 3))
	            w3.tmp[["fou.max"]] <- sapply(res.tmp, FUN = function(lst)
	                   if(inherits(lst, "try-error")) rep(NA, 3) else 
	                   c(qboot.95 = rlevel.mine(lst, k.blocks = 20*ei*40000/blocklen.tmp),
	                   qboot.99 = rlevel.mine(lst, k.blocks = 100*ei*40000/blocklen.tmp),
	                   qboot.995 = rlevel.mine(lst, k.blocks = 200*ei*40000/blocklen.tmp)))
	            # for weighted fourier maxima:
	            res.tmp <- replicate(RR, { 
	            	   vec1 <- sample(FirstStar.bootII[[ii]][[jj]][[kk]]["wfou.max", randR.ind[[ll]]], replace = TRUE);      		               try(gev.mine(vec1, init.par = growblock.gevpar[[ii]][[jj]][[kk]][,"wfou.max"]))
	                   }, simplify = FALSE)
	            v3.tmp[["wfou.max"]] <- sapply(res.tmp, FUN = function(lst)
	                   if(!inherits(lst, "try-error")) lst$par.ests else rep(NA, 3))
	            w3.tmp[["wfou.max"]] <- sapply(res.tmp, FUN = function(lst)
	                   if(inherits(lst, "try-error")) rep(NA, 3) else 
	                   c(qboot.95 = rlevel.mine(lst, k.blocks = 20*ei*40000/blocklen.tmp),
	                   qboot.99 = rlevel.mine(lst, k.blocks = 100*ei*40000/blocklen.tmp),
	                   qboot.995 = rlevel.mine(lst, k.blocks = 200*ei*40000/blocklen.tmp)))
	            v2.tmp[[kk]] <- v3.tmp
	            w2.tmp[[kk]] <- w3.tmp
	#            save.image("extrFAP7FirstStar-tmp.RData")
	            cat("ll=", ll, "ii=", ii, "jj=", jj, "kk=",kk, "ready \n")
	           }
	        v1.tmp[[jj]] <- v2.tmp
	        w1.tmp[[jj]] <- w2.tmp
	        cat("ll=", ll, "ii=",ii, "jj=", jj, "ready \n")
	       }
	    growrepboot.gevpar.tmp[[ii]] <- v1.tmp
	    growrepfullboot.quan.tmp[[ii]] <- w1.tmp
	    cat("ll=", ll, "ii=",ii, "ready \n")
	   }
    growrepboot.gevpar[[ll]] <- growrepboot.gevpar.tmp
    growrepfullboot.quan[[ll]] <- growrepfullboot.quan.tmp
    cat("ll=",ll, "ready \n")
   })


growrepboot.gevpar[[1]][[1]][[1]]   
growrepfull.quan2 <- lapply(numRep.fgev, function(lst)
# lst <- numRep.fgev[[1]] 
 lapply(lst, function(lst2)
# lst2 <- numRep.fgev[[1]][[1]] 
  lapply(lst2, function(lst1)
    {
#     lst1 <- lst2[[2]]
     res.dep <- matrix(ncol = length(lst1), nrow = 18)
     dimnames(res.dep) <- list(c("qlb.95","q.95","qub.95","qlbboot.95","qboot.95","qubboot.95",
          "qlb.99","q.99","qub.99","qlbboot.99","qboot.99","qubboot.99",
          "qlb.995","q.995","qub.995","qlbboot.995","qboot.995","qubboot.995"), names(lst1))
     res.wdep <- matrix(ncol = length(lst1), nrow = 18)
     dimnames(res.wdep) <- list(c("qlb.95","q.95","qub.95","qlbboot.95","qboot.95","qubboot.95",
          "qlb.99","q.99","qub.99","qlbboot.99","qboot.99","qubboot.99",
          "qlb.995","q.995","qub.995","qlbboot.995","qboot.995","qubboot.995"), names(lst1))
     res.fou <- matrix(ncol = length(lst1), nrow = 18)
     dimnames(res.fou) <- list(c("qlb.95","q.95","qub.95","qlbboot.95","qboot.95","qubboot.95",
          "qlb.99","q.99","qub.99","qlbboot.99","qboot.99","qubboot.99",
          "qlb.995","q.995","qub.995","qlbboot.995","qboot.995","qubboot.995"), names(lst1))
     res.wfou <- matrix(ncol = length(lst1), nrow = 18)
     dimnames(res.wfou) <- list(c("qlb.95","q.95","qub.95","qlbboot.95","qboot.95","qubboot.95",
          "qlb.99","q.99","qub.99","qlbboot.99","qboot.99","qubboot.99",
          "qlb.995","q.995","qub.995","qlbboot.995","qboot.995","qubboot.995"), names(lst1))
     for(ii in 1:length(lst1))
       {
#       	ii <- 1
        ## for the dependent testfreq selections:
   	    blocklen.tmp <- as.numeric(names(lst1))[ii] * 17
   	    gevfit.tmp <- lst1[[ii]][["dep.max"]]
     	res.dep[c("qlb.95","q.95","qub.95"), ii] <- rlevel.mine1(gevfit.tmp, k.blocks = 20*40000/blocklen.tmp)
    	res.dep[c("qlb.99","q.99","qub.99"), ii] <- rlevel.mine1(gevfit.tmp, k.blocks = 100*40000/blocklen.tmp)
    	res.dep[c("qlb.995","q.995","qub.995"), ii] <- rlevel.mine1(gevfit.tmp, k.blocks = 200*40000/blocklen.tmp)

   	    gevfit.tmp <- lst1[[ii]][["wdep.max"]]
     	res.wdep[c("qlb.95","q.95","qub.95"), ii] <- rlevel.mine1(gevfit.tmp, k.blocks = 20*40000/blocklen.tmp)
    	res.wdep[c("qlb.99","q.99","qub.99"), ii] <- rlevel.mine1(gevfit.tmp, k.blocks = 100*40000/blocklen.tmp)
    	res.wdep[c("qlb.995","q.995","qub.995"), ii] <- rlevel.mine1(gevfit.tmp, k.blocks = 200*40000/blocklen.tmp)

        ## for the simple random testfreq selections:
        ei <- 1/16
#        ei <- 1/8
   	    blocklen.tmp <- as.numeric(names(lst1))[ii]
   	    gevfit.tmp <- lst1[[ii]][["fou.max"]]
     	res.fou[c("qlb.95","q.95","qub.95"), ii] <- rlevel.mine1(gevfit.tmp, k.blocks = 20*ei*40000/blocklen.tmp)
    	res.fou[c("qlb.99","q.99","qub.99"), ii] <- rlevel.mine1(gevfit.tmp, k.blocks = 100*ei*40000/blocklen.tmp)
    	res.fou[c("qlb.995","q.995","qub.995"), ii] <- rlevel.mine1(gevfit.tmp, k.blocks = 200*ei*40000/blocklen.tmp)

   	    gevfit.tmp <- lst1[[ii]][["wfou.max"]]
     	res.wfou[c("qlb.95","q.95","qub.95"), ii] <- rlevel.mine1(gevfit.tmp, k.blocks = 20*ei*40000/blocklen.tmp)
    	res.wfou[c("qlb.99","q.99","qub.99"), ii] <- rlevel.mine1(gevfit.tmp, k.blocks = 100*ei*40000/blocklen.tmp)
    	res.wfou[c("qlb.995","q.995","qub.995"), ii] <- rlevel.mine1(gevfit.tmp, k.blocks = 200*ei*40000/blocklen.tmp)
    	}
     list(dep.max = res.dep, wdep.max = res.wdep, fou.max = res.fou, wfou.max = res.wfou)
    })))

names(growrepfull.quan2[[3]][[2]][[1]])

growrepfull.quan2[["R.200"]][["mag.g6"]][["nights25"]][["dep.max"]]
numRep.rlev[["R.200"]][["mag.g6"]][["nights25"]][["500"]]
numRep.rlev[["R.200"]][["mag.g6"]][["nights25"]][["100"]]

names(growrepboot.gevpar)       
names(growrepboot.gevpar[["R.200"]])       
names(growrepboot.gevpar[["R.200"]][["mag.g6"]])       
names(growrepboot.gevpar[["R.200"]][["mag.g6"]][["nights25"]])       
names(growrepboot.gevpar[["R.200"]][["mag.g6"]][["nights25"]][["100"]])       
names(growrepboot.gevpar[["R.200"]][["mag.g6"]][["nights25"]][["100"]][["dep.max"]])       
growrepboot.gevpar[["R.200"]][["mag.g6"]][["nights25"]][["100"]][["dep.max"]]
growrepfullboot.quan[["R.200"]][["mag.g6"]][["nights25"]][["100"]][["dep.max"]]

for(ll in 1:length(growrepfullboot.quan))
#ll <- 1
  for(ii in 1:length(growrepfullboot.quan[[ll]]))
#ii <- 3
    for(jj in 1:length(growrepfullboot.quan[[ll]][[ii]]))
#jj <- 1
      for(kk in 1:length(growrepfullboot.quan[[ll]][[ii]][[jj]]))
#kk <- 5
           {
           	nn.tmp <- names(growrepfullboot.quan[[ll]][[ii]][[jj]])[kk]
           	growrepfull.quan2[[ll]][[ii]][[jj]][["dep.max"]][c("qlbboot.95","qboot.95","qubboot.95"), nn.tmp] <- 
           	   quantile(growrepfullboot.quan[[ll]][[ii]][[jj]][[nn.tmp]][["dep.max"]]["qboot.95",], prob = 
           	   c(0.025,0.5,0.975), na.rm = TRUE)
           	growrepfull.quan2[[ll]][[ii]][[jj]][["dep.max"]][c("qlbboot.99","qboot.99","qubboot.99"), nn.tmp] <- 
           	   quantile(growrepfullboot.quan[[ll]][[ii]][[jj]][[nn.tmp]][["dep.max"]]["qboot.99",], prob = 
           	   c(0.025,0.5,0.975), na.rm = TRUE)
           	growrepfull.quan2[[ll]][[ii]][[jj]][["dep.max"]][c("qlbboot.995","qboot.995","qubboot.995"), nn.tmp] <- 
           	   quantile(growrepfullboot.quan[[ll]][[ii]][[jj]][[nn.tmp]][["dep.max"]]["qboot.995",], prob = 
           	   c(0.025,0.5,0.975), na.rm = TRUE)
           	growrepfull.quan2[[ll]][[ii]][[jj]][["wdep.max"]][c("qlbboot.95","qboot.95","qubboot.95"), nn.tmp] <- 
           	   quantile(growrepfullboot.quan[[ll]][[ii]][[jj]][[nn.tmp]][["wdep.max"]]["qboot.95",], prob = 
           	   c(0.025,0.5,0.975), na.rm = TRUE)
           	growrepfull.quan2[[ll]][[ii]][[jj]][["wdep.max"]][c("qlbboot.99","qboot.99","qubboot.99"), nn.tmp] <- 
           	   quantile(growrepfullboot.quan[[ll]][[ii]][[jj]][[nn.tmp]][["wdep.max"]]["qboot.99",], prob = 
           	   c(0.025,0.5,0.975), na.rm = TRUE)
           	growrepfull.quan2[[ll]][[ii]][[jj]][["wdep.max"]][c("qlbboot.995","qboot.995","qubboot.995"), nn.tmp] <- 
           	   quantile(growrepfullboot.quan[[ll]][[ii]][[jj]][[nn.tmp]][["wdep.max"]]["qboot.995",], prob = 
           	   c(0.025,0.5,0.975), na.rm = TRUE)
           	growrepfull.quan2[[ll]][[ii]][[jj]][["fou.max"]][c("qlbboot.95","qboot.95","qubboot.95"), nn.tmp] <- 
           	   quantile(growrepfullboot.quan[[ll]][[ii]][[jj]][[nn.tmp]][["fou.max"]]["qboot.95",], prob = 
           	   c(0.025,0.5,0.975), na.rm = TRUE)
           	growrepfull.quan2[[ll]][[ii]][[jj]][["fou.max"]][c("qlbboot.99","qboot.99","qubboot.99"), nn.tmp] <- 
           	   quantile(growrepfullboot.quan[[ll]][[ii]][[jj]][[nn.tmp]][["fou.max"]]["qboot.99",], prob = 
           	   c(0.025,0.5,0.975), na.rm = TRUE)
           	growrepfull.quan2[[ll]][[ii]][[jj]][["fou.max"]][c("qlbboot.995","qboot.995","qubboot.995"), nn.tmp] <- 
           	   quantile(growrepfullboot.quan[[ll]][[ii]][[jj]][[nn.tmp]][["fou.max"]]["qboot.995",], prob = 
           	   c(0.025,0.5,0.975), na.rm = TRUE)
           	growrepfull.quan2[[ll]][[ii]][[jj]][["wfou.max"]][c("qlbboot.95","qboot.95","qubboot.95"), nn.tmp] <- 
           	   quantile(growrepfullboot.quan[[ll]][[ii]][[jj]][[nn.tmp]][["wfou.max"]]["qboot.95",], prob = 
           	   c(0.025,0.5,0.975), na.rm = TRUE)
           	growrepfull.quan2[[ll]][[ii]][[jj]][["wfou.max"]][c("qlbboot.99","qboot.99","qubboot.99"), nn.tmp] <- 
           	   quantile(growrepfullboot.quan[[ll]][[ii]][[jj]][[nn.tmp]][["wfou.max"]]["qboot.99",], prob = 
           	   c(0.025,0.5,0.975), na.rm = TRUE)
           	growrepfull.quan2[[ll]][[ii]][[jj]][["wfou.max"]][c("qlbboot.995","qboot.995","qubboot.995"), nn.tmp] <- 
           	   quantile(growrepfullboot.quan[[ll]][[ii]][[jj]][[nn.tmp]][["wfou.max"]]["qboot.995",], prob = 
           	   c(0.025,0.5,0.975), na.rm = TRUE)
          	}

































plot(seq(0, 1, by = 0.01), pbeta(seq(0, 1, by = 0.01), 
    shape1 = 10, shape2 = 2), type = "l", 
    xlim = c(0,1), ylim = c(0,1), ylab = "", xlab = "")
lines(seq(0, 1, by = 0.01), pbeta(seq(0, 1, by = 0.01), 
    shape1 = 10, shape2 = 2)^5, col = "deepskyblue")
lines(seq(0, 1, by = 0.01), pbeta(seq(0, 1, by = 0.01), 
    shape1 = 10, shape2 = 2)^25, col = "blue")
lines(seq(0, 1, by = 0.01), pbeta(seq(0, 1, by = 0.01), 
    shape1 = 10, shape2 = 2)^100, col = "red")
abline(h=0:32)
abline(v=0:32)




plot(1:32, 1:32, type = "n", xlim = c(0,32), ylim = c(0,32), ylab = "", xlab = "")
abline(h=0:32)
abline(v=0:32)




