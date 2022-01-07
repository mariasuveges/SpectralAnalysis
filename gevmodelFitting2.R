## 
## 
## GOAL OF THE WORKSPACE:
##
##                                  CONSTRAINED OPTIMIZATION OF THE GEV:
##                                         MU = 1 + SIGMA / XI
##


WORKSPACE.NAME <- "gevmodelFitting2.RData"
save.image(WORKSPACE.NAME)

load("results_baluev1/pars_and_maxima_dfr.RObj")

## --------------------------------------------------------------------------------------------------------
## Functions
## --------------------------------------------------------------------------------------------------------


boundaryconstr.gev.fun <- function(pars, z)
   {
   	m <- length(z)
   	xi <- pars[1]
   	sig <- pars[2]
   	locsc <- 1 + xi*(z - 1 - sig/xi) / sig
   	locsc[locsc <= 0] <- 1e-120
   	- m*log(sig) - (1 + 1/xi) * sum(log(locsc)) - sum(locsc^(-1/xi))
   }
   
 
## --------------------------------------------------------------------------------------------------------
## Application
## --------------------------------------------------------------------------------------------------------

## --------------------------------------------------------------------------------------------------------
## Try out

z.tmp <- maxima.train.dfr[1001:1050, 1]
b.tmp <- optim(c(-0.01, sqrt(6*var(z.tmp))/pi), boundaryconstr.gev.fun,
   #lower = c(-0.5, 0.0000001), upper = c(0,Inf),
   z = z.tmp, control = list(fnscale = -1, maxit = 10000),
   hessian = TRUE, method = "BFGS")
n.tmp <- optim(c(-0.01, sqrt(6*var(z.tmp))/pi), boundaryconstr.gev.fun,
   #lower = c(-0.5, 0.0000001), upper = c(0,Inf),
   z = z.tmp, control = list(fnscale = -1, maxit = 10000),
   hessian = TRUE, method = "Nelder-Mead")
## Two different solutions. Which is true? $value is higher with Nelder-Mead,
## that is a bigger likelihood maximum 
sqrt(diag(-solve(n.tmp$hessian)))

n.tmp <- optim(b.tmp$par, boundaryconstr.gev.fun,
   #lower = c(-0.5, 0.0000001), upper = c(0,Inf),
   z = z.tmp, control = list(fnscale = -1, maxit = 10000),
   hessian = TRUE, method = "Nelder-Mead")
## Went back to its former values.
b.tmp <- optim(n.tmp$par, boundaryconstr.gev.fun,
   #lower = c(-0.5, 0.0000001), upper = c(0,Inf),
   z = z.tmp, control = list(fnscale = -1, maxit = 10000),
   hessian = TRUE, method = "BFGS")
## Stayed at the Nelder-Mead maximum. Why is it telling that
## In log(locsc) : NaNs produced, is this 


## --------------------------------------------------------------------------------------------------------
## Wholesale application

## Fit the gev distribution to the last 20, 50, 100, 200, 500 maxima at each location
nb.max.for.est <- c(20,30,50,100,200,500)
colnames.for.gevpars <- list(c("xi20","sig20","mu20","xise20","sigse20","muse20"),
   c("xi30","sig30","mu30","xise30","sigse30","muse30"),
   c("xi50","sig50","mu50","xise50","sigse50","muse50"), 
   c("xi100","sig100","mu100","xise100","sigse100","muse100"), 
   c("xi200","sig200","mu200","xise200","sigse200","muse200"),
   c("xi500","sig500","mu500","xise500","sigse500","muse500"))

## For the training set:

dfr.tmp <- data.frame(matrix(ncol = 37, nrow = nrow(pars.train.dfr)))
colnames(dfr.tmp) <- c("name", unlist(colnames.for.gevpars))
dfr.tmp$name <- pars.train.dfr$name

gevfap.train.dfr <- vector(length(nb.max.for.est), mode = "list")
names(gevfap.train.dfr) <- as.character(nb.max.for.est)

for(jj in 1:length(nb.max.for.est))
#system.time(
#for(jj in 1)
 {
  g.tmp <- data.frame(matrix(ncol = ncol(maxima.train.dfr),
    nrow = 1000))
  colnames(g.tmp) <- colnames(maxima.train.dfr)
  for(ii in 1:nrow(pars.train.dfr))
#  for(ii in 5)
   {
   	nn <- as.character(pars.train.dfr$name[ii])
    zass <- maxima.train.dfr[1:1000, nn]
    zest <- maxima.train.dfr[1001:(1000+nb.max.for.est[jj]),nn]
    r.tmp <- try(optim(c(-0.01, sqrt(6*var(zest))/pi), boundaryconstr.gev.fun,
	   z = zest, control = list(fnscale = -1, maxit = 10000),
	   hessian = TRUE, method = "Nelder-Mead"))
	s.tmp <- sqrt(diag(-solve(r.tmp$hessian)))
    dfr.tmp[dfr.tmp$name == nn, colnames.for.gevpars[[jj]]] <-
       if(inherits(r.tmp, "try-error")) {
          rep(NA, 6)
       } else {
          c(r.tmp$par[1], r.tmp$par[2], 1+r.tmp$par[2]/r.tmp$par[1],
            s.tmp, -s.tmp[1]*r.tmp$par[2]/r.tmp$par[1]^2 + s.tmp[2]/r.tmp$par[1])
       }
    g.tmp[, nn] <- if(inherits(r.tmp, "try-error")) {
          rep(NA, length(zass))
       } else {
          1 - sapply(zass, evd::pgev,
          shape = r.tmp$par[1], loc = 1+r.tmp$par[2]/r.tmp$par[1], scale = r.tmp$par[2])
       }
   }
  gevfap.train.dfr[[as.character(nb.max.for.est[jj])]] <- g.tmp
  cat(jj, "is ready \n")
 }
#)
dfr1 <- merge(pars.train.dfr[, c("name","ra","dec","lambda","beta","n","vart","region",
   "q95","q99","fitted.q95","fitted.q99","fr.above.thr95","fr.above.thr99",
   "fr.baluev.sign05","fr.baluev.sign01")], dfr.tmp)
pars.train.dfr <- dfr1

## How many impossible values (above or below an estimated endoint of distribution)?
sum(apply(gevfap.train.dfr[[1]], 2, function(vec) any(is.na(vec))))
## all could be computed
sum(apply(gevfap.train.dfr[[6]], 2, function(vec) any(vec == 1)))
## For 6,5,4,3, no such value (no maximum above an upper boundary),
## for 2, there are 3, for 1, 6. Better than without the constraint. 
sum(apply(gevfap.train.dfr[[6]], 2, function(vec) any(vec == 0)))
## All zero, well, ok. Even with positive fitted xi, it did not find smaller maxima than
## the lower limit.


## For the test set:

dfr.tmp <- data.frame(matrix(ncol = 37, nrow = nrow(pars.test.dfr)))
colnames(dfr.tmp) <- c("name", unlist(colnames.for.gevpars))
dfr.tmp$name <- pars.test.dfr$name

gevfap.test.dfr <- vector(length(nb.max.for.est), mode = "list")
names(gevfap.test.dfr) <- as.character(nb.max.for.est)

for(jj in 1:length(nb.max.for.est))
#system.time(
#for(jj in 1)
 {
  g.tmp <- data.frame(matrix(ncol = ncol(maxima.test.dfr),
    nrow = 1000))
  colnames(g.tmp) <- colnames(maxima.test.dfr)
  for(ii in 1:nrow(pars.test.dfr))
#  for(ii in 5)
   {
   	nn <- as.character(pars.test.dfr$name[ii])
    zass <- maxima.test.dfr[1:1000, nn]
    zest <- maxima.test.dfr[1001:(1000+nb.max.for.est[jj]),nn]
    r.tmp <- try(optim(c(-0.01, sqrt(6*var(zest))/pi), boundaryconstr.gev.fun,
	   z = zest, control = list(fnscale = -1, maxit = 10000),
	   hessian = TRUE, method = "Nelder-Mead"))
    dfr.tmp[dfr.tmp$name == nn, colnames.for.gevpars[[jj]]] <-
       if(inherits(r.tmp, "try-error")) {
          rep(NA, 6)
       } else {
          c(r.tmp$par[1], r.tmp$par[2], 1+r.tmp$par[2]/r.tmp$par[1],
            s.tmp, -s.tmp[1]*r.tmp$par[2]/r.tmp$par[1]^2 + s.tmp[2]/r.tmp$par[1])
       }
    g.tmp[, nn] <- if(inherits(r.tmp, "try-error")) {
          rep(NA, length(zass))
       } else {
          1 - sapply(zass, evd::pgev,
          shape = r.tmp$par[1], loc = 1+r.tmp$par[2]/r.tmp$par[1], scale = r.tmp$par[2])
       }
   }
  gevfap.test.dfr[[as.character(nb.max.for.est[jj])]] <- g.tmp
  cat(jj, "is ready \n")
 }
#)
dfr0 <- merge(pars.test.dfr[, c("name","ra","dec","lambda","beta","n","vart","region",
   "q95","q99","fitted.q95","fitted.q99","fr.above.thr95","fr.above.thr99",
   "fr.baluev.sign05","fr.baluev.sign01")], dfr.tmp)
   
pars.test.dfr <- dfr0

## No warnings at all about NaNs 

## How many failed optimization, not reached likelihood maximum?
sum(apply(dfr.tmp, 1, function(vec) any(is.na(vec))))   ## 6, ok, the non-existing data.

## How many impossible values (above or below an estimated endoint of distribution)?
sum(apply(gevfap.test.dfr[[6]], 2, function(vec) any(is.na(vec))))
## all could be computed, apart from the 6 missing lications.
sum(apply(gevfap.test.dfr[[6]], 2, function(vec) any(vec == 1)), na.rm = TRUE)
## For 6,5,4,3, no such value (no maximum above an upper boundary),
## for 2, there are 6, for 1, 25. Better than without the constraint. 
sum(apply(gevfap.train.dfr[[1]], 2, function(vec) any(vec == 0)), na.rm = TRUE)
## All zero, well, ok. Even with positive fitted xi, it did not find smaller maxima than
## the lower limit. All this looks pretty positive. 


## --------------------------------------------------------------------------------------------------------
## Check the number of significant values in the test set
## --------------------------------------------------------------------------------------------------------

pars.train.dfr$fr.gev20.sign05 <- NA
pars.train.dfr$fr.gev20.sign01 <- NA
pars.test.dfr$fr.gev20.sign05 <- NA
pars.test.dfr$fr.gev20.sign01 <- NA
pars.train.dfr$fr.gev30.sign05 <- NA
pars.train.dfr$fr.gev30.sign01 <- NA
pars.test.dfr$fr.gev30.sign05 <- NA
pars.test.dfr$fr.gev30.sign01 <- NA
pars.train.dfr$fr.gev50.sign05 <- NA
pars.train.dfr$fr.gev50.sign01 <- NA
pars.test.dfr$fr.gev50.sign05 <- NA
pars.test.dfr$fr.gev50.sign01 <- NA
pars.train.dfr$fr.gev100.sign05 <- NA
pars.train.dfr$fr.gev100.sign01 <- NA
pars.test.dfr$fr.gev100.sign05 <- NA
pars.test.dfr$fr.gev100.sign01 <- NA
pars.train.dfr$fr.gev200.sign05 <- NA
pars.train.dfr$fr.gev200.sign01 <- NA
pars.test.dfr$fr.gev200.sign05 <- NA
pars.test.dfr$fr.gev200.sign01 <- NA
pars.train.dfr$fr.gev500.sign05 <- NA
pars.train.dfr$fr.gev500.sign01 <- NA
pars.test.dfr$fr.gev500.sign05 <- NA
pars.test.dfr$fr.gev500.sign01 <- NA


## Compute the success rates on the first 1000 simulations 

for(ii in 1:6)
   {
	n.tmp <- as.character(pars.train.dfr$name[apply(pars.train.dfr[, (16+(ii-1)*6+1):(16+ii*6)],
	    1, function(vec) all(!is.na(vec)))])
	dfr.tmp <- gevfap.train.dfr[[ii]][, n.tmp]
	k.tmp <- nrow(dfr.tmp)
	for(nn in n.tmp)
	   {
	   	pars.train.dfr[pars.train.dfr$name == nn, 52+(ii-1)*2+1] <-
	   	   sum(dfr.tmp[, nn] < 0.05, na.rm = TRUE) / sum(!is.na(dfr.tmp[,nn]))
	   	pars.train.dfr[pars.train.dfr$name == nn, 52+(ii-1)*2+2] <-
	   	   sum(dfr.tmp[, nn] < 0.01, na.rm = TRUE) / sum(!is.na(dfr.tmp[,nn]))
       }
    cat(ii, "is ready \n")
   }

for(ii in 1:6)
   {
	n.tmp <- as.character(pars.test.dfr$name[apply(pars.test.dfr[, (16+(ii-1)*6+1):(16+ii*6)],
	    1, function(vec) all(!is.na(vec)))])
	dfr.tmp <- gevfap.test.dfr[[ii]][, n.tmp]
	k.tmp <- nrow(dfr.tmp)
	for(nn in n.tmp)
	   {
	   	pars.test.dfr[pars.test.dfr$name == nn, 52+(ii-1)*2+1] <-
	   	   sum(dfr.tmp[, nn] < 0.05, na.rm = TRUE) / sum(!is.na(dfr.tmp[,nn]))
	   	pars.test.dfr[pars.test.dfr$name == nn, 52+(ii-1)*2+2] <-
	   	   sum(dfr.tmp[, nn] < 0.01, na.rm = TRUE) / sum(!is.na(dfr.tmp[,nn]))
       }
    cat(ii, "is ready \n")
   }

quartz(height = 7.5, width = 12.5)
par(mfcol = c(3,3), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.6,0))

for(ii in c(55,59,63))
   {
	plot(pars.test.dfr$n[pars.test.dfr$region == "line"], 
	   pars.test.dfr[pars.test.dfr$region == "line", ii], 
	   type = "p", pch = "-", ylim = c(0,0.15),
	   main = "GEV FAP \n Line", ylab = "Fraction above threshold", xlab = "N")
	abline(h = 0.05, col = "violetred")
	points(pars.test.dfr$n[pars.test.dfr$region == "line"],
	   pars.test.dfr[pars.test.dfr$region == "line", ii+1], 
	   type = "p", pch = "-", col = "grey")
	abline(h = 0.01, col = "orange")
	
	plot(pars.test.dfr$n[pars.test.dfr$region == "ecl"], 
	   pars.test.dfr[pars.test.dfr$region == "ecl", ii], 
	   type = "p", pch = "-", ylim = c(0,0.15),
	   main = "Rectangle", ylab = "Fraction above threshold", xlab = "N")
	abline(h = 0.05, col = "violetred")
	points(pars.test.dfr$n[pars.test.dfr$region == "ecl"],
	   pars.test.dfr[pars.test.dfr$region == "ecl", ii+1], 
	   type = "p", pch = "-", col = "grey")
	abline(h = 0.01, col = "orange")
	
	plot(pars.test.dfr$n[pars.test.dfr$region == "randompos"], 
	   pars.test.dfr[pars.test.dfr$region == "randompos", ii], 
	   type = "p", pch = "-", ylim = c(0,0.15),
	   main = "Random positions", ylab = "Fraction above threshold", xlab = "N")
	abline(h = 0.05, col = "violetred")
	points(pars.test.dfr$n[pars.test.dfr$region == "randompos"],
	   pars.test.dfr[pars.test.dfr$region == "randompos", ii+1], 
	   type = "p", pch = "-", col = "grey")
	abline(h = 0.01, col = "orange")
   }


















