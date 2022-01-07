

# WORKSPACE.NAME <- "EMbetaFunctions.RData"
# save.image(WORKSPACE.NAME)


## ------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------
## Incomplete-data model for the beta-distributions
## ------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------

## ------------------------------------------------------------------------------------
## Functions for simple log-likelihood fit for beta parameters
## ------------------------------------------------------------------------------------

## Some simple functions for fitting individual beta parameters, as used in
## firstLookLine2.RData:
logbeta.fun <- function(par, faps)
   {
   	sum(dbeta(faps, shape1 = par[1], shape2 = par[2], log = TRUE))
   }
fitbeta.fun <- function(faps, control = list(fnscale = -1, maxit = 10000), ...)
   {
   	f.tmp <- faps
   	f.tmp[f.tmp == 1] <- 1-1e-10
   	f.tmp[f.tmp == 0] <- 1e-10
   	m.tmp <- mean(faps)
   	v.tmp <- var(faps)
  	par.init <- c(m.tmp, 1-m.tmp) * (m.tmp - m.tmp^2 - v.tmp) / v.tmp
  	optim(par = par.init, logbeta.fun, faps = f.tmp, control = control, hessian = TRUE, ...)
   }

## ------------------------------------------------------------------------------------
## Functions for EM algorithm
## ------------------------------------------------------------------------------------

## beta.logpar.fun: the (one-component) beta distribution, with a possible linear relation
## for both parameters; the number of covariates X[i] in both parameters are the
## same
##  f(p) = B( exp(c[0] + c[1]*X[1] + ... + c[K]*X[K]), exp(d[0] + d[1]*X[1] + ... + d[K]*X[K]) ) *
##         * p ^ exp(c[0] + c[1]*X[1] + ... + c[K]*X[K] - 1) * 
##         * (1-p) ^ exp(d[0] + d[1]*X[1] + ... + d[K]*X[K] - 1)
## Arguments: 
##   pars: the parameters of the linear regression; its length is always twice
##         the length of vars + 2; first half: components of alpha, second half:
##         components of beta; if it is only 2, then these are the constant terms,
##         that is, alpha and beta without dependence on any aliases or variance.
##   fap: the fap.
##   vars: the covariates; if NULL, then the one-component beta distribution does not
##         depend on any variables (constant model fro its parameters), if it is a
##         number or a vector, then covariates in the model.
##   if.log: whether I need its logarithm (for log-likelihood) or the density 
beta.logpar.fun <- function(pars, fap, vars = NULL, if.log = FALSE)
  {
   par.alpha <- pars[1:(length(pars)/2)]
   par.beta <- pars[-(1:(length(pars)/2))]
   covar <- c(1, as.numeric(vars))
   alpha <- exp(sum(par.alpha * covar))
   beta <- exp(sum(par.beta * covar))
   dbeta(fap, shape1 = alpha, shape2 = beta, log = if.log)
  }

## Checking with the uniform:  
# beta.fun(c(1,1), c(0.01), FALSE)
## 1, ok that's the uniform.
# beta.fun(c(1,1,1,0), fap = c(0.9,0.8,0.7,0.8,0.9), vars = 2, FALSE)
# beta.fun(c(1,0.5,1,0), fap = c(0.9,0.8,0.7,0.8,0.9), vars = 2, FALSE)
# beta.logpar.fun(c(0,1,0,0), fap = c(0.9,0.8,0.7,0.8,0.9), vars = 2, FALSE)
# beta.logpar.fun(c(0,0.5,0,0), fap = c(0.9,0.8,0.7,0.8,0.9), vars = 2, FALSE)


## Tsingle.logpar.fun: compute the class (=band) probability vector for one location,
## given some parameter estimates for that location, and the explanatory variables
## (aliases etc.)
## Arguments: 
##   tau: vector the estimated overall class (band) fractions; length is equal to
##        the number of bands (that is, nrow(pars))
##   pars: the matrix of parameters of the linear regressions in L bands;
##         at the moment, only one similar model formula is allowed for all the 
##         bands (no varying covariates from band to band), but with different
##         parameters;
##         nrow(pars) is L; the columns of correspond to the coefficients of the
##         covariates (vars) in the model. Their number is twice the length of vars:
##         if vars is a scalar, then ncol(pars)=2 (alpha and beta), if vars is a
##         vector, ncol(pars) = 2*length(vars); first half: components of alpha,
##         second half: components of beta
##   vars.list: a two-element list; first is the fap vector, all faps at that location;
##         the second is the covariates; if NULL, then the one-component beta
##         distribution does not depend on any variables (constant model fro its
##         parameters), if it is a number or a vector, then one or more covariates
##         in the model.
## Result: a vector of the probabilities for the FAP (given in the first element in
##         vars) to belong to each of the classes (length L)
# pars <- matrix(c(1,0,1,0, 1,0.5,1,0), ncol = 4, byrow = T)
# vars.list <- list(c(0.9,0.8,0.7,0.8,0.9), 2)
# tau <- c(0.5,0.5)
## The faps are close to 1, so the second class (second row of pars, with a larger alpha
## parameter that shifts the mean towards 1) should be more likely than the first row,
## a uniform.
Tsingle.logpar.fun <- function(tau, pars, vars.list)
   {
   	fap <- vars.list[[1]]
   	vars <- vars.list[[2]]
   	## for all FAP values and the covariates, compute all the models in the 
   	## different classes/bands (the likelihood conditional on the class/band):
   	f.tmp <- apply(pars, MAR = 1, FUN = function(vec)
         sapply(fap, beta.logpar.fun, vars = vars, pars = vec, if.log = FALSE))
   	## take the column-wise product of all the individual faps (????):
   	g.tmp <- apply(f.tmp, MAR = 2, prod)
    g.tmp[g.tmp == Inf] <- 1e+304
   	## multiply each with the prior probability of the class:
   	t.tmp <- tau * g.tmp
   	## finally, normalize with the sum to get a real probability:
   	t.tmp / sum(t.tmp)
   }
## 'sapply'ing this to a list of locations gives a matrix with L rows
## and ncol equal to the number of locations.
## Very good.


## tau.update.fun: Updating tau (the class/band probabilities)
## Arguments: 
##   Tmatrix: a matrix of L x R,
##       row index: class/band index (L), column index: locations (for the line,
##       it will be R = 900, for the ecliptic grid, R = 1538) 
tau.update.fun <- function(Tmatrix)
   {
    apply(Tmatrix, M = 1, FUN = mean, na.rm = TRUE)
   }
# tau.update.fun(l.tmp)
## Looks all right.   


## wloglik.beta.logpar.fun: the weighted contribution of one band to the
## expected log-likelihood, conditional on the current estimates of 
## the class variables; that is, the function which should be optimized
## for pars in the M-step. The optimization is done separately for each band. 
## wloglikgrad.beta.logpar.fun: the gradient of wloglik.beta.logpar.fun
## Arguments: 
##   pars: the parameters of the linear regression of the logarithm of the 
##       beta distribution parameters; first half: components of log(alpha), second
##       half: components of log(beta); if it is only 2, then these are the constant
##       terms, that is, log(alpha) and log(beta) without dependence on any aliases or
##       variance. As the optimization is done separately for all, this is only the parameter
##       vector for only one band.
##   vars.list: the covariates, a list of two-element lists. The first level corresponds
##       to the different locations. At each location, a two-element list is given, of
##       which the first element is the vector of FAPs (750 FAPs per location,
##       for the LineGrid or EclGrid). The second is the vector of explanatory
##       variables (alias peak heights) characterizing the location, in the same order
##       as in pars. Its length should be half of the length of pars minus 1.
##   weights: the probabilities for each location to belong a particular
##       class/band. It is a row of the matrix Tmatrix of dimensions
##       L x R, the output of applying Tsingle.logpar.fun, by fixing a particular
##       value for its row index. In the vector, the index indexes the locations
##       (for the line, it will be R=900, for the ecliptic grid, R=1538) 

## pars: a vector, parameters of only one band
## vars.list: a list of two-element lists
## weights: the row of Tmatrix corresponding to the same band as pars. 
wloglik.beta.logpar.fun2 <- function(pars, vars.list, weights)
   {
#   	pars <- pars.tmp2[1,]
#   	weights <- t.tmp[1,]
#   	vars.list <- vars.list.tmp
   	## Attach the  weights to all locations: this is a vector of length L,
   	## the weight with which the location contributes to the total
   	## log-likelihood of each of the bands
   	vars.modlist.tmp <- list()
   	for(ii in 1:length(vars.list))
   	   vars.modlist.tmp[[ii]] <- list(fap = vars.list[[ii]]$fap, 
   	   covar = vars.list[[ii]]$covar, weights = weights[ii])
   	r.tmp <- sapply(vars.modlist.tmp, function(lst, pars)
   	   {
#   	   	lst <- vars.modlist.tmp[[3]]
     	f.tmp <- sapply(lst[[1]], beta.logpar.fun, vars = lst[[2]], pars = pars, if.log = TRUE)
        lst[[3]]*sum(f.tmp)
   	   }, pars = pars)
   	sum(r.tmp)
   }
wloglik.beta.logpar.fun1 <- function(pars, vars.list, weights)
   {
#   	pars <- pars.tmp2[1,]
#   	weights <- t.tmp[1,]
#   	vars.list <- vars.list.tmp
   	r.tmp <- sapply(vars.list, function(lst, pars)
   	   {
#   	   	lst <- vars.modlist.tmp[[3]]
     	f.tmp <- sapply(lst[[1]], beta.logpar.fun, vars = lst[[2]], pars = pars, if.log = TRUE)
        sum(f.tmp)
   	   }, pars = pars)
   	sum(r.tmp * weights)
   }
wloglik.beta.logpar.fun <- function(pars, vars.list, weights)
   {
#   	pars <- pars.tmp2[1,]
#   	weights <- t.tmp[1,]
#   	vars.list <- vars.list.tmp
    par.alpha <- pars[1:(length(pars)/2)]
    par.beta <- pars[-(1:(length(pars)/2))]
   	r.tmp <- sapply(vars.list, function(lst, pars)
   	   {
#   	   	lst <- vars.modlist.tmp[[3]]
        covar <- c(1, as.numeric(lst[[2]]))
        alpha <- exp(sum(par.alpha * covar))
        beta <- exp(sum(par.beta * covar))
        n <- length(lst[[1]])
#     	f.tmp <- dbeta(lst[[1]], shape1 = alpha, shape2 = beta, log = TRUE)
        f.tmp <- - n * lbeta(alpha, beta) + 
                 (alpha - 1) * sum(log(lst[[1]])) +  
                 (beta - 1) * sum(log(1-lst[[1]]))
        sum(f.tmp)
   	   }, pars = pars)
   	sum(r.tmp * weights)
   }
system.time(wloglik.beta.logpar.fun2(pars, vars.list, weights))
system.time(wloglik.beta.logpar.fun1(pars, vars.list, weights))
system.time(wloglik.beta.logpar.fun(pars, vars.list, weights))


wloglikgrad.logpar.fun <- function(pars, vars.list, weights)
   {
    par.alpha <- pars[1:(length(pars)/2)]
    par.beta <- pars[-(1:(length(pars)/2))]
    g.tmp <- sapply(vars.list, function(lst)
       {
        covar <- c(1,as.numeric(lst[[2]]))
        alpha <- exp(sum(par.alpha * covar))
        beta <- exp(sum(par.beta * covar))
        n <- length(lst[[1]])
        sumlogfap <- sum(log(lst[[1]]))
        sumlogcompfap <- sum(log(1-lst[[1]]))
        dga <- n * (digamma(alpha+beta) - digamma(alpha)) + sumlogfap
        dgb <- n * (digamma(alpha+beta) - digamma(beta)) + sumlogcompfap
        # c(dga * alpha,
          # dga * alpha * lst[[2]],
          # dgb * beta,
          # dgb * beta * lst[[2]])
        c(dga * alpha * covar, dgb * beta * covar)
        })
   	matrix(weights, nrow = 1) %*% t(g.tmp)
   }
# wloglikgrad.logpar.fun(pars = res0$par[1,], vars = vars.list.tmp, 
   # weights = res0$location.band.prob[1,])
# wloglikgrad.logpar.fun(pars = res0$par[2,], vars = vars.list.tmp, 
   # weights = res0$location.band.prob[2,])
# wloglikgrad.logpar.fun(pars, vars.list, weights)


## EM.logpar.fun: the EM algorithm for the band-discriminative regression estimation
## of the beta parameters.
## Arguments: 
##   pars.init: the initial value for the matrix of parameters of the linear
##       regression in L bands; at the moment, only one single model formula 
##       is allowed for all the bands (no varying covariates from band to band),
##       but with different parameters.
##       nrow(pars) = L; the columns correspond to the coefficients of the
##       covariates (vars) in the model. Their number is twice the length of vars:
##       if vars is a scalar, then ncol(pars)=2 (alpha and beta), if vars is a
##       vector, ncol(pars) = 2*length(vars); first half: components of alpha
##       = alpha0 + alpha1*covar1 + ..., second half: components of beta = beta0 +
##       beta1*covar1 + ...
##   tau.init: initial guess for the vector of the estimated overall class (band)
##       fractions; length is equal to the number of bands (that is, nrow(pars))
##   vars.list: the covariates, a list of two-element lists. The first level corresponds
##       to the different locations. At each location, a two-element list is given, of
##       which the first element is the vector of FAPs (750 FAPs per location,
##       for the LineGrid or EclGrid). The second is the vector of explanatory
##       variables (alias peak heights) characterizing the location, in the same order
##       as in pars. Its length should be (half of the length of pars) minus 1.
##   convergence.crit: relative stopping criterion for the log-likelihood maximization.
##   control.optim: the list of desired control parameters for optim.
##   hessian: Should it compute the Hessian?
##   ...: further parameters for optim (method, box constraints, etc.) For method, 
##       default is the default of optim, that is, Nelder-Mead.
## Use of this function:
## 1. predict the class of every location on the line/ecliptic/all over, to
##    have an initial guess for class proportion (==> tau.init);
## 2. (a) fit individual beta distributions at each location;
##    (b) fit classwise linear regressions as in firstLookLine2.RData to the beta 
##    distribution parameters of the locations in each class, to have an initial
##    parameter guess in each class (==> pars.init)
## 3. Run EM.beta.fun  

EM.logpar.fun <- function(pars.init,
                          tau.init, 
                          vars.list, 
                          convergence.crit = 1e-6, 
                          tmp.filename = "partialRes_EMbetaFunctions.RObjTmp",
                          control.optim = list(maxit = 10000, fnscale = -1),
                          method = "Nelder-Mead", 
                          hessian = TRUE, ...)
   {
    ## Initialization:
    pars <- pars.init
   	loglik.value <- -1e9
   	conv.crit <- convergence.crit + 1
   	## The faps produced by the pipeline contain often 1s or 0s, due to rounding
   	## of numbers within machine-precision from 1 or 0 (10^(-70)...). These are
   	## intractable for dbeta, they must be changed for a reasonable number.
   	vars.list1 <- lapply(vars.list, function(lst)
   	   {
	   	f.tmp <- lst[[1]]
	   	f.tmp[f.tmp == 1] <- 1-1e-11
	   	f.tmp[f.tmp == 0] <- 1e-11
   	   	list(fap = f.tmp, covar = lst[[2]])
   	   }
   	   )
	## E-step: compute the individual class probabilities for all locations 
	Tmatrix <- sapply(vars.list1, function(lst)
	   Tsingle.logpar.fun(vars.list = lst, tau = tau.init, pars = pars))
    tau <- tau.update.fun(Tmatrix)
   	jj <- 1
   	loglik.store <- loglik.value
    models.store <- list()
    models.store[[jj]] <- NULL
    Tmatrix.store <- list()
    Tmatrix.store[[jj]] <- Tmatrix
    tau.store <- list()
    tau.store[[jj]] <- tau
    cat(jj, ": initialization ready \n")
   	while(conv.crit > convergence.crit)
   	   {
	    jj <- jj + 1
	   	loglik.value0 <- loglik.value
	   	## to be able to optimize Q in the M-step:
	   	## M-step: using Tmatrix as weights, optimize the classwise likelihood
	   	r.tmp <- list()
	   	for(ii in 1:nrow(Tmatrix))
	   		{
		   	 pars.tmp <- pars[ii,]
		   	 T.tmp <- Tmatrix[ii,]
	   		 r.tmp[[ii]] <- optim(pars.tmp, fn = wloglik.beta.logpar.fun,
	   		   vars.list = vars.list1, weights = T.tmp, control = control.optim, 
	   		   hessian = hessian, method = method, ...)
	   		}
	    ## Update the matrix of pars (L x added dim of the two linear models):
	    pars <- t(sapply(r.tmp, function(lst) lst$par))
	    ## Update the probabilities of each location to belong to each of the bands
	    ## (E-step):  
	   	Tmatrix <- sapply(vars.list1, function(lst) 
	   	    Tsingle.logpar.fun(vars.list = lst, tau = tau, pars = pars))
        ## Update the band probabilities:
	   	tau <- tau.update.fun(Tmatrix)
	   	## Extract the new log-likelihood maximum: 
	    loglik.value <- sum(sapply(r.tmp, function(lst) lst$value))
	    ## Compute the convergence criterion:
	    conv.crit <- abs((loglik.value - loglik.value0) / loglik.value)
	    loglik.store[jj] <- loglik.value
	    models.store[[jj]] <- r.tmp
	    Tmatrix.store[[jj]] <- Tmatrix
	    tau.store[[jj]] <- tau
	    save(loglik.store, Tmatrix.store, tau.store, models.store,
	       file = tmp.filename)
	    cat(jj, " ready; conv.crit = ", conv.crit, "\n")
   	   }
 	v.tmp <- list(parameters = pars, band.frequency = tau,
 	    location.band.prob = Tmatrix, optimization = r.tmp, 
 	    history = list(logliks = loglik.store, models = models.store,
 	    Tmatrices = Tmatrix.store, taus = tau.store))
# 	save.image(WORKSPACE.NAME)
 	return(v.tmp)
 	rm(r.tmp, v.tmp, Tmatrix, models.store, Tmatrix.store, tau.store, loglik.store)
 	gc()
   }











###############################################################################################
## Tests of the procedures
## Toy examples
############################################################################################### 



# # ## -------------------------------------------------------------------------------------
# ## A toy example with well-separated parameters:

# ## Covariates: let them be random uniform variables
# y1 <- runif(25)
# y2 <- runif(25)

# truepars.tmp <- matrix(c(0,-0.6,0,-0.2, 0,1,-0.5,0.4), ncol = 4, byrow = T)
# a1 <- exp(0 - 0.6 * y1)
# b1 <- exp(0 - 0.2 * y1)
# a2 <- exp(0 + 1 * y2)
# b2 <- exp(-0.5 + 0.4 * y2)
# # quartz()
# # par(mfrow = c(1,2), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.7,0))
# # plot(y1, a1, pch = 16, cex = 0.5, ylim = range(c(a1,a2)))
# # points(y2, a2, pch = 16, cex = 0.5, col = 2)
# # plot(y1, b1, pch = 16, cex = 0.5, ylim = range(c(b1,b2)))
# # points(y2, b2, pch = 16, cex = 0.5, col = 2)

# sw.tmp <- data.frame(band = rep(c(1,2), each = 25), y = c(y1,y2), a = c(a1,a2), b = c(b1,b2))

# vars.list.tmp <- apply(sw.tmp[,c("a","b","y")], MAR = 1, FUN = function(vec)
   # {
   	# list(fap = rbeta(400, shape1 = vec["a"], shape2 = vec["b"]),
   	     # covar = vec["y"])
   # })

# sw.tmp$supposed.band <- sw.tmp$band 
# misclass1 <- sample(1:25, size = 7)
# misclass2 <- sample(26:50, size = 12)
# sw.tmp$supposed.band2 <- sw.tmp$band 
# sw.tmp$supposed.band2[misclass1] <- 2 
# sw.tmp$supposed.band2[misclass2] <- 1 

# tau.init.tmp <- table(sw.tmp[,c("supposed.band")]) / 50

# betapar.tmp <- data.frame(t(sapply(vars.list.tmp, function(lst) 
   # {
# #   	lst <- vars.list.tmp[[1]]
   	# r.tmp <- try(fitbeta.fun(lst[[1]]))
   	# r.tmp$par
   	# })))
# colnames(betapar.tmp) <- c("a.est","b.est")
# sw.tmp$a.est <- betapar.tmp$a.est
# sw.tmp$b.est <- betapar.tmp$b.est
# # sw.tmp[1:10,]
# # sw.tmp[,c("a","a.est","b","b.est")]
# # quartz(height = 3.5, width = 7)
# # par(mfrow = c(1,2), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.7,0))
# # plot(sw.tmp[,c("a")], sw.tmp[,c("a.est")], pch = 16, cex = 0.7, ylim = c(0,6),
   # # col = c(rep("blue", 25), rep("deepskyblue", 25)))
# # abline(c(0,1))
# # plot(sw.tmp[,c("b")], sw.tmp[,c("b.est")], pch = 16, cex = 0.7, ylim = c(0,10),
   # # col = c(rep("blue", 25), rep("deepskyblue", 25)))
# # abline(c(0,1))

# pars.init.tmp <- matrix(c(lm(log(a.est) ~ y, data = sw.tmp, subset = supposed.band == 1)$coef,
    # lm(log(b.est) ~ y, data = sw.tmp, subset = supposed.band == 1)$coef,
    # lm(log(a.est) ~ y, data = sw.tmp, subset = supposed.band == 2)$coef,
    # lm(log(b.est) ~ y, data = sw.tmp, subset = supposed.band == 2)$coef), 
    # ncol = 4, byrow = TRUE)
# # pars.init.tmp


# system.time(
# res0 <- EM.logpar.fun(pars.init = pars.init.tmp, tau.init = tau.init.tmp,
   # vars.list = vars.list.tmp, method = "Nelder-Mead")
# )
# ## Time: 1200s
# system.time(
# res1 <- EM.logpar.fun(pars.init = pars.init.tmp, tau.init = tau.init.tmp,
   # vars.list = vars.list.tmp, method = "BFGS")
# )
# ## Time: 840s
# system.time(
# res2 <- EM.logpar.fun(pars.init = pars.init.tmp, tau.init = tau.init.tmp,
   # vars.list = vars.list.tmp, method = "BFGS", gr = wloglikgrad.logpar.fun)
# )
# ## Time: 198s (!)
# ## For the latter two, NaNs are produced, but where the hell can they be produced?
# ## The results are the same for the first two, so it seems either the NM produces
# ## them too but does not shout, or they do not count materially. In any case,
# ## the effect cannot be big, as the estimates are very close to the true values.
# 2*sqrt(diag(-solve(res0$optimization[[1]]$hessian)))
# 2*sqrt(diag(-solve(res0$optimization[[2]]$hessian)))
# res0$parameters
# truepars.tmp
# 2*sqrt(diag(-solve(res1$optimization[[1]]$hessian)))
# 2*sqrt(diag(-solve(res1$optimization[[2]]$hessian)))
# res1$parameters
# truepars.tmp
# 2*sqrt(diag(-solve(res2$optimization[[1]]$hessian)))
# 2*sqrt(diag(-solve(res2$optimization[[2]]$hessian)))
# res2$parameters
# truepars.tmp
# ## And all of them contain in the confidence interval the true values.
# wloglikgrad.logpar.fun(pars = res2$par[1,], vars = vars.list.tmp, 
   # weights = res2$location.band.prob[1,])
# wloglikgrad.logpar.fun(pars = res2$par[2,], vars = vars.list.tmp, 
   # weights = res2$location.band.prob[1,])
# ## not exactly zero, but close to.

# ## ------------------------------------------------------------------------------------
# ## Test case 1, more similar to the FAPs
# ## ------------------------------------------------------------------------------------

# ## 1.step:
# ## Two classes of locations: 60 locations of class 1, 140 of class 2.

# ## 2.step:
# ## One predictor; based on something like a partial spectral window as I use it,
# ## but randomly generated, not from the sky
# ##   0. to define the two classes, take twice 16 pairs of parameters defining twice
# ##      16 beta distributions in order to generate, 
# ##      from which I will 
# ##   1. generate 60 times the 16-tuples Y[i] from the first class (i = 1,...,60),
# ##      and 120 times 16-tuples Y[i] (j = 61,...,180); 
# ##      these are the simulated "spectral window peaks" at the 60 class-1 locations,
# ##      and at the 120 class-2 locations
# ##   2. take the sum of the six highest to get the covariate X[i]

# ## 3.step:
# ## For location i in class 1, take log(a[i]) = c1 + d1*X[i]
# ##                                 log(b[i]) = e1 + f1*X[i]
# ## as parameters of the beta distribution at location i (i = 1,...,60).
# ## For location i in class 2, take log(a[i]) = c2 + d2*X[i]
# ##                                 log(b[i]) = e2 + f2*X[i]
# ## as parameters of the beta distribution at location i (i = 61,...,180).

# ## 4.step:
# ## generate 100 values from Beta(a[i], b[i]) ("faps") for each location (i = 1,...,180).
# ## Make a two-element list with the 100 faps as the first component and X[i] as the second
# ## for each location i = 1,...,180 (vars.list).

# ## 5.step:
# ## Mis-label 20 randomly chosen locations of class 1 and 40 of class 2 as the other
# ## class. These will be the supposed class labels. Compute the supposed class
# ## proportions from the wrong labels (tau.init).

# ## 6.step:
# ## Estimate the beta parameters for all locations, regardless to their true or 
# ## supposed class. Fit lm to the collection of beta parameters in both supposed
# ## classes with X[i] as a covariate (pars.init).

# ## 7.step:
# ## Try EM.beta.fun on these data.


# ## -------------------------------------------------------------------------------------
# ## Well-separated parameters and good initial parameters, due to 
# ## no mislabeling of locations:


# ## 1.step:
# ## ok.

# ## 2.step (generate predictors, similar to our "real" case):
# m.testcase11 <-    c(0.6, 0.7, 0.75,0.7, 0.4, 0.25,0.2, 0.7, 0.7, 0.3, 0.25,0.3, 0.75,0.2, 0.15,0.7)
# m.testcase12 <-    c(0.3, 0.25,0.4, 0.6, 0.2, 0.8, 0.85,0.05,0.75,0.3, 0.65,0.2, 0.7, 0.2, 0.25,0.55)

# al.testcase11 <- 7*c(1.1, 0.1, 0.2, 0.2, 1.2, 0.5, 0.1, 1.1, 1.5, 1.4, 0.4, 0.4, 0.1, 1.1, 1.1, 1)
# al.testcase12 <- 6*c(0.2, 0.1, 0.5, 1.5, 1.2, 0.5, 1.1, 1.1, 0.7, 0.4, 0.4, 1.4, 0.7, 0.2, 1.2, 1)

# be.testcase11 <- (1/m.tmp1 - 1) * al.tmp1
# be.testcase12 <- (1/m.tmp2 - 1) * al.tmp2

# ## Check:
# m.testcase11 - al.testcase11/(al.testcase11 + be.testcase11)
# m.testcase12 - al.testcase12/(al.testcase12 + be.testcase12)
# ## ok within machine precision

# ## Generate random "spectral windows":
# sw.testcase11 <- apply(cbind(al.testcase11, be.testcase11), M = 1, F = function(vec)
   # rbeta(60, shape1 = vec[1], shape2 = vec[2]))
# sw.testcase12 <- apply(cbind(al.testcase12, be.testcase12), M = 1, F = function(vec)
   # rbeta(120, shape1 = vec[1], shape2 = vec[2]))
# sw.testcase1 <- data.frame(rbind(sw.testcase11, sw.testcase12))
# colnames(sw.testcase1) <- c("at4","at8","at12","at16","at24","at28","at32","at36","at40","at44",
   # "at48","at52","at56","at60","at64","at68")
# sw.testcase1$band <- c(rep(1, 60), rep(2,120)) 

# ## Generate random "spectral windows" - check the aspect:
# quartz()
# par(mfrow = c(1,1), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.7,0))
# plot(1:16, sw.tmp1[1,], type = "l", ylim = c(0,1))
# for(ii in 2:60)
   # lines(1:16, sw.tmp1[ii,])
# for(ii in 1:120)
   # lines(1:16, sw.tmp2[ii,], col = "red")
# ## Say that this is ok.

# ## Compute the sum of the six highest:
# sw.testcase1$hsix <- apply(sw.testcase1[1:16], M = 1, F = function(vec)
   # {
   	# v.tmp <- sort(vec, decreasing = TRUE)
   	# sum(v.tmp[1:6])
   # })
# # sw.testcase1[,c("hsix","band")]

# ## See their distribution:
# # quartz()
# # boxplot(sw.testcase1$hsix ~ sw.testcase1$band)
# ## Say that this is ok.

# ## 3.step (linear relationships):
# ## linear parameters:
# truepars.testcase1 <- matrix(c(-2,0.6,2.7,-0.5, 0.2,-0.2,-2,0.6), ncol = 4, byrow = TRUE)

# sw.testcase1$a[sw.testcase1$band == 1] <- sapply(sw.testcase1$hsix[sw.testcase1$band == 1],
    # function(sca) {exp(truepars.testcase1[1,2] * sca + truepars.testcase1[1,1])})
# sw.testcase1$b[sw.testcase1$band == 1] <- sapply(sw.testcase1$hsix[sw.testcase1$band == 1],
    # function(sca) {exp(truepars.testcase1[1,4] * sca + truepars.testcase1[1,3])})
# sw.testcase1$a[sw.testcase1$band == 2] <- sapply(sw.testcase1$hsix[sw.testcase1$band == 2],
    # function(sca) {exp(truepars.testcase1[2,2] * sca + truepars.testcase1[2,1])})
# sw.testcase1$b[sw.testcase1$band == 2] <- sapply(sw.testcase1$hsix[sw.testcase1$band == 2],
    # function(sca) {exp(truepars.testcase1[2,4] * sca + truepars.testcase1[2,3])})
# ## Give names to the locations:
# sw.testcase1$loc.name <- as.character(1:180)
# sw.testcase1[,c("a","b")]

# # quartz()
# # par(mfrow = c(1,2), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.7,0))
# # plot(sw.testcase1$hsix[sw.testcase1$band == 1], 
   # # sw.testcase1$a[sw.testcase1$band == 1], pch = 16, 
   # # cex = 0.5, ylim = range(sw.testcase1$a))
# # points(sw.testcase1$hsix[sw.testcase1$band == 2], 
   # # sw.testcase1$a[sw.testcase1$band == 2], pch = 16, 
   # # cex = 0.5, col = 2)
# # plot(sw.testcase1$hsix[sw.testcase1$band == 1], 
   # # sw.testcase1$b[sw.testcase1$band == 1], pch = 16, 
   # # cex = 0.5, ylim = range(sw.testcase1$b))
# # points(sw.testcase1$hsix[sw.testcase1$band == 2], 
   # # sw.testcase1$b[sw.testcase1$band == 2], pch = 16, 
   # # cex = 0.5, col = 2)

# ## 4.step (fap generation and var.list):
# vars.list.tmp <- apply(sw.testcase1[,c("a","b","hsix")], MAR = 1, FUN = function(vec)
   # {
   	# list(fap = rbeta(200, shape1 = vec["a"], shape2 = vec["b"]),
   	     # covar = vec[3])
   # })
# names(vars.list.tmp) <- sw.testcase1$loc.name

# ## 5.step (should create mis-labeling, but not in this first run):
# sw.testcase1$supposed.band <- sw.testcase1$band 
# # sw.testcase1$supposed.band[sample(1:60, size = 10)] <- 2
# # sw.testcase1$supposed.band[sample(61:180, size = 20)] <- 1
# # sw.testcase1[,c("band","supposed.band")]
# # table(sw.testcase1[,c("band","supposed.band")])
# # table(sw.testcase1[,c("supposed.band")]) / 180
# # table(sw.testcase1[,c("band")]) / 180
# ## All right.
# tau.init.tmp <- table(sw.testcase1[,c("supposed.band")]) / 180

# ## 6.step (individual beta parameters and lm under mis-labeling):
# ## first version: likelihood estimation
# betapar.tmp <- data.frame(t(sapply(vars.list.tmp, function(lst) 
   # {
# #   	lst <- vars.list.tmp[[1]]
   	# r.tmp <- try(fitbeta.fun(lst[[1]]))
   	# r.tmp$par
   	# })))
# colnames(betapar.tmp) <- c("a.est","b.est")
# sw.testcase1$a.est <- betapar.tmp$a.est
# sw.testcase1$b.est <- betapar.tmp$b.est
# # sw.testcase1[1:10,]
# sw.testcase1[,c("a","a.est","b","b.est")]
# # quartz(height = 3.5, width = 7)
# # par(mfrow = c(1,2), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.7,0))
# # plot(sw.testcase1[,c("a")], sw.testcase1[,c("a.est")], pch = 16, cex = 0.7, ylim = c(0,6),
   # # col = c(rep("blue", 60), rep("deepskyblue", 120)))
# # abline(c(0,1), col = 2)
# # plot(sw.testcase1[,c("b")], sw.testcase1[,c("b.est")], pch = 16, cex = 0.7, ylim = c(0,10),
   # # col = c(rep("blue", 60), rep("deepskyblue", 120)))
# # abline(c(0,1), col = 2)

# ## second version: moment estimation
# # betapar.tmp2 <- data.frame(t(sapply(vars.list.tmp, function(lst) 
   # # {
   	# # m.tmp <- mean(lst[[1]])
   	# # v.tmp <- var(lst[[1]])
  	# # c(m.tmp, 1-m.tmp) * (m.tmp - m.tmp^2 - v.tmp) / v.tmp
   	# # })))
# # colnames(betapar.tmp2) <- c("a.mom","b.mom")
# # sw.testcase1$a.mom <- betapar.tmp2$a.mom
# # sw.testcase1$b.mom <- betapar.tmp2$b.mom
# # # sw.testcase1[1:10,]
# # # sw.testcase1[,c("a","a.est","a.mom")]
# # # sw.testcase1[,c("b","b.est","b.mom")]
# # quartz(height = 3.5, width = 7)
# # par(mfrow = c(1,2), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.7,0))
# # plot(sw.testcase1[,c("a")], sw.testcase1[,c("a.mom")], pch = 16, cex = 0.7, ylim = c(0,4.5),
   # # col = c(rep("blue", 60), rep("deepskyblue", 120)))
# # abline(c(0,1), col = 2)
# # plot(sw.testcase1[,c("b")], sw.testcase1[,c("b.mom")], pch = 16, cex = 0.7, ylim = c(0,10),
   # # col = c(rep("blue", 60), rep("deepskyblue", 120)))
# # abline(c(0,1), col = 2)
# # ## The .est seems to be better iintial value.

# pars.init.tmp <- matrix(c(lm(log(a.est) ~ hsix, data = sw.testcase1, subset = supposed.band == 1)$coef,
    # lm(log(b.est) ~ hsix, data = sw.testcase1, subset = supposed.band == 1)$coef,
    # lm(log(a.est) ~ hsix, data = sw.testcase1, subset = supposed.band == 2)$coef,
    # lm(log(b.est) ~ hsix, data = sw.testcase1, subset = supposed.band == 2)$coef), 
    # ncol = 4, byrow = TRUE)
# pars.init.tmp

# ## Go ahead, run it:

# system.time(
# res.nm.testcase1 <- EM.logpar.fun(pars.init = pars.init.tmp, tau.init = tau.init.tmp,
   # vars.list = vars.list.tmp, method = "Nelder-Mead")
# )
# ## Time: s
# system.time(
# res.bfgs.testcase1 <- EM.logpar.fun(pars.init = pars.init.tmp, tau.init = tau.init.tmp,
   # vars.list = vars.list.tmp, method = "BFGS")
# )
# ## Time: s
# system.time(
# res.bfgsgrad.testcase1 <- EM.logpar.fun(pars.init = pars.init.tmp, tau.init = tau.init.tmp,
   # vars.list = vars.list.tmp, method = "BFGS", gr = wloglikgrad.logpar.fun)
# )
# ## Time: 
# 2*sqrt(diag(-solve(res.nm.testcase1 $optimization[[1]]$hessian)))
# 2*sqrt(diag(-solve(res.nm.testcase1 $optimization[[2]]$hessian)))
# res.nm.testcase1 $parameters
# truepars.testcase1
# 2*sqrt(diag(-solve(res.bfgs.testcase1 $optimization[[1]]$hessian)))
# 2*sqrt(diag(-solve(res.bfgs.testcase1 $optimization[[2]]$hessian)))
# res.bfgs.testcase1 $parameters
# truepars.testcase1
# 2*sqrt(diag(-solve(res.bfgsgrad.testcase1 $optimization[[1]]$hessian)))
# 2*sqrt(diag(-solve(res.bfgsgrad.testcase1 $optimization[[2]]$hessian)))
# res.bfgsgrad.testcase1 $parameters
# truepars.testcase1
# ## And all of them contain in the confidence interval the true values.
# wloglikgrad.logpar.fun(pars = res.bfgs.testcase1 $par[1,], vars = vars.list.tmp, 
   # weights = res.bfgs.testcase1 $location.band.prob[1,])
# wloglikgrad.logpar.fun(pars = res.bfgs.testcase1 $par[2,], vars = vars.list.tmp, 
   # weights = res.bfgs.testcase1 $location.band.prob[2,])
# ## The best is BFGS using the exact gradient: the gradients are 1e-4.
# ## BFGS using numeric derivatives is the second (1e-2), NM is worst (1e0 - 1e1).


# ## -------------------------------------------------------------------------------------
# ## The same case with well-separated parameters, but with worse initial parameters 
# ## due to mild mislabeling of locations:

# ## All the former parameters and simulated values are used from the first run
# ## until step 4. Then:

# misclass1.testcase1 <- lapply(c(4,8,12,16,20,25,30,35,40), sample, x = 1:60)
# misclass2.testcase1 <- lapply(c(8,16,24,32,40,50,60,70,80), sample, x = 61:180)
# misclass1.testcase1[[9]] <- sample(x = 1:60, size = 40)
# misclass1.testcase1[[9]] <- sample(x = 61:180, size = 80)
# misclass.em.testcase1 <- list()
# times.misclass.testcase1 <- numeric(9)
# for(kk in 9)
  # {
   # ## 5.step (create now mis-labeling):
   # sw.testcase1$supposed.band2 <- sw.testcase1$band 
   # sw.testcase1$supposed.band2[misclass1.testcase1[[kk]]] <- 2
   # sw.testcase1$supposed.band2[misclass2.testcase1[[kk]]] <- 1
   # tau.init.tmp <- table(sw.testcase1[,c("supposed.band2")]) / 180

   # ## 6.step (individual beta parameters and lm under mis-labeling):
   # ## The same as before, as all the simulated values are the same, apart from the linear fit.
   # pars.init.tmp <- matrix(c(
     # lm(log(a.est) ~ hsix, data = sw.testcase1, subset = supposed.band2 == 1)$coef,
     # lm(log(b.est) ~ hsix, data = sw.testcase1, subset = supposed.band2 == 1)$coef,
     # lm(log(a.est) ~ hsix, data = sw.testcase1, subset = supposed.band2 == 2)$coef,
     # lm(log(b.est) ~ hsix, data = sw.testcase1, subset = supposed.band2 == 2)$coef), 
     # ncol = 4, byrow = TRUE)
   # # pars.init.tmp
   # times.misclass.testcase1[kk] <- system.time(
      # misclass.em.testcase1[[kk]] <- EM.logpar.fun(pars.init.tmp, tau.init.tmp, vars.list.tmp, 
      # convergence.crit = 1e-6, control.optim = list(maxit = 10000, fnscale = -1), 
      # method = "BFGS", gr = wloglikgrad.logpar.fun))[3]
   # save.image(WORKSPACE.NAME)
  # }
# lapply(misclass.em.testcase1, function(lst)
    # round(lst$location.band.prob, 3))
# ## no error, perfect classes
# lapply(misclass.em.testcase1, function(lst)
    # round(lst$par, 4))
# truepars.testcase1
# lapply(misclass.em.testcase1, function(lst)
    # lst$history$taus)
# ## All find the same estimates as with the initially perfect classification, even

# ## -------------------------------------------------------------------------------------
# ## A case with parameters more similar to the real problem:

# ## 1.step:
# ## ok.

# ## 2.step (generate predictors, similar to our "real" case):
# m.tmp1 <-    c(0.5, 0.2, 0.85,0.75, 0.4, 0.8,0.2, 0.1, 0.5, 0.3, 0.25,0.3, 0.75,0.2, 0.15,0.85)
# m.tmp2 <-    c(0.4, 0.25,0.85, 0.4, 0.4, 0.8,0.15,0.05,0.7, 0.3, 0.15,0.2, 0.5, 0.2, 0.25,0.8)

# al.tmp1 <- 5*c(1.5, 1.1, 2,   1.7, 1.2, 2.2, 1.1, 1.1, 1.5, 1.4, 1.4, 1.4, 2,   1.1, 1.1, 1.9)
# al.tmp2 <- 4*c(1.2, 1.1, 2.5, 1.5, 1.2, 2,   1.1, 1.1, 1.3, 1.4, 1.4, 1.4, 1.7, 1.2, 1.2, 1.7)

# be.tmp1 <- (1/m.tmp1 - 1) * al.tmp1
# be.tmp2 <- (1/m.tmp2 - 1) * al.tmp2

# ## Check:
# m.tmp1 - al.tmp1/(al.tmp1 + be.tmp1)
# m.tmp2 - al.tmp2/(al.tmp2 + be.tmp2)
# ## ok within machine precision

# ## Generate random "spectral windows":
# sw.tmp1 <- apply(cbind(al.tmp1, be.tmp1), M = 1, F = function(vec)
   # rbeta(60, shape1 = vec[1], shape2 = vec[2]))
# sw.tmp2 <- apply(cbind(al.tmp2, be.tmp2), M = 1, F = function(vec)
   # rbeta(120, shape1 = vec[1], shape2 = vec[2]))
# sw.testcase2 <- data.frame(rbind(sw.tmp1, sw.tmp2))
# colnames(sw.testcase2) <- c("at4","at8","at12","at16","at24","at28","at32","at36","at40","at44",
   # "at48","at52","at56","at60","at64","at68")
# sw.testcase2$band <- c(rep(1, 60), rep(2,120)) 

# ## Generate random "spectral windows" - check the aspect:
# # quartz()
# # par(mfrow = c(1,1), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.7,0))
# # plot(1:16, sw.tmp1[1,], type = "l", ylim = c(0,1))
# # for(ii in 2:60)
   # # lines(1:16, sw.tmp1[ii,])
# # for(ii in 1:120)
   # # lines(1:16, sw.tmp2[ii,], col = "red")
# ## Say that this is ok.

# ## Compute the sum of the six highest:
# sw.testcase2$hsix <- apply(sw.testcase2[1:16], M = 1, F = function(vec)
   # {
   	# v.tmp <- sort(vec, decreasing = TRUE)
   	# sum(v.tmp[1:6])
   # })

# ## See their distribution, in reality, the covariate was quite overlapping:
# quartz()
# boxplot(sw.testcase2$hsix ~ sw.testcase2$band)
# ## Say that this is ok.

# ## 3.step (linear relationships):
# ## linear parameters by comparison with firstLookLine1.RData, betafits.lm3[[6]]:
# icepta.tmp1 <- -0.14
# slopea.tmp1 <- 0.0
# icepta.tmp2 <- 0.1
# slopea.tmp2 <- -0.05
# iceptb.tmp1 <- 0.87
# slopeb.tmp1 <- -0.2
# iceptb.tmp2 <- 2.3
# slopeb.tmp2 <- -0.5
# truepars.testcase2 <- matrix(c(icepta.tmp1, slopea.tmp1, iceptb.tmp1, slopeb.tmp1,
   # icepta.tmp2, slopea.tmp2, iceptb.tmp2, slopeb.tmp2), ncol = 4, byrow = T)

# sw.testcase2$a[sw.testcase2$band == 1] <- sapply(sw.testcase2$hsix[sw.testcase2$band == 1],
    # function(sca) exp(slopea.tmp1 * sca + icepta.tmp1))
# sw.testcase2$b[sw.testcase2$band == 1] <- sapply(sw.testcase2$hsix[sw.testcase2$band == 1],
    # function(sca) exp(slopeb.tmp1 * sca + iceptb.tmp1))
# sw.testcase2$a[sw.testcase2$band == 2] <- sapply(sw.testcase2$hsix[sw.testcase2$band == 2],
    # function(sca) exp(slopea.tmp2 * sca + icepta.tmp2))
# sw.testcase2$b[sw.testcase2$band == 2] <- sapply(sw.testcase2$hsix[sw.testcase2$band == 2],
    # function(sca) exp(slopeb.tmp2 * sca + iceptb.tmp2))
# ## Give names to the locations:
# sw.testcase2$loc.name <- as.character(1:180)

# #quartz()
# par(mfrow = c(1,2), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.7,0))
# plot(sw.testcase2$hsix[sw.testcase2$band == 1], 
   # sw.testcase2$a[sw.testcase2$band == 1], pch = 16, 
   # cex = 0.5, ylim = range(sw.testcase2$a), xlim = range(sw.testcase2$hsix))
# points(sw.testcase2$hsix[sw.testcase2$band == 2], 
   # sw.testcase2$a[sw.testcase2$band == 2], pch = 16, 
   # cex = 0.5, col = 2)
# plot(sw.testcase2$hsix[sw.testcase2$band == 1], 
   # sw.testcase2$b[sw.testcase2$band == 1], pch = 16, 
   # cex = 0.5, ylim = range(sw.testcase2$b), xlim = range(sw.testcase2$hsix))
# points(sw.testcase2$hsix[sw.testcase2$band == 2], 
   # sw.testcase2$b[sw.testcase2$band == 2], pch = 16, 
   # cex = 0.5, col = 2)


# ## 4.step (fap generation and var.list):
# vars.list.tmp <- apply(sw.testcase2[,c("a","b","hsix")], MAR = 1, FUN = function(vec)
   # {
   	# list(fap = rbeta(750, shape1 = vec[1], shape2 = vec[2]),
   	     # covar = vec[3])
   # })
# names(vars.list.tmp) <- sw.testcase2$loc.name

# ## 5.step (create mis-labeling):
# sw.testcase2$supposed.band <- sw.testcase2$band 
# # sw.testcase2$supposed.band[sample(1:60, size = 20)] <- 2
# # sw.testcase2$supposed.band[sample(61:180, size = 40)] <- 1
# # sw.testcase2[,c("band","supposed.band")]
# # table(sw.testcase2[,c("band","supposed.band")])
# # table(sw.testcase2[,c("supposed.band")]) / 180
# # table(sw.testcase2[,c("band")]) / 180
# ## All right.
# tau.init.tmp <- table(sw.testcase2[,c("supposed.band")]) / 180

# ## 6.step (individual beta parameters and lm under mis-labeling):
# betapar.tmp <- data.frame(t(sapply(vars.list.tmp, function(lst) 
   # {
# #   	lst <- vars.list.tmp[[1]]
   	# r.tmp <- try(fitbeta.fun(lst[[1]]))
   	# r.tmp$par
   	# })))
# colnames(betapar.tmp) <- c("a.est","b.est")
# sw.testcase2$a.est <- betapar.tmp$a.est
# sw.testcase2$b.est <- betapar.tmp$b.est
# # sw.testcase2[1:10,]
# # sw.testcase2[,c("a","a.est","b","b.est")]
# quartz(height = 5, width = 10)
# par(mfrow = c(1,2), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.7,0))
# plot(sw.testcase2[,c("a")], sw.testcase2[,c("a.est")], pch = 16, cex = 0.7)
# abline(c(0,1))
# plot(sw.testcase2[,c("b")], sw.testcase2[,c("b.est")], pch = 16, cex = 0.7)
# abline(c(0,1))

# pars.init.tmp <- matrix(c(lm(log(a.est) ~ hsix, data = sw.testcase2,
    # subset = supposed.band == 1)$coef,
    # lm(log(b.est) ~ hsix, data = sw.testcase2, subset = supposed.band == 1)$coef,
    # lm(log(a.est) ~ hsix, data = sw.testcase2, subset = supposed.band == 2)$coef,
    # lm(log(b.est) ~ hsix, data = sw.testcase2, subset = supposed.band == 2)$coef), 
    # ncol = 4, byrow = TRUE)
# truepars.testcase2
# pars.init.tmp

# #quartz()
# par(mfrow = c(1,2), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.7,0))
# plot(sw.testcase2$hsix[sw.testcase2$band == 1], 
   # log(sw.testcase2$a.est[sw.testcase2$band == 1]), pch = 16, 
   # cex = 0.5, col = "darkgrey", 
   # ylim = range(c(log(sw.testcase2$a), log(sw.testcase2$a.est))),
   # xlim = range(sw.testcase2$hsix))
# lines(sw.testcase2$hsix[sw.testcase2$band == 1], 
   # log(sw.testcase2$a[sw.testcase2$band == 1]))
# points(sw.testcase2$hsix[sw.testcase2$band == 2], 
   # log(sw.testcase2$a.est[sw.testcase2$band == 2]), pch = 16, 
   # cex = 0.5, col = "salmon")
# lines(sw.testcase2$hsix[sw.testcase2$band == 2], 
   # log(sw.testcase2$a[sw.testcase2$band == 2]), col = 2)
# abline(pars.init.tmp[1,1:2], col = 1, lty = 2)
# abline(pars.init.tmp[2,1:2], col = 2, lty = 2)
# plot(sw.testcase2$hsix[sw.testcase2$band == 1], 
   # log(sw.testcase2$b.est[sw.testcase2$band == 1]), pch = 16, 
   # cex = 0.5, col = "darkgrey", 
   # ylim = range(c(log(sw.testcase2$b), log(sw.testcase2$b.est))),
   # xlim = range(sw.testcase2$hsix))
# lines(sw.testcase2$hsix[sw.testcase2$band == 1], 
   # log(sw.testcase2$b[sw.testcase2$band == 1]))
# points(sw.testcase2$hsix[sw.testcase2$band == 2], 
   # log(sw.testcase2$b.est[sw.testcase2$band == 2]), pch = 16, 
   # cex = 0.5, col = "salmon")
# lines(sw.testcase2$hsix[sw.testcase2$band == 2], 
   # log(sw.testcase2$b[sw.testcase2$band == 2]), col = 2)
# abline(pars.init.tmp[1,3:4], col = 1, lty = 2)
# abline(pars.init.tmp[2,3:4], col = 2, lty = 2)

# convergence.crit <- 1e-6 
# control.optim <- list(maxit = 10000, fnscale = -1)

# system.time(
# nomisclass.em.testcase2 <- EM.logpar.fun(pars.init = pars.init.tmp, tau.init = tau.init.tmp,
   # vars.list = vars.list.tmp, method = "BFGS", gr = wloglikgrad.logpar.fun)
# )
# ## Conv.crit should be 1e-5?

# round(nomisclass.em.testcase2$location.band.prob, 2)
# as.numeric(which(nomisclass.em.testcase2$location.band.prob[1,] > 0.5))
# as.numeric(which(nomisclass.em.testcase2$location.band.prob[2,] > 0.5))
# ## Most of the group 1 is misclassified as class 2; only 25 remained as class 1.
# sw.testcase2$pred.nomisclass <- NA
# sw.testcase2$pred.nomisclass[nomisclass.em.testcase2$location.band.prob[1,] > 0.5] <- 1
# sw.testcase2$pred.nomisclass[nomisclass.em.testcase2$location.band.prob[2,] > 0.5] <- 2

# round(nomisclass.em.testcase2$par, 4)
# truepars.testcase2
# pars.init.tmp
# lapply(nomisclass.em.testcase2$history$models, function(lst)
    # t(sapply(lst, function(lst2) round(lst2$par, 4))))

# nomisclass.em.testcase2$history$logliks

# #quartz()
# par(mfrow = c(1,2), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.7,0))
# plot(sw.testcase2$hsix[sw.testcase2$pred.nomisclass == 1], 
   # log(sw.testcase2$a.est[sw.testcase2$pred.nomisclass == 1]), pch = 1, 
   # cex = 0.7, col = "darkgrey", 
   # ylim = range(c(log(sw.testcase2$a), log(sw.testcase2$a.est))),
   # xlim = range(sw.testcase2$hsix))
# points(sw.testcase2$hsix[sw.testcase2$pred.nomisclass == 2], 
   # log(sw.testcase2$a.est[sw.testcase2$pred.nomisclass == 2]), pch = 1, 
   # cex = 0.7, col = "salmon")
# points(sw.testcase2$hsix[sw.testcase2$pred.nomisclass == 2 & 
   # sw.testcase2$band == 1], 
   # log(sw.testcase2$a.est[sw.testcase2$pred.nomisclass == 2 & 
   # sw.testcase2$band == 1]), pch = 1, cex = 0.7, col = "blueviolet")
# points(sw.testcase2$hsix[sw.testcase2$band == 1], 
   # log(sw.testcase2$a[sw.testcase2$band == 1]), pch = 16,
   # cex = 0.4)
# points(sw.testcase2$hsix[sw.testcase2$band == 2], 
   # log(sw.testcase2$a[sw.testcase2$band == 2]), pch = 16,
   # cex = 0.4, col = 2)
# abline(pars.init.tmp[1,1:2], col = 1, lty = 2)
# abline(pars.init.tmp[2,1:2], col = 2, lty = 2)
# abline(nomisclass.em.testcase2$par[1,1:2], col = 1)
# abline(nomisclass.em.testcase2$par[2,1:2], col = 2)
# plot(sw.testcase2$hsix[sw.testcase2$pred.nomisclass == 1], 
   # log(sw.testcase2$b.est[sw.testcase2$pred.nomisclass == 1]), pch = 1, 
   # cex = 0.7, col = "darkgrey", 
   # ylim = range(c(log(sw.testcase2$b), log(sw.testcase2$b.est))),
   # xlim = range(sw.testcase2$hsix))
# points(sw.testcase2$hsix[sw.testcase2$pred.nomisclass == 2], 
   # log(sw.testcase2$b.est[sw.testcase2$pred.nomisclass == 2]), pch = 1, 
   # cex = 0.7, col = "salmon")
# points(sw.testcase2$hsix[sw.testcase2$pred.nomisclass == 2 & 
   # sw.testcase2$band == 1], 
   # log(sw.testcase2$b.est[sw.testcase2$pred.nomisclass == 2 & 
   # sw.testcase2$band == 1]), pch = 1, cex = 0.7, col = "blueviolet")
# points(sw.testcase2$hsix[sw.testcase2$band == 1], 
   # log(sw.testcase2$b[sw.testcase2$band == 1]), pch = 16,
   # cex = 0.4)
# points(sw.testcase2$hsix[sw.testcase2$band == 2], 
   # log(sw.testcase2$b[sw.testcase2$band == 2]), pch = 16, 
   # col = 2, cex = 0.4)
# abline(pars.init.tmp[1,3:4], col = 1, lty = 2)
# abline(pars.init.tmp[2,3:4], col = 2, lty = 2)
# abline(nomisclass.em.testcase2$par[1,3:4], col = 1)
# abline(nomisclass.em.testcase2$par[2,3:4], col = 2)
# ## The locations in the intersection were all allocated to the larger (red)
# ## group. Ok, this is completely reasonable: the EM algorithm wanted to 
# ## find two relatively well-separated groups, and it found such groups. It is
# ## not the truth, but it's because the truth is fuzzy.


# ## -------------------------------------------------------------------------------------
# ## The same case with worse initial parameters 
# ## due to mild mislabeling of locations:

# ## All the former parameters and simulated values are used from the first run
# ## until step 4. Then:

# misclass1.testcase2 <- lapply(c(12,30), sample, x = 1:60)
# misclass2.testcase2 <- lapply(c(24,60), sample, x = 61:180)
# # misclass1.testcase2[[9]] <- sample(x = 1:60, size = 40)
# # misclass1.testcase2[[9]] <- sample(x = 61:180, size = 80)
# misclass.em.testcase2 <- list()
# times.misclass.testcase2 <- numeric(2)
# for(kk in 1:2)
  # {
   # ## 5.step (create now mis-labeling):
   # sw.testcase2$supposed.band2 <- sw.testcase2$band 
   # sw.testcase2$supposed.band2[misclass1.testcase2[[kk]]] <- 2
   # sw.testcase2$supposed.band2[misclass2.testcase2[[kk]]] <- 1
   # tau.init.tmp <- table(sw.testcase2[,c("supposed.band2")]) / 180

   # ## 6.step (individual beta parameters and lm under mis-labeling):
   # ## The same as before, as all the simulated values are the same, apart from the linear fit.
   # pars.init.tmp <- matrix(c(
     # lm(log(a.est) ~ hsix, data = sw.testcase2, subset = supposed.band2 == 1)$coef,
     # lm(log(b.est) ~ hsix, data = sw.testcase2, subset = supposed.band2 == 1)$coef,
     # lm(log(a.est) ~ hsix, data = sw.testcase2, subset = supposed.band2 == 2)$coef,
     # lm(log(b.est) ~ hsix, data = sw.testcase2, subset = supposed.band2 == 2)$coef), 
     # ncol = 4, byrow = TRUE)
   # # pars.init.tmp
   # times.misclass.testcase2[kk] <- system.time(
      # misclass.em.testcase2[[kk]] <- EM.logpar.fun(pars.init.tmp, tau.init.tmp, vars.list.tmp, 
      # convergence.crit = 1e-6, control.optim = list(maxit = 10000, fnscale = -1), 
      # method = "BFGS", gr = wloglikgrad.logpar.fun))[3]
   # save.image(WORKSPACE.NAME)
  # }
# lapply(misclass.em.testcase2, function(lst)
    # round(lst$location.band.prob, 3))
# ## no error, perfect classes
# lapply(misclass.em.testcase2, function(lst)
    # round(lst$par, 4))
# truepars.testcase2
# lapply(misclass.em.testcase2, function(lst)
    # lst$history$taus)
# ## All find the same estimates as with the initially perfect classification, even













###############################################################################################
## The way I found it out
############################################################################################### -------------------------------------------------------------------------------------

## A simple case with well-separated parameters:

## Covariates: let them be random uniform variables
# y1 <- runif(25)
# y2 <- runif(25)

# truepars.tmp <- matrix(c(0,-0.6,0,-0.2, 0,1,-0.5,0.4), ncol = 4, byrow = T)
# a1 <- exp(0 - 0.6 * y1)
# b1 <- exp(0 - 0.2 * y1)
# a2 <- exp(0 + 1 * y2)
# b2 <- exp(-0.5 + 0.4 * y2)
# quartz()
# par(mfrow = c(1,2), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.7,0))
# plot(y1, a1, pch = 16, cex = 0.5, ylim = range(c(a1,a2)))
# points(y2, a2, pch = 16, cex = 0.5, col = 2)
# plot(y1, b1, pch = 16, cex = 0.5, ylim = range(c(b1,b2)))
# points(y2, b2, pch = 16, cex = 0.5, col = 2)

# # sw.tmp <- data.frame(band = rep(c(1,2), each = 25), y = c(y1,y2), a = c(a1,a2), b = c(b1,b2))

# vars.list.tmp <- apply(sw.tmp[,c("a","b","y")], MAR = 1, FUN = function(vec)
   # {
   	# list(fap = rbeta(400, shape1 = vec["a"], shape2 = vec["b"]),
   	     # covar = vec["y"])
   # })

# sw.tmp$supposed.band <- sw.tmp$band 
# misclass1 <- sample(1:25, size = 7)
# misclass2 <- sample(26:50, size = 12)
# sw.tmp$supposed.band2 <- sw.tmp$band 
# sw.tmp$supposed.band2[misclass1] <- 2 
# sw.tmp$supposed.band2[misclass2] <- 1 

# tau.init.tmp <- table(sw.tmp[,c("supposed.band")]) / 50

# betapar.tmp <- data.frame(t(sapply(vars.list.tmp, function(lst) 
   # {
# #   	lst <- vars.list.tmp[[1]]
   	# r.tmp <- try(fitbeta.fun(lst[[1]]))
   	# r.tmp$par
   	# })))
# colnames(betapar.tmp) <- c("a.est","b.est")
# sw.tmp$a.est <- betapar.tmp$a.est
# sw.tmp$b.est <- betapar.tmp$b.est
# sw.tmp[1:10,]
# sw.tmp[,c("a","a.est","b","b.est")]
# quartz(height = 3.5, width = 7)
# par(mfrow = c(1,2), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.7,0))
# plot(sw.tmp[,c("a")], sw.tmp[,c("a.est")], pch = 16, cex = 0.7, ylim = c(0,6),
   # col = c(rep("blue", 25), rep("deepskyblue", 25)))
# abline(c(0,1))
# plot(sw.tmp[,c("b")], sw.tmp[,c("b.est")], pch = 16, cex = 0.7, ylim = c(0,10),
   # col = c(rep("blue", 25), rep("deepskyblue", 25)))
# abline(c(0,1))

# # pars.init.tmp <- matrix(c(lm(log(a.est) ~ y, data = sw.tmp, subset = supposed.band == 1)$coef,
    # lm(log(b.est) ~ y, data = sw.tmp, subset = supposed.band == 1)$coef,
    # lm(log(a.est) ~ y, data = sw.tmp, subset = supposed.band == 2)$coef,
    # lm(log(b.est) ~ y, data = sw.tmp, subset = supposed.band == 2)$coef), 
    # ncol = 4, byrow = TRUE)
# pars.init.tmp

# ##  See the fit:

# quartz()
# par(mfrow = c(1,2), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.7,0))
# plot(sw.tmp$y[sw.tmp$band == 1], log(sw.tmp$a.est[sw.tmp$band == 1]), pch = 16,
   # cex = 0.6, col = "darkgrey", ylim = range(log(sw.tmp$a.est)))
# abline(c(pars.init.tmp[1,1], pars.init.tmp[1,2]))
# points(sw.tmp$y[sw.tmp$band == 2], log(sw.tmp$a.est[sw.tmp$band == 2]), pch = 16,
   # cex = 0.6, col = "salmon")
# abline(c(pars.init.tmp[2,1], pars.init.tmp[2,2]), col = "red")
# plot(sw.tmp$y[sw.tmp$band == 1], log(sw.tmp$b.est[sw.tmp$band == 1]), pch = 16,
   # cex = 0.6, col = "darkgrey", ylim = range(log(sw.tmp$b.est)))
# abline(c(pars.init.tmp[1,3], pars.init.tmp[1,4]))
# points(sw.tmp$y[sw.tmp$band == 2], log(sw.tmp$b.est[sw.tmp$band == 2]), pch = 16,
   # cex = 0.6, col = "salmon")
# abline(c(pars.init.tmp[2,3], pars.init.tmp[2,4]), col = "red")
# ## All right. A bit scattered, but okay.

# ## Try out all the functions, one by one:

# ## 1. beta.logpar.fun

# pars <- pars.init.tmp[1,]
# vars <- vars.list.tmp[[1]][[2]]
# ## step-by-step
# # beta.logpar.fun <- function(pars, fap, vars = NULL, if.log = FALSE)
  # # {
   # par.alpha <- pars[1:(length(pars)/2)]
   # par.beta <- pars[-(1:(length(pars)/2))]
   # covar <- c(1,vars)
   # alpha <- exp(sum(par.alpha * covar))
   # beta <- exp(sum(par.beta * covar))
   # dbeta(fap, shape1 = alpha, shape2 = beta, log = if.log)
  # # }
# ## see the computation of alpha and beta:
# alpha.init1 <- sapply(vars.list.tmp[1:25], function(lst, pars)
  # {
   # par.alpha <- pars[1:(length(pars)/2)]
   # covar <- c(1,lst[[2]])
   # exp(sum(par.alpha * covar))
   # }, pars = pars.init.tmp[1,]) 
# beta.init1 <- sapply(vars.list.tmp[1:25], function(lst, pars)
  # {
   # par.beta <- pars[-(1:(length(pars)/2))]
   # covar <- c(1,lst[[2]])
   # exp(sum(par.beta * covar))
   # }, pars = pars.init.tmp[1,]) 
# alpha.init2 <- sapply(vars.list.tmp[26:50], function(lst, pars)
  # {
   # par.alpha <- pars[1:(length(pars)/2)]
   # covar <- c(1,lst[[2]])
   # alpha <- exp(sum(par.alpha * covar))
   # }, pars = pars.init.tmp[2,]) 
# beta.init2 <- sapply(vars.list.tmp[26:50], function(lst, pars)
  # {
   # par.beta <- pars[-(1:(length(pars)/2))]
   # covar <- c(1,lst[[2]])
   # exp(sum(par.beta * covar))
   # }, pars = pars.init.tmp[2,]) 
# sw.tmp$a.lmpred <- c(alpha.init1, alpha.init2)
# sw.tmp$b.lmpred <- c(beta.init1, beta.init2)

# ## Check the linear regression for the logarithm of hte parameters:
# quartz()
# par(mfrow = c(1,2), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.7,0))
# plot(y1, a1, pch = 16, cex = 0.5, ylim = range(c(a1,a2)))
# points(sapply(vars.list.tmp[1:25], function(lst) lst[[2]]), alpha.init1, col = "darkgrey")
# points(sapply(vars.list.tmp[26:50], function(lst) lst[[2]]), alpha.init2, col = "salmon")
# points(y2, a2, pch = 16, cex = 0.5, col = 2)
# plot(y1, b1, pch = 16, cex = 0.5, ylim = range(c(b1,b2)))
# points(y2, b2, pch = 16, cex = 0.5, col = 2)
# points(sapply(vars.list.tmp[1:25], function(lst) lst[[2]]), beta.init1, col = "darkgrey")
# points(sapply(vars.list.tmp[26:50], function(lst) lst[[2]]), beta.init2, col = "salmon")

# ## another view, log(paramters) vs hte covariate as above, but check the predicted
# ## parameters too:
# quartz()
# par(mfrow = c(1,2), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.7,0))
# plot(sw.tmp$y[sw.tmp$band == 1], log(sw.tmp$a.est[sw.tmp$band == 1]), pch = 16,
   # cex = 0.6, col = "darkgrey", ylim = range(log(sw.tmp$a.est)))
# abline(c(pars.init.tmp[1,1], pars.init.tmp[1,2]))
# points(sw.tmp$y[sw.tmp$band == 1], log(sw.tmp$a.lmpred[sw.tmp$band == 1]), pch = 1)
# points(sw.tmp$y[sw.tmp$band == 2], log(sw.tmp$a.est[sw.tmp$band == 2]), pch = 16,
   # cex = 0.6, col = "salmon")
# abline(c(pars.init.tmp[2,1], pars.init.tmp[2,2]), col = "red")
# points(sw.tmp$y[sw.tmp$band == 2], log(sw.tmp$a.lmpred[sw.tmp$band == 2]), pch = 1, 
   # col = "red")
# plot(sw.tmp$y[sw.tmp$band == 1], log(sw.tmp$b.est[sw.tmp$band == 1]), pch = 16,
   # cex = 0.6, col = "darkgrey", ylim = range(log(sw.tmp$b.est)))
# abline(c(pars.init.tmp[1,3], pars.init.tmp[1,4]))
# points(sw.tmp$y[sw.tmp$band == 1], log(sw.tmp$b.lmpred[sw.tmp$band == 1]), pch = 1)
# points(sw.tmp$y[sw.tmp$band == 2], log(sw.tmp$b.est[sw.tmp$band == 2]), pch = 16,
   # cex = 0.6, col = "salmon")
# abline(c(pars.init.tmp[2,3], pars.init.tmp[2,4]), col = "red")
# points(sw.tmp$y[sw.tmp$band == 2], log(sw.tmp$b.lmpred[sw.tmp$band == 2]), pch = 1, 
   # col = "red")
# ## All right. A bit scattered, but okay.


# ## the look of the fitted density (beta.logpar.fun) and the true density:
# probs1 <- sapply(vars.list.tmp[1:25], function(lst)
  # sapply(lst[[1]], beta.logpar.fun, pars = pars.init.tmp[1,], vars = lst[[2]], if.log = FALSE))
# probs2 <- sapply(vars.list.tmp[26:50], function(lst)
  # sapply(lst[[1]], beta.logpar.fun, pars = pars.init.tmp[2,], vars = lst[[2]], if.log = FALSE))
# quartz(height = 7.5, width = 12.5)
# par(mfrow = c(5,5), mar = c(1,1,1,1), mgp = c(1.5,0.5,0))
# for(ii in 1:25)
  # {
   # plot(vars.list.tmp[[ii]][[1]], probs1[,ii], pch = 16, cex = 0.5) 
   # lines(seq(0,1,by=0.001), 
      # dbeta(seq(0,1,by=0.001), shape1=alpha.init1[ii], shape2=beta.init1[ii]), col = 2)   
   # lines(sort(vars.list.tmp[[ii]][[1]]), dbeta(sort(vars.list.tmp[[ii]][[1]]),
      # shape1=sw.tmp$a[ii], shape2=sw.tmp$b[ii]), col = "deepskyblue")   
  # }
# quartz(height = 7.5, width = 12.5)
# par(mfrow = c(5,5), mar = c(1,1,1,1), mgp = c(1.5,0.5,0))
# for(ii in 1:25)
  # {
   # plot(vars.list.tmp[[ii+25]][[1]], probs2[,ii], pch = 16, cex = 0.5) 
   # lines(seq(0,1,by=0.001), 
      # dbeta(seq(0,1,by=0.001), shape1=alpha.init2[ii], shape2=beta.init2[ii]), col = 2)   
   # lines(sort(vars.list.tmp[[ii+25]][[1]]), dbeta(sort(vars.list.tmp[[ii+25]][[1]]),
      # shape1=sw.tmp$a[ii+25], shape2=sw.tmp$b[ii+25]), col = "deepskyblue")   
  # }
# ## big difference from the true density: 14, 15, 18 in the first group
# ## no big differnece in the second.

# ## qq-plots:
# quartz(height = 7.5, width = 12.5)
# par(mfrow = c(5,5), mar = c(1,1,1,1), mgp = c(1.5,0.5,0))
# for(ii in 1:25)
  # {
   # plot(qbeta(ppoints(400), shape1=sw.tmp$a.lmpred[ii], shape2=sw.tmp$b.lmpred[ii]), 
      # sort(vars.list.tmp[[ii]][[1]]), pch = 16, cex = 0.5) 
   # abline(c(0,1), col = 2)
  # }
# quartz(height = 7.5, width = 12.5)
# par(mfrow = c(5,5), mar = c(1,1,1,1), mgp = c(1.5,0.5,0))
# for(ii in 26:50)
  # {
   # plot(qbeta(ppoints(400), shape1=sw.tmp$a.lmpred[ii], shape2=sw.tmp$b.lmpred[ii]), 
      # sort(vars.list.tmp[[ii]][[1]]), pch = 16, cex = 0.5) 
   # abline(c(0,1), col = 2)
  # }
# ## all look all right, even those few that looked ugly on the density plots.
# ## A bit suspicious that I cannot identify those that are so bad fits on the density
# ## plots.

# ## 2. Tsingle.logpar.fun
# misclass1 <- sample(1:25, size = 7)
# misclass2 <- sample(26:50, size = 12)
# sw.tmp$supposed.band2 <- sw.tmp$band 
# sw.tmp$supposed.band2[misclass1] <- 2 
# sw.tmp$supposed.band2[misclass2] <- 1 
# tau.init.tmp <- table(sw.tmp[,c("supposed.band2")]) / 50
# pars.init.tmp <- matrix(c(lm(log(a.est) ~ y, data = sw.tmp, subset = supposed.band2 == 1)$coef,
    # lm(log(b.est) ~ y, data = sw.tmp, subset = supposed.band2 == 1)$coef,
    # lm(log(a.est) ~ y, data = sw.tmp, subset = supposed.band2 == 2)$coef,
    # lm(log(b.est) ~ y, data = sw.tmp, subset = supposed.band2 == 2)$coef), 
    # ncol = 4, byrow = TRUE)
# pars.init.tmp
# truepars.tmp
# ## look at the points:
# quartz()
# par(mfrow = c(1,2), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.7,0))
# plot(sw.tmp$y[sw.tmp$supposed.band2 == 1],
   # log(sw.tmp$a.est[sw.tmp$supposed.band2 == 1]), pch = 16,
   # cex = 0.6, col = "darkgrey", ylim = range(log(sw.tmp$a.est)))
# abline(c(pars.init.tmp[1,1], pars.init.tmp[1,2]))
# points(sw.tmp$y[sw.tmp$supposed.band2 == 2],
   # log(sw.tmp$a.est[sw.tmp$supposed.band2 == 2]), pch = 16,
   # cex = 0.6, col = "salmon")
# abline(c(pars.init.tmp[2,1], pars.init.tmp[2,2]), col = "red")
# plot(sw.tmp$y[sw.tmp$supposed.band2 == 1],
   # log(sw.tmp$b.est[sw.tmp$supposed.band2 == 1]), pch = 16,
   # cex = 0.6, col = "darkgrey", ylim = range(log(sw.tmp$b.est)))
# abline(c(pars.init.tmp[1,3], pars.init.tmp[1,4]))
# points(sw.tmp$y[sw.tmp$supposed.band2 == 2],
   # log(sw.tmp$b.est[sw.tmp$supposed.band2 == 2]), pch = 16,
   # cex = 0.6, col = "salmon")
# abline(c(pars.init.tmp[2,3], pars.init.tmp[2,4]), col = "red")
# # Tsingle.logpar.fun <- function(tau, pars, vars.list)
   # # {
   	# tau <- tau.init.tmp
   	# pars <- pars.init.tmp
   	# vars.list <- vars.list.tmp[[42]]
   	# fap <- vars.list[[1]]
   	# vars <- vars.list[[2]]
   	# ## for all FAP values and the covariates, compute all the models in the 
   	# ## different classes/bands (the likelihood conditional on the class/band):
   	# f.tmp <- apply(pars, MAR = 1, FUN = function(vec)
         # sapply(fap, beta.logpar.fun, vars = vars, pars = vec, if.log = FALSE))
   	# ## take the column-wise product of all the individual faps:
   	# g.tmp <- apply(f.tmp, MAR = 2, prod)
    # g.tmp[g.tmp == Inf] <- 1e+304
   	# ## multiply each with the prior probability of the class:
   	# t.tmp <- tau * g.tmp
   	# ## finally, normalize with the sum:
   	# t.tmp / sum(t.tmp)
   # # }
# ## Looks all right, pretty much what you would expect.

# ## group-wise application:
# t.tmp <- round(sapply(vars.list.tmp, Tsingle.logpar.fun, tau = tau.init.tmp,
   	# pars = pars.init.tmp), 3)
# ## in the first group, 6 misclassification, in the second, 7.
# which(t.tmp[2,1:25] >= 0.5)
# which(t.tmp[1,26:50] >= 0.5) + 25
# misclass1
# misclass2
# ## only some of them intersects with the  misclassified data 

# ##  A detailed EM algorithm: 
# ## use these as weights in an lm, get initial values for par this way, optim from
# ## there

# sw.tmp$band2 <- NA
# sw.tmp$band2[which(t.tmp[2,] >= 0.5)] <- 2
# sw.tmp$band2[which(t.tmp[1,] > 0.5)] <- 1
# ## lm fits:
# pars.tmp2 <- matrix(c(
  # lm(log(a.est) ~ y, data = sw.tmp, weights = t.tmp[1,])$coef,
  # lm(log(b.est) ~ y, data = sw.tmp, weights = t.tmp[1,])$coef,
  # lm(log(a.est) ~ y, data = sw.tmp, weights = t.tmp[2,])$coef,
  # lm(log(b.est) ~ y, data = sw.tmp, weights = t.tmp[2,])$coef), 
  # ncol = 4, byrow = TRUE)
# ## look at the points:
# quartz()
# par(mfrow = c(1,2), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.7,0))
# plot(sw.tmp$y[sw.tmp$band2 == 1],
   # log(sw.tmp$a.est[sw.tmp$band2 == 1]), pch = 16,
   # cex = 0.6, col = "darkgrey", ylim = range(log(sw.tmp$a.est)))
# abline(c(pars.tmp2[1,1], pars.tmp2[1,2]))
# points(sw.tmp$y[sw.tmp$band2 == 2],
   # log(sw.tmp$a.est[sw.tmp$band2 == 2]), pch = 16,
   # cex = 0.6, col = "salmon")
# abline(c(pars.tmp2[2,1], pars.tmp2[2,2]), col = "red")
# plot(sw.tmp$y[sw.tmp$band2 == 1],
   # log(sw.tmp$b.est[sw.tmp$band2 == 1]), pch = 16,
   # cex = 0.6, col = "darkgrey", ylim = range(log(sw.tmp$b.est)))
# abline(c(pars.tmp2[1,3], pars.tmp2[1,4]))
# points(sw.tmp$y[sw.tmp$band2 == 2],
   # log(sw.tmp$b.est[sw.tmp$band2 == 2]), pch = 16,
   # cex = 0.6, col = "salmon")
# abline(c(pars.tmp2[2,3], pars.tmp2[2,4]), col = "red")
# ## Hm. The b of grey changed a lot.
 
## pars: a vector, parameters of only one band
## vars.list: a list of two-element list
## weights: the row of Tmatrix corresponding to the same band as pars. 
# wloglik.beta.logpar.fun <- function(pars, vars.list, weights)
   # {
# #   	pars <- pars.tmp2[1,]
# #   	weights <- t.tmp[1,]
# #   	vars.list <- vars.list.tmp
   	# ## Attach the  weights to all locations: this is a vector of length L,
   	# ## the weight with which the location contributes to the total
   	# ## log-likelihood of each of the bands
   	# vars.modlist.tmp <- list()
   	# for(ii in 1:length(vars.list))
   	   # vars.modlist.tmp[[ii]] <- list(fap = vars.list.tmp[[ii]]$fap, 
   	   # covar = vars.list.tmp[[ii]]$covar, weigths = weights[ii])
   	# r.tmp <- sapply(vars.modlist.tmp, function(lst, pars)
   	   # {
# #   	   	lst <- vars.modlist.tmp[[3]]
     	# f.tmp <- sapply(lst[[1]], beta.logpar.fun, vars = lst[[2]], pars = pars, if.log = TRUE)
        # lst[[3]]*sum(f.tmp)
   	   # }, pars = pars)
   	# sum(r.tmp)
   # }
# o1.tmp <- optim(par = pars.tmp2[1,], fn = wloglik.beta.logpar.fun, vars.list = vars.list.tmp,
   # weights = t.tmp[1,], control = list(fnscale = -1))
# o2.tmp <- optim(par = pars.tmp2[2,], fn = wloglik.beta.logpar.fun, vars.list = vars.list.tmp,
   # weights = t.tmp[2,], control = list(fnscale = -1))
# ## Ok, it seems to improve.
# tau.tmp2 <- apply(t.tmp, M = 1, FUN = mean)
     # # 1      2 
# # 0.5294 0.4706 
# ## Gets better, too.
# pars.tmp3 <-  matrix(c(o1.tmp$par, o2.tmp$par), ncol = 4, byrow = TRUE)
# t.tmp3 <- sapply(vars.list.tmp, Tsingle.logpar.fun, tau = tau.tmp2,
   	# pars = pars.tmp3)
# ## in the first group, 4 misclassifications, in the second, 7 (the same).

# o1.tmp3 <- optim(par = pars.tmp3[1,], fn = wloglik.beta.logpar.fun, vars.list = vars.list.tmp,
   # weights = t.tmp3[1,], control = list(fnscale = -1))
# o2.tmp3 <- optim(par = pars.tmp3[2,], fn = wloglik.beta.logpar.fun, vars.list = vars.list.tmp,
   # weights = t.tmp3[2,], control = list(fnscale = -1))
# pars.tmp4 <-  matrix(c(o1.tmp3$par, o2.tmp3$par), ncol = 4, byrow = TRUE)
# tau.tmp3 <- apply(t.tmp3, M = 1, FUN = mean)
     # # 1      2 
# # 0.5568363 0.4431637 
# ## Worse now, both pars[2,] and tau.
# t.tmp4 <- sapply(vars.list.tmp, Tsingle.logpar.fun, tau = tau.tmp3,
   	# pars = pars.tmp4)
# round(t.tmp4, 2)
# round(t.tmp3, 2)
# which(t.tmp4[2,1:25] >= 0.5)
# which(t.tmp4[1,26:50] >= 0.5) + 25
# ## Four misclassifications in the first band (the same), 5 in the second: now
# ## the second class improved!

# o1.tmp4 <- optim(par = pars.tmp4[1,], fn = wloglik.beta.logpar.fun, vars.list = vars.list.tmp,
   # weights = t.tmp4[1,], control = list(fnscale = -1))
# o2.tmp4 <- optim(par = pars.tmp4[2,], fn = wloglik.beta.logpar.fun, vars.list = vars.list.tmp,
   # weights = t.tmp4[2,], control = list(fnscale = -1))
# pars.tmp5 <-  matrix(c(o1.tmp4$par, o2.tmp4$par), ncol = 4, byrow = TRUE)
# tau.tmp4 <- apply(t.tmp4, M = 1, FUN = mean)
     # # 1      2 
# # 0.5269925 0.4730075 
# ## Better tau, a little improvement in the b parameters of both bands.
# t.tmp5 <- sapply(vars.list.tmp, Tsingle.logpar.fun, tau = tau.tmp4,
   	# pars = pars.tmp5)
# round(t.tmp4, 2)
# round(t.tmp5, 2)
# which(t.tmp5[2,1:25] >= 0.5)
# which(t.tmp5[1,26:50] >= 0.5) + 25
# ## 3 misclassifications in the first band (the same), 2 only in the second: now
# ## it's improved! Even for these, the probabilities are improving.

# o1.tmp5 <- optim(par = pars.tmp5[1,], fn = wloglik.beta.logpar.fun, vars.list = vars.list.tmp,
   # weights = t.tmp5[1,], control = list(fnscale = -1))
# o2.tmp5 <- optim(par = pars.tmp5[2,], fn = wloglik.beta.logpar.fun, vars.list = vars.list.tmp,
   # weights = t.tmp5[2,], control = list(fnscale = -1))
# pars.tmp6 <-  matrix(c(o1.tmp5$par, o2.tmp5$par), ncol = 4, byrow = TRUE)
# tau.tmp5 <- apply(t.tmp5, M = 1, FUN = mean)
     # # 1      2 
# # 0.4875918 0.5124082 
# ## Better tau, some improvement in the b parameters of both bands.
# t.tmp6 <- sapply(vars.list.tmp, Tsingle.logpar.fun, tau = tau.tmp5,
   	# pars = pars.tmp6)
# round(t.tmp6, 2)
# round(t.tmp5, 2)
# which(t.tmp6[2,1:25] >= 0.5)
# which(t.tmp6[1,26:50] >= 0.5) + 25
# ## Huh, this is now perfect!!! Almost.

# o1.tmp6 <- optim(par = pars.tmp6[1,], fn = wloglik.beta.logpar.fun, vars.list = vars.list.tmp,
   # weights = t.tmp6[1,], control = list(fnscale = -1))
# o2.tmp6 <- optim(par = pars.tmp6[2,], fn = wloglik.beta.logpar.fun, vars.list = vars.list.tmp,
   # weights = t.tmp6[2,], control = list(fnscale = -1))
# pars.tmp7 <-  matrix(c(o1.tmp6$par, o2.tmp6$par), ncol = 4, byrow = TRUE)
# tau.tmp6 <- apply(t.tmp6, M = 1, FUN = mean)
# ## Hey, it's done! Not a bad estimate. Now try to put it into a loop.

# # wloglikgrad.logpar.fun <- function(pars, vars.list, weights)
   # {
    # par.alpha <- pars[1:(length(pars)/2)]
    # par.beta <- pars[-(1:(length(pars)/2))]
    # g.tmp <- sapply(vars.list, function(lst)
       # {
        # covar <- c(1,lst[[2]])
        # alpha <- exp(sum(par.alpha * covar))
        # beta <- exp(sum(par.beta * covar))
        # n <- length(lst[[1]])
        # sumlogfap <- sum(log(lst[[1]]))
        # sumlogcompfap <- sum(log(1-lst[[1]]))
        # dga <- n * (digamma(alpha+beta) - digamma(alpha)) + sumlogfap
        # dgb <- n * (digamma(alpha+beta) - digamma(beta)) + sumlogcompfap
        # c(dga * alpha,
          # dga * alpha * lst[[2]],
          # dgb * beta,
          # dgb * beta * lst[[2]])
        # })
   	# matrix(weights, nrow = 1) %*% t(g.tmp)
   # }
# wloglikgrad.logpar.fun(pars = res0$par[1,], vars = vars.list.tmp, 
   # weights = res0$location.band.prob[1,])
# wloglikgrad.logpar.fun(pars = res0$par[2,], vars = vars.list.tmp, 
   # weights = res0$location.band.prob[2,])

# EM.logpar.fun <- function(pars.init,
                          # tau.init, 
                          # vars.list, 
                          # convergence.crit = 1e-6, 
                          # control.optim = list(maxit = 10000, fnscale = -1),
                          # method = "Nelder-Mead", 
                          # hessian = TRUE, ...)
   # {
    # # pars.init <- pars.init.tmp
    # # tau.init <- tau.init.tmp
    # # vars.list <- vars.list.tmp
    # # convergence.crit <- 1e-6
    # # control.optim <- list(maxit = 10000, fnscale = -1)
    # # method <- "Nelder-Mead" 
    # # hessian <- TRUE   	
    # ## Initialization:
    # pars <- pars.init
   	# loglik.value <- -1e9
   	# conv.crit <- convergence.crit + 1
   	# ## The faps produced by the pipeline contain often 1s or 0s, due to rounding
   	# ## of numbers within machine-precision from 1 or 0 (10^(-70)...). These are
   	# ## intractable for dbeta, they must be changed for a reasonable number.
   	# vars.list1 <- lapply(vars.list, function(lst)
   	   # {
	   	# f.tmp <- lst[[1]]
	   	# f.tmp[f.tmp == 1] <- 1-1e-11
	   	# f.tmp[f.tmp == 0] <- 1e-11
   	   	# list(fap = f.tmp, covar = lst[[2]])
   	   # }
   	   # )
	# Tmatrix <- sapply(vars.list1, function(lst)
	   # Tsingle.logpar.fun(vars.list = lst, tau = tau.init, pars = pars))
    # tau <- tau.update.fun(Tmatrix)
   	# jj <- 1
   	# loglik.store <- loglik.value
    # models.store <- list()
    # models.store[[jj]] <- NULL
    # Tmatrix.store <- list()
    # Tmatrix.store[[jj]] <- Tmatrix
    # tau.store <- list()
    # tau.store[[jj]] <- tau
    # cat(jj, ": initialization ready \n")
   	# while(conv.crit > convergence.crit)
   	   # {
	    # jj <- jj + 1
	   	# loglik.value0 <- loglik.value
	   	# ## E-step: compute the individual class probabilities for all faps, 
	   	# ## to be able to optimize Q in the M-step:
	   	# ## M-step: using Tmatrix, optimize the classwise likelihood, that is, Q.fun
	   	# r.tmp <- list()
	   	# for(ii in 1:nrow(Tmatrix))
	   		# {
		   	 # pars.tmp <- pars[ii,]
		   	 # T.tmp <- Tmatrix[ii,]
	   		 # r.tmp[[ii]] <- optim(pars.tmp, fn = wloglik.beta.logpar.fun,
	   		   # vars.list = vars.list1, weights = T.tmp, control = control.optim, 
	   		   # hessian = hessian, method = method, ...)
	   		# }
	    # ## Update the matrix of pars (L x added dim of the two linear models):
	    # pars <- t(sapply(r.tmp, function(lst) lst$par))
	    # ## Update the probabilities of each location to belong to each of the bands:
	   	# Tmatrix <- sapply(vars.list1, function(lst) 
	   	    # Tsingle.logpar.fun(vars.list = lst, tau = tau, pars = pars))
        # ## Update the band probabilities:
	   	# tau <- tau.update.fun(Tmatrix)
	   	# ## Extract the new log-likelihood maximum: 
	    # loglik.value <- sum(sapply(r.tmp, function(lst) lst$value))
	    # ## Compute the convergence criterion:
	    # conv.crit <- abs((loglik.value - loglik.value0) / loglik.value)
	    # loglik.store[jj] <- loglik.value
	    # models.store[[jj]] <- r.tmp
	    # Tmatrix.store[[jj]] <- Tmatrix
	    # tau.store[[jj]] <- tau
	    # save(loglik.store, Tmatrix.store, tau.store, models.store,
	       # file = "partialRes_EMbetaFunctions.RObjTmp")
	    # cat(jj, " ready; conv.crit = ", conv.crit, "\n")
   	   # }
 	# v.tmp <- list(parameters = pars, band.frequency = tau,
 	    # location.band.prob = Tmatrix, optimization = r.tmp, 
 	    # history = list(logliks = loglik.store, models = models.store,
 	    # Tmatrices = Tmatrix.store, taus = tau.store))
# # 	save.image(WORKSPACE.NAME)
 	# return(v.tmp)
 	# rm(r.tmp, v.tmp, Tmatrix, models.store, Tmatrix.store, tau.store, loglik.store)
 	# gc()
   # }

# system.time(
# res0 <- EM.logpar.fun(pars.init = pars.init.tmp, tau.init = tau.init.tmp,
   # vars.list = vars.list.tmp, method = "Nelder-Mead")
# )
# ## Time: 1200s
# system.time(
# res1 <- EM.logpar.fun(pars.init = pars.init.tmp, tau.init = tau.init.tmp,
   # vars.list = vars.list.tmp, method = "BFGS")
# )
# ## Time: 840s
# system.time(
# res2 <- EM.logpar.fun(pars.init = pars.init.tmp, tau.init = tau.init.tmp,
   # vars.list = vars.list.tmp, method = "BFGS", gr = wloglikgrad.logpar.fun)
# )
# ## Time: 198s (!)
# ## For the latter two, NaNs are produced, but where the hell can they be produced?
# ## The results are the same for the first two, so it seems either the NM produces
# ## them too but does not shout, or they do not count materially. In any case,
# ## the effect is not big, as the estimates are very close to the true values.
# 2*sqrt(diag(-solve(res0$optimization[[1]]$hessian)))
# 2*sqrt(diag(-solve(res0$optimization[[2]]$hessian)))
# res0$parameters
# truepars.tmp
# 2*sqrt(diag(-solve(res1$optimization[[1]]$hessian)))
# 2*sqrt(diag(-solve(res1$optimization[[2]]$hessian)))
# res1$parameters
# truepars.tmp
# 2*sqrt(diag(-solve(res2$optimization[[1]]$hessian)))
# 2*sqrt(diag(-solve(res2$optimization[[2]]$hessian)))
# res2$parameters
# truepars.tmp
# ## And all of them contain in the confidence interval the true values.
# wloglikgrad.logpar.fun(pars = res2$par[1,], vars = vars.list.tmp, 
   # weights = res2$location.band.prob[1,])
# wloglikgrad.logpar.fun(pars = res2$par[2,], vars = vars.list.tmp, 
   # weights = res2$location.band.prob[1,])


###############################################################################################
## Unused but potentially useful functions
############################################################################################### 

## beta.fun: the (one-component) beta distribution, with a possible linear relation
## for both parameters; the number of covariates X[i] in both parameters are the
## same
##  f(p) = B(c[0] + c[1]*X[1] + ... + c[K]*X[K], d[0] + d[1]*X[1] + ... + d[K]*X[K]) *
##         * p ^ (c[0] + c[1]*X[1] + ... + c[K]*X[K]) * 
##         * (1-p) ^ (d[0] + d[1]*X[1] + ... + d[K]*X[K])
## Arguments: 
##   pars: the parameters of the linear regression; its length is always twice
##         the length of vars + 2; first half: components of alpha, second half:
##         components of beta; if it is only 2, then these are the constant terms,
##         that is, alpha and beta without dependence on any aliases or variance.
##   fap: the fap.
##   vars: the covariates; if NULL, then the one-component beta distribution does not
##         depend on any variables (constant model fro its parameters), if it is a
##         number or a vector, then covariates in the model.
##   if.log: whether I need its logarithm (for log-likelihood) or the density 
## Possible problem: what if the linear relationship gives negative 
## beta-distribution parameters? Try to fit linearly the log of the 
## parameters? Visually, it was pretty linear on the original scale.
# beta.fun <- function(pars, fap, vars = NULL, if.log = FALSE)
  # {
   # par.alpha <- pars[1:(length(pars)/2)]
   # par.beta <- pars[-(1:(length(pars)/2))]
   # covar <- c(1,vars)
   # alpha <- sum(par.alpha * covar)
   # beta <- sum(par.beta * covar)
   # dbeta(fap, shape1 = alpha, shape2 = beta, log = if.log)
  # }


## A version of T matrix computation, for linear regression of the alpha and
## beta (not the log of them).
# Tsingle.fun <- function(tau, pars, vars.list)
   # {
   	# fap <- vars.list[[1]]
   	# vars <- vars.list[[2]]
   	# ## for all FAP values and the covariates, compute all the models in the 
   	# ## different classes/bands (the likelihood conditional on the class/band):
   	# f.tmp <- apply(pars, MAR = 1, FUN = function(vec)
         # sapply(fap, beta.fun, vars = vars, pars = vec, if.log = FALSE))
   	# ## take the column-wise product of all the individual faps (????):
   	# g.tmp <- apply(f.tmp, MAR = 2, prod)
    # g.tmp[g.tmp == Inf] <- 1e+304
   	# ## multiply each with the prior probability of the class:
   	# t.tmp <- tau * g.tmp
   	# ## finally, normalize with the sum to get a real probability:
   	# t.tmp / sum(t.tmp)
   # }


# unique(matrix(c(0,0,0,0,1,1,1,1), ncol = 4, byrow = TRUE))
# unique(matrix(c(0,0,0,0,1,1,1,1), ncol = 4))
## An interesting info: unique is automatically applied column-wise.





















