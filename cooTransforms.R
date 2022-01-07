

## --------------------------------------------------------------------------------------------------------
## Transforming functions between equatorial and ecliptic coordinate systems:
## --------------------------------------------------------------------------------------------------------

eq2ecl.fun <- function(ra, dec)
   {
   	eps <- (23 + 26/60 + 21.448/3600)*pi/180
   	c(atan((sin(pi*ra/180) * cos(eps) + tan(pi*dec/180) * sin(eps)) / cos(pi*ra/180))*180/pi,
   	  asin(sin(pi*dec/180) * cos(eps) - cos(pi*dec/180) * sin(eps) * sin(pi*ra/180))*180/pi)
   }

# eq2ecl.fun(180,0)
# sapply(-89:89, eq2ecl.fun, ra = 89.9)
ecl2eq.fun <- function(lam, beta)
   {
   	eps <- (23 + 26/60 + 21.448/3600)*pi/180
   	c(atan((sin(pi*lam/180) * cos(eps) + tan(pi*beta/180) * sin(eps)) / cos(pi*lam/180))*180/pi + 180,
   	  asin(sin(pi*beta/180) * cos(eps) - cos(pi*beta/180) * sin(eps) * sin(pi*lam/180))*180/pi)
   }
