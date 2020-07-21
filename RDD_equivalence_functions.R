############################################
############################################
###
### Functions for RDD equivalence testing
###
############################################
############################################


############################################
### Two-one sided test
##
## est: point estimate
## se: estimated standard error
## eps: positive equivalence range value
##      on scale of estimate
##      (assumed symmetric equivalence range)
## alpha: level of the test
############################################

###
rdd.tost = function(est, se, eps, alpha = 0.05) {
  if(eps < 0) { 
    stop("Epsilon must be > 0 for TOST equivalence test.")
  }
  p = max(pnorm((est - eps)/se, lower.tail = TRUE), pnorm((est + eps)/se, lower.tail = FALSE))
  rej = p < alpha
  inverted <- max(abs(est - abs(qnorm(alpha)) * se), abs(est + abs(qnorm(alpha)) * se))
  return(list(rej = rej, p = p, inverted = inverted))
}
##

############################################
### Equivlaence T-test
##
## est: point estimate
## se: estimated standard error
## eps: positive equivalence range value
##      in standard deviations
##      (assumed symmetric equivalence range)
## alpha: level of the test
############################################

###
rdd.equiv <- function(est, se, eps, alpha = 0.05, inv_int_search_tol = 0.001, max_search_grid = 3) {
  if(eps < 0) { 
    stop("Epsilon must be > 0 for equivalence t-test.")
  }
  
  p = pchisq(est^2/se^2, 1, eps^2/se^2)
  rej = p < alpha
  
  inverted <- tryCatch(uniroot(function(x) {
    pchisq(est^2/se^2, 1, x^2/se^2) - alpha}, 
    c(0.00001, max(10, max_search_grid * abs(est) * se)), tol = inv_int_search_tol)$root, 
    silent = TRUE, error = function(e) NA)
  
  ## if noncentrality parameter estimate is so close to zero that
  ## p < alpha with at ncp = 0, return NA
  if(pchisq(est^2/se^2, 1) < alpha) {
    inverted = NA
  } else {
    inverted <- tryCatch(uniroot(function(x) {
      pchisq(est^2/se^2, 1, x^2/se^2) - alpha}, 
      c(0.00001, max(10, max_search_grid * abs(est) * se)), tol = inv_int_search_tol)$root, 
      silent = TRUE, error = function(e) NA)
    
    if(is.na(inverted)) warning("No interval found, try increasing max_search_grid.")
  }
  
  return(list(rej = rej, p = p, inverted = inverted))
}
##

############################################
### Equivalence density test
### (using two-one-sided tests)
##
## estL: point estimate from left
## estR: point estimate from right
## seL: standard error estimate from the left
## seR: standard error estimate from the right
## eps: positive equivalence range value
##      in ratio units
##      eps_l will be 1/eps_U
## alpha: level of the test
############################################

##
rdd.tost.ratio = function(estL, estR, seL, seR, eps = 1.5, alpha = 0.05) {
  if(eps < 1) { 
    stop("Epsilon must be > 1 for density equivalence test.")
  }
  
  T1 = (estL - 1/eps * estR) / sqrt(seL^2 + (1/eps)^2*seR^2)
  T2 = (estL - eps * estR) / sqrt(seL^2 + (eps)^2*seR^2)
  
  inverted <- tryCatch(uniroot(function(x) {
    max(pnorm((estL - 1/x * estR) / sqrt(seL^2 + (1/x)^2*seR^2), lower.tail = FALSE), 
        pnorm((estL - x * estR) / sqrt(seL^2 + (x)^2*seR^2))) - alpha}, 
                               c(1, 100), tol = 0.0001)$root, silent = TRUE, error = function(e) NA)
  
  p = max(pnorm(T1, lower.tail = FALSE), pnorm(T2))
  return(list(p = p, inverted = inverted))
}
##
