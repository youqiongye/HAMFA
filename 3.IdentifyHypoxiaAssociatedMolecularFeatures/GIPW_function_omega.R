###################################################################################################################################
##
## ---------------------------------------------------------------------------------------------------------------------------------
## Part 1: Direct comparison
## 1. A general function to calculate weight for propensity score analysis.
## 2.Included weigth schemes are: (1) inverse probability weight, (2) average treatment effect
##    for the treated (ATT) weight, (3) average treament effect for the control weight, 
##    (3) matching weight, (4) overlap weight (Crump et al, 2006), (5) trapezoidal weight.
## 3. Output includes weighting method (weight), point estimation (est), standard deviation of point estimation (std),
##    and 95% confidence interval.
## 4. Code is module based, where omega function is the only variation part across different weighting 
##    methods.
##
## ---------------------------------------------------------------------------------------------------------------------------------
## Part 2: Double-robust estimation (DR)
## 1. DR point estimation.
## 2. DR sandwich type variance estimation.
## ---------------------------------------------------------------------------------------------------------------------------------
## Part 3: Balance checking

## ---------------------------------------------------------------------------------------------------------------------------------
## Latest update: February 09, 2015
####################################################################################################################################

# For funtions logit, inv.logit;
library(gdata);
library(gtools);
library(gmodels);
library(gplots);
library(formula.tools);  # For function lhs.vars;
library(Hmisc);  # For weighted mean and weighted variance;
library(Matrix);  # for function bdiag;

GIPW.std.omega <- function(dat, form.ps, weight, trt.est=T, delta=0.002, K=4) {
  # Obtain std estimator from sandwich method
  #
  # Args
  #   dat: data on which analysis is to be performed
  #   form.ps: formula for propensity score model
  #   weight: weighting method to be used
  #   trt.est: perfomr treatment effect analys? Default is TRUE, otherwise output only weights and ps.
  #   delta: closeness to non-differential point in omega function, delta = 0.002 by default
  #   K: coefficient of trapezoidal weight, K is the slope of trapezoidal edge, K = 4 by default
  #
  # Return
  #   A list composed of point estimation (est), standard deviation (std), and 95% confidence interval
  
  n <- nrow(dat);  # number of sujects
  
  out.ps <- ps.model(dat, as.formula(form.ps));
  ps.hat <- out.ps$ps.hat;  # estimated ps
  dat$ps.hat <- ps.hat;
  beta.hat <- as.numeric(coef(out.ps$fm));
  
  Q <- dat$Z*dat$ps.hat + (1 - dat$Z)*(1 - dat$ps.hat);  # denominator of generic weight
  omega <- sapply(dat$ps, calc.omega, weight = weight, delta = delta, K = K);
  W <- omega/Q;  # generic weight
  dat$W <- W;
  
  if ( !(trt.est) ) {
    ans <- list(weight = weight, W=W, ps=ps.hat);
    return(ans);
  } else {
    
    ## direct comparison for point estimation
    mu1.hat <- sum(dat$W * dat$Z * dat$Y)/sum(dat$W * dat$Z);
    mu0.hat <- sum(dat$W * (1-dat$Z)*dat$Y)/sum(dat$W * (1-dat$Z));
    est <- mu1.hat - mu0.hat;
    
    ## standard deviation estimation
    Amat <- Bmat <- 0;  # An and Bn matrix for variance calculation
    for (i in 1 : n) {
      Xi <- as.numeric(dat[i, c("X0", names(coef(out.ps$fm))[-1])]);  # X0 is the key word for intercept
      Zi <- dat[i, "Z"];
      Yi <- dat[i, "Y"];
      ei <- calc.ps.Xbeta(Xi, beta.hat);
      ei.deriv1 <- calc.ps.deriv1(Xi, beta.hat);
      ei.deriv2 <- calc.ps.deriv2(Xi, beta.hat);
      omega.ei <- omega[i];
      omegaei.deriv <- omega.derive.ei(ei, weight, delta, K);
      Qi <- Q[i];
      Qi.deriv <- 2*Zi - 1;
      
      phi <- c(Zi*(Yi-mu1.hat)*omega.ei/Qi,
               (1-Zi)*(Yi-mu0.hat)*omega.ei/Qi,
               (Zi-ei)/(ei*(1-ei))*ei.deriv1);
      
      Bmat <- Bmat + outer(phi, phi);  # Bn matrix
      
      # first row of phi's first derivative w.r.t theta
      first.row <- c( -Zi*omega.ei/Qi, 0, Zi*(Yi-mu1.hat)*ei.deriv1*(Qi*omegaei.deriv - omega.ei*Qi.deriv)/Qi^2);
      
      # second row of phi's first derivative w.r.t theta
      second.row <- c( 0, -(1-Zi)*omega.ei/Qi, (1-Zi)*(Yi-mu0.hat)*ei.deriv1*(Qi*omegaei.deriv - omega.ei*Qi.deriv)/Qi^2);
      
      # third row of phi's first derivative w.r.t theta
      tmp0 <- matrix(0, nrow=length(beta.hat), ncol=2);
      tmp1 <- -ei*(1-ei)*colVec(Xi)%*%Xi;
      third.row <- cbind(tmp0, tmp1);
      
      phi.deriv <- rbind(first.row, second.row, third.row);
      Amat <- Amat + phi.deriv;
    }
    
    Amat <- Amat/n;
    Bmat <- Bmat/n;
    Amat.inv <- solve(Amat);
    var.mat <- ( Amat.inv %*% Bmat %*% t(Amat.inv))/n;
    tmp1 <- c(1, -1, rep(0, length(beta.hat)));
    var.est <- rowVec(tmp1) %*% var.mat %*% colVec(tmp1);
    std <- sqrt(as.numeric(var.est));
    CI.lower <- est - 1.96*std;
    CI.upper <- est + 1.96*std;
    
    ans <- list(weight = weight, est = est, std = std,
                CI.lower = CI.lower, CI.upper = CI.upper, W=W, ps=ps.hat);
    return(ans);
  }
}

######################################################################################
##
## Common supporting functions
##
######################################################################################
rowVec <- function(x) {
  t(x)
}

colVec <- function(x) {
  t(t(x))
}

ps.model <- function(dat, form) {
  # ps.model to calculate propensity score
  # 
  # Args:
  #   dat: the data from which estimated ps to be calculated
  #   form: the formula used in propensity score model
  #
  # Return:
  #   Estimated propensity score and glm fitting
  
  fm <- glm(form, data = dat, family = binomial(link = "logit"))
  ps.hat <- as.numeric(predict(fm, newdata = dat, type = "response"))
  ps.hat <- pmin(pmax(0.000001, ps.hat), 0.999999)  # ps.hat cannot be exactly 0 or 1
  return(list(ps.hat = ps.hat, fm = fm))  
}

calc.omega <- function(ps, weight, delta, K) {
  # Calculate omega for each weighting method
  #
  # Args:
  #   ps: estimated proopensity score,
  #   weight: weighting method.
  #   K: wighting coefficient for trapezoidal weighting
  #
  # Return:
  #   A vector of length(ps)
  
  ans <- 0
  if (weight == "IPW") {
    ans <- 1
  } else if (weight == "MW") {
    ans <- calc.omega.MW(ps, delta)
  } else if (weight == "ATT") {
    ans <- ps
  } else if (weight == "ATC") {
    ans <- 1-ps
  } else if (weight == "OVERLAP") {
    ans <- 4*ps*(1-ps)
  } else if (weight == "TRAPEZOIDAL") {
    ans <- calc.omega.trapzd(ps=ps, delta=delta, K=K)
  } else {
    stop("Error in calc.omega: weight method does not exist!")
    # ans <- user defined omega function
  }
  
  ans
}

calc.ps.Xbeta <- function(Xmat, beta) {
  # Calculate the propensity score, e(X, beta)
  #
  # Args:
  #   X: a vector or a matrix with each column being X_i
  #   beta: the coefficient, of the same length as the nrow(X)
  #
  # Return:
  #   A scalar of matching weight
  
  Xmat <- as.matrix(Xmat)
  tmp <- as.numeric(rowVec(beta) %*% Xmat)
  tmp <- exp(tmp)
  names(tmp) <- NULL
  return(tmp/(1 + tmp))
}

calc.ps.deriv1 <- function(Xmat, beta) {
  # Calculate the derivative of propensity score w.r.t. beta.
  #
  # Args:
  #   X: a matrix with each column being X_i,
  #   beta: the coefficient, of the same length as the nrow(X).
  #
  # Return:
  #   A vector of the same dimension as beta.
  
  Xmat <- as.matrix(Xmat)
  tmp.ps <- calc.ps.Xbeta(Xmat, beta)
  ans <- tmp.ps*(1 - tmp.ps)*t(Xmat)  # Xmat is row vector from a single line of data
  names(ans) <- rownames(ans) <- NULL
  return(t(ans))
}

calc.ps.deriv2 <- function( Xi, beta ) {
  # Calculate the second derivative of propensity score w.r.t beta.
  #
  # Args:
  #   Xi: a vector (not a matrix!)
  #   beta:  a vector of coefficients, of the same length as Xi
  # 
  # Return:
  #   A square matrix of length(Xi) by length(beta)
  
  Xi <- colVec(Xi)
  tmp.ps <- calc.ps.Xbeta(Xi, beta)
  tmp.deriv1 <- calc.ps.deriv1(Xi, beta)
  ans <- Xi %*% rowVec(tmp.deriv1)
  names(ans) <- rownames(ans) <- NULL
  ans <- (1 - 2*tmp.ps)*ans
  return(ans)
}

calc.omega.MW <- function (ps, delta) {
  # Calculate omega value for MW method
  #
  # Args:
  # ps: propensity score, scalar
  # delta: closeness to non-differential point
  #
  # Return:
  # Value of omega, scalar
  
  ans <- 0
  if (ps <= 0.5 - delta) {
    ans <- 2*ps
  } else if (ps >= 0.5 + delta) {
    ans <- 2*(1 - ps)
  } else {
    ans <- approx.omega.MW(ps, delta)
  }
  
  ans
}

approx.omega.MW <- function (ps, delta) {
  # Approximate omega function of MW at non-differential point, eta(ps)
  #
  # Args:
  # ps: propensity score, scalar
  # delta: closeness to non-differential point
  #
  # Return:
  # Approximated omega at non-differential point, scalar
  
  A <- solve.A.MW(delta)
  ans <- rowVec(A) %*% c(1, ps, ps^2, ps^3)
  
  ans
}

solve.A.MW <- function(delta) {
  # Get coefficients for approximated cubic polynomial for omega (eta(e)) function of MW
  #
  # Args:
  #   delta: pre-defined closeness to midpiece of ps.
  #
  # Return
  #   A vecotor of solved coefficients
  
  if ( delta < 0.00001 ) { 
    stop("*** ERROR in solve.a: delta too small ***")
  }
  
  tmp1 <- 0.5 - delta
  tmp2 <- 0.5 + delta
  
  D <- matrix(c(1, tmp1, tmp1^2, tmp1^3,
                0, 1, 2*tmp1, 3*tmp1^2,
                1, tmp2, tmp2^2, tmp2^3,
                0, 1, 2*tmp2, 3*tmp2^2),
              ncol = 4, nrow = 4, byrow = TRUE)
  C <- 2*c(tmp1, 1, tmp1, -1)
  A <- solve(D) %*% C  # coefficients of cubic polynomial
  
  A
}

calc.omega.trapzd <- function (ps, delta, K) {
  # Calculate omega value for trapezoidal weighting method
  #
  # Args:
  # ps: propensity score, scalar
  # delta: closeness to non-differential point
  # K: trapezoidal weighting coffecient
  #
  # Return:
  # Value of omega, scalar
  
  ans <- 0
  if ( (0 < ps) & (ps <= 1/K - delta) ) {
    ans <- K*ps
  } else if ( (1/K + delta <= ps) & (ps <= 1 - 1/K - delta) ) {
    ans <- 1
  } else if ( (1 - 1/K + delta <= ps) & (ps < 1) ) {
    ans <- K*(1 - ps)
  } else {
    ans <- approx.omega.trapzd(ps, delta, K)
  }

  ans
}

approx.omega.trapzd <- function (ps, delta, K) {
  # Approximate omega function of trapezoidal weight at non-differential point, eta(ps)
  #
  # Args:
  # ps: propensity score, scalar
  # delta: closeness to non-differential point
  # K: trapezoidal weight coefficient
  #
  # Return:
  # Approximated omega at non-differential point, scalar
  
  A <- 0
  if ( (1/K - delta < ps) & (ps < 1/K + delta) ) {
    A <- solve.A.trapzd1st(delta, K)
  } else {
    A <- solve.A.trapzd2nd(delta, K)
  }  
  ans <- rowVec(A) %*% c(1, ps, ps^2, ps^3)
  
  ans
}

solve.A.trapzd1st <- function (delta = delta, K = K) {
  # Get coefficients for approximated cubic polynomial for omega (eta(e)) function of 
  # trapezoidal weighting at the first non-differential pivot
  #
  # Args:
  #   delta: pre-defined closeness to midpiece of ps.
  #   K: coefficient of trapezoidal weight
  #
  # Return
  #   A vecotor solved coefficients
  
  if ( delta < 0.00001 ) { 
    stop("*** ERROR in solve.a: delta too small ***")
  }
  
  tmp1 <- 1/K - delta
  tmp2 <- 1/K + delta
  
  D <- matrix(c(1, tmp1, tmp1^2, tmp1^3,
                0, 1, 2*tmp1, 3*tmp1^2,
                1, tmp2, tmp2^2, tmp2^3,
                0, 1, 2*tmp2, 3*tmp2^2),
              ncol = 4, nrow = 4, byrow = TRUE)
  C <- 2*c(K*tmp1, K, 1, 0)
  A <- solve(D) %*% C  # coefficients of cubic polynomial

  A
}

solve.A.trapzd2nd <- function (delta = delta, K = K) {
  # Get coefficients for approximated cubic polynomial for omega (eta(e)) function of 
  # trapezoidal weighting at the second non-differential pivot
  #
  # Args:
  #   delta: pre-defined closeness to midpiece of ps.
  #   K: coefficient of trapezoidal weight
  #
  # Return
  #   A vecotor solved coefficients
  
  if ( delta < 0.00001 ) { 
    stop("*** ERROR in solve.a: delta too small ***")
  }
  
  tmp1 <- 1 - 1/K - delta
  tmp2 <- 1 - 1/K + delta
  
  D <- matrix(c(1, tmp1, tmp1^2, tmp1^3,
                0, 1, 2*tmp1, 3*tmp1^2,
                1, tmp2, tmp2^2, tmp2^3,
                0, 1, 2*tmp2, 3*tmp2^2),
              ncol = 4, nrow = 4, byrow = TRUE)
  C <- 2*c(1, 0, K*(1/K - delta), -K)
  A <- solve(D) %*% C  # coefficients of cubic polynomial
  
  A
}

omega.derive.ei <- function(ps = ei, weight = weight, delta = delta, K = K) {
  # Calculate the first derivative of omega function w.r.t propensity score (ei)
  #
  # Args:
  #    ps: estimated propensity score
  #    weight: selected weighting method
  #    delta: pre-defined closeness to non-differential points
  #    K: weighting coefficient for trapezoidal weighting
  #
  # Returns:
  #    The first derivative of omega w.r.t to propensity score
  
  ans <- 0
  if (weight == "IPW") {
    ans <- 0
  } else if (weight == "ATT") {
    ans <- 1
  } else if (weight == "ATC") {
    ans <- -1
  } else if (weight == "OVERLAP") {
    ans <- 4*(1 - 2*ps)
  } else if (weight == "MW") {
    ans <- omega.derive.ei.MW (ps, delta)      
  } else if (weight == "TRAPEZOIDAL") {
    ans <- omega.derive.ei.trapzd (ps, delta, K)    
  } else {
    stop("User defined first-order derivative of omega function is not provided!")
  }
}

omega.derive.ei.MW <- function (ps, delta) {
  # calculate of the first derivative of omega function w.r.t propensity score with approximation
  # in MW method
  #
  # Args:
  #    ps: propensity score
  #    delta: pre-defined closeness to non-differential points
  #
  # Returns:
  #    The first derivative of omega (MW method) w.r.t to propensity score,
  #    approximation at non-differential point is included.
  
  ans <- 0
  if ( (0 < ps) & (ps <= 0.5 - delta) ) {
    ans <- 2*ps
  } else if ( (0.5 + delta <= ps) & (ps <1) ) {
    ans <- 2*(1 - ps)
  } else {
    A <- solve.A.MW(delta)
    ans <- A[2] + 2*A[3]*ps + 3*A[4]*ps^2
  }
  
  ans
}

omega.derive.ei.trapzd <- function (ps, delta, K) {
  # calculate of the first derivative of omega function w.r.t propensity score with approximation
  # in Trapezoidal weighting method
  #
  # Args:
  #    ps: propensity score
  #    delta: pre-defined closeness to non-differential points
  #    K: coefficient of trapezoidal weighting
  #
  # Returns:
  #    The first derivative of omega (trapezoidal weighting method) w.r.t to propensity score,
  #    approximation at non-differential point is included.
  
  ans <- 0
  
  if ( (0 < ps) & (ps <= 1/K - delta) ) {
    ans <- K
  } else if ( (1/K - delta < ps) & (ps < 1/K + delta) ) {
    A <- solve.A.trapzd1st(delta = delta, K = K)
    ans <- A[2] + 2*A[3]*ps + 3*A[4]*ps^2
  } else if ( (1/K + delta <= ps) & (ps <= 1 - 1/K - delta) ) {
    ans <- 0
  } else if ( (1 - 1/K + delta <= ps) & (ps < 1) ) {
    ans <- -K
  } else {
    A <- solve.A.trapzd2nd(delta = delta, K = K)
    ans <- A[2] + 2*A[3]*ps + 3*A[4]*ps^2
  }

  ans
}


