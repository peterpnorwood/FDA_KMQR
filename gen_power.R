## Loading libraries
library(refund)
library(fda)
library(fda.usc)
library(quantreg)
library(parallel)

## Sourcing dependencies
setwd("~/KMQR")
source("KMQR_FINAL.R")

## complete_process ##
## Author: Peter Norwood and Jacob Alfieri
## Purpose: Generates Power Levels at Different Signal Levels

## Inputs:
## d  : the delta (signal) vector
## n  : a number specifying the number of observations (curves)
## k  : a number specifying the power of the exponent of the integral
##      for y generation
## x2 : if TRUE then we will use x2 as the value for generating y
##      if FALSE then we will use the average of x1,x2,x3 for generating y

## Output:
## output : a list containing the power levels for each delta level
## [[i]]  : The power level at the ith delta level
##          left number is the QR 


complete_process <- function(d, n, k, x2) {
  
  curves_d <- lapply(seq(1,1000),function(i){return(gen_data(d,n,k,x2))})
  pvalues_d <- sapply(curves_d,list_interface)
  perc_reject <- apply(pvalues_d<0.05,1,mean,na.rm=TRUE)
  return(perc_reject)
}

######################################################################################################
######################################################################################################
######################################################################################################

## Interface function ##
## Takes the list of simulated functional data and passes it to the calc_pvalues function
list_interface <- function(functional_data,pow_0) {
  p_values <- calc_pvalues(functional_data[[1]],functional_data[[2]],functional_data[[3]],pow_0=1)
  return(p_values)
}

######################################################################################################
######################################################################################################
######################################################################################################

## calc_pvalues ##
## Author: Peter Norwood and Jacob Alfieri
## Purpose: Generates pvalues for the QR and KMQR  tests

## Inputs:
## func_curves: a matrix of functional curves, where each row is a curve
##              and each column is a point defining it
## y: a vector of y values 
## x_cov: a vector or dataframe of covariates, where each row is a curve and each 
##        column is a seperate covariate
## pow_0: a number, containing power for KMQR kernel
## tau_0: a number,  containing the quantile for quantile regression

## Output:
## output: a vector of two pvalues

###### ---------------------------------- #################

calc_pvalues <- function(func_curves,y,x_cov,pow_0,tau_0 = .5) {
  
  ## Calculates FPCA scores
  fpca_scores <- fpca.sc(func_curves, var = T, pve = .99)$scores
  
  ## Quantile regression fit with parametric and nonparametric covariates
  fpca_qr <- rq(y~x_cov+fpca_scores, tau = tau_0)
  ## Quantile regression fit only with parametric covariates
  reg_qr <- rq(y~x_cov, tau = tau_0)
  
  ## Returns p-values for each delta level
  p_values <- as.matrix(c(anova(fpca_qr,reg_qr)$table$pvalue, 
                          final_function(y,fpca_scores,x_cov,tau_0,B)))
  
  row.names(p_values) <- c("QR", "KMQR")
  return(p_values)
  
}

######################################################################################################
######################################################################################################
######################################################################################################

## gen_data ###
## Author: Peter Norwood and Jacob Alfieri
## Purpose: Simulates a functional data set

## Inputs:
## delta: a number specifying the strength of the association between y and
##        functional curves
## n: a number specifying the number of observations (curves)
## p: a number specifying the power of the exponent of the integral
##    for y generation, default 1
## x2true : if TRUE then we will use x2 as the value for generating y
##          if FALSE then we will use the average of x1,x2,x3 for generating y

## Output:
## output: a list containing the simulated data set
##         [[1]] = functional curves (matrix) each row is an observation (curve),
##         each column is point defining the curve
##         [[2]] = y values (vector)
##         [[3]] = x covariates (vector)

###### ---------------------------------- #################


gen_data <- function(delta, n, p, x2true) {
  
  ## Randomly samples all coefficients outside of sapply, so that
  ## they can also be used in generating yi's
  
  x1 <- 2*rnorm(n)
  x2 <- sqrt(2)*rnorm(n)
  x3 <- rnorm(n)
  coef <- data.frame(x1,x2,x3)
  ## Calculates constants outside of sapply
  t = seq(0,1, len = 51)
  basis <- cbind(1, sin(2*pi*t), cos(2*pi*t))
  
  ## Generates functional curves
  func_curves <- t(apply(coef, 1, function(coef_0){
    ## Multiplies functional values by coefficients, adds noise and stores to dataframe func_curves
    curve <- (basis %*% c(coef_0[1], coef_0[2], coef_0[3]) + rnorm(length(t)))
  }))
  
  ## Adds row names to func_curves
  rownames(func_curves) <- paste0("Xt", 1:n)
  
  ## Randomly samples covariates outside of sapply, so that
  ## they can be stored for use in the analysis
  x_covar <- rnorm(n)
  
  ## Calculates y value for each curve
  ## Can use x2 or and average of x1, x2, x3 (with smaller delta values)
  avg <- apply(coef, 1, mean)
  
  if(x2true == TRUE)  {
    y <- 1 + delta*(x2*.5)^p + rnorm(n)
  } else {
    y <- 1 + delta*(avg*.5)^p + rnorm(n)
  }
  
  ## Combines func_curves and y into a single list
  output <- list(func_curves,y,x_covar)
  return(output)
  
}








