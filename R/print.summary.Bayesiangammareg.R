print.summary.Bayesiangammareg <-
function(x, ...){
  
  cat (" \n         #################################################################
       ###                  Bayesian Gamma Regression                ###
       ################################################################# \n")
  
  cat("\n Call: \n")
  print(x$call)
  cat("\n")
  
  printCoefmat(x$coefficients, digits=4)
  
  cat("\n deviance: \n")
  print(x$deviance)
  
  cat("\n AIC: \n")
  print(x$AIC)
  
  cat("\n BIC: \n")
  print(x$BIC)
}
