print.Bayesiangammareg <-
function(x, ...)
{
  
  cat (" \n         ################################################################
       ###                 Bayesian Gamma Regression                ###
       ################################################################ \n")
  
  cat("\n Call: \n")
  
  print(x$call)
  cat("\n Coefficients: \n")
  print(x$coefficients)
  
}
