criteria <-
function(X,gammaresiduals){
  if(is.null(gammaresiduals)){
    stop("Residual data is not included")
  }
  if(is.null(X)){
    stop("Variable data is not included")
  }
  deviance = sum(gammaresiduals$deviance^2)
  AIC <- deviance + (2*(dim(X)[2]))  
  BIC <- deviance + (log(dim(X)[1])*(dim(X)[2]))
  
  criteria <- list()
  
  criteria$deviance <- deviance
  criteria$AIC <- AIC
  criteria$BIC <- BIC
  
  return(criteria)
}
