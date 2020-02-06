summary.Bayesiangammareg <-
function(object, ...){
  
  se <- object$desv
  
  TAB <- cbind( Coefficient = coef(object),
                Desv. = se,
                L.CredIntv = object$interv[,1],
                U.CredIntv = object$interv[,2] 
  )
  
  colnames(TAB) <- c("Estimate", "Est.Error", "L.CredIntv",  "U.CredIntv")
  
  gamma.resid <- gammaresiduals(object$Y,object$X,object)
  criteria <- criteria(object$X,gamma.resid)
  
  res <- list(call=object$call, coefficients=TAB, deviance=criteria$deviance, AIC=criteria$AIC, BIC=criteria$BIC)  
  
  class(res) <- "summary.Bayesiangammareg"
  res  
}
