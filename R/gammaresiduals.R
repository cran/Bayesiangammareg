gammaresiduals <-
function(Y,X,model){
  Y <- as.matrix(Y)
  residuals <- model$residuals
  variance <- model$variance
  phi <- model$precision
  yestimado <-  model$fitted.values
  
  #Absolute residuals  
  rabs<-abs(residuals)
  
  #Standardized Weighted Residual 1
  rp<-residuals/sqrt(variance)
  
  # Res Deviance
  rd = -2*(log(Y/yestimado) - (Y-yestimado)/yestimado)
  
  #Residuals astesrisc
  rast= (log(Y) + log(phi/yestimado) - digamma(phi))/sqrt(trigamma(phi)) 
  
  gammaresiduals<- list()
  gammaresiduals$abs <- rabs
  gammaresiduals$pearson <-rp
  gammaresiduals$deviance <- rd
  gammaresiduals$rgamma<- rast
  
  return(gammaresiduals)
  
}
