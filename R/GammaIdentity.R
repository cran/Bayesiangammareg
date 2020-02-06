GammaIdentity <-
function(Y,X,Z,nsim,bpri,Bpri,gpri,Gpri,burn,jump,graph1,graph2){
  
  dmvnorm1<-function(p,x,mu,sigma){
    deno<-(2*pi)^(p/2)*abs(det(sigma))^(1/2)
    num<-exp((-1/2)*(t(x-mu)%*%solve(sigma)%*%(x-mu)))
    densidad<-num/deno
    return(densidad)
  }
  
  muproposal <- function(Y,X,Z,betas,gammas,bpri,Bpri){
    eta1 <- X%*%betas  
    mu  <- eta1        
    eta2 <- Z%*%gammas 
    phi <- exp(eta2)   
    Y.mono <- Y           
    sigmab <- (mu^2)/phi  
    
    if (det(diag(as.vector(sigmab)))==Inf | det(diag(as.vector(sigmab)))==0){
      B.pos <- (solve(solve(Bpri) + t(X)%*%X))
      b.pos <- B.pos%*%(solve(Bpri)%*%bpri + t(X)%*%Y.mono)
    }
    else{
      B.pos <- (solve(solve(Bpri)+ t(X)%*%solve(diag(as.vector(sigmab)))%*%X))
      b.pos <- B.pos%*%(solve(Bpri)%*%bpri + t(X)%*%solve(diag(as.vector(sigmab)))%*%Y.mono)
    }
    value1=drop(rmvnorm(1,b.pos,B.pos))
    return(t(t(value1)))
  }
  
  ##propuesta gamma normal##  
  gammaproposal <- function(Y,X,Z,betas,gammas,gpri,Gpri){
    eta1 <- X%*%betas
    mu <- eta1 
    eta2 <- Z%*%gammas
    phi <- exp(eta2)
    sigma <- (1/(phi^2))*((trigamma(phi)-(1/phi))^(-1))
    Y.tilde <- (eta2 - (1/phi)*((trigamma(phi) - (1/phi))^(-1))*(digamma(phi) - log(phi)-log(Y)+log(mu)-1+ Y/mu)) 
    
    if (det(diag(as.vector(sigma)))==Inf | det(diag(as.vector(sigma)))==0){
      G.pos <- (solve(solve(Gpri) + t(Z)%*%Z))
      g.pos <- G.pos%*%(solve(Gpri)%*%gpri + t(Z)%*%Y.tilde)
    }
    else{
      
      G.pos <- solve(solve(Gpri)+ t(Z)%*%solve(diag(as.vector(sigma)))%*%Z)
      g.pos <- G.pos%*%(solve(Gpri)%*%gpri + t(Z)%*%solve(diag(as.vector(sigma)))%*%Y.tilde)
    }
    value2=drop(rmvnorm(1,g.pos,G.pos))
    return(t(t(value2)))
  }
  
  mukernel <- function(X,Z,Y,betas.n,betas.v,gammas.v,bpri,Bpri){
    eta1 <- X%*%betas.v
    mu <- eta1 
    eta2 <- Z%*%gammas.v
    phi <- exp(eta2)
    Y.mono <- Y
    sigmab <- (mu^2)/phi
    if (det(diag(as.vector(sigmab)))==Inf | det(diag(as.vector(sigmab)))==0){
      B.pos <- (solve(solve(Bpri) + t(X)%*%X))
      b.pos <- B.pos%*%(solve(Bpri)%*%bpri + t(X)%*%Y.mono)
    }
    else{
      B.pos <- (solve(solve(Bpri)+ t(X)%*%solve(diag(as.vector(sigmab)))%*%X))
      b.pos <- B.pos%*%(solve(Bpri)%*%bpri + t(X)%*%solve(diag(as.vector(sigmab)))%*%Y.mono)
    }
    value3=drop(dmvnorm(t(betas.n),b.pos,B.pos))
    return(t(t(value3)))
  }  
  
  gammakernel <- function(X,Z,Y,gammas.n,betas.v,gammas.v,gpri,Gpri){
    eta1 <- X%*%betas.v
    mu <- eta1 
    eta2 <- Z%*%gammas.v
    phi <- exp(eta2)
    sigma <- (1/(phi^2))*((trigamma(phi)-(1/phi))^(-1))
    Y.tilde <- (eta2 - (1/phi)*((trigamma(phi) - (1/phi))^(-1))*(digamma(phi) - log(phi)-log(Y)+log(mu)-1+ Y/mu)) 
    if (det(diag(as.vector(sigma)))==Inf | det(diag(as.vector(sigma)))==0){
      G.pos <- (solve(solve(Gpri) + t(Z)%*%Z))
      g.pos <- G.pos%*%(solve(Gpri)%*%gpri + t(Z)%*%Y.tilde)
    }
    else{
      
      G.pos <- solve(solve(Gpri)+ t(Z)%*%solve(diag(as.vector(sigma)))%*%Z)
      g.pos <- G.pos%*%(solve(Gpri)%*%gpri + t(Z)%*%solve(diag(as.vector(sigma)))%*%Y.tilde)
    }
    value4=drop(dmvnorm(t(gammas.n),g.pos,G.pos))
    return(t(t(value4)))
  }
  
  
  dpostb<-function(X,Z,Y,betas,gammas,bpri,Bpri){
    eta1<-X%*%betas
    mu<-eta1
    eta2<-Z%*%gammas
    alfa<-exp(eta2)
    L<-prod(dgamma(Y,shape =alfa ,scale =mu/alfa ))##verosimilitud
    P<-dmvnorm(t(betas),bpri,Bpri)  #priori
    value<-L*P
    return(value)
  }
  
  dpostg<-function(X,Z,Y,betas,gammas,gpri,Gpri){
    eta1<-X%*%betas
    mu<-eta1
    eta2<-Z%*%gammas
    alfa<-exp(eta2)
    L<-prod(dgamma(Y,shape=alfa,scale=mu/alfa)) ##verosimilitud
    P<-dmvnorm(t(gammas),gpri,Gpri)
    value<-P*L
    return(value)
  }
  
  ### Cepeda - Metropolis - Hastings
  Y = as.matrix(Y)
  if (is.null(X) | is.null(Z) | is.null(Y)) {
    stop("There is no data")
  }
  if (burn > 1 | burn < 0) {
    stop("Burn must be a proportion between 0 and 1")
  }
  if (nsim <= 0) {
    stop("the number of simulations must be greater than 0")
  }
  if (jump < 0 | jump > nsim) {
    stop("Jumper must be a positive number lesser than nsim")
  }
  ind1<-rep(0,nsim)
  ind2<-rep(0,nsim)
  betas.ind <- matrix(bpri,nrow=ncol(X))
  gammas.ind <- matrix(gpri,nrow=ncol(Z))
  
  beta.mcmc<-matrix(NA,nrow=nsim,ncol=ncol(X))
  gamma.mcmc<-matrix(NA,nrow=nsim,ncol=ncol(Z))  
  
  for(i in 1:nsim) {
    #Betas
    betas.sim <- matrix(muproposal(Y, X, Z, betas.ind, gammas.ind, bpri, Bpri), nrow = ncol(X))
    
    q1.mu <- mukernel(X, Z, Y, betas.sim, betas.ind, gammas.ind,   bpri, Bpri)
    q2.mu <- mukernel(X, Z, Y, betas.ind, betas.sim, gammas.ind,  bpri, Bpri)
    p1.mu <- dpostb(X, Z, Y, betas.sim, gammas.ind, bpri, Bpri)
    p2.mu <- dpostb(X, Z, Y, betas.ind, gammas.ind, bpri,  Bpri)
    
    
    ###gammas###
    gammas.sim <- matrix(gammaproposal(Y, X, Z, betas.sim,  gammas.ind, gpri, Gpri), nrow = ncol(Z))
    q1.gamma <- gammakernel(X, Z, Y, gammas.sim, betas.sim, gammas.ind, gpri, Gpri)
    q2.gamma <- gammakernel(X, Z, Y, gammas.ind, betas.sim, gammas.sim, gpri, Gpri)
    p1.gamma <- dpostg(X, Z, Y, betas.sim, gammas.sim, gpri,   Gpri)
    p2.gamma <- dpostg(X, Z, Y, betas.sim, gammas.ind, gpri,   Gpri)
    
    if(p1.mu==0 |q2.mu == 0){
      alfa1=10
    }
    else
    {
      alfa1<-((p1.mu/p2.mu)*(q1.mu/q2.mu))
    }
    
    Mu.val<-min(1,alfa1)
    u<-runif(1)
    if (u <=Mu.val) {
      betas.ind <- betas.sim
      ind1[i] = 1
    }
    
    beta.mcmc[i,]<-betas.ind
    beta.mcmc <- as.ts(beta.mcmc)
    
    if(p1.gamma==0 |q2.gamma == 0){
      alfa2=10
    }
    else
    {
      alfa2<-((p1.gamma/p2.gamma)*(q1.gamma/q2.gamma))
    }
    
    Gamma.val<-min(1,alfa2)
    u<-runif(1)
    if (u <=Gamma.val) {
      gammas.ind <- gammas.sim
      ind2[i] = 1
    }
    gamma.mcmc[i,]<-gammas.ind
    gamma.mcmc <- as.ts(gamma.mcmc)
    
##    if (i%%1000 == 0)##
  ##    cat("Burn-in iteration : ", i, "\n")##
  }
  ###extraccion###
  tburn <- nsim*burn
  extr <- seq(0,(nsim-tburn),jump)
  
  betas.burn <-as.matrix(beta.mcmc[(tburn+1):nrow(beta.mcmc),])
  gammas.burn <-as.matrix(gamma.mcmc[(tburn+1):nrow(gamma.mcmc),])
  
  beta.mcmc.auto <- as.matrix(betas.burn[extr,])
  beta.mcmc.auto <- as.ts(beta.mcmc.auto)
  gamma.mcmc.auto <- as.matrix(gammas.burn[extr,])
  gamma.mcmc.auto <- as.ts(gamma.mcmc.auto)
  
  #Beta y Gamma estimations
  Bestimado <- colMeans(beta.mcmc.auto)
  Gammaest <- colMeans(gamma.mcmc.auto)
  
  #estandar errors of beta and gamma
  DesvBeta <- matrix(apply(beta.mcmc.auto,2,sd))
  DesvGamma <- matrix(apply(gamma.mcmc.auto,2,sd))
  
  #Forma
  phi <- exp(Z%*%Gammaest)
  
  #estimate values of the dependent variable
  yestimado = X%*%Bestimado
  
  #estimate variance of the dependent variable
  variance<- (yestimado^2)/phi
  
  #Residuals 
  residuals = as.matrix(Y) - yestimado
  
  B1 <- matrix(0, ncol(X),1)
  B2 <- matrix(0, ncol(X),1)
  
  
  #Credibility intervals for beta
  for(i in 1:ncol(X)){
    B1[i,]<-quantile(beta.mcmc.auto[,i],0.025)
    B2[i,]<-quantile(beta.mcmc.auto[,i],0.975)
    B <- cbind(B1,B2)
  }
  
  # Credibility intervals for gamma
  
  G1 <- matrix(0, ncol(Z),1)
  G2 <- matrix(0, ncol(Z),1)
  
  for(i in 1:ncol(Z)){
    G1[i,]<-quantile(gamma.mcmc.auto[,i],0.025)
    G2[i,]<-quantile(gamma.mcmc.auto[,i],0.975)
    G <- cbind(G1,G2)
  }
  
  #Absolute residuals  
  rabs<-abs(residuals)
  #Standardized Weighted Residual 1
  rp<-residuals/sqrt(variance)
  # Res Deviance
  rd = -2*(log(Y/yestimado) - (Y-yestimado)/yestimado)
  #Residuals astesrisc
  rast= (log(Y) + log(phi/yestimado) - digamma(phi))/sqrt(trigamma(phi)) 
  
  deviance = sum(rd^2)
  AIC <- deviance + (2*(dim(X)[2]))  
  BIC <- deviance+(log(dim(X)[1])*(dim(X)[2]))
  
  
  if (graph1==TRUE) {
    
    for(i in 1:ncol(X)){
      dev.new()
      ts.plot(beta.mcmc[,i], main=paste("Complete chain for beta",i), xlab="number of iterations", ylab=paste("parameter beta",i))
    }
    
    for(i in 1:ncol(Z)){
      dev.new()
      ts.plot(gamma.mcmc[,i], main=paste("Complete chain for gamma",i), xlab="number of iterations", ylab=paste("parameter gamma",i))
      
    }
    
  } else{
  }
  
  
  if (graph2==TRUE) {
    
    for(i in 1:ncol(X)){
      dev.new()
      ts.plot(beta.mcmc.auto[,i], main=paste("Burn chain for beta",i), xlab="number of iterations", ylab=paste("parameter beta",i))
      
    }
    
    for(i in 1:ncol(Z)){
      dev.new()
      ts.plot(gamma.mcmc.auto[,i], main=paste("Burn chain for gamma",i), xlab="number of iterations", ylab=paste("parameter gamma",i))
      
    }
    
  } else{
  }
  

return(list(Bestimado=Bestimado,Gammaest=Gammaest,DesvBeta=DesvBeta, DesvGamma=DesvGamma, yestimado=yestimado,  phi=phi, 
beta.mcmc=beta.mcmc, gamma.mcmc=gamma.mcmc, beta.mcmc.auto=beta.mcmc.auto, gamma.mcmc.auto=gamma.mcmc.auto, 
variance=variance,residuals=residuals,B=B,G=G,rabs=rabs, rp=rp,rd=rd,  rast=rast, deviance=deviance, 
AIC=AIC, BIC=BIC, Y = Y,X=X,Z=Z))
}
