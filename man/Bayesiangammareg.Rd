\name{Bayesiangammareg}
\alias{Bayesiangammareg}
\title{Bayesian Gamma Regression: Joint Mean and Shape Modeling}
\description{ Function to do Bayesian Gamma Regression: Joint
Mean and Shape Modeling}
\usage{
Bayesiangammareg(Y, X, Z, nsim, bpri, Bpri, gpri, Gpri, burn, jump,
graph1, graph2, meanlink = "log")
}
\arguments{
  \item{Y}{object of class matrix, with the dependent variable.}
  \item{X}{object of class matrix, with the variables for modeling the mean.}
  \item{Z}{object of class matrix, with the variables for modeling the shape.}
  \item{nsim}{a number that indicate the number of iterations.}
  \item{bpri}{a vector with the initial values of beta.}
  \item{Bpri}{a matrix with the initial values of the variance of beta.}
  \item{gpri}{a vector with the initial values of gamma.}
  \item{Gpri}{a matrix with the initial values of the variance of gamma.}
  \item{burn}{a proportion that indicate the number of iterations to be burn 
  at the beginning of the chain.}
  \item{jump}{a number that indicate the distance between samples of the 
  autocorrelated the
chain, to be excluded from the final chain.}
  \item{graph1}{if it is TRUE present the graph of the chains without jump and burn.}
  \item{graph2}{if it is TRUE present the graph of the chains with jump and burn.}
  \item{meanlink}{represent the link function, logarithm or identity.}
}
\details{The Bayesian Gamma regression allows the joint modeling of the mean 
and the shape of a gamma distributed variable, using a Bayesian estimation 
algorithm proposed by Cepeda-Cuervo (2001).
}
\value{object of class bayesiangammareg with:
  \item{coefficients}{object of class matrix with the estimated coefficients of beta and gamma.}
  \item{desv}{object of class matrix with the estimated desviations of beta and gamma.}
  \item{interv}{object of class matrix with the estimated confidence intervals of
  beta and gamma.}
  \item{fitted.values}{object of class matrix with the fitted values of y.}
  \item{residuals}{object of class matrix with the residuals of the regression.}
  \item{beta.mcmc}{object of class matrix with the complete chains for beta.}
  \item{gamma.mcmc}{object of class matrix with the complete chains for gamma.}
  \item{beta.mcmc.short}{object of class matrix with the chains for beta after the burned process.}
  \item{gamma.mcmc.short}{object of class matrix with the chains for gamma after the burned process.}
  \item{call}{Call.}
}
\references{1. Cepeda-Cuervo E. (2001) Modelagem da variabilidade em modelos lineares generalizados. Ph.D. tesis. Instituto de Matematicas. Universidade Federal do Rio do Janeiro.
2. Cepeda-Cuervo E. and Gamerman D. (2005). Bayesian Methodology for modeling parameters in the two-parameter exponential family. Estadistica 57, 93 105.
}
\author{
Arturo Camargo Lozano   \email{bacamargol@unal.edu.co},
Edilberto Cepeda-Cuervo \email{ecepedac@unal.edu.co}
}
\examples{
X1 <- rep(1,50)
X2 <- runif(50,0,30)
X3 <- runif(50,0,20)
X4 <- runif(50,10,20)
mui <- 15 + 3*X2 + 2*X3
alphai <- exp(3 + 0.15*X2 + 0.15*X4)
Y <- rgamma(50,shape=alphai,scale=mui/alphai)
X <- cbind(X1,X2,X3)
Z <- cbind(X1,X2,X4)
bpri <- c(1,1,1)
Bpri <- diag(10^(3),nrow=ncol(X),ncol=ncol(X))
gpri <- c(0,0,0)
Gpri <- diag(10^(3),nrow=ncol(Z),ncol=ncol(Z))
burn <- 0
jump <- 1
nsim <- 300
graph1=FALSE
graph2=FALSE
Bayesiangammareg(Y,X,Z,nsim,bpri,Bpri,gpri,Gpri,burn,jump,graph1,graph2,"ide")
}
\keyword{Bayesian}
\keyword{Metropolis Hastings}
\keyword{Mean-Covariance modeling}
\keyword{Gamma}
