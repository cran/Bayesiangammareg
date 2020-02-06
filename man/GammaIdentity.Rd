\name{GammaIdentity}
\alias{GammaIdentity}
\title{Bayesian Gamma Regression with link Identity for the Model of Mean.}
\description{Function to do Bayesian Gamma Regression link Identity: Joint Mean and Shape
modeling with Identity link for Mean.}
\usage{
GammaIdentity(Y, X, Z, nsim, bpri, Bpri, gpri, Gpri, burn, jump, graph1, graph2)
}
\arguments{
  \item{Y}{Object of class matrix, with the dependent variable.}
  \item{X}{Object of class matrix, with the variables for modeling the mean.}
  \item{Z}{Object of class matrix, with the variables for modeling the shape.}
  \item{nsim}{a number that indicate the number of iterations.}
  \item{bpri}{a vector with the initial values of beta.}
  \item{Bpri}{a matrix with the initial values of the variance of beta.}
  \item{gpri}{a vector with the initial values of gamma.}
  \item{Gpri}{a matrix with the initial values of the variance of gamma.}
  \item{burn}{a proportion that indicate the number of iterations to be burn at
  the beginning of the chain.}
  \item{jump}{a number that indicate the distance between samples of 
  the autocorrelated the chain, to be excluded from the final chain.}
  \item{graph1}{if it is TRUE present the graph of the chains without jump and burn.}
  \item{graph2}{if it is TRUE present the graph of the chains with jump and burn.}
}
%%\details{}
\value{
object of class \code{bayesiangammareg} with the following:
  \item{Bestimado}{object of class matrix with the estimated coefficients of beta}
  \item{Gammaest}{object of class matrix with the estimated coefficients of gamma}
  \item{X}{object of class matrix, with the variables for modelling the mean}
  \item{Z}{object of class matrix, with the variables for modelling the precision}
  \item{DesvBeta}{object of class matrix with the estimated desviations of beta}
  \item{DesvGamma}{object of class matrix with the estimated desviations of gamma}
  \item{B}{object of class matrix with the B values}
  \item{G}{object of class matrix with the G values}
  \item{yestimado}{object of class matrix with the fitted values of y}
  \item{residuals}{object of class matrix with the residuals of the regression}
  \item{phi}{object of class matrix with the precision terms of the regression}
  \item{variance}{object of class matrix with the variance terms of the regression}
  \item{beta.mcmc}{object of class matrix with the complete chains for beta}
  \item{gamma.mcmc}{object of class matrix with the complete chains for gamma}
  \item{beta.mcmc.auto}{object of class matrix with the chains for beta after the burned process}
  \item{gamma.mcmc.auto}{object of class matrix with the chains for gamma after the burned process}
  }
\references{1. Cepeda-Cuervo E. (2001) Modelagem da variabilidade em modelos lineares generalizados. Ph.D. tesis. Instituto de Matematicas. Universidade Federal do Rio do Janeiro.
2. Cepeda-Cuervo E. and Gamerman D. (2005). Bayesian Methodology for modeling parameters in the two-parameter exponential family. Estadistica 57, 93 105.
3. Cepeda Cuervo E. and Gamerman D. (2001). Bayesian Modeling of Variance Heterogeneity in Normal Regression Models. Brazilian Journal of Probability and Statistics. 14, 207-221.
}
\author{
Arturo Camargo Lozano   \email{bacamargol@unal.edu.co},
Edilberto Cepeda-Cuervo \email{ecepedac@unal.edu.co} }
%%\note{}
%%\seealso{}
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