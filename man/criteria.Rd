\name{criteria}
\alias{criteria}
\title{Criteria for Comparison the Bayesian Gamma Regression.}
\description{Performs the comparison criterias for the Bayesian Gamma regression}
\usage{
criteria(X, gammaresiduals)
}
\arguments{
  \item{X}{object of class matrix, with the independent variable for the mean.}
  \item{gammaresiduals}{object of class bayesiangammareg, with the residuals of the 
  Bayesian Gamma regression, that can be calculated by the function gammaresiduals}
}
\details{This function calculate the residuals of a Bayesian Gamma regression.}
\value{
\item{deviance}{the deviance criteria}
\item{AIC}{the AIC criteria}
\item{BIC}{the BIC criteria}
}
\references{1. Cepeda-Cuervo E. (2001) Modelagem da variabilidade em modelos lineares generalizados. Ph.D. tesis. Instituto de Matematicas. Universidade Federal do Rio do Janeiro.
2. Cepeda-Cuervo E. and Gamerman D. (2005). Bayesian Methodology for modeling parameters in the two-parameter exponential family. Estadistica 57, 93 105.
3. Espinheira, P. L., Ferrari, S. L., and Cribari-Neto, F. On beta regression residuals. Journal of Applied Statistics 
4. Cepeda-Cuervo E., Corrales, M., Cifuentes, M. V., and Zarate, H. (2016). On Gamma Regression Residuals.
}
\author{
Arturo Camargo Lozano   \email{bacamargol@unal.edu.co},
Edilberto Cepeda-Cuervo \email{ecepedac@unal.edu.co}}
\keyword{Bayesian}
\keyword{Metropolis Hastings}
\keyword{Mean-Covariance modeling}
\keyword{Gamma}
