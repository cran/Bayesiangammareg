\name{gammaresiduals}
\alias{gammaresiduals}
\title{Residuals of the Gamma Regression}
\description{This function calculates the Gamma regression residuals}
\usage{
gammaresiduals(Y, X, model)
}
\arguments{
  \item{Y}{ object of class matrix, with the dependent variable.}
  \item{X}{ object of class matrix, with the independent variable.}
  \item{model}{ object of class Bayesiangammareg.}
}
%%\details{}
\value{
  \item{rabs}{ Pearson absolute residuals}
  \item{rp}{ Pearson residuals}
  \item{rd}{ deviance residuals}
  \item{rast}{ Asteric residuals}
}
\references{1. Cepeda-Cuervo E. (2001) Modelagem da variabilidade em modelos lineares generalizados. Ph.D. tesis. Instituto de Matematicas. Universidade Federal do Rio do Janeiro.
2. Cepeda-Cuervo E. and Gamerman D. (2005). Bayesian Methodology for modeling parameters in the two-parameter exponential family. Estadistica 57, 93 105.
3. Cepeda Cuervo E., Corrales, M., Cifuentes, M. V., and Zarate, H. (2016). On Gamma Regression Residuals.}
\author{
Arturo Camargo Lozano   \email{bacamargol@unal.edu.co},
Edilberto Cepeda Cuervo \email{ecepedac@unal.edu.co}}
\keyword{Bayesian}
\keyword{Metropolis Hastings}
\keyword{Mean-Covariance modeling}
\keyword{Gamma}
