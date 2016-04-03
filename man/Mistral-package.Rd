\name{mistral-package}
\alias{mistral-package}
\alias{mistral}
\docType{package}
\title{Methods In Structural Reliability Analysis}

\description{
Provide tools for structural reliability analysis (failure probability, quantile).
}

\details{
 \tabular{lll}{
   Package: \tab mistral\cr
   Type: \tab Package\cr
   Version: \tab 2.1.0
   Date: \tab 2016-04-03\cr
   License: \tab CeCILL
  }
This package provides tools for structural reliability analysis:
\itemize{
    \item Calculate failure probability with FORM method and importance sampling,
    \item Calculate failure probability with crude Monte Carlo method,
    \item Calculate failure probability with Subset Simulation algorithm,
    \item Calculate failure probability with Monotonic Reliability Methods (MRM),
    \item Calculate failure probability with metamodel based algorithms : AKMCS, SMART and MetaIS,
    \item Calculate failure probability with a metamodel based Subset Simulation : S2MART,
    \item Wilks formula: Compute a quantile (or tolerance interval) with a given confidence level from a i.i.d. sample,
    \item Wilks formula: Compute the minimal sample size to estimate a quantile with a given confidence level,
    \item Calculate a quantile under monotonicity constraints.
}
}

\author{Clement Walter, Gilles Defaux, Bertrand Iooss, Vincent Moutoussamy, with contributions from Nicolas Bousquet, Claire Cannamela and Paul Lemaitre
(maintainer: Bertrand Iooss \email{biooss@yahoo.fr})}

\references{
O. Ditlevsen and H.O. Madsen. Structural reliability methods, Wiley, 1996.

M. Lemaire, A. Chateauneuf and J. Mitteau. Structural reliability, Wiley Online Library, 2009.

J. Morio and M. Balesdent. Estimation of rare event probabilities in complex aerospace and other systems, WP, 2016.

S.S. Wilks. Determination of Sample Sizes for Setting Tolerance Limits. Annals Mathematical Statistics, 12:91-96, 1941.

}

\examples{

########## FORM ###########

 distribution = list()
 distribution[[1]] = list("gamma",c(2,1))
 distribution[[2]] = list("gamma",c(3,1))

 f <- function(X){
   X[1]/sum(X) - qbeta((1e-5),2,3)
 }

 res <- FORM(f, u.dep = c(0,0.1), inputDist = distribution, 
     N.calls = 1000, eps = 1e-7, Method = "HLRF", IS = "TRUE", 
     q = 0.1, copula = "unif")

########### Wilks ##########

N <- WilksFormula(0.95,0.95,order=1)
print(N)

}
\keyword{ package }
