\name{mistral-package}
\alias{mistral-package}
\alias{mistral}
\docType{package}
\title{Methods in structural reliability}

\description{
Provide tools for structural reliability analysis (failure probability, quantile).
}

\details{
 \tabular{lll}{
   Package: \tab mistral\cr
   Type: \tab Package\cr
   Version: \tab 1.0-1
   Date: \tab 2013-03-15\cr
   License: \tab CeCILL
  }
This package provides tools for structural reliability analysis:
\itemize{
    \item Calculate failure probability by FORM method and important sampling, compute the design point, give the importance factors.
    \item Compute Wilks formula (minimal size of a sample to estimate quantile with a confidence level).
}
}

\author{Vincent Moutoussamy, Nicolas Bousquet, Bertrand Iooss
(maintainer: Bertrand Iooss \email{biooss@yahoo.fr})}

\references{
H.A. David and H.N. Nagaraja. Order statistics, Wiley, 2003.

O. Ditlevsen and H.O. Madsen. Structural reliability methods, Wiley, 1996

M. Lemaire, A. Chateauneuf and J. Mitteau. Structural reliability, Wiley Online Library, 2009.

W.T. Nutt and G.B. Wallis. Evaluation of nuclear safety from the outputs of computer codes in the presence of uncertainties. Reliability Engineering and System Safety, 83:57-77, 2004.

S.S. Wilks. Determination of Sample Sizes for Setting Tolerance Limits. Annals Mathematical Statistics, 12:91-96, 1941.
}

\examples{

########## FORM ###########

 distribution = list()
 distribution[[1]] = list("gamma",c(2,1))
 distribution[[2]] = list("gamma",c(3,1))

 func <- function(X){
   X[1]/sum(X) - qbeta((1e-5),2,3)
 }

 res.list <- FORM(func, u.dep = c(0,0.1), choice.law = distribution, N.calls = 1000, eps = 1e-7, Method = "HLRF", IS = "TRUE", q = 0.1, copula = "unif")

########### Wilks ##########

N <- WilksFormula(0.95,0.95,order=1)
print(N)

}
\keyword{ package }
