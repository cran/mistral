\name{FORM}
\alias{FORM}
\title{ FORM method}

\description{
  Calculate failure probability by FORM method and important sampling.
}

\usage{
FORM(f, u.dep, choice.law, N.calls = 1e4, eps = 1e-7,
     Method = "HLRF", IS = FALSE, q = 0.5, copula = "unif")
}

\arguments{
  \item{f}{
     the failure fonction
   }
  \item{u.dep}{
     a vector, starting point to the research of the design point  
   }
  \item{choice.law}{
     a list which contains the name of the densityof each input and their parameters
   }
  \item{N.calls}{
     Number of calls to f allowed
   }
  \item{eps}{
     stop criterion : distance of two points between two iterations
   }
  \item{Method}{
     Choice of the method to the research of the design point: "AR" for Abdo-Rackwitz and "HLRF" for Hasofer-Lindt-Rackwitz-Fiessler
   }
  \item{IS}{
    "TRUE" for using importance Sampling method (applied after FORM which gives the importance density). Default ="FALSE".
   }
  \item{q}{
     ratio of N.calls for the reserach of the design point by FORM. Default = 0.5. 1-q = the remaining ratio to use importance sampling. 
   }
  \item{copula}{
     choice of the copula. Default = "unif" (uniform copula)
   }
}

\value{
  \item{proba.def}{Failure probability}
  \item{indice.reliab}{reliability index}
  \item{compt.f}{Number of calls of f}
  \item{Design.Point}{Coordinates of the design point}
  \item{fact.imp}{Importance factors}
  \item{variance}{Standard error of the probability estimator (if IS = TRUE)}
  \item{Interval.conf}{Confidence interval of the estimator at 0.95 (if IS = TRUE)}
  \item{DOE}{List which contains the design of experiments}
}

\author{
  Vincent Moutoussamy
}

\examples{ 
 distribution = list()
 distribution[[1]] = list("gamma",c(2,1))
 distribution[[2]] = list("gamma",c(3,1))

 func <- function(X){
    X[1]/sum(X) - qbeta((1e-5),2,3)
 }

 res.list <- FORM(func,
                  u.dep = c(0,0.1),
                  choice.law = distribution,
                  N.calls = 1000,
                  eps = 1e-7,
                  Method = "HLRF",
                  IS = "TRUE",
                  q = 0.1,
                  copula = "unif")
}

\references{
O. Ditlevsen and H.O. Madsen. Structural reliability methods, Wiley, 1996

M. Lemaire, A. Chateauneuf and J. Mitteau. Structural reliability, Wiley Online Library, 2009.
}
