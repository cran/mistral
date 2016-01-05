\name{FORM}
\alias{FORM}
\title{ FORM method}

\description{
  Calculate failure probability by FORM method and important sampling.
}

\usage{
FORM(f, u.dep, inputDistribution, N.calls, eps = 1e-7,
     Method = "HLRF", IS = FALSE, q = 0.5, copula = "unif")
}

\arguments{
  \item{f}{
     A failure fonction
   }
  \item{u.dep}{
     A vector, starting point to the research of the design point  
   }
  \item{inputDistribution}{
     A list which contains the name of the input distribution and their parameters. For the input "i", inputDistribution[[i]] = list("name_law",c(parameters1,..., parametersN))
   }
  \item{N.calls}{
     Number of calls to f allowed
   }
  \item{eps}{
     Stop criterion : distance of two points between two iterations
   }
  \item{Method}{
     Choice of the method to research the design point: "AR" for Abdo-Rackwitz and "HLRF" for Hasofer-Lindt-Rackwitz-Fiessler
   }
  \item{IS}{
    "TRUE" for using importance Sampling method (applied after FORM which provides the importance density). Default = "FALSE".
   }
  \item{q}{
     Ratio of N.calls for the reserach of the design point by FORM. Default = 0.5. 1-q = the remaining ratio to use importance sampling. 
   }
  \item{copula}{
     Choice of the copula. Default = "unif" (uniform copula)
   }
}

\value{
  \item{proba.def}{Failure probability}
  \item{indice.reliab}{Reliability index}
  \item{compt.f}{Number of calls to f}
  \item{Design.Point}{Coordinates of the design point}
  \item{fact.imp}{Importance factors}
  \item{variance}{Standard error of the probability estimator (if IS = TRUE)}
  \item{Interval.conf}{Confidence interval of the estimator at 0.95 (if IS = TRUE)}
  \item{DOE}{List which contains the design of experiments}
}

\author{
  Vincent Moutoussamy
}

\details{
    This function estimate the probability that the output of the failure function is negative using FORM algorithm.
    The importance sampling procedure estimate a probability using a gaussian distribution centered in the design point with a covariance matrix equal to the indendity.
}

\references{
O. Ditlevsen and H.O. Madsen. Structural reliability methods, Wiley, 1996

M. Lemaire, A. Chateauneuf and J. Mitteau. Structural reliability, Wiley Online Library, 2009.
}

\examples{ 
\dontrun{
 distribution = list()
 distribution[[1]] = list("gamma",c(2,1))
 distribution[[2]] = list("gamma",c(3,1))
 
 f <- function(X){
    X[1]/sum(X) - qbeta((1e-5),2,3)
 }

 res.list <- FORM(f,
                  u.dep = c(0, 0.1),
                  inputDistribution = distribution,
                  N.calls = 1000,
                  eps = 1e-7,
                  Method = "HLRF",
                  IS = "TRUE",
                  q = 0.1,
                  copula = "unif")
}
}


