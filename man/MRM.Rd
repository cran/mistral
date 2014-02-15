\name{MRM}
\alias{MRM}
\title{MRM method}
\description{
  Calculate failure probability by MRM method.
}
\usage{
MRM(f, ndim, choice.law, dir.monot, N.calls, Method, ordre.p = 0, silent = TRUE)
}
\arguments{
 \item{f}{
    the failure fonction
  }
  \item{ndim}{
    dimension of the inputs
  }
  \item{choice.law}{
    a list of length \sQuote{ndim} which contain name of input distribution and their parameters. For the input \sQuote{i}, choice.law[[i]] = list(\sQuote{name_law},c(parameters1,..., parametersN) )
  }
  \item{dir.monot}{
    vector of size \code{ndim} which represent the monotonicity of the failure function. \code{dir.monot[i]} = -1 (resp. 1) if the code is decreasing (resp. increasing) in direction i.
  }
  \item{N.calls}{
     Number of calls to f allowed
  }
  \item{Method}{
    there are two methods available. "MC_monotone"" is an adapation of the Monte Carlo method under constraints of monotony. "MRM" is a sequential sampling method.
  }
  \item{ordre.p}{
    order of magnitude of the search probability (for \code{ndim} >= 3)
  }
  \item{silent}{
    if silent = TRUE, print curent number of call to f. Default: FALSE.
  }
}

\details{
    These methods compute the probability that the output of the failure function is negative

   \code{ordre.p} is useful to compute the bounds in dimension >= 3, its computation is a santard Monte Carlo method. For dimension 2, there is an exat computation without sampling.

}

\value{
    \item{Um}{Exact lower bounds of the failure probability }
    \item{UM}{Exact upper bounds of the failure probability }
    \item{MLE}{Maximum likelihood estimator of the failure probability }
    \item{IC.inf}{Lower bound of the confidence interval of the failure probability based on MLE}
    \item{IC.sup}{Upper bound of the confidence interval of the failure probability based on MLE}
    \item{CV.MLE}{Coefficient of variation of the MLE}
    \item{N.tot}{Total number of simulation (just for "MC_monotone")}
}

\author{
  Vincent Moutoussamy and Nicolas Bousquet
}

\references{
    Bousquet, N. (2012) Accelerated monte carlo estimation of exceedance probabilities under monotonicity constraints. Annales de la Faculte des Sciences de Toulouse. XXI(3), 557-592.
}


\examples{
\dontrun{
 choice.law <- list()
 choice.law[[1]] <- list("norm",c(4,1))
 choice.law[[2]] <- list("norm",c(0,1))

 ndim <- 2
 dir.monot <- c(1, -1)
 N.calls <- 80

  f <- function(Input){
    return(Input[1] - Input[2])
  }
  
  res.MRM <- MRM(f, ndim, choice.law, dir.monot, N.calls, Method = "MRM",
                  ordre.p = 0, silent = FALSE)

  N <- 1:dim(res.MRM)[1]
  
  plot(N, res.MRM[, 1],
        col = "blue", lwd=2, type='l', ylim=c(0,1e-2),
        xlab="Number of call to the failure function",
        ylab="")
  lines(N, res.MRM[, 2], col = "blue", lwd = 2)
  lines(N, res.MRM[, 3], col = "red", lwd = 2)
  lines(N, res.MRM[, 4], col = "orange", lwd = 2, lty = 2)
  lines(N, res.MRM[, 5], col = "orange", lwd = 2, lty = 2)
  legend("topright",
          c("Exact Bounds", "MLE", "Confidence interval"), 
          col = c("blue", "red", "orange"),
          text.col = c("blue", "red", "orange"),
          lty = c(1, 1, 2),
          merge = TRUE)


  #It can be long:
  res.MC_monotone <- MRM(f, ndim, choice.law, dir.monot, N.calls, Method = "MC_monotone",
                          ordre.p = 0, silent = FALSE)
}
}