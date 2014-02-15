\name{MonteCarlo}
\alias{MonteCarlo}
\title{
Crude Monte Carlo method
}
\description{
Estimate a failure probability using crude Monte Carlo method
}
\usage{
MonteCarlo(dimension,
          LimitStateFunction,
          N_max         = 5e+05,
          N_batch       = 1000,
          failure       = 0,
          precision     = 0.05,
          plot          = FALSE,
          output_dir    = NULL,
          verbose       = 0)
}

\arguments{
  \item{dimension}{an integer giving the dimension of the input space.}
  \item{LimitStateFunction}{the failure function.}
  \item{N_max}{Maximum number of calls to the \code{LimiStateFunction}.}
  \item{N_batch}{Size of batch for Monte-Carlo population increase.}
  \item{failure}{the value defining the failure domain F = \{ x | \code{LimitStateFunction}(x) < \code{failure} \}.}
  \item{precision}{A target maximum value of the coefficient of variation}
  \item{plot}{a boolean parameter specifying if function and samples should be plotted. Note that this option is only to be used when working on \sQuote{light} limit state functions as it requires the calculus of this function on a grid of size 161x161 (plot is done a -8:8 x -8:8 grid with 161 meshes.}
  \item{output_dir}{optional. If plots are to be saved in .jpeg in a given directory. This variable will be pasted with "_Monte_Carlo_brut.jpeg" to get the full output directory.}
  \item{verbose}{Eiher 0 for an almost no output message, or 1 for medium size or 2 for full size}
}

\value{
  An object of class \code{list} containing the failure probability and some more outputs as described below:

  \item{proba}{The estimated failure probability.}
  \item{cov}{The coefficient of variation of the Monte-Carlo probability estimate.}
  \item{Ncall}{The total number of calls to the \code{LimitStateFunction}.}
}

\author{
  Clement Walter\cr
  \email{clement.walter@cea.fr}
}

\note{
Problem is supposed to be defined in the standard space. If not, use \code{\link{UtoX}} to do so.
Furthermore, each time a set of vector is defined as a matrix, \sQuote{nrow} = \code{dimension} and \sQuote{ncol} = number of vector.
}

\seealso{
  \code{\link{SubsetSimulation}}
}

\examples{
#Try Naive Monte Carlo on a given function with different failure level
#Limit state function defined by Kiureghian & Dakessian :
kiureghian = function(x, b=5, kappa=0.5, e=0.1) {b - x[2] - kappa*(x[1]-e)^2}

\dontrun{
res = list()
res[[1]] = MonteCarlo(2,kiureghian,failure = 0,plot=TRUE)
res[[2]] = MonteCarlo(2,kiureghian,failure = 1,plot=TRUE)
res[[3]] = MonteCarlo(2,kiureghian,failure = -1,plot=TRUE)

}


#Try Naive Monte Carlo on a given function and change number of points.
#Limit state function defined by Waarts :
waarts = function(u) { min(
  	(3+(u[1]-u[2])^2/10 - (u[1]+u[2])/sqrt(2)),
		(3+(u[1]-u[2])^2/10 + (u[1]+u[2])/sqrt(2)),
		u[1]-u[2]+7/sqrt(2),
		u[2]-u[1]+7/sqrt(2))
}
\dontrun{
res = list()
res[[1]] = MonteCarlo(2,waarts,N_max = 10000)
res[[2]] = MonteCarlo(2,waarts,N_max = 100000)
res[[3]] = MonteCarlo(2,waarts,N_max = 500000)
}

}
