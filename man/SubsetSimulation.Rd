\name{SubsetSimulation}
\alias{SubsetSimulation}
\title{
Subset Simulation Monte-Carlo
}
\description{
Estimate a failure probability with a Monte-Carlo method applied on nested subsets.
}
\usage{
SubsetSimulation(dimension,
                 limit_state_function,
                 proposal_pdf,
                 pdf            = dnorm,
                 rpdf           = rnorm,
                 cutoff_prob    = 0.1,
                 n_init_samples = 10000,
				         burnin 	     	= 20,
				         thinning 		  = 4,
                 plot           = FALSE,
                 output_dir     = NULL,
			        	 verbose 		    = 0)
}

\arguments{
  \item{dimension}{an integer giving the dimension of the input space.}
  \item{limit_state_function}{the failure function.}
  \item{proposal_pdf}{proposal PDF for Metropolis-Hastings algorithm. Default random walk is uniform on a radius of 2.}
  \item{pdf}{PDF of the input space to be used in Metropolis-Hastings algorithm.}
  \item{rpdf}{a random generator for the input space PDF.}
  \item{cutoff_prob}{the cut-off probability for each subset.}
  \item{n_init_samples}{number of samples to be used for each subset.}
  \item{burnin}{burnin parameter  for Metropolis-Hastings algorithm.}
  \item{thinning}{thinning parameter for Metropolis-Hastings algorithm.}
  \item{plot}{a boolean parameter specifying if function and samples should be plotted. The plot is refreshed at each iteration with the new data. Note that this option is only to be used when working on \sQuote{light} limit state functions as it requires the calculus of this function on a grid of size 161x161 (plot is done a -8:8 x -8:8 grid with 161 meshes.}
  \item{output_dir}{optional. If plots are to be saved in .jpeg in a given directory. This variable will be pasted with "_Subset_Simulation.jpeg" to get the full output directory.}
  \item{verbose}{Eiher 0 for an almost no output message, or 1 for medium size or 2 for full size}
}

\details{
This algorithm uses the property of conditional probabilities on nested subsets to calculate a given probability defined by a limit state function.

It operates iteratively on \sQuote{populations} to estimate the quantile corresponding to a probability of \code{cutoff_prob}. Then, it generates samples conditionnaly to this threshold, until found threshold be lower than 0.

Finally, the estimate is the product of the conditional probabilities.
}

\value{
  An object of class \code{list} containing the failure probability and some more outputs as described below:

  \item{proba}{The estimated failure probability.}
  \item{cov}{The coefficient of variation of the Monte-Carlo probability estimate.}
  \item{Ncall}{The total number of calls to the \code{limit_state_function}.}
}

\author{
  Clement Walter\cr
  \email{clement.walter@cea.fr}
}

\note{
Problem is supposed to be defined in the standard space. If not, use \code{\link{UtoX}} to do so.
Furthermore, each time a set of vector is defined as a matrix, \sQuote{nrow} = \code{dimension} and \sQuote{ncol} = number of vector.
}

\references{
  \itemize{
    \item
      S.-K. Au, J. L. Beck:\cr
      \emph{Estimation of small failure probabilities in high dimensions by Subset Simulation} \cr
      Probabilistic Engineering Mechanics (2001)
  }
}

\seealso{
  \code{\link{MonteCarlo}}
  \code{\link{S2MART}}
}

\examples{
#Try Subset Simulation Monte Carlo on a given function and change number of points.
#Limit state function defined by Kiureghian & Dakessian :
kiureghian = function(x, b=5, kappa=0.5, e=0.1) {b - x[2] - kappa*(x[1]-e)^2}

\dontrun{
res = list()
res[[1]] = SubsetSimulation(2,kiureghian,n_init_samples=10000)
res[[2]] = SubsetSimulation(2,kiureghian,n_init_samples=100000)
res[[3]] = SubsetSimulation(2,kiureghian,n_init_samples=500000)
}

#Try Subset Simulation Monte Carlo on a given function with different failure level
#Limit state function defined by Waarts :
waarts = function(u) { min(
    (3+(u[1]-u[2])^2/10 - (u[1]+u[2])/sqrt(2)),
		(3+(u[1]-u[2])^2/10 + (u[1]+u[2])/sqrt(2)),
		u[1]-u[2]+7/sqrt(2),
		u[2]-u[1]+7/sqrt(2))
}
\dontrun{
res = list()
res[[1]] = SubsetSimulation(2,waarts,failure=0,plot=TRUE)
res[[2]] = SubsetSimulation(2,waarts,failure=1,plot=TRUE)
res[[3]] = SubsetSimulation(2,waarts,failure=-1,plot=TRUE)
}

}