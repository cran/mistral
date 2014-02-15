\name{S2MART}
\alias{S2MART}
\title{Subset by Support vector Margin Algorithm for Reliability esTimation}

\description{
\code{S2MART} lets easily combine metamodel based reliability algorithm SMART and subset simulation algorithm. Basically, it introduces a metamodeling step at each subset simulation threshold, making number of necessary samples lower and the probability estimation better according to subset simulation by itself.
}

\usage{
S2MART(dimension,
      limit_state_function,
      Nn 			        = 100,
      alpha_quantile 	= 0.1,
      ...,
      plot			      = FALSE,
      limited_plot		= FALSE,
      output_dir	  	= NULL,
      verbose 			  = 0)
}

\arguments{
  \item{dimension}{an integer giving the dimension of the input space.}
  \item{limit_state_function}{the failure fonction.}
  \item{Nn}{the number of samples to evaluate the quantiles in the subset step.}
  \item{alpha_quantile}{cutoff probability for the subsets.}
	\item{...}{parameters of SMART algorithm.}
  \item{plot}{a boolean parameter specifying if function and samples should be plotted. The plot is refreshed at each iteration with the new data. Note that this option is only to be used when working on \sQuote{light} limit state functions as it requires the calculus of this function on a grid of size 161x161 (plot is done a -8:8 x -8:8 grid with 161 meshes.}
  \item{limited_plot}{only a final plot with \code{limit_state_function}, final DOE and metamodels. Should be used with \code{plot}==FALSE. As for \code{plot} it requires the calculus of the \code{limit_state_function} on a grid of size 161x161.}
  \item{output_dir}{optional. If plots are to be saved in .jpeg in a given directory. This variable will be pasted with "_Subset.jpeg" to get the full output directory.}
  \item{verbose}{Eiher 0 for an almost no output message, or 1 for medium size or 2 for full size}
}

\value{
  An object of class \code{list} containing the failure probability and some more outputs as described below:

  \item{proba}{The estimated failure probability.}
  \item{cov}{The coefficient of variation of the Monte-Carlo probability estimate.}
  \item{Ncall}{The total number of calls to the \code{limit_state_function}.}
  \item{learn_db}{The final learning database, ie. all points where \code{limit_state_function} has been calculated.}
  \item{lsf_value}{The value of the \code{limit_state_function} on the learning database.}
  \item{meta_model}{The final metamodel. An S4 object from \pkg{DiceKriging}.}
}

\details{
  S2MART algorithm is based on the idea that subset simulations conditional probabilities are estimated with a relatively poor precision as it requires calls to the expensive-to-evaluate limit state function and does not take benefit from its numerous calls to the limit state function in the Metropolis-Hastings algorithm. In this scope, the key concept is to reduce the subset simulation population to its minimum and use it only to estimate crudely the next quantile. Then the use of a metamodel-based algorithm lets refine the border and calculate an accurate estimation of the conditional probability by the mean of a crude Monte-Carlo.

  In this scope, a compromise has to be found between the two sources of calls to the limit state function as total number of calls = (\code{Nn} + number of calls to refine the metamodel) x (number of subsets) :
  \itemize{
     \item{\code{Nn} calls to find the next threshold value : the bigger \code{Nn}, the more accurate the \sQuote{decreasing speed} specified by the \code{alpha_quantile} value and so the smaller the number of subsets}
     \item{total number of calls to refine the metamodel at each threshold}
  }
}

\note{
Problem is supposed to be defined in the standard space. If not, use \code{\link{UtoX}} to do so.
}

\references{
  \itemize{
    \item
      J.-M. Bourinet, F. Deheeger, M. Lemaire:\cr
      \emph{Assessing small failure probabilities by combined Subset Simulation and Support Vector Machines}\cr
      Structural Safety (2011)

    \item
      F. Deheeger:\cr
      \emph{Couplage m?cano-fiabiliste : 2SMART - m?thodologie d'apprentissage stochastique en fiabilit?}\cr
      PhD. Thesis, Universit? Blaise Pascal - Clermont II, 2008

    \item
      S.-K. Au, J. L. Beck:\cr
      \emph{Estimation of small failure probabilities in high dimensions by Subset Simulation} \cr
      Probabilistic Engineering Mechanics (2001)

    \item
      A. Der Kiureghian, T. Dakessian:\cr
      \emph{Multiple design points in first and second-order reliability}\cr
      Structural Safety, vol.20 (1998)

    \item
      P.-H. Waarts:\cr
      \emph{Structural reliability using finite element methods: an appraisal of DARS:\cr Directional Adaptive Response Surface Sampling}\cr
      PhD. Thesis, Technical University of Delft, The Netherlands, 2000
  }
}

\author{
  Clement Walter\cr
  \email{clement.walter@cea.fr}
}

\seealso{
  \code{\link{SMART}}
  \code{\link{SubsetSimulation}}
  \code{\link{MonteCarlo}}
  \code{\link[DiceKriging]{km}} (in package \pkg{DiceKriging})
  \code{\link[e1071]{svm}} (in package \pkg{e1071})
}

\examples{
#Limit state function defined by Kiureghian & Dakessian :
kiureghian = function(x, b=5, kappa=0.5, e=0.1) {b - x[2] - kappa*(x[1]-e)^2}

\dontrun{
res = S2MART(dimension = 2,
            limit_state_function = kiureghian,
            N1 = 1000, N2 = 5000, N3 = 10000,
            plot = TRUE)

#Compare with crude Monte-Carlo reference value
reference = MonteCarlo(2, kiureghian, N_max = 500000)
}

#See impact of metamodel-based subset simulation with Waarts function :
waarts = function(u) { min(
		(3+(u[1]-u[2])^2/10 - (u[1]+u[2])/sqrt(2)),
		(3+(u[1]-u[2])^2/10 + (u[1]+u[2])/sqrt(2)),
		u[1]-u[2]+7/sqrt(2),
		u[2]-u[1]+7/sqrt(2))
	}

\dontrun{
res = list()
res$SMART = SMART(dimension  = 2, limit_state_function = waarts, plot=TRUE)
res$S2MART = S2MART(dimension = 2,
                    limit_state_function = waarts,
                    N1 = 1000, N2 = 5000, N3 = 10000,
                    plot=TRUE)
res$SS = SubsetSimulation(dimension = 2, waarts, n_init_samples = 10000)
res$MC = MonteCarlo(2, waarts, N_max = 500000)
}
}
