% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/S2MART.R
\name{S2MART}
\alias{S2MART}
\title{Subset by Support vector Margin Algorithm for Reliability esTimation}
\usage{
S2MART(dimension, lsf, Nn = 100, alpha_quantile = 0.1, failure = 0, ...,
  plot = FALSE, output_dir = NULL, verbose = 0)
}
\arguments{
\item{dimension}{the dimension of the input space}

\item{lsf}{the function defining the failure domain. Failure is lsf(X) < \code{failure}}

\item{Nn}{number of samples to evaluate the quantiles in the subset step}

\item{alpha_quantile}{cutoff probability for the subsets}

\item{failure}{the failure threshold}

\item{...}{All others parameters of the metamodel based algorithm}

\item{plot}{to produce a plot of the failure and safety domain. Note that this requires a lot of
calls to the \code{lsf} and is thus only for training purpose}

\item{output_dir}{to save the plot into the given directory. This will be pasted with "_S2MART.pdf"}

\item{verbose}{either 0 for almost no output, 1 for medium size output and 2 for all outputs}
}
\value{
An object of class \code{list} containing the failure probability
and some more outputs as described below:
\item{p}{The estimated failure probability.}
\item{cov}{The coefficient of variation of the Monte-Carlo probability
estimate.}
\item{Ncall}{The total number of calls to the \code{lsf}.}
\item{learn_db}{The final learning database, ie. all points where \code{lsf}
has been calculated.}
\item{lsf_value}{The value of the \code{lsf} on the learning database.}
\item{meta_model}{The final metamodel. An object from \pkg{e1071}.}
}
\description{
\code{S2MART} introduces a metamodeling step at each subset simulation
threshold, making number of necessary samples lower and the probability estimation
better according to subset simulation by itself.
}
\details{
S2MART algorithm is based on the idea that subset simulations conditional
  probabilities are estimated with a relatively poor precision as it
  requires calls to the expensive-to-evaluate limit state function and
  does not take benefit from its numerous calls to the limit state function
  in the Metropolis-Hastings algorithm. In this scope, the key concept is
  to reduce the subset simulation population to its minimum and use it only
  to estimate crudely the next quantile. Then the use of a metamodel-based
  algorithm lets refine the border and calculate an accurate estimation of
  the conditional probability by the mean of a crude Monte-Carlo.

In this scope, a compromise has to be found between the two sources of
calls to the limit state function as total number of calls = (\code{Nn} +
number of calls to refine the metamodel) x (number of subsets) :
 \itemize{
   \item{\code{Nn} calls to find the next threshold value : the bigger \code{Nn},
   the more accurate the \sQuote{decreasing speed} specified by the
   \code{alpha_quantile} value and so the smaller the number of subsets}
   \item{total number of calls to refine the metamodel at each threshold}
 }
}
\note{
Problem is supposed to be defined in the standard space. If not,
use \code{\link{UtoX}} to do so. Furthermore, each time a set of vector
is defined as a matrix, \sQuote{nrow} = \code{dimension} and
\sQuote{ncol} = number of vector to be consistent with \code{as.matrix}
transformation of a vector.

Algorithm calls lsf(X) (where X is a matrix as defined previously) and
expects a vector in return. This allows the user to optimise the computation
of a batch of points, either by vectorial computation, or by the use of
external codes (optimised C or C++ codes for example) and/or parallel
computation; see examples in \link{MonteCarlo}.
}
\examples{
\dontrun{
  res = S2MART(dimension = 2,
               lsf = kiureghian,
               N1 = 1000, N2 = 5000, N3 = 10000,
               plot = TRUE)
  
  #Compare with crude Monte-Carlo reference value
  reference = MonteCarlo(2, kiureghian, N_max = 500000)
}

#See impact of metamodel-based subset simulation with Waarts function :
\dontrun{
  res = list()
  # SMART stands for the pure metamodel based algorithm targeting directly the
  # failure domain. This is not recommended by its authors which for this purpose
  # designed S2MART : Subset-SMART
  res$SMART = mistral:::SMART(dimension  = 2, lsf = waarts, plot=TRUE)
  res$S2MART = S2MART(dimension = 2,
                      lsf = waarts,
                      N1 = 1000, N2 = 5000, N3 = 10000,
                      plot=TRUE)
  res$SS = SubsetSimulation(dimension = 2, waarts, n_init_samples = 10000)
 res$MC = MonteCarlo(2, waarts, N_max = 500000)
}

}
\author{
Clement WALTER \email{clement.walter@cea.fr}
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
\seealso{
\code{\link{SMART}}
\code{\link{SubsetSimulation}}
\code{\link{MonteCarlo}}
\code{\link[DiceKriging]{km}} (in package \pkg{DiceKriging})
\code{\link[e1071]{svm}} (in package \pkg{e1071})
}

