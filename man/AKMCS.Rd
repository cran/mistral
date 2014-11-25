\name{AKMCS}
\alias{AKMCS}
\title{Active learning reliability method combining Kriging and Monte Carlo Simulation}

\description{
Estimate a failure probability with the AKMCS method.
}

\usage{
AKMCS(dimension,
      limit_state_function,
      N                 = 500000,
      N1                = 10*dimension,
      Nmax              = 200,
      learn_db          = NULL,
      lsf_value         = NULL,
      failure           = 0.0,
      precision         = 0.05,
      meta_model        = NULL,
      kernel            = "matern5_2",
      learn_each_train  = FALSE,
      crit_min          = 2,
      limit_fun_MH      = NULL,
      sampling_strategy = "MH",
      first_DOE         = "Gaussian",
      seeds             = NULL,
      seeds_eval        = NULL,
      burnin            = 30,
      thinning          = 4,
      plot              = FALSE,
      limited_plot      = FALSE,
      add               = FALSE,
      output_dir        = NULL,
      z_MH              = NULL,
      z_lsf             = NULL,
      verbose 		      = 0)
}

\arguments{
  \item{dimension}{an integer giving the dimension of the input space.}
  \item{limit_state_function}{the failure function.}
  \item{N}{an integer defining the Monte-Carlo population size for probability estimation.}
  \item{N1}{an integer defining the size of the first Design Of Experiment got by clustering of the \code{N} standard gaussian samples.}
  \item{Nmax}{an integer defining the maximum number of calls to the limit state function during refinement steps. This means total number of call will be \code{N1} + \code{Nmax}.}
  \item{learn_db}{optional. A matrix of already known points, with dim : \code{dimension} x number_of_vector.}
  \item{lsf_value}{values of the \code{limit_state_function} on the vectors given in \code{learn_db}.}
  \item{failure}{the value defining the failure domain F = \{ x | \code{limit_state_function}(x) < \code{failure} \}.}
  \item{precision}{the maximum value of the coefficient of variation for \code{proba}. If the first run with \code{N} gives a too large \code{cov}, then approximate necessary \code{N} is derived from \code{cov} and Monte-Carlo estimate is run again.}
  \item{meta_model}{optional. If a kriging based metamodel has already been fitted to the data (from \pkg{DiceKriging} package) it can be given as an input to keep the same parameters.}
  \item{kernel}{a specified kernel to be used for the metamodel. See \pkg{DiceKriging} for available options.}
  \item{learn_each_train}{specify if hyperparameters of the model should be evaluated each time points are added to the learning database ("TRUE") or only the first time ("FALSE").}
  \item{crit_min}{the minimum value of the criterion to be used for refinement step.}
  \item{limit_fun_MH}{optional. If the working space is to be reduced to some subset defining by a function, eg. in case of use in a Subset Simulation algorithm. As for the \code{limit_state_function}, failure domain is defined by points whom values of \code{limit_fun_MH} are negative.}
  \item{sampling_strategy}{either "AR" or "MH", to specify which sampling strategy is to be used when generating Monte-Carlo population in a case of subset simulation : "AR" stands for \sQuote{accept-reject} while "MH" stands for Metropolis-Hastings.}
  \item{first_DOE}{Either "Gaussian" or "Uniform", to specify the population on which clustering if done}
  \item{seeds}{optional. If \code{sampling_strategy}=="MH", seeds from which MH algorithm starts. This should be a matrix with \sQuote{nrow} = \code{dimension} and \sQuote{ncol} = number of vector.}
  \item{seeds_eval}{optional. The value of the \code{limit_fun_MH} on the \code{seeds}.}
  \item{burnin}{a burnin parameter for Metropolis-Hastings algorithm.}
  \item{thinning}{a thinning parameter for Metropolis-Hastings algorithm. \code{thinning} = 0 means no thinning.}
  \item{plot}{a boolean parameter specifying if function and samples should be plotted. The plot is refreshed at each iteration with the new data. Note that this option is only to be used when working on \sQuote{light} limit state functions as it requires the calculus of this function on a grid of size 161x161 (plot are done a -8:8 x -8:8 grid with 161 meshes.}
  \item{limited_plot}{only a final plot with \code{limit_state_function}, final DOE and metamodel. Should be used with \code{plot}==FALSE.}
  \item{add}{optional. "TRUE" if plots are to be added to the current active device.}
  \item{output_dir}{optional. If plots are to be saved in .jpeg in a given directory. This variable will be pasted with "_AKMCS.jpeg" to get the full output directory.}
  \item{z_MH}{optional. For plots, if metamodel has already been evaluated on the grid then \code{z_MH} (from outer function) can be provided to avoid extra computational time.}
  \item{z_lsf}{optional. For plots, if LSF has already been evaluated on the grid then \code{z_lsf} (from outer function) can be provided to avoid extra computational time.}
  \item{verbose}{Eiher 0 for an almost no output message, or 1 for medium size or 2 for full size}
}

\value{
  An object of class \code{list} containing the failure probability and some more outputs as described below:

  \item{proba}{The estimated failure probability.}
  \item{cov}{The coefficient of variation of the Monte-Carlo probability estimate.}
  \item{Ncall}{The total number of calls to the \code{limit_state_function}.}
  \item{learn_db}{The final learning database, ie. all points where \code{limit_state_function} has been calculated.}
  \item{lsf_value}{The value of the \code{limit_state_function} on the learning database.}
  \item{meta_fun}{The metamodel approximation of the \code{limit_state_function}. A call output is a list containing the value and the standard deviation.}
  \item{meta_model}{The final metamodel. An S4 object from \pkg{DiceKriging}.}
  \item{points}{Points in the failure domain according to the metamodel.}
  \item{meta_eval}{Evaluation of the metamodel on these points.}
  \item{z_meta}{If \code{plot}==TRUE, the evaluation of the metamodel on the plot grid.}
}

\details{
  AKMCS strategy is based on a original Monte-Carlo population which is classified with a kriging-based metamodel. This means that no sampling is done during refinements steps. Indeed, it tries to classify this Monte-Carlo population with a confidence greater than a given value, for instance \sQuote{distance} to the failure should be greater than \code{crit_min} standard deviation.

  Thus, while this criterion is not verified, the point minimizing it is added to the learning database and then evaluated.

  Finally, once all points are classified or when the maximum number of calls has been reached, crude Monte-Carlo is performed. A final test controlling the size of this population regarding the targeted coefficient of variation is done; if it is too small then a new population of sufficient size (considering ordre of magnitude of found probability) is generated, and algorithm run again.
}

\note{
Problem is supposed to be defined in the standard space. If not, use \code{\link{UtoX}} to do so.
Furthermore, each time a set of vector is defined as a matrix, \sQuote{nrow} = \code{dimension} and \sQuote{ncol} = number of vector.
}

\references{
  \itemize{
    \item
      B. Echard, N. Gayton, M. Lemaire:\cr
      \emph{AK-MCS : an Active learning reliability method combining Kriging and Monte Carlo Simulation}\cr
      Structural Safety, Elsevier, 2011.\cr

    \item
      B. Echard, N. Gayton, M. Lemaire and N. Relun:\cr
      \emph{A combined Importance Sampling and Kriging reliability method for small failure probabilities with time-demanding numerical models}\cr
      Reliability Engineering \& System Safety,2012\cr

    \item
      B. Echard, N. Gayton and A. Bignonnet:\cr
      \emph{A reliability analysis method for fatigue design}\cr
      International Journal of Fatigue, 2014\cr
}
}

\author{
  Clement Walter\cr
  \email{clement.walter@cea.fr}
}

\seealso{
  \code{\link{SubsetSimulation}}
  \code{\link{MonteCarlo}}
  \code{\link[DiceKriging]{km}} (in package \pkg{DiceKriging})
}

\examples{
#Limit state function defined by Kiureghian & Dakessian :
kiureghian = function(x, b=5, kappa=0.5, e=0.1) {
x = as.matrix(x)
b - x[2,] - kappa*(x[1,]-e)^2
}

\dontrun{res = AKMCS(dimension=2,limit_state_function=kiureghian,plot=TRUE)

#Compare with crude Monte-Carlo reference value
N = 500000
U = matrix(rnorm(dimension*N),dimension,N)
G = apply(U,2,kiureghian)
P = mean(G<0)
cov = sqrt((1-P)/(N*P))
}

#See impact of kernel choice with Waarts function :
waarts = function(x) {
x = as.matrix(x)
apply(x, 2, function(u) {min(
		(3+(u[1]-u[2])^2/10 - (u[1]+u[2])/sqrt(2)),
		(3+(u[1]-u[2])^2/10 + (u[1]+u[2])/sqrt(2)),
		u[1]-u[2]+7/sqrt(2),
		u[2]-u[1]+7/sqrt(2))})
}

\dontrun{
res = list()
res$matern5_2 = AKMCS(2, waarts, plot=TRUE)
res$matern3_2 = AKMCS(2, waarts, kernel="matern3_2", plot=TRUE)
res$gaussian  = AKMCS(2, waarts, kernel="gauss", plot=TRUE)
res$exp       = AKMCS(2, waarts, kernel="exp", plot=TRUE)

#Compare with crude Monte-Carlo reference value
N = 500000
U = matrix(rnorm(dimension*N),dimension,N)
G = apply(U,2,waarts)
P = mean(G<0)
cov = sqrt((1-P)/(N*P))
}

}


