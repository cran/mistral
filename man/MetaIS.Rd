\name{MetaIS}
\alias{MetaIS}
\title{Metamodel based Impotance Sampling}

\description{ Estimate failure probability by MetaIS method.}

\usage{
MetaIS(dimension,
      limit_state_function,
      N                 = 500000,
      N_alpha           = 100,
      N_DOE             = 2*dimension,
      N1                = N_DOE*30,  
      Ru                = 8,
      Nmin              = 30,
      Nmax              = 200,
      Ncall_max         = 1000,
      precision         = 0.05,
      N_seeds           = 1,
      Niter_seed        = 10000,
      N_alphaLOO        = 5000,
      K_alphaLOO        = 2*dimension,
      alpha_int         = c(0.1,10),
      k_margin          = 1.96,
      learn_db          = NULL,
      lsf_value         = NULL,
      failure           = 0,
      meta_model        = NULL,
      kernel            = "matern5_2",
      learn_each_train  = FALSE,
      limit_fun_MH      = NULL,
      sampling_strategy = "MH",
      seeds             = NULL,
      seeds_eval        = NULL,
      burnin            = 20,
      thinning          = 4,
      plot              = FALSE,
      limited_plot      = FALSE,
      add               = FALSE,
      output_dir        = NULL,
      z_MH              = NULL,
      z_lsf             = NULL,
      verbose		        = 0)
}

\arguments{
  \item{dimension}{an integer giving the dimension of the input space.}
  \item{limit_state_function}{the failure fonction.}
  \item{N}{an integer defining the size of Monte-Carlo population for augmented failure probability estimation.}
  \item{N_alpha}{an integer defining the initial size of Monte-Carlo population for alpha estimate.}
  \item{N_DOE}{an integer defining the size of initial DOE got by clustering of the \code{N1} samples.}
  \item{N1}{an integer defining the size of the initial uniform population sampled in a hypersphere of radius Ru.}
  \item{Ru}{a real defining the radius of the hypersphere for the initial uniform sampling.}
  \item{Nmin}{an integer defining the minimum number of calls to the \code{limit_state_function} for the construction step.}
  \item{Nmax}{an integer defining the maximum number of calls to the \code{limit_state_function} for the construction step.}
  \item{Ncall_max}{an integer defining the maximum number of calls to the \code{limit_state_function} for the whole algorithm.}
  \item{precision}{a real defining the targeted maximal value of the coefficient of variation.}
  \item{N_seeds}{an integer defining the number of seeds for Metropolis-Hastings algorithm while generating into the margin (according to MarginProbability*gauss density).}
  \item{Niter_seed}{maximum number of iterations for the research of a seed for alphaLOO refinement sampling.}
  \item{N_alphaLOO}{an integer defining the number of points to sample at each refinement step.}
  \item{K_alphaLOO}{an integer defining the number of clusters at each refinement step.}
  \item{alpha_int}{a 2-d real vector defining the range for alphaLOO criterion : \code{alpha_int}[1] < alphaLOO < \code{alpha_int}[2].}
  \item{k_margin}{a real defining the margin width according to standard gaussian law; default is 1.96 which means points outside of the margin are classified with more than 97,5\%.}
  \item{learn_db}{optional. A matrix of already known points, with dim : \code{dimension} x number_of_vector.}
  \item{lsf_value}{values of the limit state function on the vectors given in \code{learn_db}.}
  \item{failure}{the value defining the failure domain F = \{ x | \code{limit_state_function}(x) < \code{failure} \}.}
  \item{meta_model}{optional. If a kriging based metamodel has already been fitted to the data (from \pkg{DiceKriging} package) it can be given as an input to keep the same parameters.}
  \item{kernel}{a specified kernel to be used for the metamodel. See \pkg{DiceKriging} for available options.}
  \item{learn_each_train}{specify if hyperparameters of the model should be evaluated each time points are added to the learning database ("TRUE") or only the first time ("FALSE").}
  \item{limit_fun_MH}{optional. If the working space is to be reduced to some subset defining by a function, eg. in case of use in a Subset Simulation algorithm. As for the \code{limit_state_function}, failure domain is defined by points whom values of \code{limit_fun_MH} are negative.}
  \item{sampling_strategy}{either "AR" or "MH", to specify which sampling strategy is to be used when generating Monte-Carlo population in a case of subset simulation : "AR" stands for \sQuote{accept-reject} while "MH" stands for Metropolis-Hastings.}
  \item{seeds}{optional. If \code{sampling_strategy}=="MH", seeds from which MH algorithm starts. This should be a matrix with \sQuote{nrow} = \code{dimension} and \sQuote{ncol} = number of vector.}
  \item{seeds_eval}{optional. The value of the \code{limit_fun_MH} on the \code{seeds}.}
  \item{burnin}{a burnin parameter for Metropolis-Hastings algorithm. To be used only for Monte-Carlo population sampling and set to 0 elsewhere.}
  \item{thinning}{a thinning parameter for Metropolis-Hastings algorithm. \code{thinning} = 0 means no thinning. To be used only for Monte-Carlo population sampling and set to 0 elsewhere.}
  \item{plot}{a boolean parameter specifying if function and samples should be plotted. The plot is refreshed at each iteration with the new data. Note that this option is only to be used when working on \sQuote{light} limit state functions as it requires the calculus of this function on a grid of size 161x161 (plot is done a -8:8 x -8:8 grid with 161 meshes.}
  \item{limited_plot}{only a final plot with \code{limit_state_function}, final DOE and metamodel. Should be used with \code{plot}==FALSE.}
  \item{add}{optional. "TRUE" if plots are to be added to the current active device.}
  \item{output_dir}{optional. If plots are to be saved in .jpeg in a given directory. This variable will be pasted with "_MetaIS.jpeg" to get the full output directory.}
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

  MetaIS is an Important Sampling based probability estimator. It makes use of a kriging surogate to approximate the optimal density function, replacing the indicatrice by its kriging pendant, the probability of being in the failure domain. In this context, the normallizing constant of this quasi-optimal PDF is called the \sQuote{augmented failure probability} and the modified probability \sQuote{alpha}.
  
  After a first uniform Design of Experiments, MetaIS uses an alpha-Leave-One-Out criterion combined with a margin sampling strategy to refine a kriging-based metamodel. Samples are generated according to the weighted margin probability with Metropolis-Hastings algorithm and some are selected by clustering; the \code{N_seeds} are got from an accept-reject strategy on a standard population.

  Once criterion is reached or maximum number of call done, the augmented failure probability is estimated with a crude Monte-Carlo. Then, a new population is generated according to the quasi-optimal instrumenal PDF; \code{burnin} and \code{thinning} are used here and alpha is evaluated. While the coefficient of variation of alpha estimate is greater than a given threshold and some computation spots still available (defined by \code{Ncall_max}) the estimate is refined with extra calculus.

  The final probability is the product of p_epsilon and alpha, and final squared coefficient of variation is the sum of p_epsilon and alpha one's.
}

\note{
Problem is supposed to be defined in the standard space. If not, use \code{\link{UtoX}} to do so.
Furthermore, each time a set of vector is defined as a matrix, \sQuote{nrow} = \code{dimension} and \sQuote{ncol} = number of vector.
}

\references{
  \itemize{
    \item
      V. Dubourg:\cr
      Meta-modeles adaptatifs pour l'analyse de fiabilite et l'optimisation sous containte fiabiliste\cr
      PhD Thesis, Universite Blaise Pascal - Clermont II,2011\cr

    \item
      V. Dubourg, B. Sudret, F. Deheeger:\cr
      Metamodel-based importance sampling for structural reliability analysis Original Research Article\cr
      Probabilistic Engineering Mechanics, Volume 33, July 2013, Pages 47-57\cr

    \item
      V. Dubourg, B. Sudret:\cr
      Metamodel-based importance sampling for reliability sensitivity analysis.\cr
      Accepted for publication in Structural Safety, special issue in the honor of Prof. Wilson Tang.(2013)\cr

    \item
      V. Dubourg, B. Sudret and J.-M. Bourinet:\cr
      Reliability-based design optimization using kriging surrogates and subset simulation.\cr
      Struct. Multidisc. Optim.(2011)\cr
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

\dontrun{res = MetaIS(dimension=2,limit_state_function=kiureghian,plot=TRUE)

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
res$matern5_2 = MetaIS(2,waarts,plot=TRUE)
res$matern3_2 = MetaIS(2,waarts,kernel="matern3_2",plot=TRUE)
res$gaussian = MetaIS(2,waarts,kernel="gauss",plot=TRUE)
res$exp = MetaIS(2,waarts,kernel="exp",plot=TRUE)

#Compare with crude Monte-Carlo reference value
N = 500000
U = matrix(rnorm(dimension*N),dimension,N)
G = apply(U,2,waarts)
P = mean(G<0)
cov = sqrt((1-P)/(N*P))
}

}