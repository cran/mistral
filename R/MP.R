#' @title Moving Particles
#'
#' @description This function runs the Moving Particles algorithm for estimating extreme probability
#' and quantile.
#'
#' @author Clement WALTER \email{clementwalter@icloud.com}
#'
#' @aliases LPA
#'
#' @details \code{MP} is a wrap up of \code{\link{IRW}} for probability and quantile estimation.
#'  By construction, the several calls to \code{\link{IRW}} are parallel (\pkg{foreach})
#' and so is the algorithm. Especially, with \code{N.batch=1}, this is the Last Particle
#' Algorithm, which is a specific version of \code{\link{SubsetSimulation}} with \code{p_0 = 1-1/N}.
#' However, note that this algorithm not only gives a quantile or a probability estimate
#' but also an estimate of the whole cdf until the given threshold \code{q}.
#'
#' The probability estimator only requires to generate several random walks as it is the estimation
#' of the parameter of a Poisson random variable. The quantile estimator is a little bit more complicated
#' and requires a 2-passes algorithm. It is thus not exactly fully parallel as cluster/cores have to
#' communicate after the first pass. During the first pass, particles are moved a given number of
#' times, during the second pass particles are moved until the farthest event reach during the first
#' pass. Hence, the random process is completely simulated until this given state.
#'
#' For an easy user experiment, all the parameters are defined by default with the optimised values
#' as described in the reference paper (see References below) and a typical use will only specify
#' \code{N} and \code{N.batch}.
#'
#' @note The \code{alpha} parameter is set to 0.05 by default. Indeed it should not be
#' set too small as it is defined approximating the Poisson distribution with the Gaussian one.
#' However if no estimate is produce then the algorithm can be restarted for the few missing events.
#' In any cases, setting \code{Niter_1fold = -N/N.batch*log(p)} gives 100\% chances to produces
#' a quantile estimator.
#'
#' @seealso
#' \code{\link{SubsetSimulation}}
#' \code{\link{MonteCarlo}}
#' \code{\link{IRW}}
#'
#' @references
#'   \itemize{
#'    \item A. Guyader, N. Hengartner and E. Matzner-Lober:\cr
#'     \emph{Simulation and estimation of extreme quantiles and extreme
#'     probabilities}\cr
#'     Applied Mathematics and Optimization, 64(2), 171-196.\cr
#'    \item C. Walter:\cr
#'    \emph{Moving Particles: a parallel optimal Multilevel Splitting
#'    method with application in quantiles estimation and meta-model
#'    based algorithms}\cr
#'    Structural Safety, 55, 10-25.\cr
#'    \item E. Simonnet:\cr
#'    \emph{Combinatorial analysis of the adaptive last particle method}\cr
#'    Statistics and Computing, 1-20.\cr
#'  }
#'
#' @examples
#' \dontrun{
#' # Estimate some probability and quantile with the parabolic lsf
#' p.est <- MP(2, kiureghian, N = 100, q = 0) # estimate P(lsf(X) < 0)
#' p.est <- MP(2, kiureghian, N = 100, q = 7.8, lower.tail = FALSE) # estimate P(lsf(X) > 7.8)
#'
#' q.est <- MP(2, kiureghian, N = 100, p = 1e-3) # estimate q such that P(lsf(X) < q) = 1e-3
#' q.est <- MP(2, kiureghian, N = 100, p = 1e-3, lower.tail = FALSE) # estimate q such
#' # that P(lsf(X) > q) = 1e-3
#'
#'
#' # plot the empirical cdf
#' plot(xplot <- seq(-3, p.est$L_max, l = 100), sapply(xplot, p.est$ecdf_MP))
#'
#' # check validity range
#' p.est$ecdf_MP(p.est$L_max - 1)
#' # this example will fail because the quantile is greater than the limit
#' tryCatch({
#'    p.est$ecdf_MP(p.est$L_max + 0.1)},
#'    error = function(cond) message(cond))
#'
#' # Run in parallel
#' library(doParallel)
#' registerDoParallel()
#' p.est <- MP(2, kiureghian, N = 100, q = 0, N.batch = getDoParWorkers())
#' }
#'
#' @import ggplot2
#' @import foreach
#' @import iterators
#' @import doParallel
#' @export
MP = function(dimension,
              #' @param dimension the dimension of the input space.
              lsf,
              #' @param lsf the function defining the RV of interest Y = lsf(X).
              N = 100,
              #' @param N the total number of particles,
              N.batch = foreach::getDoParWorkers(),
              #' @param N.batch the number of parallel batches for the algorithm. Each batch will then
              #' have \code{N/N.batch} particles. Typically this could be \code{detectCores()} or some
              #' other machine-derived parameters. Note that \code{N/N.batch} has to be an integer.
              p,
              #' @param p a given probability to estimate the corresponding quantile (as in qxxxx functions).
              q,
              #' @param q a given quantile to estimate the corresponding probability (as in pxxxx functions).
              lower.tail = TRUE,
              #' @param lower.tail as for pxxxx functions, TRUE for estimating P(lsf(X) < q), FALSE
              #' for P(lsf(X) > q).
              Niter_1fold,
              #' @param Niter_1fold a function = fun(N) giving the deterministic number of iterations
              #' for the first pass.
              alpha = 0.05,
              #' @param alpha when using default \code{Niter_1fold} function, this is the risk not to have
              #' simulated enough samples to produce a quantile estimator.
              compute_confidence = FALSE,
              #' @param compute_confidence if \code{TRUE}, the algorithm runs a little bit longer to produces
              #' a 95\% interval on the quantile estimator.
              verbose = 0,
              #' @param verbose to control level of print (either 0, or 1, or 2).
              chi2 = FALSE,
              #' @param chi2 for a chi2 test on the number of events.
              breaks = N.batch/5,
              #' @param breaks for the final histogram is \code{chi2 == TRUE}.
              ...) {
  #' @param ... further arguments past to \code{\link{IRW}}.
  
  cat('========================================================================\n')
  cat('                       Beginning of MP algorithm\n')
  cat('========================================================================\n')
  ## STEP 0 : INITIALISATION
  
  # Fix NOTE issue with R CMD check
  k <- NULL
  
  # Define transpose for list of lists with same fields, eg output of foreach
  t.list <- function(l){
    lapply(split(do.call("c", l), names(l[[1]])), unname)
  }
  arg <- list(...)
  arg$dimension <- dimension
  arg$lsf <- lsf
  arg$N <- N/N.batch
  if(lower.tail==TRUE){
    lsf_dec = arg$lsf
    arg$lsf = function(x) -1*lsf_dec(x)
  }
  
  estim.proba <- missing(p)
  estim.quantile <- !estim.proba
  
  if(estim.quantile){ # quantile estimation
    
    if(missing(Niter_1fold)){
      Niter_1fold = function(N,n = N.batch) {
        t = -log(p)
        t = t + 1.96*sqrt(t/n/N)*compute_confidence
        if(n>1){
          bn = sqrt(2*log(n)) - 0.5*( log(log(n)) + log(4*pi) )/sqrt(2*log(n))
          beta = bn - log(log(1/alpha))/sqrt(2*log(n))
          delta = beta^2 + 4*N*t
        }
        else{
          bn = 0
          beta = 0
          delta = 0
        }
        res = ceiling(0.5*(beta^2 + 2*N*t - beta*sqrt(delta) ))
        return(res)
      }
    }
    moves = c()
    
    cat(" =============================================== \n")
    cat(" STEP 1 : MOVE PARTICLES A GIVEN NUMBER OF TIMES \n")
    cat(" =============================================== \n\n")
    
    cat(" * Number of deterministic event per algorithm =",Niter_1fold(N),"\n\n")
    
    cat(" ### PARALLEL PART ###\n")
    cat(" * backend:", ifelse(is.null(foreach::getDoParName()), "no registered backend", foreach::getDoParName()), "\n")
    cat(" * N.batch =", N.batch, "\n\n")
    mp.1 <- foreach(k=icount(N.batch)) %dopar% {
      arg$Nevent = Niter_1fold(N)
      if(verbose<2){
        capture.output(mp <- do.call(IRW,arg))
      } else {
        mp <- do.call(IRW, arg)
      }
      mp[c('L', 'Ncall', 'X', 'y')]
    }
    
    t.mp.1 <- t(mp.1)
    Ncall <- unlist(t.mp.1$Ncall)
    L <- unlist(t.mp.1$L)
    L_max <- max(L)
    # cat("\n ### END OF PARALLEL PART; furthest point =", (-1)^lower.tail*L_max,"\n\n")
    
    cat(" ===================================================== \n")
    cat(" STEP 2 : RESTART ALGORITHM UNTIL",(-1)^lower.tail*L_max," \n")
    cat(" ===================================================== \n\n")
    
    cat(" ### PARALLEL PART ###\n")
    cat(" * backend:", foreach::getDoParName(), "\n")
    cat(" * N.batch =", N.batch, "\n\n")
    mp.2 <- foreach(i=iter(mp.1)) %dopar% {
      arg = c(list(
        X = i$X,
        y = i$y,
        q = L_max,
        last.return = FALSE),
        arg)
      if(verbose<2){
        capture.output(mp <- do.call(IRW,arg))
      } else {
        mp <- do.call(IRW,arg)
      }
      tryCatch(
        if(max(mp$L) > L_max){stop("Last time of algorithms should not be greater than furthest time of the first step \n")},
        warning = function(cond){}
      )
      mp$L <- mp$L[-1] # First time equals last one from first pass
      mp[c('L', 'Ncall', 'X', 'y')]
    }
    t.mp.2 <- t(mp.2)
    moves <- sapply(t.mp.2$L, length) + Niter_1fold(N)
    X = do.call(cbind, t.mp.2$X)
    y = unlist(t.mp.2$y)
    Ncall <- Ncall + unlist(t.mp.2$Ncall)
    L = sort(c(L, unlist(t.mp.2$L)))
    # cat("\n ### END OF PARALLEL PART \n\n")
    
  }
  else{ # probability estimation
    arg$q <- (-1)^lower.tail*q
    
    # cat(" ===================================== \n")
    # cat(" STEP 1 : MOVE PARTICLES UNTIL FAILURE \n")
    # cat(" ===================================== \n\n")
    # 
    cat(" ### PARALLEL PART ###\n")
    cat(" * backend:", foreach::getDoParName(), "\n")
    cat(" * N.batch =", N.batch, "\n\n")
    mp <- foreach(k=icount(N.batch)) %dopar% {
      #cat(" * TRAIN N0",k,"in",N.batch,"\n")
      if(verbose<2){capture.output(mp <- do.call(IRW, arg))}
      else{mp <- do.call(IRW, arg)}
      mp[c('Ncall', 'L', 'X', 'y')]
    }
    t.mp <- t(mp)
    Ncall <- L <- X <- y <- NULL
    lapply(mp, function(l) {
      Ncall <<- c(Ncall, l$Ncall)
      L <<- c(L, l$L)
      X <<- cbind(X, l$X)
      y <<- c(y, l$y)
    })
    moves <- sapply(t.mp$L, length)
    # cat(" ### END OF PARALLEL PART \n\n")
  }
  
  cat("========================================================================\n")
  cat("                        End of MP algorithm \n")
  cat("========================================================================\n\n")
  
  Ntot = length(y)
  if(estim.proba){
    p = (1 - 1/Ntot)^sum(L<(q*(-1)^lower.tail))
    p_int = c( p*exp(-1.96*sqrt(-log(p)/Ntot)) , p*exp(+1.96*sqrt(-log(p)/Ntot)) )
    cov = sqrt(-log(p)/Ntot)
    L_max <- q
    q_int <- c(q, q)
  }
  else{
    m = ceiling(-Ntot*log(p))
    q = 0.5*(L[m-1] + L[m])*(-1)^lower.tail
    q_int = c( L[floor(m-1.96*sqrt(m))] , L[ceiling(m+1.96*sqrt(m))] )*(-1)^lower.tail
    cov = NA
    p_int <- c(p, p)
    L_max <- (-1)^lower.tail*L_max
  }
  L <- (-1)^lower.tail*L
  ecdf_MP <- local({
    lower.tail <- lower.tail
    L_max <- L_max
    L <- L;
    Ntot <- Ntot
    function(q) {
      sapply(q, function(q) {
        if(lower.tail==TRUE){
          Nevent <- sum(L>q)
          if(q<L_max) {
            message(paste('Empirical cdf valid only for q >=', L_max))
            Nevent <- NA
          }
        }
        else{
          Nevent <- sum(L<q)
          if(q>L_max) {
            message(paste('Empirical cdf valid only for q =<', L_max))
            Nevent <- NA
          }
        }
        (1-1/Ntot)^Nevent
      })
    }
  })
  
  if(chi2 == TRUE) {
    res = hist(moves, breaks = breaks, freq = FALSE)
    cont = res$counts; l = length(cont)
    chi2_p = 1:l;
    br = res$breaks
    br[1] = -Inf; br[length(br)] = Inf;
    for(i in 1:l){
      chi2_p[i] = ppois(br[i+1],mean(moves)) - ppois(br[i], mean(moves))
    }
    capture.output(chi2.test <- chisq.test(x = cont, p = chi2_p))
    p.val = 1 - pchisq(chi2.test$statistic, df = (l-2))
  }
  
  cat("   - p =",p,"\n")
  cat("   - q =",q,"\n")
  if(estim.proba){
    cat("   - 95% confidence intervalle :",p_int[1],"< p <",p_int[2],"\n")
  }
  else{
    if(lower.tail==FALSE){
      cat("   - 95% confidence intervalle :",q_int[1],"< q <",q_int[2],"\n")
    }
    else{
      cat("   - 95% confidence intervalle :",q_int[2],"< q <",q_int[1],"\n")
    }
    cat("   - Maximum number of moves/algorithm =",max(moves),"\n")
    cat("   - L_max =",L_max,"\n")
    cat("   - Number of moves =",length(L),"\n")
    cat("   - Targeted number of moves =",ceiling(m+compute_confidence*1.96*sqrt(m)),"\n")
  }
  cat("   - Total number of calls =",sum(Ncall),"\n")
  if(chi2 == TRUE) {cat("   - Chi-2 test =", chi2.test$statistic,"; p-value =", p.val,"\n")}
  
  #' @return An object of class \code{list} containing the outputs described below:
  res = list(p = p,
             #' \item{p}{the estimated probability or the reference for the quantile estimate.}
             q = q,
             #' \item{q}{the estimated quantile or the reference for the probability estimate.}
             cv = ifelse(estim.proba, sqrt(sum(moves))/Ntot, 0),
             #' \item{cv}{the coefficient of variation of the probability estimator.}
             ecdf = ecdf_MP,
             #' \item{ecdf}{the empirical cdf.}
             L = L,
             #' \item{L}{the states of the random walk.}
             L_max = L_max,
             #' \item{L_max}{the farthest state reached by the random process. Validity range
             #' for the \code{ecdf} is then (-Inf, L_max] or [L_max, Inf).}
             times = L,
             #' \item{times}{the \emph{times} of the random process.}
             Ncall = Ncall,
             #' \item{Ncall}{the total number of calls to the \code{lsf}.}
             X = X,
             #' \item{X}{the \code{N} particles in their final state.}
             y = y,
             #' \item{y}{the value of the \code{lsf(X)}.}
             moves = moves)
  #' \item{moves}{a vector containing the number of moves for each batch.}
  
  if(estim.proba){
    res <- c(res, list(
      p_int = p_int,
      #' \item{p_int}{a 95\% confidence intervalle on the probability estimate.}
      cov = cov
      #' \item{cov}{the coefficient of variation of the estimator}
    ))
  }
  else{
    res <- c(res, list(
      q_int = q_int
      #' \item{q_int}{a 95\% confidence intervall on the quantile estimate.}
    ))
  }
  if(chi2 == TRUE) {res = c(res, list(chi2 = chi2.test))}
  #' \item{chi2}{the output of the chisq.test function.}
  return(res)
}
