#---------------------------------------------------------------------------------#
#--------------------------------Algorithm FORM :--------------------------------#
#---------------------------------------------------------------------------------#
#' @export
FORM <- function(f, u.dep, inputDistribution, N.calls, eps = 1e-7, 
                 Method = "HLRF", IS = FALSE, q = 0.5, copula = "unif"){
  
  # Few error messages : 
  if (mode(f) != "function") stop("f is not a function")
  if( mode(inputDistribution) != "list") stop("inputDistribution is not a list")
  if (missing(inputDistribution)) stop("inputDistribution is missing")
  if (missing(N.calls)) stop("N.calls is missing")
  
  ndim <- length(inputDistribution)
  
  #-------------------------------------------------------------------------------#
  #--------------- Create the list whith paramters of the input ------------------#
  transformtionToInputSpace <- function(inputDistribution){
    
    InputDist <- list()
    InputDist <- inputDistribution
    
    for(i in 1:ndim){
      
      nparam <- length(inputDistribution[[i]][[2]])
      
      for(j in 1:nparam){
        InputDist[[i]]$q <- paste("q",InputDist[[i]][[1]], sep = "");
        InputDist[[i]]$p <- paste("p",InputDist[[i]][[1]], sep = "");
        InputDist[[i]]$d <- paste("d",InputDist[[i]][[1]], sep = "");
        InputDist[[i]]$r <- paste("r",InputDist[[i]][[1]], sep = "");
      }
    }
    
    return(InputDist)
  }
  
  InputDist <- transformtionToInputSpace(inputDistribution)
  
  #------ Transformation to the gaussian space ------#

  G <- function(X){
    ndim <- length(InputDist)
    XU   <- pnorm(X)
    
    for(i in 1:ndim){
      XU[i] <- do.call(InputDist[[i]]$q, c(list(XU[i, drop = FALSE]), InputDist[[i]][[2]]))
    }
    return(f(XU))
  }

  #------ Gradient------#
  GRAD <- function(func,X){
    GG <- 0
    df <- 0
    RES <- list()

    RES[[1]] <- 0
    RES.2 <- NULL
    RES.3 <- NULL

    dfunc <- 0;
    
    for (i in 1:ndim) {
      dX <- X
      dX[i] <- dX[i] + 1e-6 
      f.dX <- func(dX)
      f.X <- func(X)
      dfunc[i] <- (f.dX-f.X)/1e-6
      RES.2 <- rbind(RES.2, dX)
      RES.3 <- rbind(RES.3, f.dX)
      RES.2 <- rbind(RES.2, X)
      RES.3 <- rbind(RES.3, f.X)     
    }
    RES[[2]] <- RES.2
    RES[[3]] <- RES.3
    RES[[1]] <- dfunc

    return(RES)
  }
  #----------------------------------------------------------------------#
  DOE      <- list()                       # Design of experiments

  DOE[[1]] <- NULL                         # DOE[[1]] contains x
  DOE[[2]] <- NULL                         # DOE[[2]] contains f(x)

  List.res <- list()				 # List of the output
  u.new    <- 0
  cp       <- 0   			       #Numbers of calls to f
  fact.imp <- 0					 #importance factor
  N.FORM   <- floor(N.calls*q) 		 #Number of calls reserved to FORM
  
  u.ans <- G(u.dep); cp <- cp + 1

  DOE[[1]] <- u.dep
  DOE[[2]] <- u.ans

  res     <- GRAD(G, u.dep); cp <- cp + 2*ndim

  g <- res[[1]]

  DOE[[1]] <- rbind(DOE[[1]],res[[2]])
  DOE[[2]] <- rbind(DOE[[2]],res[[3]])

  err   <- 1					       #distance of two tested points

  #----------------------------------Algorithme HLRF :----------------- #
  if(Method == "HLRF"){
    
    a <- g/(sqrt(g%*%g))
     
    while( (err > eps) & (cp <= (N.FORM - (1 + 2*ndim)) ) ){ # modif BIS
      
      u.new <- (t(u.dep)%*%a)*a - u.ans*a/sqrt(g%*%g)
      err   <- sqrt(t(u.dep-u.new)%*%(u.dep-u.new))
      u.dep <- u.new
      
      u.ans <- G(u.new); cp <- cp + 1
      DOE[[1]] <- rbind(DOE[[1]],u.new)
      DOE[[2]] <- rbind(DOE[[2]],u.ans)    
 
      res   <- GRAD(G, u.dep); cp <- cp + 2*ndim

      g <- res[[1]]

      DOE[[1]] <- rbind(DOE[[1]],res[[2]])
      DOE[[2]] <- rbind(DOE[[2]],res[[3]])      

      a <- ifelse( g == 0,c(0,0), g/(sqrt(g%*%g)) )
    }
    
    B 	 <- -a%*%u.dep 				 #indices of reliability
    P        <- pnorm(-B)   				 #failure probability
    fact.imp <- a^2				       
    List.res <- list(P, B, cp, u.dep, fact.imp, DOE);
    
    names(List.res) <- c("proba.def", "indice.reliab", "compt.f", "Design.Point", "fact.imp", "DOE")
    if(IS == FALSE){return(List.res)}
    }
  #-----------------------------End of the algorithm HLRF :--------------------------#
  
  #-----------------------------Algorithm Abdo-Rackwitz :--------------------------#
  if(Method == "AR"){
    
    lambda <- 2*(u.ans - (t(g)%*%u.dep))/(g%*%g)
    
    while( (err > eps) & (cp <= (N.FORM - (1 + 2*ndim)) ) ){ # modif BIS
        
      u.new <- -lambda*g/2
      err   <- sqrt(t(u.dep-u.new)%*%(u.dep-u.new))
      u.dep <- u.new
      
      u.ans <- G(u.new); cp <- cp + 1
      DOE[[1]] <- rbind(DOE[[1]],u.new)
      DOE[[2]] <- rbind(DOE[[2]],u.ans)    
 
      res   <- GRAD(G, u.dep); cp <- cp + 2*ndim

      g <- res[[1]]

      DOE[[1]] <- rbind(DOE[[1]], res[[2]])
      DOE[[2]] <- rbind(DOE[[2]], res[[3]])
      
      lambda <- 2*(u.ans - (t(g)%*%u.dep))/(g%*%g)
      if(cp >= N.FORM){break}
    }
    
    a <- g/(sqrt(g%*%g))
    B <- -a%*%u.dep  			
    P <- pnorm(-B)  				

    fact.imp <- a^2				
    List.res <- list(P, B, cp, u.dep, fact.imp, DOE);
    
    names(List.res) <- c("proba.def", "indice.reliab", "compt.f", "Design.Point", "fact.imp", "DOE")
    if(IS == FALSE){return(List.res)}
    }
  #------------------------End of the algorithm of Abdo-Rackwitz----------------------#
  
  #---------------------------Importance Sampling :---------------------------------#
  if(IS == TRUE){
    alpha  <- 0.05				
    N.IS   <- N.calls - cp			                 #Number of call reserved to the importance sampling
    norm.u <- u.dep%*%u.dep
    P.temp <- 0
    P 	 <- 0
    VAR 	 <- 0
    s	 <- 0
    Z 	 <- numeric(N.IS)

    v.temp <- matrix(rnorm(N.IS*ndim,0,1),ncol=ndim)
    v.temp <- v.temp  +  matrix(rep(u.dep,N.IS),ncol=ndim,byrow=TRUE)     #One simulate points around u.dep

    zz <- v.temp -  matrix(rep(u.dep,N.IS), ncol = ndim, byrow = TRUE)
    s  <- exp( u.dep%*%u.dep/2 - as.numeric(v.temp%*%u.dep) )

    test.v.temp <- apply(v.temp, MARGIN = 1, G)
    DOE[[1]] <- rbind(DOE[[1]], v.temp)
    DOE[[2]] <- rbind(DOE[[2]], matrix(test.v.temp,ncol=1))
    indic.Def <- 1*(test.v.temp <=0) 
    P.temp <- indic.Def*s

    VAR <- sum(P.temp^2) 
    cp  <- cp + N.IS					            #Numbers of call to the failure fonction inside the "sapply"
    P   <- sum(P.temp)/N.IS    			            #failure probability
    B   <- -qnorm(P)  	        			            #indices of reliability
    
    VAR <- (sum(P.temp^2)/N.IS - P^2)/(N.IS-1)              #standard error of the estimator
    
    IC.inf     <- P - qnorm(1-alpha/2)*sqrt(VAR) 		#low bound of confidence interval
    IC.sup     <- P + qnorm(1-alpha/2)*sqrt(VAR)		#high bound of confidence interval

    Inter.conf <- c(IC.inf,IC.sup)
    fact.imp   <- a^2					          
    List.res   <- list(P, B, cp, u.dep, fact.imp, VAR, Inter.conf, DOE)
    
    names(List.res) <-c("proba.def", "indice.reliab", "compt.f", "Design.Point", "fact.imp", "variance", "Interval.conf", "DOE");
    return(List.res)
  } 
  #----------------------------End Importance Sampling----------------------------#
  
  }
#---------------------------------------------------------------------------------#
#----------------------------- End of algorithm  ---------------------------------#
#------------------------------------FORM-----------------------------------------#
#---------------------------------------------------------------------------------#
