#-----------------------------------------------------------------------------------------------------------------#
#
#
#
#                                                  Algorithme MRM :
#
#
#
#-----------------------------------------------------------------------------------------------------------------#


#-----------------------------------------------------------------------------------------------------------------#
# Input :											    	    
#-----------------------------------------------------------------------------------------------------------------#

# f : Failure function.								    
# ndim : Dimension of input
# choice.law : a list of length "ndim" which contain name of input distribution and their parameters. For the input "i", choice.law[[i]] = list("name_law",c(parameters1,..., parametersN)  
# dir.monot : vector of size "ndim" which represent the monotonicity of the failure function. dir.monot[i] = -1 (resp. 1) if the code is decreasing (resp. increasing) in direction i.
# N.calls : Number of calls to f allowed
# Method = there is two methods available. "MC_monotone" is an adapation of the Monte Carlo method under constraints of monotony. "MRM" is a sequential sampling method.
# ordre.p : order of magnitude of the search probability (for ndim >= 3)
# silent : if silent = TRUE, print curent number of call to f. Default: FALSE.
#-----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#
# Output :	
#-----------------------------------------------------------------------------------------------------------------#
#Um: Exact lower bounds of the failure probability 
#UM: Exact upper bounds of the failure probability 
#MLE: Maximum likelihood estimator of the failure probability 
#IC.inf: Lower bound of the confidence interval of the failure probability based on MLE
#IC.sup: Upper bound of the confidence interval of the failure probability based on MLE
#CV.MLE: Coefficient of variation of the MLE
#N.tot: Total number of simulation (just for "MC_monotone")
#-----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#


MRM <- function(f, ndim, choice.law, dir.monot, N.calls, Method,  ordre.p = 0, silent = TRUE){


  InputDistribution <- function(choice.law){

    InputDist <- list()
    InputDist <- choice.law

    for(i in 1:ndim){
      nparam <- length(choice.law[[i]][[2]])
        for(j in 1:nparam){
          InputDist[[i]]$q <- paste("q", InputDist[[i]][[1]], sep = "");
          InputDist[[i]]$p <- paste("p", InputDist[[i]][[1]], sep = "");
          InputDist[[i]]$d <- paste("d", InputDist[[i]][[1]], sep = "");
          InputDist[[i]]$r <- paste("r", InputDist[[i]][[1]], sep = "");
        }
    }
    InputDist
  }

  InputDist <- InputDistribution(choice.law)

#-----------------------------------------------------------------------------------------------------------------#
#	                           Transformation in the uniform space  
#-----------------------------------------------------------------------------------------------------------------#
  G <- function(X){
    XU <- numeric()
    for(i in 1:ndim){
      if(dir.monot[i] == -1){X[i] <- 1 - X[i]}
        XU[i] <- do.call(InputDist[[i]]$q,c(list(X[i, drop = FALSE]), InputDist[[i]][[2]]))
      }
    return(f(XU))
  }
#-----------------------------------------------------------------------------------------------------------------#
#                Initialisation of the method : a dichotomie on the Uniforme space
#-----------------------------------------------------------------------------------------------------------------#

  Intersect <- function(ndim, FUNC){

    a     <- 2
    k     <- 2
    res   <- list()
    u.new <- 0
    temp  <- 0
    u.dep <- list()
    out   <- list()
    comp  <- 2

    u.dep[[1]]  <- rep(1/2, ndim)
    temp        <- FUNC(u.dep[[1]])
    u.dep[[2]]  <- sign(temp)
    cp 	        <- 1 						# compteur d'appel à G
    u.other     <- u.dep
    u.new       <- u.dep[[1]]
    LIST        <- list()
    LIST[[1]]   <- u.dep
    list.set    <- 0
    list.set[1] <- LIST[[1]][2]

    if(temp > 0){
      u.new <- u.dep[[1]] - 1/(a^k)
    }else{
      u.new <- u.dep[[1]] + 1/(a^k)
    }

    eps <- ( u.new - u.dep[[1]] )%*%( u.new - u.dep[[1]] )

    if( ( u.other[[2]] != sign(temp)) & (eps > 1e-7) ){
      exist   <- exist + 1
      u.other <- list( u.dep[[1]], sign(temp) )
    }

    k          <- k + 1
    sign.0     <- sign(temp)
    sign.other <- -sign.0

    while(sign(temp) != sign.other){         #N.dicho définie en début de programme
 
      u.dep[[1]] <- u.new 
      temp       <- FUNC(u.dep[[1]])
      u.dep[[2]] <- sign(temp)
      cp 	 <- cp + 1

      if(temp > 0){
        u.new <- u.dep[[1]] - 1/(a^k)
      }else{
        u.new <- u.dep[[1]] + 1/(a^k)
      }

      eps <- ( u.new - u.dep[[1]] )%*%( u.new - u.dep[[1]] )
      k              <- k + 1
      LIST[[comp]]   <- u.dep
      list.set[comp] <- LIST[[comp]][[2]]
      comp           <- comp + 1
    }

    return(LIST)
    list.set <- as.numeric(list.set)

    if( abs(sum(list.set)) == length(LIST)){

      res[[1]] <- LIST[[length(LIST)]][[1]]
      res[[2]] <- LIST[[length(LIST)]][[2]]
      out <- list(res, cp)
      return(out)

    }else{

      u.dep[[1]] <- LIST[[max(which(list.set == -1))]][[1]]
      u.dep[[2]] <- LIST[[max(which(list.set == -1))]][[2]]

      u.other[[1]] <- LIST[[max(which(list.set == 1))]][[1]]
      u.other[[2]] <- LIST[[max(which(list.set == 1))]][[2]]

      res[[1]] <- rbind(u.dep[[1]], u.other[[1]])
      res[[2]] <- c(u.dep[[2]], u.other[[2]])

      out <- list(res, cp)

      return(out)
    }

  }

#-----------------------------------------------------------------------------------------------------------------#
#
#                         		 Test if the points of the set "x" are
#				                  "smaller" or "greater"  than "y".
#					             If set == 1 : Return TRUE if x[i, ] <= y				   
#					             If set == 2 : Return TRUE if x[i, ] => y		
#		   
#-----------------------------------------------------------------------------------------------------------------#

  is.dominant <- function(x, y, ndim, set){

     if((set != 1)&(set != 2)){
	 print("ERROR : set must to be equal to 1 or 2.")
	 break()
	} 

    dominant <- NULL;

# If x is a vector
    if( is.null(dim(x)) ){
	if(set == 2){
          if ( sum(x >= y) == ndim ){
            return(TRUE)
          }else{
            return(FALSE)
          } 
        }else{
          if( sum(x <= y) == ndim ){
            return(TRUE)
          }else{
            return(FALSE)
          }
	}
    }
	
    y.1 <- NULL
    Y.2 <- NULL
    y.1 <- rep(y,dim(x)[1])
    y.2 <- matrix(y.1, ncol = ndim, byrow = TRUE)

    if(set == 2){
      dominant <- apply(x >= y.2, 1, sum) == ndim
    }else{
      dominant <- apply(x <= y.2, 1, sum) == ndim
    }

    return(dominant)
  }



#-----------------------------------------------------------------------------------------------------------------#
#		              Function which compute the exact bounds
#-----------------------------------------------------------------------------------------------------------------#

  Volume.bounds <- function(X.MC, S, set){ 

   if(set == 1){S <- 1 - S}
 
    if(is.null(dim(S))){
      if(set == 1){
        #return(prod(1 - S))
        return(1 - prod(S))
      }
      if(set == 2){
        return(prod(S))
      }
    }

    DS <- dim(S)[1] 
    if(ndim == 2){
      S    <- S[order(S[,1]),]
      res  <- diag(outer(S[,1], c(0,S[1:(DS - 1), 1]), "-"))
      res1 <- res%*%S[,2]
      res1 <- ifelse(set == 1, 1 - res1, res1)
      return(res1) 
    }


    MC.VOL <- 0
    RES.VOL <- 0
    MM <- dim(X.MC)[1]

    for(i in 1:(DS-1)){
      u     <- apply(S, MARGIN = 1, prod)
      u.max <- which.max(u)
      uu    <- S[u.max, ]

      S     <- S[-u.max, ]
      ss.1  <- is.dominant(X.MC, uu, ndim, 1)
      ss.2  <- is.dominant(X.MC, uu, ndim, 2)

      RES.VOL <- RES.VOL + sum(ss.1)

      X.MC <- X.MC[which((ss.1 == 0)&(ss.2 == 0) ) , ]
    }

    uu <- S
    if(set == 1){
      tt <- is.dominant(X.MC, uu, ndim, 2)
      RES.VOL <- RES.VOL + sum(tt)
      RES.VOL <- 1 - RES.VOL/MM
    }else{
      tt <- is.dominant(X.MC, uu, ndim, 1)
      RES.VOL <- RES.VOL + sum(tt)
      RES.VOL <- RES.VOL/MM
    }
    return(RES.VOL) 

  }

#-----------------------------------------------------------------------------------------------------------------
#				                     Gives the frontier of a set
#-----------------------------------------------------------------------------------------------------------------

  Frontier <- function(S, set){
    if(is.null(dim(S)) |( dim(S)[1] == 1)){
      return(S)
    }
    R <- NULL
    if(set == 1){
      S <- 1 - S
    }
    while(!is.null(dim(S))){
      aa  <- apply(S, MARGIN = 1, prod)
      temp <- S[which.max(aa), ]
      R <- rbind(R, temp)
      S <- S[-which.max(aa), ]
      ss <- is.dominant(S, temp, ndim, -1)
      if(!is.null(dim(S))){
        S <- S[which(ss == FALSE),]
      }else{
        S <- matrix(S, ncol= ndim)
        S <- S[which(ss == FALSE), ]
        R = rbind(R, S)
        if(set == 1){
          return(1 - R)
        }else{
          return(R)
        }
      }
    }
   
    if(set == 1){
      return(1 - R)
    }else{
      return(R)
    }
  }


#-----------------------------------------------------------------------------------------------------------------
#
#
#
#     				Monte Carlo method under monotonicity constraints
#
#
#-----------------------------------------------------------------------------------------------------------------
  Method.MC.monotone <- function(N.calls){ 

    if(ndim > 2){
      UU <- runif(ndim*10^(ordre.p + 2))
      UU <- matrix(UU, ncol=ndim, byrow = TRUE)
      D.UU <- dim(UU)[1]
    }

    NN    <- 0        #Number of call to f
    N.tot <- 0
    res   <- NULL

    Z.safe <- NULL
    Z.fail <- NULL

    is.Call <- 0    
    while(NN < N.calls){
      if(silent == FALSE){
        if(N.tot%%100 == 0){print(NN);flush.console();}
      }

      U <- runif(ndim)

      N.tot <- N.tot + 1

      if( is.null(Z.safe)& is.null(Z.fail) ){
        t.u <- G(U)
        NN  <- NN + 1

        is.Call[NN] <- N.tot
      }

      if(is.null(Z.safe)&( !is.null(Z.fail)) ){
        ttf <- is.dominant(Z.fail, U, ndim, set = 2)
        if( sum(ttf) == 0 ){
          t.u <- G(U)
          NN  <- NN + 1

          is.Call[NN] <- N.tot          
        }else{
          t.u <- -1
        }
      }
      
      if(!is.null(Z.safe) & is.null(Z.fail) ){       
        tts <- is.dominant(Z.safe, U, ndim, set = 1)
        if(sum(tts) == 0){
          t.u <- G(U)
          NN  <- NN + 1
          is.Call[NN] <- N.tot          
        }else{
          t.u <- 1
        }
      }
    
      if((!is.null(Z.safe)) &( !is.null(Z.fail)) ){      
        ttf <- is.dominant(Z.fail, U, ndim, set = 2)       
        tts <- is.dominant(Z.safe, U, ndim, set = 1)
        if( (sum(tts) == 0) & (sum(ttf) == 0) ){
          t.u <- G(U)
          NN  <- NN + 1
          is.Call[NN] <- N.tot
        }
        if( (sum(tts) == 0)& (sum(ttf) != 0) ){
          t.u <- -1
        }
        if( (sum(tts) != 0)& (sum(ttf) == 0) ){
          t.u  <- 1
        }
        
      }
      
      if(t.u <= 0){
        Z.fail <- rbind(Z.fail, U)
        res    <- c(res, 1)
      }else{        
        Z.safe <- rbind(Z.safe, U)
        res    <- c(res, 0)
      }

    }

    I <- 1:N.tot

    alpha <- 0.05
 
    cum.res <- cumsum(res)

    estimation_MC <- cum.res/I                                        #Monte Carlo Estimator
    Var_MC        <- (estimation_MC)*(1 - estimation_MC)/I            #Variance of the estimator
    IC.inf        <- estimation_MC - qnorm(1 - alpha/2)*sqrt(Var_MC)  #Confidence Interval
    IC.sup        <- estimation_MC + qnorm(1 - alpha/2)*sqrt(Var_MC)  #Confidence Interval
    CV_MC         <- 100*sqrt(Var_MC)/estimation_MC

    if(is.null(Z.fail)){
      Um <- 0
    }else{
      ZF <- Frontier(Z.fail, -1)        
      Um <- Volume.bounds(UU, ZF, -1)    
    }

    if(is.null(Z.safe)){
      UM <-1
    }else{ 
      ZS <- Frontier(Z.safe, set = 1) 
      UM <- Volume.bounds(UU, ZS, 1)  
    }
    
    return(list(cbind(IC.inf, IC.inf, estimation_MC, CV_MC, Var_MC)[is.Call, ], Um, UM, N.tot))
  }


#-----------------------------------------------------------------------------------------------------------------
#	               	              Maximum Likelihood estimator
#-----------------------------------------------------------------------------------------------------------------

    # p.k = (p.k^-, p.k^+)
    # p : unknow
    # signature[k] = 1 si H(Y[k,]) <= 0, 0 otherwise

    log.likehood <- function(p , p.k, signature){
      gamma <- (p - p.k[,1])/(p.k[,2] - p.k[,1])
      u     <- (gamma^signature)*((1 - gamma)^(1 - signature))
      return(prod(u))
    }

#-----------------------------------------------------------------------------------------------------------------
    as.binary <- function (x) { 
      base <- 2;
      r <- numeric(ndim)
      for (i in ndim:1){ 
        r[i] <- x%%base 
	  x  <- x%/%base
      } 
      return(r) 
   }


#-----------------------------------------------------------------------------------------------------------------
#		                             To simulate 1 point uniformy aroud "x"
#-----------------------------------------------------------------------------------------------------------------
    
    SIM <- function(x, W){

      B <- 0
      # One split the space around "x" in (2^d - 1) set, and one compute the volume of each set
      B <- apply( matrix(W, ncol = 1), 
                  MARGIN = 1, 
                  function(v){
                    Z <- as.binary(v)
                    v <- 0
                    u <- 0
                    for(j in 1:ndim){
                      u[j] <- ifelse(Z[j] == 0, 1 - x[j], x[j])
                    }     
                    return(prod(u))  
                  }
                 )

      B <- cumsum(B)
      U <- runif(1, 0, max(B))

      #One set is choose randomly
      pos <- ifelse(U < B[1], 1, which.max(B[B <= U]) + 1)

      Z <- as.binary(W[pos])
      A <- 0
      
      #One sample in that set
      for(i in 1:ndim){
        A[i] = ifelse(Z[i] == 0, runif(1, x[i], 1), runif(1, 0, x[i]))  
      }
      return(A)      
    }

#-----------------------------------------------------------------------------------------------------------------
#                   Sample CP points in the non dominated space
#-----------------------------------------------------------------------------------------------------------------

    Sim.non.dominated.space <- function(CP, Z.safe, Z.fail, W){
      CP1 <- 0;
      Y   <- NULL
      Y.temp  <- apply(1 - Z.safe, MARGIN = 1, prod)
      Y.temp1 <- Z.safe[which.max(Y.temp),]
      while(CP1 < CP){
        Y.temp2 <- SIM(Y.temp1, W)
        tts1 <- is.dominant(Z.safe, Y.temp2, ndim, 1)
        ttf1 <- is.dominant(Z.fail, Y.temp2, ndim, 2)
        if( (sum(tts1) == 0 ) & ( sum(ttf1) == 0) ){
          Y <- rbind(Y,Y.temp2)
          CP1 <- CP1 + 1
        } 
      }
      return(Y)
    }

#-----------------------------------------------------------------------------------------------------------------
#
#
#
# 						            METHOD : MRM
#
#
#
#-----------------------------------------------------------------------------------------------------------------

  Choix.method.1.1 <- function(N.calls, H){


    #Create of the sample to compute the bounds
    if(ndim > 2){
      UU   <- runif(ndim*10^(ordre.p + 2))
      UU   <- matrix(UU, ncol=ndim, byrow = TRUE)
      D.UU <- dim(UU)[1]
    }

    if(ndim == 2){UU <- 0; D.UU <- 0}
    
    V <- list()

    V <- Intersect(ndim, H)    

    list.set <- 0
    for(i in 1:length(V)){
      list.set[i] <- V[[i]][[2]]
    }

    u.dep   <- list()
    u.other <- list()

    u.dep[[1]] <- V[[max(which(list.set == -1))]][[1]]
    u.dep[[2]] <- V[[max(which(list.set == -1))]][[2]]

    u.other[[1]] <- V[[max(which(list.set == 1))]][[1]]
    u.other[[2]] <- V[[max(which(list.set == 1))]][[2]]

    Z.fail <- t(as.matrix(u.dep[[1]]))
    Z.safe <- t(as.matrix(u.other[[1]]))

    cp     <- length(V)  
    
    um <- 0
    uM <- 1

    Um  <- 0
    UM  <- 1

    eps   <- 1e-7
    alpha <- 0.05

    SIGN   <- 0
    ICinf  <- 0
    ICsup  <- 0
    VAR    <- 0
    CV.MLE <- 0
    MLE    <- 0

    ZS <- NULL
    ZF <- NULL

    um <- prod(V[[cp]][[1]])
    uM <- 1 - prod(1 - V[[cp]][[1]])

    j  <- 1
    Um <- um
    UM <- uM
    W  <- 1:(2^(ndim) - 1)

    while(cp < N.calls){
      if(silent == FALSE){
        print(paste("compteur =",cp));flush.console()
      }
      uu <-  Sim.non.dominated.space (1, Z.safe, Z.fail, W)
#########################################################################
#       Appel au code
#########################################################################
      H.u <- H(uu)
      SIGN[j] <- (1-sign(H.u))/2

      if(H.u > 0){
        Z.safe.old <- rbind(uu, Z.safe)
        ss     <- is.dominant(Z.safe, uu, ndim, 2)
        Z.safe <- Z.safe[which(ss == FALSE), ]
        Z.safe <- rbind(uu, Z.safe)

        vol <- Volume.bounds(UU,Z.safe, 1)

         Um[j+1] <- Um[j]
         if(vol >= UM[j]){
           UM[j + 1] <- UM[j]
         }else{
           UM[j + 1] <- vol
         }
         CC <- ifelse(cp == N.calls, 1, 0)
         ZS <- rbind(ZS, uu)

      }else{
        Z.fail.old <- rbind(uu, Z.fail)

        ff     <- is.dominant(Z.fail, uu, ndim, 1)
        Z.fail <- Z.fail[which(ff == FALSE), ]
        Z.fail <- rbind(uu, Z.fail)

        vol <- Volume.bounds(UU, Z.fail, 2)

        UM[j+1] <- UM[j]
        if(vol <= Um[j]){
           Um[j+1] <- Um[j]
        }else{
          Um[j+1] <- vol
        }

        ZF <- rbind(ZF, uu)
      }

      cp <- cp + 1

      MLE.test <- optimize(f = log.likehood,
                               interval = c(Um[j],UM[j]),
                               maximum = TRUE,
                               signature = SIGN,
                               p.k = cbind(Um[1:j], UM[1:j])
                               )   
      MLE[j] <- as.numeric(MLE.test[1])

      VAR <- sum( 1/((MLE - Um[1:j])*(UM[1:j]- MLE)))

      bn  <- 1/VAR
      an  <- eps*VAR^(5/2)/abs( sum( 1/((MLE + eps - Um[1:j])*(UM[1:j] - MLE - eps))) - sum( 1/((MLE - Um[1:j])*(UM[1:j]- MLE)))  )

      ICinf[j]  <- MLE[j] - qnorm(1 - alpha/2)/sqrt(VAR - alpha/an)
      ICsup[j]  <- MLE[j] + qnorm(1 - alpha/2)/sqrt(VAR + alpha/an) 

      CV.MLE[j] <- 100/(sqrt(VAR)*MLE[j])   

      j <- j + 1 

    } #end of "while cp < N.calls"

    RR <- cbind(Um[1:(j-1)], UM[1:(j-1)], MLE[1:j-1], ICinf[1:(j-1)] , ICsup[1:(j-1)], CV.MLE[1:(j-1)])
    return(RR)
  }

#-----------------------------------------------------------------------------------------------------------------
#
#
#
#
# 					               FIN	METHODE 1.1
#
#
#
#-----------------------------------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------------------------------
#
#
#
#
# 					            	FIN DE LA PARTIE 2
#
#
#
#-----------------------------------------------------------------------------------------------------------------  
  if(Method == "MRM"){
    RESULT <- Choix.method.1.1(N.calls, G)
  }


  if(Method == "MC_monotone"){
    RESULT <- Method.MC.monotone(N.calls)
  }

  return(RESULT)

}
#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------
#						Fin de MRM	
#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------


