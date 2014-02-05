## -----------------------------------------------------------------------------
## Fonction AK-MCS
## -----------------------------------------------------------------------------
##    Copyright (C) 2013
##    Developpement initial : C. WALTER
##    Modifications         : G. DEFAUX
##    CEA
## -----------------------------------------------------------------------------

AKMCS = function(dimension,
			limit_state_function,
			N    = 500000,          # Monte-Carlo population size
			N1   = 10*dimension,    # Size of the first DOE
			Nmax = 200,             # Maximum number of calls to the LSF
			learn_db  = NULL,       # Coordinates of alredy known points
			lsf_value = NULL,       # Value of the LSF on these points
			failure   = 0,          # Failure threshold
			precision = 0.05,       # Maximum desired cov on the Monte-Carlo estimate
			meta_model = NULL,      # Provide here a kriging metamodel from km if wanted
			kernel = "matern5_2",   # Specify the kernel to use for km
			learn_each_train = FALSE,# Specify if kernel parameters are re-estimated at each train
      crit_min = 2,           # Minimum value of the criteria to be used for refinement
			limit_fun_MH = NULL,    # Define an area of exclusion with a limit function, eg in metaSS
			sampling_strategy = "MH",# Either MH for Metropolis-Hastings of AR for accept-reject
			first_DOE = "Gaussian", # Either Gaussian or Uniform, to specify the population on which clustering if done
			seeds = NULL,           # If some points are already known to be in the appropriate subdomain, eg in metaSS
			seeds_eval = NULL,      # Value of the metamodel on these points
			burnin = 30,            # Burnin parameter for MH
			thinning = 4,           # Thinning parameter for MH
			plot = FALSE,           # Set to TRUE for a full plot, ie refresh at each iteration
			limited_plot = FALSE,   # Set to TRUE for a final plot with final DOE, metamodel and LSF
			add = FALSE,            # If plots are to be added to a current device
			output_dir = NULL,      # If plots are to be saved in jpeg in a given directory
			z_MH = NULL,            # For plots, if metamodel has already been evaluated on the grid
			z_lsf = NULL,           # For plots, if LSF has already been evaluated on the grid
      verbose = 0) {          # Either 0 for almost no output, 1 for medium size output and 2 for all outputs
  

cat("==========================================================================================\n")
cat("                              Beginning of AK-MCS algorithm\n")
cat("==========================================================================================\n\n")


## STEP 0 : INITIALISATION

Ncall  = 0
cov    = Inf
U      = list(N=NULL,  N1=NULL)
G      = list(N1=NULL, g=lsf_value)
z_meta = NA
STOP   = 0
Nfailure = 0

# plotting part
if(plot == TRUE){
  if(verbose>0){cat(" * 2D PLOT : SET UP \n\n")}
	x_plot = c(-80:80)/10
	y_plot = c(-80:80)/10
	z = rbind(rep(x_plot,length(y_plot)),sort(rep(y_plot,length(x_plot))))
	if(is.null(z_lsf)) {
		z_lsf = outer(x_plot,y_plot,function(x,y){z = cbind(x,y); apply(z,1,limit_state_function)})
	}
	if(add == FALSE){
		if(is.null(output_dir)) {dev.new("x11",title="AKMCS")}
		else{
			fileDir = paste(output_dir,"_AKMCS.png",sep="")
			png(fileDir)
		}
	}
	par(pty="s")
	plot(x_plot, y_plot, xlab="x", ylab="y", type="n")

	tryCatch(points(learn_db[1,],learn_db[2,],col=2,pch=3))
	if(!is.null(limit_fun_MH)){
		if(is.null(z_MH)) { z_MH = outer(x_plot,y_plot,function(x,y){z = cbind(x,y); z=t(z);limit_fun_MH(z)$mean})}
		contour(x_plot,y_plot,z_MH,level=0,labels="Subset limit",method="edge",add=TRUE)
  }
	contour(x_plot, y_plot, z_lsf, level=failure, labels=paste("LSF=",failure,sep=""), method="edge", add=TRUE)
}


while(cov>precision){

if(Nfailure==0){
	cat(" ============================================= \n")
	cat(" STEP 1 : GENERATION OF THE WORKING POPULATION \n")
	cat(" ============================================= \n\n")
	

	if(is.null(limit_fun_MH)){
		if(verbose>0){cat(" * Generate N =",N,"standard Gaussian samples\n\n")}
		U$N = matrix(rnorm(dimension*N, mean=0, sd=1),dimension,N)
	}
	else{
		if(sampling_strategy=="MH"){
		  if(verbose>0){cat(" * Generate N =",N,"points from",dim(seeds)[2],"seeds with Metropolis-Hastings algorithm\n")}
			gen_pop = generateWithlrmM(seeds=seeds,
						  seeds_eval=seeds_eval,
						  N=N,
						  limit_f=limit_fun_MH,
						  burnin=burnin,
						  thinning=thinning)
			U$N = gen_pop$points
		}
		else{
		  if(verbose>0){cat(" * Generate the N =",N,"Monte-Carlo population with an accept-reject strategy on a standard Gaussian sampling\n")}
			rand = function(dimension,N) {matrix(rnorm(dimension*N),dimension,N)}
			U$N = generateWithAR(dimension=dimension,
					    N=N,
					    limit_f=limit_fun_MH,
					    rand=rand)
		}
	}


	cat(" ================== \n")
	cat(" STEP 2 : FIRST DoE \n")
	cat(" ================== \n\n")

	switch(first_DOE,
		Gaussian = {
		  if(verbose>0){cat(" * Get N1 =",N1,"points by clustering of the N =",N,"points\n")}
			U$N1 = t(kmeans(t(U$N), centers=N1,iter.max=20)$centers)
		},
		Uniform = {
		  if(verbose>0){cat(" * Get N1 =",N1,"points with a uniform sampling in a hyper sphere of radius max(radius(N points)) \n")}
			radius = max(apply(U$N,2,function(x){sqrt(x%*%x)}));
			U$N1 = t( kmeans( t(runifSphere(dimension,10000*dimension,radius)), centers=N1, iter.max=20)$centers )
		},
		stop("Wrong first DOE sampling strategy\n")
	)
	
	if(verbose>0){cat(" * Assessment of g on these points\n")}
	G$N1 = apply(U$N1,2,limit_state_function);Ncall = Ncall + N1
	    
	if(plot==TRUE){
	  if(verbose>0){cat(" * 2D PLOT : First DoE (Red +) \n")}
			if(is.null(seeds)) {
				radius = apply(U$N,2,function(x) {sqrt(sum(x^2))})
				Ru = max(radius)
				rm(radius)
				symbols(0,0, circles=Ru, inches=FALSE, add=TRUE)
			}
			else {contour(x_plot,
				y_plot,
				z_MH,
				level=0,
				labels="Subset limit",
				method="edge",
				add=TRUE)
			}
			points(U$N1[1,], U$N1[2,], pch=3, col='red')
	}

	if(verbose>0){cat(" * Add points to the learning database\n")}
	if(is.null(learn_db)){
		learn_db = U$N1
		G$g = G$N1
	}
	else{
		learn_db = cbind(learn_db,U$N1)
		G$g = c(G$g,G$N1)
	}
}
if(verbose>0){cat(" * Train the model :\n")}
	if(is.null(meta_model) || learn_each_train==TRUE || Nfailure>0) {
	  if(verbose>1){cat("    - Learn hyperparameters !!! \n")}
		meta = trainModel(design   = learn_db,
						  response = (G$g-failure),
						  kernel   = kernel,
						  type="Kriging")
    Nfailure = 0
	}
	else {
	  if(verbose>1){cat("    - Use previous hyperparameters !!! \n")}
		meta = trainModel(meta_model,
						  updesign=U$N1,
						  upresponse=(G$N1-failure),
						  type="Kriging")
	}
	
	meta_model = meta$model
	meta_fun   = meta$fun
	
	#plotting part
	if(plot == TRUE){
	  if(verbose>0){cat(" * 2D PLOT : FIRST APPROXIMATED LSF USING KRIGING \n")}
		z_meta = meta_fun(z)
		z_meta$mean = matrix(z_meta$mean,length(x_plot),length(y_plot))
		z_meta$sd   = matrix(z_meta$sd,length(x_plot),length(y_plot))
		z_crit      = abs(z_meta$mean)/z_meta$sd
		contour(x_plot, y_plot,
			z_meta$mean,
			level=0,
			labels="Metamodel",
			method="edge",
			col = 4, add=TRUE)
		contour(x_plot,y_plot,
			z_crit,
			level=crit_min,
			labels=paste("U =",crit_min),
			method="edge",
			lty="dashed",add=TRUE,col=4)
	}

if(verbose>0){cat(" * Evaluate the meta-model on the work population\n")}
	meta_pred = meta_fun(U$N)
	G_meta    = meta_pred$mean
	sd        = meta_pred$sd
	
if(verbose>0){cat(" * Evaluate criterion on the work population\n")}
	criterion = abs(G_meta)/sd
	minC = min(criterion)
if(verbose>1){cat("    - minimum value of the criterion =",minC,"\n\n")}


	cat(" ======================= \n")
	cat(" STEP 3 : UPDATE THE DoE \n")
	cat(" ======================= \n")
	
	k = 0;
	while(minC<crit_min & k<Nmax) {
				
		k = k+1;
		if(verbose>0){cat("\n * ITERATION ",k,"\n")}
		if(verbose>1){cat("   - min < 2 & k =",k,"< Nmax =",Nmax,"=> improve the model\n")}
		
		candidate = matrix(U$N[,which.min(criterion)],ncol=1)
		eval      = limit_state_function(candidate);Ncall = Ncall + 1
		learn_db  = cbind(learn_db,candidate)
		G$g       = c(G$g,eval)

		#plotting part
		if(plot==TRUE){
		  if(verbose>1){cat("   - 2D PLOT : UPDATE \n")}
			points(candidate[1], candidate[2], pch=21, col='blue3')
		}

		if(verbose>1){cat("   - Train the model\n")}
		if(learn_each_train==TRUE) {
		  if(verbose>1){cat("   - Learn hyperparameters !!! \n")}
			meta = trainModel(design=learn_db,
							  response=(G$g-failure),
							  kernel=kernel,
							  type="Kriging")
      Nfailure = 0
		}
		else {
		  if(verbose>1){cat("     + Use previous hyperparameters !!! \n")}
			meta = trainModel(meta_model,
							  updesign=candidate,
							  upresponse=(eval-failure),
							  type="Kriging")
		}
		meta_model = meta$model
		meta_fun   = meta$fun
	
		if(verbose>1){cat("   - Evaluate the meta-model on the work population\n")}
		meta_pred = meta_fun(U$N)
		G_meta    = meta_pred$mean
		sd        = meta_pred$sd
		
		if(verbose>1){cat("   - Evaluate criterion on the work population\n")}
		criterion = abs(G_meta)/sd
		minC      = min(criterion)
		cat("     + minimum value of the criterion =",minC,"\n")

		#plotting part
		if(plot==TRUE){
		  if(verbose>0){cat("   - 2D PLOT : END ITERATION",k,"\n")}
			z_meta = meta_fun(z)
			z_meta$mean = matrix(z_meta$mean,length(x_plot),length(y_plot))
			z_meta$sd = matrix(z_meta$sd,length(x_plot),length(y_plot))
			z_crit = abs(z_meta$mean)/z_meta$sd
			plot(x_plot,y_plot,xlab="x",ylab="y",type="n")
			contour(x_plot,y_plot,
				z_lsf,
				level=failure,
				labels=paste("LSF=",failure,sep=""),
				method="edge",add=TRUE)
			contour(x_plot,y_plot,
				z_meta$mean,
				level=0,
				labels="Metamodel",
				method="edge",
				add=TRUE,
				col=4)
			contour(x_plot,y_plot,
				z_crit,
				level=crit_min,
				labels=paste("U =",crit_min),
				method="edge",
				lty="dashed",
				add=TRUE,
				col=4)
			if(is.null(seeds)) {symbols(0,0,circles=Ru,inches=F,add=TRUE)}
			else {contour(x_plot,y_plot,z_MH,level=0,labels="Subset limit",method="edge",add=TRUE)}
			points(learn_db[1,],learn_db[2,],pch=3,col=2)
		}
	}

	cat(" ======================================================================================= \n")
	cat(" STEP 4 : EVALUATE FAILURE PROBABILITY WITH A MONTE-CARLO ESTIMATOR USING THE META-MODEL\n")
	cat(" ======================================================================================= \n")
	
	P = 1/N*sum(G_meta<0)
if(verbose>0){cat(" * Pf =",P,"\n")}
	cov = sqrt((1-P)/(N*P))
if(verbose>0){cat(" * cov =",cov,"\n")}
	
	if( cov > precision) {
		if(P>0){
			N = ceiling((1-P)/(precision^2*P))
			cat("   => cov too large ; this order of magnitude for the probability brings N =",N,"\n")
		}
		else {
			Nfailure = sum(G$g<failure)
			cat("   => cov too large, only",Nfailure,"failings points in the learn_db\n")
			ind   = which.min(G$g)
			seeds_DOE = learn_db[,ind]
			eval_seeds_DOE = G$g[ind]
			cat(" * Trying to find failing points with remaining calls to the LSF \n")
			gen_pop  = getSeedBelowThreshold(seeds_DOE,
											 eval_seeds_DOE,
											 threshold=failure,
											 limit_f=limit_state_function,
											 niter_max=Nmax-Ncall)
			learn_db = cbind(learn_db,gen_pop$samples)
			G$g = c(G$g,gen_pop$eval_samples)
			Ncall = Ncall + gen_pop$Ncall
			Nfailure = sum(G$g<failure)
			if(Nfailure==0) {stop("   - No failing points found and Ncall =",Ncall,"= Nmax =",Nmax,"\n")}
			else {
        cat("   -",Nfailure,"failings points found in",gen_pop$Ncall,"call(s) ; return to the beginning of the algorithm \n")
        if(plot==TRUE) {
          cat("   - 2D PLOT : UPDATE\n")
          points(learn_db[1,],learn_db[2,],pch=3,col=2)
        }
			}
		}
	}
  else{
    cat(" * cov < precision =",precision,"; End of AKMCS algorithm")
    cat("   - Pf =",P,"\n")
    cat("   - cov =",cov,"\n")
    cat("   - Ncall =",Ncall,"\n")
  }
}

#plotting part
if(limited_plot==TRUE){
cat("\n * 2D PLOT : LSF, FINAL DATABASE AND METAMODEL")

	x_plot = c(-80:80)/10
	y_plot = c(-80:80)/10
	if(is.null(z_lsf)) {z_lsf = outer(x_plot,y_plot,function(x,y){z = cbind(x,y); apply(z,1,limit_state_function)})}
	if(add==FALSE){
		if(is.null(output_dir)) {dev.new("x11",title="AKMCS")}
		else{
			fileDir = paste(output_dir,"_AKMCS.jpeg",sep="")
			jpeg(fileDir)
		}
		par(pty="s")
		plot(x_plot,y_plot,type="n")
	}

	z = rbind(rep(x_plot,length(y_plot)),sort(rep(y_plot,length(x_plot))))
	z_meta = meta_fun(z)
	z_meta$mean = matrix(z_meta$mean,length(x_plot),length(y_plot))
	z_meta$sd = matrix(z_meta$sd,length(x_plot),length(y_plot))
	z_crit = abs(z_meta$mean)/z_meta$sd
	contour(x_plot,y_plot,z_meta$mean,level=0,labels="Metamodel",method="edge",add=TRUE,col=4)
	contour(x_plot,y_plot,z_crit,level=crit_min,labels=paste("U =",crit_min),method="edge",lty="dashed",add=TRUE,col=4)
	contour(x_plot,y_plot,z_lsf,level=failure,labels=paste("LSF=",failure,sep=""),method="edge",add=TRUE)
	points(learn_db[1,],learn_db[2,],col=2,pch=3)
}

if( (plot || limited_plot) & (add==FALSE) & !is.null(output_dir) ) { dev.off() }

points = U$N[,(G_meta<0)]
meta_eval = G_meta[G_meta<0]

if(plot+limited_plot) {
	res = list(proba=P, 
			   cov=cov,
			   Ncall=Ncall,
			   learn_db=learn_db,
			   lsf_value=G$g,
			   meta_fun=meta_fun,
			   meta_model=meta_model,
			   points=points,
			   meta_eval=meta_eval,
			   z_meta=z_meta$mean); }
else {res = list(proba=P,
				 cov=cov,
				 Ncall=Ncall,
				 learn_db=learn_db,
				 lsf_value=G$g,
				 meta_fun=meta_fun,
				 meta_model=meta_model,
				 points=points,
				 meta_eval=meta_eval); }

return(res)

}
