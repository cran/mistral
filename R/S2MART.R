## -----------------------------------------------------------------------------
## Fonction S2MART
## -----------------------------------------------------------------------------
##    Copyright (C) 2013
##    Developpement : C. WALTER
##    CEA
## -----------------------------------------------------------------------------

S2MART = function(dimension,limit_state_function,
			## Algorithm parameters
			Nn = 100,		          # Number of samples to evaluate the quantiles in the subset step
			alpha_quantile = 0.1,	# Cutoff probability for the subsets
			## Meta-model choice

			...,			            # All others parameters of the metamodel based algorithm
			## Plot information
			plot=FALSE,		
			limited_plot=FALSE,
			output_dir=NULL,
      verbose = 0) {        # Either 0 for almost no output, 1 for medium size output and 2 for all outputs

cat("==========================================================================================\n")
cat("                 Beginning of Metamodelisation on Subset Simulation algorithm \n")
cat("==========================================================================================\n")

	i = 0;			# Subset number
	y0 = Inf;		# first level
	y = NA;			# will be a vector containing the levels
	meta_fun = list(NA);	# initialise meta_fun list, that will contain limit_state_function approximation at each iteration
	z_meta = list(NA)	# initialise z_meta list, that will contain evaluation of the surrogate model on the grid
	P = 1;			# probability
	Ncall = 0;		# number of calls to the limit_state_function
	delta2 = 0;		# theoretical cov
	z_lsf = NULL		# For plotting part

	if(plot==TRUE && !is.null(output_dir)) {output_dir = paste(output_dir,"_S2MART",sep="")}

	#Define the list variables used in the core loop
	U = list(Nn=matrix(nrow=dimension,ncol=Nn))
	#G stands for the value of the limit state function on these points
	G = list(g=NA,#value on learn_db points, ie all the value already calculated at a given iteration
		Nn=NA*c(1:Nn))
	
	#beginning of the core loop
	while (y0>0) {

		i = i+1;
		cat("\n SUBSET NUMBER ",i,"\n")
		cat(" -------------\n\n")

		if(plot){
			if(verbose>0){cat(" * 2D PLOT : SET-UP \n")}
			if(is.null(output_dir)){dev.new("x11",title=paste("S2MART subset number",i))}
 			else{fileDir = paste(output_dir,"_Subset",i,".jpeg",sep="");jpeg(fileDir)}
			par(pty="s")
			x_plot = c(-80:80)/10
			y_plot = c(-80:80)/10
			z_lsf = outer(x_plot,y_plot,function(x,y){z = cbind(x,y); apply(z,1,limit_state_function)})
			plot(x_plot,y_plot,xlab="x",ylab="y",type="n")
		}

		cat("===========================================================================\n")
		cat(" STEP 1 : Subset Simulation part \n")
		cat("===========================================================================\n\n")
	
		if(i==1){
		  if(verbose>0){cat(" * Generate Nn =",Nn," standard gaussian samples\n")}
			U$Nn = matrix(rnorm(dimension*Nn, mean=0, sd=1),dimension,Nn)
		}
		else {
		  if(verbose>0){cat(" * Generate Nn = ",Nn," points from the alpha*Nn points lying in F(",i-1,") with MH algorithm\n",sep="")}
			U$Nn = generateWithlrmM(seeds=U$Nn[,G$Nn<y0],seeds_eval=G$Nn[G$Nn<y0],N=Nn,limit_f=meta_fun[[i-1]])$points
		}
	
		#assessment of g on these points
		if(verbose>0){cat(" * Assessment of the LSF on these points\n")}
		G$Nn = apply(U$Nn,2,limit_state_function);Ncall = Ncall + Nn
	
		#Determination of y[i] as alpha-quantile of Nn points Un=U$Nn
		if(verbose>0){cat(" * Determination of y[",i,"] as alpha-quantile of these samples\n",sep="")}
		y0 <- y[i] <- getQuantile(data=G$Nn,alpha=alpha_quantile)

	
		#Add points U$Nn to the learning database
		if(verbose>0){cat(" * Add points U$Nn to the learning database\n\n")}
		if(i==1) {
			learn_db = cbind(seq(0,0,l=dimension),U$Nn)
			g0 = limit_state_function(seq(0,0,l=dimension));Ncall = Ncall + 1;
			G$g = c(g0,G$Nn)
		}
		else {
			learn_db = cbind(learn_db,U$Nn)
			G$g = c(G$g,G$Nn)
		}



		cat("===========================================================================\n")
		cat(" STEP 2 : Metamodel algorithm part \n")
		cat("===========================================================================\n\n")
		if(i==1){
			arg = list(dimension=dimension,
					limit_state_function=limit_state_function,
					failure = y0,
					learn_db = learn_db,
					lsf_value = G$g,
					plot = plot,
					z_lsf = z_lsf,
					add = TRUE,
					output_dir = output_dir,...)
			    meta_step = do.call(SMART,arg)
		}
		else{
			arg = list(dimension=dimension,
					limit_state_function = limit_state_function,
					failure = y0,
					learn_db = learn_db,
					lsf_value = G$g,
					seeds = seeds,
					seeds_eval = seeds_meta,
					limit_fun_MH = meta_fun[[i-1]],
					z_MH = z_meta[[i-1]],
					z_lsf = z_lsf,
					plot = plot,
					add = TRUE,
					output_dir = output_dir,...)
				  meta_step = do.call(SMART,arg)
		}

		if(!is.null(output_dir)){
		  if(verbose>0){cat("\n * 2D PLOT : CLOSE DEVICE \n")}
			dev.off()
		}

		P = P*meta_step$proba
		delta2 = tryCatch(delta2 + (meta_step$cov)^2,error = function(cond) {return(NA)})
		Ncall = Ncall + meta_step$Ncall
		learn_db = meta_step$learn_db
		G$g = meta_step$lsf_value
		meta_fun[[i]] = meta_step$meta_fun
		meta_model = meta_step$meta_model
		seeds = meta_step$points
		seeds_meta = meta_step$meta_eval
		if(plot) {z_meta[[i]] = meta_step$z_meta}

		if(y0>0) {cat("\n * Current threshold =",y0,"> 0 => start a new subset\n")
			  cat("   - Current probability =",P,"\n")
			  cat("   - Current number of call =",Ncall,"\n")}
		else {cat("\n * Current threshold =",y0,"=> end of the algorithm\n")
		      cat("   - Final probability =",P,"\n")
		      cat("   - Total number of call =",Ncall,"\n")}
	}
	
	if(plot+limited_plot) {
	  if(verbose>0){cat("\n 2D PLOT : FINAL PLOT \n")}
		if(is.null(output_dir)){dev.new("x11",title="Models comparison")}
		else{fileDir = paste(output_dir,"_Models_comparison.jpeg",sep="");jpeg(fileDir)}
		plot(x_plot,y_plot,xlab="x",ylab="y",type="n")
		points(learn_db[1,],learn_db[2,],col=2,pch=3)

		if(limited_plot==TRUE){
			x_plot = c(-80:80)/10
			y_plot = c(-80:80)/10
			z = rbind(rep(x_plot,length(y_plot)),sort(rep(y_plot,length(x_plot))))
			z_lsf = outer(x_plot,y_plot,function(x,y){z = cbind(x,y); apply(z,1,limit_state_function)})
			for(level in c(1:i)){
				z_meta[[level]] = meta_fun[[level]](z)$mean
			}
		}
		for(j in c(1:i)) {
			contour(x_plot,y_plot,z_meta[[j]],level=0,labels=paste("Meta[",j,"] classifier",sep=""),col=j,method="edge",add=TRUE)
			contour(x_plot,y_plot,z_lsf,level=y[j],col=j,method="edge",add=TRUE)
		}

		if(!is.null(output_dir)){dev.off()}
	}
	
	res = list(proba=P,
             cov=sqrt(delta2),
             Ncall=Ncall,
             learn_db=learn_db,
             lsf_value=G$g,
             meta_model=meta_model
             );
	
	return(res)
}
