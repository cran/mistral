## -----------------------------------------------------------------------------
## Fonction MetaIS
## -----------------------------------------------------------------------------
##    Copyright (C) 2013
##    Developpement : C. WALTER
##    CEA
## -----------------------------------------------------------------------------

MetaIS = function(dimension, limit_state_function,
			## Algorithm parameters
			N = 500000,             # size of Monte-Carlo population for P_epsilon estimate
			N_alpha = 100,          # initial size of Monte-Carlo population for alpha estimate
			N_DOE = 2*dimension,    # size of initial DOE got by clustering of the N1 samples
			N1 = N_DOE*30,          # size of the initial uniform population sampled in a hypersphere of radius Ru
			Ru = 8,                 # radius of the hypersphere for the initial sampling
			Nmin = 30,              # minimum number of call for the construction step
			Nmax = 200,             # maximum number of call for the construction step
			Ncall_max = 1000,       # maximum number of call for the whole algorithm
			precision = 0.05,       # desired maximal value of cov
			N_seeds = 1,            # number of seeds for MH algoritm while generating into the margin (according to MP*gauss)
			Niter_seed = 10000,     # maximum number of iteration for the research of a seed for alphaLOO refinement sampling
			N_alphaLOO = 5000,      # number of points to sample at each refinement step
			K_alphaLOO = 2*dimension, #number of clusters at each refinement step
			alpha_int = c(0.1,10),  # range for alpha to stop construction step
			k_margin = 1.96,        # margin width ; this means points are classified with more than 97,5%
			## Subset parameters
      learn_db  = NULL,       # Coordinates of alredy known points
      lsf_value = NULL,       # Value of the LSF on these points
      failure   = 0,          # Failure threshold
      meta_model = NULL,      # Provide here a kriging metamodel from km if wanted     
      kernel = "matern5_2",   # Specify the kernel to use for km
      learn_each_train = FALSE,# Specify if kernel parameters are re-estimated at each train
      limit_fun_MH = NULL,    # Define an area of exclusion with a limit function, eg in metaSS
      sampling_strategy = "MH",# Either MH for Metropolis-Hastings of AR for accept-reject
      seeds = NULL,           # If some points are already known to be in the appropriate subdomain, eg in metaSS
      seeds_eval = NULL,      # Value of the metamodel on these points
      burnin = 20,            # Burnin parameter for MH
      thinning = 4,           # Thinning parameter for MH
      ## plot parameter
      plot = FALSE,           # Set to TRUE for a full plot, ie refresh at each iteration
      limited_plot = FALSE,   # Set to TRUE for a final plot with final DOE, metamodel and LSF
      add = FALSE,            # If plots are to be added to a current device
      output_dir = NULL,      # If plots are to be saved in jpeg in a given directory
      z_MH = NULL,            # For plots, if metamodel has already been evaluated on the grid
      z_lsf = NULL,           # For plots, if LSF has already been evaluated on the grid
      verbose = 0) {          # Either 0 for almost no output, 1 for medium size output and 2 for all outputs


cat("==========================================================================================\n")
cat("                              Beginning of Meta-IS algorithm \n")
cat("==========================================================================================\n\n")

cat("===========================================================================\n")
cat(" STEP 1 : Adaptative construction of h the approximated optimal density \n")
cat("===========================================================================\n\n")

## Init
Ncall = 0
cov_epsilon = Inf
U = list(N=NULL,N1=NULL,N_DOE=NULL)
G = list(N=NA*c(1:N),N_DOE=NULL,g=lsf_value)
G_meta = list(N=NULL,N1=NULL)
z_meta = NA;
ITER = 0;


## plotting part
if(plot == TRUE || limited_plot == TRUE){
	if(verbose>0){cat(" * 2D PLOT : SET-UP \n")}
	x_plot = seq(-Ru,Ru,l=20*Ru+1)
	y_plot = seq(-Ru,Ru,l=20*Ru+1)
	z = rbind(rep(x_plot,length(y_plot)),sort(rep(y_plot,length(x_plot))))
	z_lsf = outer(x_plot,y_plot,function(x,y){z = cbind(x,y); 
                                              apply(z,1,limit_state_function)})
	if(add==FALSE) {
		if(is.null(output_dir)) {dev.new("x11",title="MetaIS")}
		else{
			fileDir = paste(output_dir,"_MetaIS.jpeg",sep="")
			jpeg(fileDir)
		}
		par(pty="s")
		plot(x_plot,y_plot,type="n")
	}
}
	
while(cov_epsilon > precision){

	if(plot == TRUE){
	  if(verbose>0){cat(" * 2D PLOT \n")}
		if(is.null(limit_fun_MH)) {
			symbols(0,0, circles=Ru, inches=F, add=TRUE, lty=2)
		}
		else {
			if(is.null(z_MH)) {z_MH = outer(x_plot,y_plot,function(x,y){z = cbind(x,y); 
								    z=t(z);limit_fun_MH(z)$mean})}
			contour(x_plot,y_plot,z_MH,
		level=0,labels="Subset limit",method="edge",add=TRUE)
		}
		tryCatch(points(learn_db[1,],learn_db[2,],col=2,pch=3))
		contour(x_plot,y_plot,z_lsf,
	    level=failure,labels=paste("LSF=",failure,sep=""),method="edge",add=TRUE)
	}


	cat("\n A- REFINEMENT OF PROBABILISTIC CLASSIFICATION FUNCTION PI \n")
	cat("    ======================================================== \n\n")

	ITER = ITER + 1
	cat(" ITERATION ",ITER,"\n")
	cat(" -------------\n\n")

	if(N_DOE>0){
	  if(verbose>0){cat(" * Generate N1 =",N1,"samples uniformly distributed in a hypersphere of radius Ru =",Ru,"\n")}
		U$N1 = runifSphere(dimension,N1,radius=Ru)

	  if(verbose>0){cat(" * Get N_DOE =",N_DOE,"points by clustering of the N1 =",N1,"points\n")}
		if(!is.null(limit_fun_MH)) {
			ind = limit_fun_MH(U$N1)$mean #in a Subset algorithm, select points in Fi-1
			prop = sum(ind<0)/N1 #get accept-reject proportion
			U$N1 = runifSphere(dimension,ceiling(N1/prop),radius=Ru) #generate more points to take AR proportion into account
			ind = limit_fun_MH(U$N1)$mean
			U$N1 = U$N1[,ind<0] #finally get ~N1 points uniformly distributed in the subdomain Fi-1
		}
		U$N_DOE = t(kmeans(t(U$N1), centers=N_DOE,iter.max=20)$centers)
		
	  if(verbose>0){cat(" * Assessment of performance function G on these points\n")}
		G$N_DOE = limit_state_function(U$N_DOE);Ncall = Ncall + N_DOE
	
	  if(verbose>0){cat(" * Add points to the learning database\n")}
		if(is.null(learn_db)){
			learn_db = U$N_DOE
			G$g = G$N_DOE
		}
		else{
			learn_db = cbind(learn_db,U$N_DOE)
			G$g = c(G$g,G$N_DOE)
		}

		if(plot==TRUE){
		  if(verbose>0){cat(" * 2D PLOT : FIRST DoE \n")}
			points(U$N1[1,],U$N1[2,],col=4)
			points(U$N_DOE[1,],U$N_DOE[2,],pch=3,col=2)
		}
	}

	if(verbose>0){cat(" * Train the model\n")}
	if(is.null(meta_model) || learn_each_train==TRUE) {
	  if(verbose>1){cat("   - Learn hyperparameters !!! \n")}
		meta = trainModel(design=learn_db,
		      response=(G$g-failure),
		      kernel=kernel,
		      type="Kriging")
	}
	else {
	  if(verbose>1){cat("   - Use previous hyperparameters !!! \n")}
		meta = trainModel(meta_model,
		      updesign=U$N_DOE,
		      upresponse=(G$N_DOE-failure),
		      type="Kriging")
	}
	
	if(verbose>0){cat("\n * UPDATE quantities based on kriging surrogate model : MP, wMP, pi\n")}
	meta_model = meta$model
	meta_fun = meta$fun
	MP = function(x,k=k_margin) {
		x = as.matrix(x)
		G_meta = meta_fun(x)
		res = pnorm((k*G_meta$sd - G_meta$mean)/G_meta$sd) - pnorm(-(k*G_meta$sd + G_meta$mean)/G_meta$sd)
		return(res)
	}
	wMP = function(x) {
		x = as.matrix(x)
		MP(x)*apply(x,2,function(u) {exp(-0.5*t(u)%*%u)})
	}
	pi = function(x) {
		x = as.matrix(x)
		G_meta = meta_fun(x)
		pnorm(-G_meta$mean/G_meta$sd)
	}

	#plotting part
	if(plot==TRUE){
	  if(verbose>0){cat(" * 2D PLOT \n")}
		z_meta = meta_fun(z)
		z_meta$mean = matrix(z_meta$mean,length(x_plot),length(y_plot))
		z_meta$sd = matrix(z_meta$sd,length(x_plot),length(y_plot))
		z_crit = abs(z_meta$mean)/z_meta$sd
		contour(x_plot,y_plot,z_meta$mean,
	    level=0,labels="Metamodel",method="edge",add=TRUE,col=4)
		contour(x_plot,y_plot,z_crit,
	    level=k_margin,labels=paste("-",k_margin,sep=""),
	    method="edge",lty="dashed",add=TRUE,col=4)
		contour(x_plot,y_plot,z_crit,
	    level=-k_margin,labels=paste("+",k_margin,sep=""),
	    method="edge",lty="dashed",add=TRUE,col=4)
	}

	k = 0;
	Nmax = min(Nmax,Ncall_max - N_alpha);
	if(verbose>0){cat(" * Calculate alphaLOO \n")}
	LOO = leaveOneOut.km(meta_model, type="UK", trend.reestim=FALSE)
	piLOO = pnorm(-LOO$mean/LOO$sd)
	notNull = piLOO>(10^(-16))
	if(verbose>1){cat("   -",sum(piLOO<10^(-16)),"samples not considered as pi<10^-16\n")}
	alphaLOO = mean(1*(G$g[notNull]<failure)/piLOO[notNull])
	if(verbose>1){cat("   - alphaLOO =",alphaLOO,"\n")}
	criterion = (alpha_int[1]<alphaLOO)*(alphaLOO<alpha_int[2])*(k>=Nmin) + (k>Nmax)

	while(!criterion) {
	  if(verbose>0){cat(" * Criterion not reached :\n")}
	  if(verbose>1){cat("   - alphaLOO :",alpha_int[1],"<",alphaLOO,"<",alpha_int[2],"\n")}
	  if(verbose>1){cat("   - k :",Nmin,"<",k,"<",Nmax,"\n")}
		ITER = ITER + 1
		cat("\n\n ITERATION ",ITER,"\n")
		cat(" -------------\n\n")
	  if(verbose>0){cat(" * Find seeds using accept-reject strategy on a standard population\n")}
		n     = 0;
		niter = 0;
		candidate = matrix(NA,dimension,N_seeds)
		while(n<N_seeds && niter<=Niter_seed) {
			niter = niter + 1
			missing_seeds = N_seeds - n
			tmp = matrix(rnorm(dimension*missing_seeds),dimension,missing_seeds)
			is_margin = wMP(tmp)
			is_margin = is_margin>(10^(-16))
			found_seeds = sum(is_margin)
			if(found_seeds>0) candidate[,n + 1:found_seeds] = tmp[,is_margin]
			n = sum(!is.na(candidate[1,])) #number of NAs stands for number of missing vectors
		}
		if(niter>Niter_seed && n<N_seeds){
			criterion = TRUE;
			if(verbose>1){cat("   - No seed found in",niter,"iterations, end of refinement step\n")}
		}
		else{
		  if(verbose>1){cat("   -",n,"seed(s) founded after",niter,"iterations\n")}

			#plotting part
			if(plot==TRUE){
			  if(verbose>0){cat(" * 2D PLOT \n")}
				points(candidate[1,],candidate[2,],col=4)
			}

		  if(verbose>0){cat(" * Generate N_alphaLOO =",N_alphaLOO,"samples with the weighted margin probability\n")}
			gen_pop = generateWithlrmM(seeds=candidate,
				    N=N_alphaLOO,
				    p=wMP,
				    modified=FALSE,
				    burnin=0,thinning=0)
			candidate = gen_pop$points

			#plotting part
			if(plot==TRUE){
			  if(verbose>0){cat(" * 2D PLOT \n")}
				points(candidate[1,],candidate[2,],col=4)
			}

		  if(verbose>0){cat(" * Get K_alphaLOO =",K_alphaLOO,"points by clustering\n")}
			candidate = tryCatch(as.matrix(t(kmeans(t(candidate), centers=K_alphaLOO, iter.max=30)$centers)),
						error = function(cond){
							message(cond);
							r = rankMatrix(candidate);
							res = as.matrix(t(kmeans(t(candidate), centers=r, iter.max=30)$centers))
							return(res)
						})

		  if(verbose>0){cat(" * Calculate performance function G on candidates\n")}
			eval = limit_state_function(candidate);Ncall = Ncall + K_alphaLOO

		  if(verbose>0){cat(" * Add points to he learning database\n")}
			learn_db = cbind(learn_db,candidate)
			G$g = c(G$g,eval)

		  if(verbose>0){cat(" * Train the model\n")}
			if(learn_each_train==TRUE) {
			  if(verbose>1){cat("   - Learn hyperparameters\n")}
				meta = trainModel(design=learn_db,
						  response=(G$g-failure),
						  kernel=kernel,type="Kriging")
			}
			else {
			  if(verbose>1){cat("   - Use previous hyperparameters\n")}
				meta = trainModel(meta_model,
						  updesign=candidate,
						  upresponse=(eval-failure),type="Kriging")
			}

		  if(verbose>0){cat("\n * UPDATE quantities based on kriging surrogate model : MP, wMP, pi\n")}
			meta_model = meta$model
			meta_fun = meta$fun
			MP = function(x,k=k_margin) {
				x = as.matrix(x)
				G_meta = meta_fun(x)
				res = pnorm((k*G_meta$sd - G_meta$mean)/G_meta$sd) - pnorm(-(k*G_meta$sd + G_meta$mean)/G_meta$sd)
				return(res)
			}
			wMP = function(x) {
				x = as.matrix(x)
				MP(x)*apply(x,2,function(u){exp(-0.5*t(u)%*%u)})
			}
			pi = function(x) {
				x = as.matrix(x)
				G_meta = meta_fun(x)
				pnorm(-G_meta$mean/G_meta$sd)
			}

			k = k + K_alphaLOO
		  if(verbose>0){cat(" * Calculate alphaLOO\n")}
			LOO = leaveOneOut.km(meta_model, type="UK", trend.reestim=FALSE)
			piLOO = pnorm(-LOO$mean/LOO$sd)
			notNull = piLOO>(10^(-16))
		  if(verbose>1){cat("   -",sum(piLOO<10^(-16)),"samples not considered as pi<10^-16\n")}
			alphaLOO = mean(1*(G$g[notNull]<failure)/piLOO[notNull])
		  if(verbose>1){cat("   - alphaLOO =",alphaLOO,"\n")}

			criterion = (alpha_int[1]<alphaLOO)*(alphaLOO<alpha_int[2])*(k>=Nmin) + (k>Nmax)
	
			#plotting part
			if(plot==TRUE | (limited_plot && criterion) ){
			  if(verbose>0){cat(" * 2D PLOT \n")}
				z_meta = meta_fun(z)
				z_meta$mean = matrix(z_meta$mean,length(x_plot),length(y_plot))
				z_meta$sd = matrix(z_meta$sd,length(x_plot),length(y_plot))
				z_crit = abs(z_meta$mean)/z_meta$sd
				contour(x_plot,y_plot,z_meta$mean,
					level=0,labels="Metamodel",
					method="edge",add=FALSE,col=4)
				contour(x_plot,y_plot,z_crit,
					level=k_margin,labels=paste("-",k_margin,sep=""),
					method="edge",lty="dashed",add=TRUE,col=4)
				contour(x_plot,y_plot,z_crit,
					level=-k_margin,labels=paste("+",k_margin,sep=""),
					method="edge",lty="dashed",add=TRUE,col=4)
				contour(x_plot,y_plot,z_lsf,
					level=failure,labels=paste("LSF=",failure,sep=""),
					method="edge",add=TRUE)
				if(is.null(seeds)) {
					symbols(0,0,circles=Ru,inches=F,add=TRUE, lty=3)
				}
				else {contour(x_plot,y_plot,z_MH,
					      level=0,labels="Subset limit",
					      method="edge",add=TRUE)}
				points(learn_db[1,],learn_db[2,],col=2,pch=3)
			}
		}
	}

	cat("\n B- ESTIMATE AUGMENTED FAILURE PROBABILITY USING MC ESTIMATOR  \n")
	cat("    =========================================================== \n\n")

	if(!is.null(limit_fun_MH)){
		if(sampling_strategy=="MH"){
		  if(verbose>0){cat(" * Generate Monte-Carlo population with Metropolis-Hastings algorithm\n")}
			gen_pop = generateWithlrmM(seeds=seeds,N=N,limit_f=limit_fun_MH,burnin=burnin,thinning=thinning,VA_function=pi)
			U$N = gen_pop$points
			Ind = gen_pop$VA_values
			P_epsilon = gen_pop$est
			MC_var = gen_pop$var
			cov_epsilon = gen_pop$delta
			MC_gamma = gen_pop$gamma
		}
		else{
		  if(verbose>0){cat(" * Generate Monte-Carlo population with an accept-reject strategy on a standard gaussian sampling\n")}
			rand = function(dimension,N) {matrix(rnorm(dimension*N),dimension,N)}
			U$N = generateWithAR(N,dimension,limit_f=limit_fun_MH,rand=rand)

		  if(verbose>0){cat(" * Calculate Monte-Carlo estimate\n")}
			Ind = pi(U$N)
			P_epsilon <- MC_est <- mean(Ind)
			VA_var = var(Ind)
			MC_var = VA_var/N
			cov_epsilon = sqrt(MC_var)/MC_est
			MC_gamma = 0
		}
	}
	else{
	  if(verbose>0){cat(" * Generate standard gaussian samples\n")}
		U$N = matrix(rnorm(dimension*N),dimension,N)
	  if(verbose>0){cat(" * Calculate Monte-Carlo estimate\n")}
		Ind = pi(U$N)
		P_epsilon <- MC_est <- mean(Ind)
		VA_var = var(Ind)
		MC_var = VA_var/N
		cov_epsilon = sqrt(MC_var)/MC_est
		MC_gamma = 0
	}

	G_meta$N = meta_fun(U$N)$mean
	fail_points = (G_meta$N<0)
	points=U$N[,fail_points]
	meta_eval=G_meta$N[fail_points]
	cat(" P_epsilon =",P_epsilon,"\n")
	cat(" cov_epsilon =",cov_epsilon,"\n")

	if(cov_epsilon>precision) {
		N = ceiling((1+MC_gamma)*VA_var/(precision^2*MC_est^2))
		cat(" * cov_epsilon too large ; this order of magnitude for the probabilty brings N =",N,"\n")
	}
}


#Adaptative importance sampling scheme
cat("\n===========================================================================\n")
cat("\n STEP 2 : Adaptative importance sampling scheme \n")
cat("\n===========================================================================\n\n")

if(verbose>0){cat(" * Define h PDF \n")}
	h = function(x) {
		x = as.matrix(x)
		res = pi(x)*apply(x,2,function(u) {prod(dnorm(u))})/P_epsilon
		return(res)
	}

if(verbose>0){cat(" * Generate samples according to h with Metropolis-Hastings and calculate MC estimator\n")}
if(verbose>1){cat("   - Calculate approximated optimal density on the learning database\n")}
	h_learn_db = h(learn_db)
	
if(verbose>1){cat("   - Find samples whose value is > 0\n")}
	if(dim(learn_db[,h_learn_db>0])[2]>0) {
	  if(verbose>1){cat("     ? seeds =",dim(as.matrix(learn_db[,h_learn_db>0]))[2],"samples in learn_db | h > 0\n")}
		seeds_alpha = learn_db[,h_learn_db>0]
	}
	else{
	  if(verbose>1){cat("     ? No samples in learn_db | h > 0 \n")}
	  if(verbose>0){cat("       => Calculate quasi-optimal density fonction on previous Monte-Carlo population\n")}
		h_U = h(U$N)
		cat("     ? Select points | h > 0\n")
		seeds_alpha = U$N[,h_U>0]
	}
	N_alpha_max = Ncall_max - Ncall

if(verbose>1){cat("   - Generate N_alpha_max =",N_alpha_max,"points from",dim(seeds_alpha)[2],"seeds with r-mM algorithm\n")}
	gen_pop = generateWithlrmM(seeds=seeds_alpha,N=N_alpha_max,p=h,modified=FALSE)
	U$N_alpha_max = gen_pop$points
	inDB = FindInDatabase(mat=U$N_alpha_max,db=learn_db)
	G$N_alpha_max = c(1:N_alpha_max)*NA
	G$N_alpha_max[!is.na(inDB)] = G$g[inDB[!is.na(inDB)]]

	if(plot==TRUE) {
	  if(verbose>1){cat("   - 2D PLOT \n")}
		points(U$N_alpha_max[1,],U$N_alpha_max[2,],pch=8,col=3)
	}

	cov_alpha = Inf
	while((cov_alpha>precision)*(N_alpha<N_alpha_max)) {
	  if(verbose>1){cat("   - Monte-Carlo estimator of alpha\n")}
	  if(verbose>1){cat("     ? Select randomly N_alpha = ",N_alpha," in the working population\n")}
		indices = tryCatch(c(indices,sample(c(1:N_alpha_max)[-indices],(N_alpha-length(indices)),replace=F)),
				error=function(cond) {return(sample(c(1:N_alpha_max),N_alpha,replace=F))})
		U$Nalpha = U$N_alpha_max[,indices]
		G$Nalpha = G$N_alpha_max[indices]

	  if(verbose>1){cat("     ? Evaluate limit_state_function on these samples if necessary\n")}
		isNAinG = is.na(G$Nalpha);
	  if(verbose>1){cat("      +",sum(isNAinG),"samples to evaluate in N_alpha =",N_alpha," samples\n")}
		G$Nalpha[isNAinG] = limit_state_function(U$Nalpha[,isNAinG]); Ncall = Ncall + length(isNAinG)
		learn_db = cbind(learn_db,U$Nalpha[,isNAinG])
		G$g = c(G$g,G$Nalpha[isNAinG])

	  if(verbose>1){cat("     ? Evaluate meta-model derived indicatrice pi on these samples\n")}
		Ind_meta = pi(U$Nalpha)
		Ind_lsf = (G$Nalpha<failure)

	  if(verbose>0){cat("#Evaluate kriging indicatrice likeness\n")}
	  if(verbose>0){cat(" mean(Ind_meta) =",mean(Ind_meta),"\n")}

	  if(verbose>0){cat("#Evaluate alpha estimator\n")}
		Ind = Ind_lsf/Ind_meta
		alpha = mean(Ind)
		VA_var = var(Ind)
		MC_var = VA_var/N_alpha
		VA_cov = sd(Ind)/alpha
		cov_alpha <- MC_cov <- sqrt(1/N_alpha)*VA_cov
		cat(" alpha =",alpha,"\n")
		cat(" cov_alpha =",cov_alpha,"\n")

		if(cov_alpha>precision) {
			N_alpha = ceiling(VA_var/(precision*alpha)^2)
			cat("#cov_alpha too large ; this order of magnitude for alpha brings N_alpha =",N_alpha,"\n")
			if(N_alpha>N_alpha_max) {
				cat("#N_alpha =",N_alpha,"> N_alpha_max =",N_alpha_max," => N_alpha = N_alpha_max\n")
				N_alpha = N_alpha_max;
			}
		}
	}
	if(N_alpha==N_alpha_max){
		U$Nalpha = U$N_alpha_max
		G$Nalpha = G$N_alpha_max

		if(verbose>0){cat("#Evaluate limit_state_function on these samples if necessary\n")}
		isNAinG = is.na(G$Nalpha);
		if(verbose>1){cat(sum(isNAinG),"samples to evaluate in N_alpha =",N_alpha," samples\n")}
		G$Nalpha[isNAinG] = limit_state_function(U$Nalpha[,isNAinG]); Ncall = Ncall + length(isNAinG)
		learn_db = cbind(learn_db,U$Nalpha[,isNAinG])
		G$g = c(G$g,G$Nalpha[isNAinG])

		if(verbose>0){cat("#Evaluate meta-model derived indicatrice pi on these samples\n")}
		Ind_meta = pi(U$Nalpha)
		Ind_lsf = (G$Nalpha<failure)

		if(verbose>0){cat("#Evaluate alpha estimator\n")}
		Ind = Ind_lsf/Ind_meta
		alpha = mean(Ind)
		VA_var = var(Ind)
		MC_var = VA_var/N_alpha
		VA_cov = sd(Ind)/alpha
		cov_alpha <- MC_cov <- sqrt(1/N_alpha)*VA_cov
		cat(" alpha =",alpha,"\n")
		cat(" cov_alpha =",cov_alpha,"\n")
	}

#Results
cat("==========================================================================================",
"                              End of Meta-IS algorithm",
"==========================================================================================",sep="\n")
P = P_epsilon*alpha
cov = sqrt(cov_epsilon^2 + cov_alpha^2 + cov_epsilon^2*cov_alpha^2)

cat("P_epsilon =",P_epsilon,"\n
alpha =",alpha,"\n
P =",P,"\n
cov_epsilon =",cov_epsilon,"\n
cov_alpha =",cov_alpha,"\n
cov =",cov,"\n")

if(plot + limited_plot){
    res = list(proba=P,
               cov=cov,
               Ncall=Ncall,
               learn_db=learn_db,
               lsf_value=G$g,
               meta_fun=meta_fun,
               meta_model=meta_model,
               points=points,
               meta_eval=meta_eval,
               z_meta=z_meta$mean);
}
else {res = list(proba=P,
		cov=cov,
		Ncall=Ncall,
		learn_db=learn_db,
		lsf_value=G$g,
		meta_fun=meta_fun,
		meta_model=meta_model,
		points=points,
		meta_eval=meta_eval)}



return(res)

}
 
