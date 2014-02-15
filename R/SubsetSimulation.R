## -----------------------------------------------------------------------------
## Fonction SubsetSimulation
## -----------------------------------------------------------------------------
##    Copyright (C) 2013
##    Developpement : C. WALTER
##    CEA
## -----------------------------------------------------------------------------

SubsetSimulation = function(dimension,
          limit_state_function,
			    proposal_pdf,             # Proposal PDF for MH
			    pdf = dnorm,              # PDF of the input space
			    rpdf = rnorm,             # Random generator for the PDF
			    cutoff_prob = 0.1,        # Cutoff probability for the subsets
			    n_init_samples = 10000,   # Number of samples per subset
          burnin = 20,
          thinning = 4,
			    plot = FALSE,
			    output_dir = NULL,
          verbose = 0) {            # Either 0 for almost no output, 1 for medium size output and 2 for all outputs

cat("==========================================================================================\n")
cat("                              Beginning of Subset Simulation algorithm \n")
cat("==========================================================================================\n\n")
  
#general meaning of variables used here
# U stores the coordinates of the n_init_samples=N samples ; U = matrix(dimension x N)
# G stores the value of the limit state function on these N samlpes ; G = vector(N)
# G_sorted the list got from sort(G) : G[[1]] = vector N, value of G sorted ; G[[2]] = vector N, indices of G
# P the final probability, updated at each stage
# seeds stores the coordinates of the N*p0 samples initializing the MC at each stage ; seeds = matrix(dimension x floor(p0*N)
# I stands for the indicatrice function at each stage and stores the value for each MC sample ; I = matrix(n(seeds) x floor(1/p0)+1)
# R is the covariance between samples ; R = vector floor(1/p0)+1 | R[k] = covariance of samples at length k
# gamma : see definition in Au & Beck's paper, implemented as R list
# delta list of c.o.v for each conditional probability

if(missing(proposal_pdf)){
  proposal_pdf = function(x,y=NA,width=2) {
    if(is.na(y[1])){
      runif(length(x),min=-width/2,max=width/2)+x
    }
    else{
      1/width^length(x)*(max(abs(x-y))<=2)
    }
  }
}

cat("  * n_init_samples :",n_init_samples,"\n\n");

cat(" ================== \n")
cat(" STEP 1 : FIRST DoE \n")
cat(" ================== \n\n")

if(verbose>0){cat("  * generate the first N =",n_init_samples,"samples by crude MC \n")}
#generate the first N samples by crude MC
U = matrix(NA, nrow=dimension, ncol=n_init_samples);
for (i in c(1:n_init_samples)){
	U[,i] = rpdf(dimension);
}

if(verbose>0){cat("  * evaluate the limit state function on these points \n")}
#calcul the limit state function on these points
G = apply(U, 2, limit_state_function);
Ncall = n_init_samples
G_sorted = sort(G,index.return=TRUE);

if(verbose>0){cat("  * find the quantile q0 verifying P[g(U) < q0] = p0 \n")}
#find the quantile q0 verifying P(g(U)<q0) = p0
q0_rank = floor(cutoff_prob*n_init_samples)
q0 = G_sorted[[1]][q0_rank];
cat("      q0 =",q0,"\n")

if(verbose>0){cat("  * Evaluate actual probability \n")}
#calculate the actual probability
if(q0>=0){
	P = q0_rank/n_init_samples
	subset_prob = cutoff_prob
}
else{
	#set q0 = 0
	q0 = 0

	#calculate how many negative values in the last N g(U)
	s = sum((G-abs(G))/(2*G))

	#and so get the probability
	subset_prob = s/n_init_samples
	P = subset_prob
}
n_subset = 1
cat("      P =",P,"\n")
cat("      q0 =",q0,"\n")

if(verbose>0){cat("  * Evaluate COV for the first subset \n")}
#calculate delta the cov for the first subset
delta2 = (1-subset_prob)/subset_prob/n_init_samples

#plotting part
if(plot==TRUE){
  if(verbose>0){cat("  * 2D plot : FIRST STEP \n")}
	x = c(-80:80)/10
	y = c(-80:80)/10
	z = outer(x,y,function(x,y){z = cbind(x,y); apply(z,1,limit_state_function)})
	if(is.null(output_dir)) {dev.new("x11",title="Subset Simulation")}
	else{
		fileDir = paste(output_dir,"_Subset_Simulation.jpeg",sep="")
		jpeg(fileDir)
	}
	par(pty="s")
	plot(x,y,type="n")
}

cat("\n ================ \n")
cat(" SUBSET ALGORITHM \n")
cat(" ================ \n\n")


#beginning of the subset algorithm
while(q0>0){

    
    n_subset = n_subset + 1
    cat("\n  * STEP ",n_subset," :\n")
    
	#plotting part
	if(plot==TRUE){
	  if(verbose>1){cat("    - 2D PLOT \n")}
		points(U[1,], U[2,], col=n_subset, pch=n_subset)
		contour(x,y,z, level=q0, method="edge", add=TRUE)
	}

	
    if(verbose>1){cat("    - get the floor(p*N)=",q0_rank,"seeds for next step Metropolis-Hastings algorithm \n")}
	#get the floor(p*N)=q0_rank seeds for next step Metropolis-Hastings algorithm
	seeds_ind = G_sorted[[2]][1:q0_rank];
	seeds = U[,seeds_ind];
	n_seeds = q0_rank

    if(verbose>1){cat("    - Get N =",n_init_samples,"samples from theses seeds with Metropolis-Hastings algorithm\n")}
	gen_pop = generateWithlrmM(seeds = seeds,
                               seeds_eval = G[seeds_ind]-q0,
                               N = n_init_samples,
                               limit_f = function(x) {limit_state_function(x) - q0},
                               p = pdf,
                               q = proposal_pdf,
                               modified = TRUE,
                               burnin = burnin, thinning = thinning)
	U = gen_pop$points
	G = gen_pop$eval+q0
	chain_length = gen_pop$chain_length
	Ncall = Ncall + gen_pop$Ncall

    if(verbose>1){cat("    - find the quantile q verifying P(g(U)<q) = p \n")}
	#find the quantile q verifying P(g(U)<q) = p
	G_sorted = sort(G,index.return=TRUE);
	q0 = G_sorted[[1]][floor(cutoff_prob*n_init_samples)];
    if(verbose>1){cat("      q0 =",q0,"\n")}

    if(verbose>1){cat("    - calculate the actual probability \n")}
	#calculate the actual probability
	if(q0 >= 0){
		P = P*floor(cutoff_prob*n_init_samples)/n_init_samples
		subset_prob = cutoff_prob
	}
	else{
		#set q0 = 0
		q0 = 0

		#calculate how many negative values in the last N g(U)
		s = sum((G-abs(G))/(2*G))

		#and so get the probability
		subset_prob = s/n_init_samples
		P = P*subset_prob
	}
	cat("      P =",P,"\n")
	cat("      q0 =",q0,"\n")

    if(verbose>1){cat("    - calculate the covariance between samples \n")}
	#calculate the covariance between samples
	Ind = eval(G<=q0)*1
	if(verbose>1){MC_stat = MCMCcovariance(n_seeds = n_seeds, 
                             chain_length = chain_length, 
                             VA_values = Ind, 
                             VA_esp = subset_prob, 
                             VA_var = subset_prob*(1-subset_prob))}
    else{capture.output(MC_stat <- MCMCcovariance(n_seeds = n_seeds, 
                                                 chain_length = chain_length, 
                                                 VA_values = Ind, 
                                                 VA_esp = subset_prob, 
                                                 VA_var = subset_prob*(1-subset_prob)))}

    if(verbose>1){cat("    - Update delta2 \n")}
	#calculate delta
	delta2 = delta2 + (MC_stat$cov)^2

}

#final plot
if(plot==TRUE){
	points(U[1,],U[2,],col=n_subset+1*(G<0),pch=n_subset)
	contour(x,y,z,level=0,method="edge",add=TRUE)
}
if(!is.null(output_dir)) {dev.off()}

#get the global cov estimated with sum of squared delta[[i]]
delta = sqrt(delta2)

#return the result
result = list(proba=P,cov=delta,Ncall=Ncall)

cat("==========================================================================================\n")
cat("                              Beginning of Subset Simulation algorithm \n")
cat("==========================================================================================\n\n")

cat("   - proba =",P,"\n")
cat("   - cov =",delta,"\n")
cat("   - Ncall =",Ncall,"\n")

return(result)

}
