## -----------------------------------------------------------------------------
## Fonction Monte Carlo
## -----------------------------------------------------------------------------
##    Copyright (C) 2013
##    Developpement : C. WALTER
##    CEA
## -----------------------------------------------------------------------------

MonteCarlo = function(dimension, 
                      LimitStateFunction,
                      N_max 	= 500000,
		      N_batch 	= 1000,
                      failure	= 0,
		      precision	= 0.05,
                      plot 	= FALSE, 
                      output_dir = NULL,
		      verbose 	= 0){

cat("========================================================================\n")
cat("                 Beginning of Monte-Carlo algorithm\n")
cat("========================================================================\n\n")
  
# Step 0 : Initialization
cov = Inf;
Ncall = 0;

if(verbose>0){cat(" * STEP 1 : FIRST SAMPLING AND ESTIMATION \n")}

if(verbose>1){cat("   - Generate N_batch = ",N_batch," standard samples\n")}
U = matrix(rnorm(dimension*N_batch),dimension,N_batch)

if(verbose>1){cat("   - Evaluate LSF on these samples\n")}
G = apply(U, 2, LimitStateFunction); Ncall = Ncall + N_batch
indG <- which(G < failure)

if(verbose>1){cat("   - Evaluate failure probability and corresponding CoV \n")}
P = mean(G < failure)
cov = sqrt((1-P)/(Ncall*P))

Nrem = min(N_max - Ncall, N_batch)

if( (verbose>0) && (cov > precision) && (Nrem > 0) ) {cat(" * STEP 2 : LOOP UNTIL COV < PRECISION \n")}

while( (cov > precision) && (Nrem > 0) ) {

	if(verbose>0){cat(" * cov =",cov,">",precision,"and",N_max - Ncall,"remaining calls to the LSF \n")}
	if(verbose>1){cat("   - Generate N =",Nrem,"standard samples\n")}
	U_new = matrix(rnorm(dimension*Nrem),dimension,Nrem)
	U = cbind(U,U_new)

	if(verbose>1){cat("   - Evaluate LSF on these samples\n")}
	G = c(G,apply(U_new, 2, LimitStateFunction)); Ncall = Ncall + Nrem
	indG <- which(G < failure)
	
	if(verbose>1){cat("   - Evaluate failure probability and corresponding CoV \n")}
	P = mean(G < failure)
	cov = sqrt((1-P)/(Ncall*P))

	if(verbose>1){cat("   - P =",P,"\n")}
	if(verbose>1){cat("   - cov =",cov,"\n")}

	Nrem = min(N_max - Ncall, N_batch)
}

if(plot == TRUE) {
	#plotting part
	if(verbose>0){cat(" * 2D plot : LSF and samples \n")}
	x = c(-80:80)/10
	y = c(-80:80)/10
	z = outer(x,y,function(x,y){z = cbind(x,y); apply(z,1,LimitStateFunction)})
	if(is.null(output_dir)) {dev.new("x11",title="Monte Carlo estimation")}
	else{
		fileDir = paste(output_dir,"_Monte_Carlo_brut.jpeg",sep="")
		jpeg(fileDir)
	}
	par(pty="s")
	plot(x, y, type="n")
	points(U[1,], U[2,], col='#00000022',pch=19, cex=0.8)
	points(U[1,indG], U[2,indG], col='#FF000022',pch=19, cex=0.8)
	contour(x, y, z, level=failure, labcex = 0.8, method="edge", add=TRUE, col='blue3', lwd=2)
	if(!is.null(output_dir)) {dev.off()}
}

cat("========================================================================\n")
cat("                    End of Monte-Carlo algorithm\n")
cat("========================================================================\n\n")

cat("   - proba =",P,"\n")
cat("   - cov =",cov,"\n")
cat("   - Ncall =",Ncall,"\n")

res = list(proba=P, cov=cov, Ncall=Ncall)
return(res)
}