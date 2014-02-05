## -----------------------------------------------------------------------------
## Fonction generateWithAR
## -----------------------------------------------------------------------------
##    Copyright (C) 2013
##    Developpement : C. WALTER
##    CEA
## -----------------------------------------------------------------------------

generateWithAR = function(N,dimension,limit_f,rand) {

	tic = proc.time()[3]
	if(missing(rand)) {rand = function(dimension,N) {matrix(rnorm(dimension*N),dimension,N)}}
	limit_fun = function(x) {tryCatch(limit_f(x)$mean,error = function(cond) {return(limit_f(x))})}

	cat("#Find accept-reject ratio\n")
	n_inF = 0
	Ninit = N
	while(n_inF==0) {
		U = rand(dimension=dimension,N=Ninit)
		ind = limit_fun(U) #domain is defined by limit_f<0
		inF = (ind<0)
		n_inF = sum(inF)
		Ninit = 10*N
	}
	prop = n_inF/Ninit #get accept-reject proportion
	cat(" ratio =",prop,"\n")

	cat("#Generate",(ceiling(N/prop)-n_inF),"new points\n")
	tmp = rand(dimension=dimension,N=(ceiling(N/prop)-n_inF)) #generate more points to take AR proportion into account

	cat("#Test new samples\n")
	ind = limit_fun(tmp)

	cat("#Get final samples by concatenating first and second run\n")
	U = cbind(U[,inF],tmp[,ind<0])
	U = U[,1:N]

	toc = proc.time()[3] - tic
	cat(" targeted population size =",N,"\n")
	cat("",N+ceiling(N/prop)-n_inF,"samples generated in",toc,"sec,",dim(U)[2],"points kept after accept-reject test\n")
	return(U)
}