## -----------------------------------------------------------------------------
## Fonction getSeedBelowThreshold
## -----------------------------------------------------------------------------
##    Copyright (C) 2013
##    Developpement : C. WALTER
##    CEA
## -----------------------------------------------------------------------------

getSeedBelowThreshold = function(seeds,eval_seeds,threshold,limit_f,niter_max=20) {

#this function draws Markov chain from seeds to make them being below a given threshold regarding to a given limit state function
#each sample is generated as for a classic Markov chain but is kept only if is value is lower than the seed one's
#then iteratively this proc produces points including into a domain defined by limit_f<threshold

	seeds = as.matrix(seeds)
	Nseeds = dim(seeds)[2]
	samples = seeds
	eval_samples = eval_seeds
	Ncall = 0

	for(i in 1:Nseeds){
		y = eval_seeds[i]
		limit_MH = function(x) {tryCatch(limit_f(x)$mean-y,
						  error=function(cond){limit_f(x)-y}
					)}
		x0 = as.matrix(seeds[,i])
		candidate = MetropolisHastings(x0=x0,eval_x0=0,chain_length=1,modified=FALSE,limit_fun = limit_MH,burnin=0,thinning=0)
		if(candidate$eval[2]==0 && !is.na(candidate$eval_samples[2])){ # this means transition was rejected
		  samples = cbind(samples,candidate$samples[,2])
		  eval_samples = c(eval_samples,candidate$eval_samples[2]) 
		}
		else{
		  samples = cbind(samples,candidate$points[,2])
		  eval_samples = c(eval_samples,candidate$eval[2]+y)
		}
		seeds[,i] = candidate$points[,2]
		eval_seeds[i] = candidate$eval[2]+y
		Ncall = Ncall+1

		niter = 1;
		while(eval_seeds[i]>threshold && niter<=niter_max) {
			y = eval_seeds[i]
			limit_MH = function(x) {tryCatch(limit_f(x)$mean-y,error=function(cond){limit_f(x)-y})}
			candidate = MetropolisHastings(x0=seeds[,i],eval_x0=0,chain_length=1,modified=FALSE,limit_fun = limit_MH,burnin=0,thinning=0)
			if(candidate$eval[2]==0 && !is.na(candidate$eval_samples[2])){ # this means transition was rejected
			  samples = cbind(samples,candidate$samples[,2])
			  eval_samples = c(eval_samples,candidate$eval_samples[2]) 
			}
      else{
        samples = cbind(samples,candidate$points[,2])
        eval_samples = c(eval_samples,candidate$eval[2]+y)
      }
			seeds[,i] = candidate$points[,2]
			eval_seeds[i] = candidate$eval[2]+y
			Ncall = Ncall+1
			niter = niter + 1;
		}
		if(niter>niter_max){
			cat("#For seed number",i,"no new seed found after",niter,"iterations, NAs inserted\n")
			seeds[,i] = NA
			eval_seeds[i] = NA
		}

	}

	res = list(seeds=seeds,eval=eval_seeds,samples=samples,eval_samples=eval_samples,Ncall=Ncall)

	return(res)
}