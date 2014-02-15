## -----------------------------------------------------------------------------
## Fonction runifSphere
## -----------------------------------------------------------------------------
##    Copyright (C) 2013
##    Developpement : C. WALTER
##    CEA
## -----------------------------------------------------------------------------

runifSphere = function(dimension,N,radius) {

    tmp = rnorm(dimension*N, mean=0, sd=1)
    dim(tmp) = c(dimension,N)
    tmp = apply(tmp, 2, function(x) {dim(x) = c(1,length(x));x/norm(x,"F")*(runif(1,min=0,max=1))^(1/dimension)*radius})

    return(tmp)
}

