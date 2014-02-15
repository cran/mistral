## -----------------------------------------------------------------------------
## Fonction FindInDatabase
## -----------------------------------------------------------------------------
##    Copyright (C) 2013
##    Developpement : C. WALTER
##    CEA
## -----------------------------------------------------------------------------

FindInDatabase = function(mat,db,orientation=2) {

#considering mat and db matrices are two sets of vector of a same vector space, this function allows for seeking the first ones into the second ones
#orientation sets if vector are the lines =1 or columns =2 of the matrices

#Check dimension consistency
if(orientation==1) {
	if(!dim(mat)[2]==dim(db)[2]) {stop(paste("orientation =",orientation,"while matrices doesn't have the same number of columns"))}
	else{mat=t(mat);db=t(db);}
}
if(orientation==2) {
	if(!dim(mat)[1]==dim(db)[1]) stop(paste("orientation =",orientation,"while matrices doesn't have the same number of lines"))
}

is_in_db = apply(mat,2,function(x) {
					test = apply(db,2,function(y) {all(x==y)});
					ind = ifelse(sum(test)>0,which.max(test),NA)
					return(ind)
				}
	)

return(is_in_db)

}