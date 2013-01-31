########################################################################
# compute drawdown probabilities
########################################################################

# make nper draws from a normal with mean sr
# and unit variance. compute the cumulative
# sum, then return the max drawdown. repeat
# n times. returns in units of standard deviations.
#
# generator should be a function which takes
# an integer and mu and generates that many
# i.i.d. draws from a distribution with
# mean mu and sd 1. defaults to a normal.
rdd <- function(n,nper,sr=0,generator=NULL) {
	# generator
	if (is.null(generator)) {
		generator <- function(ngen,sr) { rnorm(ngen,mean=sr,sd=1) }
	}
	# to avoid memory bonks, do this instead:
	MAXEL <- 5e7
	ncut <- ceiling(n*nper/MAXEL)
	.rdd_single <- function(n,nper,sr,generator) {
		rets <- generator(n*nper,sr)
		rets <- matrix(rets,nrow=nper)
		crets <- apply(rets,2,cumsum)
		rm(rets)
		ddown <- apply(crets,2,function(x){ cummax(x) - x })
		retval <- apply(ddown,2,max)
		return(retval)
	}
	retval <- replicate(ncut,.rdd_single(ceiling(n/ncut),nper,sr,generator),
											simplify="array")
	dim(retval) <- c(length(retval),1)
	retval <- retval[1:n]
 #par(mfrow=c(3,1))
 #plot(crets) 
 #plot(apply(crets,2,cummax)) 
 #plot(ddown) 
	return(retval)
}

tgen <- function(n,sr,df=4) { 
	sr + sqrt((df-2)/df) * rt(n,df) 
}

#for (sr in seq(0.1,0.4,by=0.1)) {
	#zzz <- rdd(1024,253,sr)
	#print(quantile(zzz,probs=0.975,type=7))
#}

#for (sr in seq(0.1,0.4,by=0.1)) {
	#zzz <- rdd(1024,253,sr,tgen)
	#print(quantile(zzz,probs=0.975,type=7))
#}


## compute the 97.5% upper quantile of the prediction
## interval on max drawdown, as a factor of # sds
## versus sr

#srx <- seq(0,2.5/sqrt(253),length.out=32)
#sry <- sapply(srx,function(sr) { 
	#zzz <- rdd(2^18,253,sr,tgen)
	#retval <- quantile(zzz,probs=0.975,type=7)[[1]]
#})

#plot(srx,1/sry)
#yfit <- lm(1/sry ~ srx)
#summary(yfit)
#abline(yfit)

#srg <- sapply(srx,function(sr) { 
	#zzz <- rdd(2^14,253,sr)
	#retval <- quantile(zzz,probs=0.975,type=7)[[1]]
#})
#summary(lm(1/srg ~ srx))

#plot(srx,sry-srg)



#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
