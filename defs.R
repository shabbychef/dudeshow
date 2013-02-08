# 
# * Fri Dec 28 2012 05:26:33 PM Steven E. Pav <steven@cerebellumcapital.com>
#
# base R definitions for dude paper

# compiler flags!

FINAL.VERSION <- TRUE
#FINAL.VERSION <- FALSE


safe.install <- function(pkg.name,do.load=FALSE) {
	if (! length(which(.packages(all.available=TRUE) %in% pkg.name))) {
		install.packages(pkg.name,repos="http://cran.cnr.berkeley.edu/")
	}
	if (do.load) {
		library(pkg.name,character.only=TRUE)
	}
}

dummy <- lapply(c("LambertW","MASS","quantmod","knitr","formatR"),safe.install,do.load=TRUE)
dummy <- lapply(c("fPortfolio","TTR","xtable"),safe.install,do.load=FALSE)

# set the knitr options ... for everyone!
opts_knit$set(echo=FALSE)
opts_knit$set(eps=TRUE)
opts_knit$set(warning=FALSE)
opts_knit$set(message=FALSE)
opts_knit$set(cache=TRUE)
opts_knit$set(cache.path="cache/")
opts_knit$set(results="asis")
opts_knit$set(fig.ext="eps")

# see http://biasedestimators.blogspot.com/2011/02/pretty-printing-r-functions-in-sweave.html
prettyPrint <- function(f,w=60) {
	tidy.source(text=capture.output(print.function(f, source=FALSE)),
							width.cutoff=w)
}

# some constants
dpy <- 253

# make this small to reduce the run time for compiles based on
# monte carlo studies; the default value shall be 1000 for
# the 'normal' high-resolution compile of this document.
# a value less than 100 shall be ignored;
mc.resolution <- ifelse(FINAL.VERSION,1000,200)

########################################################################
# utils:
#make this like cci matlab?
rel.returns <- function(x,k=1) {
	retval <- diff(log(x),lag=k)
	return(retval)
}
rr2lr <- function(x) {
	retval <- log(1+x)
	return(retval)
}
lr2rr <- function(x) { 
	retval <- exp(x) - 1
	return(retval)
}
lr2mtm <- function(x) { 
	retval <- exp(cumsum(x))
	return(retval)
}
rr2mtm <- function(x) { 
	retval <- lr2mtm(rr2lr(x))
	return(retval)
}

########################################################################
## some data#FOLDUP

getSymbols(c("^GSPC","^VIX"),from="1970-01-01")
SPX.adjclose <- GSPC[,4]
SPX.start <- time(SPX.adjclose[1])
SPX.end <- time(SPX.adjclose[length(SPX.adjclose)])
SPX.rets <- rel.returns(GSPC[,4])
SPX.rets <- as.vector(SPX.rets[2:length(SPX.rets),1])

VIX.spot <- VIX[,4]

#UNFOLD
########################################################################
## expectation of the t:#FOLDUP
# the geometric bias of Sharpe ratio
f_tbias <- function(n) { 
	sqrt((n-1) / 2) * exp(lgamma((n-2)/2) - lgamma((n-1)/2))
}

#approximate tbias
f_apx_tbias1 <- function(n) { 
	1 + (0.75 / n)
}
f_apx_tbias2 <- function(n) { 
#1 + (0.75 / n) + (2 / n**2)
	1 + (0.75 / (n - 1)) + (32 / (25 * ((n-1) ** 2)))
}
#UNFOLD
########################################################################
# power #FOLDUP

# the power of the univariate t-test;
f_tpower <- function(n,rho = 0,alpha = 0.05) {
	f_tpower_ncp(ncp = sqrt(n) * rho,n = n,alpha = alpha)
}

# the power of the univariate t-test as function of the ncp
f_tpower_ncp <- function(ncp,n,alpha = 0.05) {
	pt(qt(1-alpha,n-1),df = n-1,ncp = ncp,lower.tail = FALSE)
}

# the power of an f-test 
f_fpower <- function(df1,df2,ncp,alpha = 0.05) {
	pf(qf(alpha,df1=df1,df2=df2,ncp=0,lower.tail=FALSE),df1 = df1,df2 = df2,ncp = ncp,lower.tail = FALSE)
}

# the power of the Hotelling test
f_hpower <- function(n,p,rhosq,alpha = 0.05) {
	f_fpower(df1 = p,df2 = n - p,ncp = n * rhosq,alpha = alpha)
}
#UNFOLD
########################################################################
# sample size computations#FOLDUP

# find the sample size for a given power of the univariate hotelling test
f_hreqsize <- function(rhosq,p,powr = 0.80,alpha = 0.05) {
	#2FIX: get a sane upper bound!
	zz <- uniroot(function(n,p,rhosq,alpha,powr)(f_hpower(n,p,rhosq,alpha) - powr),
								c(max(8,p+1),2000000 * 10 / (rhosq)),p = p,rhosq = rhosq,powr = powr,alpha = alpha)
	return(zz$root)
}

#estimate the required sample size using a fit...
#this was all done at alpha = 0.05, so that is not a variable.
f_est_hreqsize <- function(rhosq,p,powr = 0.80) {
	estv <- exp(2.259 - log(rhosq) + 0.363 * log(p) + 1.31 * log(powr) - 0.0757 * log(powr) * log(p))
}

# find the sample size for a given power of the univariate t-test
f_treqsize <- function(rho,powr = 0.80,alpha = 0.05) {
	zz <- uniroot(function(n,rho,alpha,powr)(f_tpower(n,rho,alpha) - powr),
								c(8,20 * 10 / (rho*rho)),rho = rho,powr = powr,alpha = alpha)
	return(zz$root)
}
#UNFOLD
########################################################################
#inversions#FOLDUP

#find the non-centrality parameter; that is, find
#ncp such that pt(t,df,ncp) = alpha
f_nct_cdf_ncp <- function(t,df,alpha) {
	#find approximate endpoints;
	flo <- min(-1,t - qnorm((1 + alpha)/2))
	fhi <- max(1,t - qnorm((alpha)/2))
	#expand them until pt(t,df,flo) > alpha and pt(t,df,fhi) < alpha
	while ((pt(t,df,flo) < alpha) && (flo > -100000)) { flo <- 2 * flo }
	while ((pt(t,df,fhi) > alpha) && (fhi < 100000)) { fhi <- 2 * fhi }
	ncp <- uniroot(function(ncp,t,df,alpha)(pt(t,df,ncp) - alpha),
								 c(flo,fhi),t = t,df = df,alpha = alpha)
	return(ncp$root)
}

#find some t such that
# phi(t-ncp) = f(t;df,ncp) 
#where phi is the density function of the normal
#and f(x;df,ncp) is the density function of a noncentral
#t-distribution with n d.o.f. and noncentrality parameter ncp

f_eq_t_pdf_disc <- function(df,ncp = 0) {
	lims <- c(-abs(ncp) - 20,ncp)
	foo <- uniroot(function(t,ncp,df)(dnorm(t,mean = ncp,sd=1) -
																		dt(t,df=df,ncp=ncp)),
								 lims,df = df,ncp = ncp)
	return(foo$root)
}

#find 
# max_t | Phi(t-ncp) - F(t;df,ncp) |
#where Phi is the distribution function of the normal
#and F(x;df,ncp) is the distribution function of a noncentral
#t-distribution with n d.o.f. and noncentrality parameter ncp

f_max_t_cdf_disc <- function(df,ncp = 0) {
	opt_t <- f_eq_t_pdf_disc(df=df,ncp=ncp)
	disc <- abs(pnorm(opt_t,mean = ncp,sd=1) -
							pt(opt_t,df=df,ncp=ncp))
	return(disc)
}
#UNFOLD

#compute the asymptotic mean and variance of the sqrt of a
#non-central F distribution

f_sqrt_ncf_asym_mu <- function(df1,df2,ncp = 0) {
	return(sqrt((df2 / (df2 - 2)) * (df1 + ncp) / df1))
}
f_sqrt_ncf_asym_var <- function(df1,df2,ncp = 0) {
	return((1 / (2 * df1)) * 
				 (((df1 + ncp) / (df2 - 2)) + (2 * ncp + df1) / (df1 + ncp)))
}
f_sqrt_ncf_apx_pow <- function(df1,df2,ncp,alpha = 0.05) {
	zalp <- qnorm(1 - alpha)
	numr <- 1 - f_sqrt_ncf_asym_mu(df1,df2,ncp = ncp) + zalp / sqrt(2 * df1)
	deno <- sqrt(f_sqrt_ncf_asym_var(df1,df2,ncp = ncp))
	return(1 - pnorm(numr / deno))
}


# vanilla Sharpe ratio in terms of whatever input units
f_vsharpe <- function(rets) {
	return(mean(rets) / sd(rets))
}

# annualized Sharpe ratio
f_asharpe <- function(rets, dpy = 253) {
	return(sqrt(dpy) * f_vsharpe(rets))
}

# the t statistic
f_tstat <- function(rets) { 
	return(sqrt(length(rets)) * f_vsharpe(rets))
}

########################################################################
# confidence intervals

## confidence intervals on the Sharpe Ratio#FOLDUP

#the sample.sr should *not* be annualized
f_sr_se_shab <- function(sample.sr,n) {
	cn <- f_tbias(n)
	dn <- (n-1) / ((n-3) * cn * cn)
	W  <- (sample.sr / cn) ** 2
	se <- sqrt((dn/n) + W * (dn - 1))
}

f_sr_se_walck <- function(sample.sr,n) {
	se <- sqrt((1/n) + 0.5 * sample.sr ** 2 / (n - 1))
}

f_sr_se_lo <- function(sample.sr,n) {
	se <- sqrt((1 + 0.5 * sample.sr ** 2) / (n - 1))
}

f_sr_ci_shab <- function(sample.sr,n,alpha = 0.05) {
	cn <- f_tbias(n)
	medv <- sample.sr / cn
	se <- f_sr_se_shab(sample.sr,n)
	zalp <- qnorm(1 - alpha / 2)
	cilo <- medv - zalp * se
	cihi <- medv + zalp * se
	return(list('lo' = cilo,'hi' = cihi))
}

f_sr_ci_lo <- function(sample.sr,n,alpha = 0.05) {
	se <- f_sr_se_lo(sample.sr,n)
	zalp <- qnorm(1 - alpha / 2)
	cilo <- sample.sr - zalp * se
	cihi <- sample.sr + zalp * se
	return(list('lo' = cilo,'hi' = cihi))
}

# Walck gives this normal approximation
f_sr_ci_walck <- function(sample.sr,n,alpha = 0.05) {
	se <- f_sr_se_walck(sample.sr,n)
	zalp <- qnorm(1 - alpha / 2)
	midp <- sample.sr * (1 - 1 / (4 * (n - 1)))
	cilo <- midp - zalp * se
	cihi <- midp + zalp * se
	return(list('lo' = cilo,'hi' = cihi))
}

# these are the 'exact' symmetric CI, which I first saw in 
# Scholz' paper. I thought they were novel at that time.:w
f_sr_ci_scholz <- function(sample.sr,n,alpha = 0.05) {
	sn <- sqrt(n)
	t <- sample.sr * sn
	cilo <- (1 / sn) * f_nct_cdf_ncp(t,df = n-1,alpha = 1 - alpha/2)
	cihi <- (1 / sn) * f_nct_cdf_ncp(t,df = n-1,alpha = alpha/2)
	return(list('lo' = cilo,'hi' = cihi))
}
#UNFOLD

# inference on F's ncp#FOLDUP

# confidence distribution, gives CIs
qcofncp <- function(p,Fs,df1,df2,ub=NULL) {
	# return max{0 <= ncp <= ub | pf(Fs,df1,df2,ncp) >= 1 - p}
	# or 0 if none exist
	f.zer <- pf(Fs,df1,df2,0)
	if (f.zer < (1-p)) {
		return(0)
	} else {
		if (is.null(ub)) {
			ub <- 1.0
			fpf <- pf(Fs,df1,df2,ub)
			while (fpf >= (1-p)) {
				ub <- 2.0 * ub
				fpf <- pf(Fs,df1,df2,ub)
			}
			lb <- 0.5 * ub
		} else {
			lb <- 0
			fpf <- pf(Fs,df1,df2,ub)
		}
		# now call uniroot
		zerf <- function(z,n,xv,limv) { pf(xv,p,n-p,n*z^2) - limv }
		ncp <- uniroot(function(ncp,Fs,df1,df2,tgt){pf(Fs,df1,df2,ncp) - tgt},
									 c(lb,ub),Fs=Fs,df1=df1,df2=df2,tgt=1-p)$root
		return(ncp)
	}
}

# use same to construct confidence intervals
fncp.ci <- function(F,df1,df2,alpha.lo=0.025,alpha.hi=1-alpha.lo) {
	if (alpha.hi >= 1) {
		cihi <- Inf
	} else {
		cihi <- qcofncp(alpha.hi,F,df1,df2)
	}

	if (alpha.lo <= 0) {
		cilo <- 0
	} else {
		if (is.finite(cihi)) {
			cilo <- qcofncp(alpha.lo,F,df1,df2,ub=cihi)
		} else {
			cilo <- qcofncp(alpha.lo,F,df1,df2)
		}
	}
	return(list('lo' = cilo,'hi' = cihi))
}

#MLE of the ncp based on a single F-stat
fncp.mle <- function(Fs,df1,df2,ub=NULL,lb=0) {
	if (Fs <= 1) { return(0.0) }  # Spruill's Thm 3.1, eqn 8
	if (is.null(ub)) {
		prevdpf <- -Inf
		ub <- 1
		dpf <- df(Fs,df1,df2,ncp=ub,log=TRUE)
		while (prevdpf < dpf) {
			prevdpf <- dpf
			ub <- 2 * ub
			dpf <- df(Fs,df1,df2,ncp=ub,log=TRUE)
		}
		lb <- ifelse(ub > 2,ub/4,lb)
	}
	ncp.mle <- optimize(function(ncp){df(Fs,df1,df2,ncp=ncp,log=TRUE)},
											c(lb,ub),maximum=TRUE)$maximum;
	return(ncp.mle)
}
#UNFOLD

# Make hotelling distribution stuff like df,pf,qf,rf
f_hot2F <- function(T2,p,n) {
	return(T2 * (n - p) / (p * (n-1)))
}
f_F2hot <- function(F,p,n) {
	return(F * (p * (n-1)) / (n - p))
}

# inference on Hotelling's ncp, by extension#FOLDUP

# confidence distribution, gives CIs
qcoT2ncp <- function(plev,T2,p,n,ub=NULL) {
	# convert to F
	Fs <- f_hot2F(T2=T2,p=p,n=n)
	if (!is.null(ub)) {
		ub <- f_hot2F(T2=ub,p=p,n=n)
	}
	# delegate
	F.ncp <- qcofncp(plev,Fs,df1=p,df2=n-p,ub=ub)
	# they have the same ncp; no back conversion
	return(F.ncp)
}

# use same to construct confidence intervals
T2ncp.ci <- function(T2,p,n,alpha.lo=0.025,alpha.hi=1-alpha.lo) {
	# convert to F
	Fs <- f_hot2F(T2=T2,p=p,n=n)
	# delegate
	F.ci <- fncp.ci(Fs,df1=p,df2=n-p,alpha.lo=alpha.lo,alpha.hi=alpha.hi)
	# they have the same ncp; no back conversion
	return(F.ci)
}

#MLE of the ncp
T2ncp.mle <- function(T2,p,n,ub=NULL) {
	# convert to F
	Fs <- f_hot2F(T2=T2,p=p,n=n)
	if (!is.null(ub)) {
		ub <- f_hot2F(T2=ub,p=p,n=n)
	}
	# delegate
	F.mle <- fncp.mle(Fs,df1=p,df2=n-p,ub=ub)
	# they have the same ncp; no back conversion
	return(F.mle)
}
#UNFOLD

# SR^* is a Hotelling, basically. 

# convert SR^* <-> T2
f_srstar2hot <- function(sample.sr,n) {
	return(n * sample.sr^2)
}
f_hot2srstar <- function(T2,n) {
	return(sqrt(T2 / n))
}

# inference on SR^*'s ncp, by extension#FOLDUP

# confidence distribution, gives CIs
qcosrstarncp <- function(plev,srs,p,n,ub=NULL) {
	# convert to T2
	T2 <- f_srstar2hot(sample.sr=srs,n=n)
	if (!is.null(ub)) {
		ub <- f_srstar2hot(sample.sr=ub,n=n)
	}
	# delegate
	T2.ncp <- qcoT2ncp(plev,T2,p=p,n=n,ub=ub)
	# convert back
	return(f_hot2srstar(T2.ncp,n=n))
}

# use same to construct confidence intervals
srstarncp.ci <- function(srs,p,n,alpha.lo=0.025,alpha.hi=1-alpha.lo) {
	# convert to T2
	T2 <- f_srstar2hot(sample.sr=srs,n=n)
	# delegate
	T2.ci <- T2ncp.ci(T2,p=p,n=n,alpha.lo=alpha.lo,alpha.hi=alpha.hi)
	# convert back
	return(list('lo' = f_hot2srstar(T2.ci$lo,n=n),'hi' = f_hot2srstar(T2.ci$hi,n=n)))
}

#MLE of the ncp
srstarncp.mle <- function(srs,p,n,ub=NULL) {
	# convert to T2
	T2 <- f_srstar2hot(sample.sr=srs,n=n)
	if (!is.null(ub)) {
		ub <- f_srstar2hot(sample.sr=ub,n=n)
	}
	# delegate
	T2.mle <- T2ncp.mle(T2,p=p,n=n,ub=ub)
	# convert back
	return(f_hot2srstar(T2.mle,n=n))
}
#UNFOLD

## confidence intervals on optimal Sharpe#FOLDUP
#2FIX: export the guts of this...
#(1 - alpha) confidence interval on optimal SNR under assumption that portfolio
#optimizes in-sample Sharpe ratio, with n observations on p assets.
f_srstar_ci <- function(sample.sr,n,p,alpha = 0.05) {
	xval <- (n * (n-p) / (p * (n-1))) * sample.sr^2
	zerf <- function(z,n,xv,limv) { pf(xv,p,n-p,n*z^2) - limv }

	pfzero <- pf(xval, p, n-p, 0)
	if (pfzero < alpha / 2) {
		cihi <- 0
		cilo <- 0
	} else {
		#approximate upper bound 
		up <- 2 * (sample.sr + pnorm(1 - alpha / 4) / sqrt(2 * n))
		fup <- zerf(up,n,xv = xval,limv = 1 - alpha/2)
		while (fup < 0) {
			up <- 2 * up;
			fup <- zerf(up,n,xv = xval,limv = 1 - alpha/2)
		}

		if (pfzero < 1 - alpha / 2) {
			cilo <- 0
		} else {
			cilo <- uniroot(zerf,
											c(0,up),n = n,xv = xval,limv = 1 - alpha/2)$root
		}
		cihi <- uniroot(zerf,
										c(cilo,up),n = n,xv = xval,limv = alpha/2)$root
	}
	return(list('lo' = cilo,'hi' = cihi))
}
#UNFOLD

########################################################################
#utilities#FOLDUP
# logspace function
lseq <- function(from,to,length.out) { 
	exp(seq(log(from),log(to),length.out = length.out))
}
#UNFOLD

# convert SR^* <-> T2
f_srstar2hot <- function(sample.sr,n) {
	return(n * sample.sr^2)
}
f_hot2srstar <- function(T2,n) {
	return(sqrt(T2 / n))
}

########################################################################
# distributions#FOLDUP

# density, distribution, quantile, and generator for (noncentral) Hotelling
# distribution; this is just a rescaling of the (noncentral) F distribution.
dhot <- function(x, p, n, ncp = 0, log = FALSE) {
	z <- df(f_hot2F(x,p = p,n = n),df1 = p,df2 = n - p,ncp = ncp,log = log)
	if (log) {
		return(log(f_F2hot(1)) + z)
	} else {
		return(f_F2hot(z))
	}
}

phot <- function(q, p, n, ncp = 0, lower.tail = TRUE, log.p = FALSE) {
	return(pf(f_hot2F(q,p = p,n = n),df1 = p,df2 = n - p,ncp = ncp, 
						lower.tail = lower.tail, log.p = log.p))
}
qhot <- function(pct, p, n, ncp = 0, lower.tail = TRUE, log.p = FALSE) {
	return(f_F2hot(qf(pct,df1 = p,df2 = n - p,ncp = ncp, 
						lower.tail = lower.tail, log.p = log.p),p = p,n = n))
}
rhot <- function(ngen, p, n, ncp = 0) {
	return(f_F2hot(rf(ngen, df1 = p,df2 = n - p,ncp = ncp),p = p,n = n))
}

#Hotelling 
gen_hot_T2 <- function(n,p = 1,df = 10,mean = 0,sd = 1) {
	#2FIX: this was just cut and paste from gen_t and is not
  #correct .
	return(mean + sqrt((df-2)/df) * sd * rt(n,df = df))
}
#UNFOLD
########################################################################
# AR time series #FOLDUP

#generate an AR(1) series:
AR1.filter <- function(phi,epsilons) {
	n <- length(epsilons)
	series <- numeric(n)
	series[1] <- epsilons[1] / (1 - phi)
	for(iii in 2:n) {
		series[iii] <- phi * series[iii-1] + epsilons[iii]
	}
	return(series)
}

AR1.series <- function(n,phi,mu,sd,burn.in=100) {
	series <- mu + AR1.filter(phi,rnorm(n + burn.in,mean=0,sd = sd))
	return(series[(burn.in + 1):(burn.in + n)])
}
#UNFOLD
########################################################################
# a bunch of rv generators of fixed mean and sd;#FOLDUP

#gaussian
gen_norm <- rnorm

#t(10)
gen_t <- function(n,df = 10,mean = 0,sd = 1) {
	return(mean + sqrt((df-2)/df) * sd * rt(n,df = df))
}

#Tukey h
gen_tukey_h <- function(n,h = 0.1,mean = 0,sd = 1) {
	Gauss_input = create_LambertW_input("normal", beta=c(0,1))
	params = list(delta = c(h))
	LW.Gauss = create_LambertW_output(input = Gauss_input, theta = params)
	#get the moments of this distribution
	moms <- mLambertW(beta=c(0,1),distname=c("normal"),delta = h,gamma = 0, alpha = 1)
	samp <- LW.Gauss$rY(params)(n=n)
	samp <- mean  + (sd/moms$sd) * (samp - moms$mean)
}

#Lambert x Gaussian
gen_lambert_w <- function(n,dl = 0.1,mean = 0,sd = 1) {
	Gauss_input = create_LambertW_input("normal", beta=c(0,1))
	params = list(delta = c(0), gamma=c(dl), alpha = 1)
	LW.Gauss = create_LambertW_output(input = Gauss_input, theta = params)
	#get the moments of this distribution
	moms <- mLambertW(beta=c(0,1),distname=c("normal"),delta = 0,gamma = dl, alpha = 1)
	samp <- LW.Gauss$rY(params)(n=n)
	samp <- mean  + (sd/moms$sd) * (samp - moms$mean)
}

#sample from SP500
gen_SP500 <- function(n,mean = 0,sd = 1) {
	Z <- SPX.rets
	Z <- mean + (sd/sd(Z)) * (Z - mean(Z))
	samp <- sample(Z,n,replace=TRUE)
}

#sample from symmetric SP500
gen_sym_SP500 <- function(n,mean = 0,sd = 1) {
	Z <- c(SPX.rets,-SPX.rets)
	Z <- mean + (sd/sd(Z)) * (Z - mean(Z))
	samp <- sample(Z,n,replace=TRUE)
}


#generate the t-statistic of one of these things; for the normal, it is just
#a (noncentral) t; otherwise we do it via actual sampling

gen_1t_of_norm <- function(n,mean = 0,sd = 1) {
	ncp = sqrt(n) * mean / sd
	return(rt(1, df=n-1, ncp=ncp))
}

gen_1t_of_t <- function(n,df = 10,mean = 0,sd = 1) {
	return(f_tstat(gen_t(n,df = df,mean = mean,sd = sd)))
}

#etc...


#now the moments of these things
moms_norm <- function(mean = 0,sd = 1) {
	moms <- NULL
	moms$skew <- 0
	moms$kurtosis <- 0
	moms$mean <- mean
	moms$sd <- sd
	return(moms)
}

moms_t <- function(df = 10,mean = 0,sd = 1) {
	moms <- NULL
	moms$skew <- 0
	moms$kurtosis <- 6 / (df - 4)
	moms$mean <- mean
	moms$sd <- sd
	return(moms)
}

moms_tukey_h <- function(h = 0.1,mean = 0,sd = 1) {
	Gauss_input = create_LambertW_input("normal", beta=c(0,1))
	params = list(delta = c(h))
	LW.Gauss = create_LambertW_output(input = Gauss_input, theta = params)
	#get the moments of this distribution
	moms <- mLambertW(beta=c(0,1),distname=c("normal"),delta = h,gamma = 0, alpha = 1)
	moms$mean <- mean
	moms$sd <- sd
	return(moms)
}
moms_lambert_w <- function(dl = 0.1,mean = 0,sd = 1) {
	Gauss_input = create_LambertW_input("normal", beta=c(0,1))
	params = list(delta = c(0), gamma=c(dl), alpha = 1)
	LW.Gauss = create_LambertW_output(input = Gauss_input, theta = params)
	#get the moments of this distribution
	moms <- mLambertW(beta=c(0,1),distname=c("normal"),delta = 0,gamma = dl, alpha = 1)
	moms$mean <- mean
	moms$sd <- sd
	return(moms)
}
moms_SP500 <- function(mean = 0,sd = 1) {
	Z <- SPX.rets
	moms <- NULL
	moms$skew <- skewness(Z)
	moms$kurtosis <- kurtosis(Z) - 3 
	moms$mean <- mean
	moms$sd <- sd
	return(moms)
}
moms_sym_SP500 <- function(mean = 0,sd = 1) {
	Z <- c(SPX.rets,-SPX.rets)
	moms <- NULL
	moms$skew <- 0
	moms$kurtosis <- kurtosis(Z) - 3 
	moms$mean <- mean
	moms$sd <- sd
	return(moms)
}
#UNFOLD

# I told you
# "a value less than 100 shall be ignored"
mc.resolution <- max(mc.resolution,100)

## generate the bibliography#FOLDUP
##see also
##http://r.789695.n4.nabble.com/Automating-citations-in-Sweave-td872079.html

#FID <- file("rauto.bib", "w")  # open an output file connection

#cite.by.name <- function(x){ 
	#res <- toBibtex(citation(x)) 
	#if (is.list(res)) res <- res[[1]] 
	#res[1] <- sub("{",paste("{",x,sep=''),res[1],fixed=TRUE) 
	##2FIX: multiple citations; bleah;
	#cat(res,file = FID, sep = "\n")
	#return(NULL)
#} 
#z <- sapply( .packages(TRUE), function(x) try( cite.by.name(x) ) )
#close(FID)
##UNFOLD

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
