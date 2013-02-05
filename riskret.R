# * Mon Feb 04 2013 03:18:12 PM Steven E. Pav <steven@cerebellumcapital.com>
# generate risk/return points to show the risk return space?

safe.install <- function(pkg.name,do.load=FALSE) {
	if (! length(which(.packages(all.available=TRUE) %in% pkg.name))) {
		install.packages(pkg.name,repos=c("http://cran.cnr.berkeley.edu/",
																			"http://R-Forge.R-project.org"))
	}
	if (do.load) {
		library(pkg.name,character.only=TRUE)
	}
}

dummy <- lapply(c("tikzDevice","grDevices"),safe.install,do.load=TRUE)

require(tikzDevice)
require(grDevices)

# set up population#FOLDUP
set.seed(118423)

#npnt <- 64
npnt <- 512
epy <- 253

#sr <- -0.6 + abs(rnorm(npnt,mean=0.5,sd=0.666))
sr <- rnorm(npnt,mean=0.3,sd=0.600)
xdf <- 5
sg <- 0.002 + 0.0090 * sqrt((1/xdf) * rchisq(npnt, xdf))
mu <- sg * sr / sqrt(epy)

# *very* optimistic this:
rfr <- 2.00
#UNFOLD

# plot func#FOLDUP
plot.all <- function(mu,sg,rfr,epy=253) {
	par(new=FALSE,usr=c(0,30,-20,40))

	xval <- 100 * sqrt(epy) * sg
	yval <-100 * epy * mu

	xat <- seq(0,40,by=10)
	yat <- seq(-30,50,by=15)


	# use xaxp and yaxp to control the x/y tick points?
	plot(100 * sqrt(epy) * sg,100 * epy * mu,
			 xlim=xat[c(1,length(xat))],
			 ylim=yat[c(1,length(yat))],
			 xaxt="n",yaxt="n",
			 xlab='annualized volatility',ylab='annualized return')
	abline(h=0,col='green',lty='dashed')

	tlab <- sapply(xat,function(dd) { sprintf('%g \\%%',dd) })
	axis(1,at=xat,label=tlab)

	tlab <- sapply(yat,function(dd) { sprintf('%g \\%%',dd) })
	axis(2,at=yat,label=tlab)

	points(0,rfr,type='p',col='red')

	x0 <- quantile(xval,probs=0.25)
	y0 <- quantile(yval,probs=0.76)
	x1 <- x0 + 2.1 * (quantile(xval,probs=0.01) - x0)
	#y1 <- y0 + 0.8 * (quantile(yval,probs=0.99) - y0)
	y1 <- y0 - (x1 - x0)

	# convex hull AKA efficient frontier#FOLDUP
	foohull <- chull(c(xval,max(xval),max(yval)),c(yval,max(yval),min(yval)))
	foohull <- foohull[-length(foohull)]
	foohull <- foohull[-length(foohull)]

	foook <- diff(xval[foohull]) >= 0
	tophull <- foohull[c(foook,TRUE)]

	#plot(xval,yval)
	lines(xval[tophull],yval[tophull],lwd=2,col='cyan')
	#UNFOLD
	
	arrows(x0,y0,x1,y1,lwd=2,col='blue',length=0.075)

	text(x0 + 0.55 * (x1-x0),
			 y0 + 0.9 * (y1-y0),
			 "preference",col='blue')

	#title(main="risk reward")
}#UNFOLD

#plot.all(mu,sg,rfr,epy)

#require(tikzDevice)
# The following wwill create riskret.tex in the working
# directory the first time this is run it may take a long time because the
# process of calulating string widths for proper placement is
# computationally intensive, the results will get cached for the current R
# session or will get permenantly cached if you set
# options( tikzMetricsDictionary='/path/to/dictionary' ) which will be
# created if it does not exist.  Also if the flag standAlone is not set to
# TRUE then a file is created which can be included with \include{}
#tikz('riskret.tex', standAlone = TRUE, width=5, height=5)
tikz('riskret.tex', width=5, height=5)

plot.all(mu,sg,rfr,epy)

#Close the device
dev.off()

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
