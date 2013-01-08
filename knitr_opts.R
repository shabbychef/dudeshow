# 
# * Mon Dec 31 2012 11:50:47 AM Steven E. Pav <steven@cerebellumcapital.com>
#
# base knitr definitions for dude paper

require(knitr)

# set the knitr options ... for everyone!
opts_knit$set(echo=FALSE)
opts_knit$set(eps=TRUE)
opts_knit$set(warning=FALSE)
opts_knit$set(message=FALSE)
opts_knit$set(cache=TRUE)
opts_knit$set(cache.path="cache/")
opts_knit$set(results="asis")
opts_knit$set(fig.ext="eps")
opts_knit$set(fig.path="figure/")
opts_knit$set(dev=c("pdf","cairo_ps"))

