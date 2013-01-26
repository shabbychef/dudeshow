# 
# * Mon Dec 31 2012 11:50:47 AM Steven E. Pav <steven@cerebellumcapital.com>
#
# base knitr definitions for dude paper

require(knitr)

# set the knitr options ... for everyone!
opts_knit$set(progress=TRUE)

opts_chunk$set(echo=FALSE,warning=FALSE,message=FALSE,print=FALSE)
opts_chunk$set(results="asis")
opts_chunk$set(cache=TRUE,cache.path="cache/")
opts_chunk$set(fig.path="figure/",dev=c("pdf","cairo_ps"))
#opts_knit$set(eps=TRUE)


