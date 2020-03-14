#-----------------------------------------------------------------
# Comparison between Sentinel-2 and WorldView-2 PSRI estimations
# Test of collinearity between sensors with emphasis in spatial resolutions
# Gerardo E. Soto - March 2017

#------------------------------------------
# Input

 library(spatstat)
 library(maptools)
 library(raster)
 library(rgdal)
 library(MASS)

 setwd("C:/Users/Gerardoku-PC/Documents/Gerardo/Trabajos/Manuscritos/Carpintero/Nidos")
 setwd("C:/Users/Gerardo/Documents/Gerardo/Trabajos/Manuscritos/Carpintero/Nest clustering")
 list.files()

 table=read.table("nidos totales.csv", header=T, sep=",")
  
 psri.wv=readGDAL("C:/Users/Gerardo/Documents/Gerardo/Trabajos/Manuscritos/Carpintero/Nest clustering/shps/psri.tif")
 psri.s2=readGDAL("C:/Users/Gerardo/Documents/Gerardo/Trabajos/Manuscritos/Carpintero/Nest clustering/shps/psri_s2a_align.tif")

 ## rasterize
 psri.wv=raster(psri.wv)
 psri.s2=raster(psri.s2)

 psri.wv[psri.wv < -10] <- NA
 psri.s2[psri.s2 < -0.75] <- NA
 psri.s2[psri.s2 > -0.35] <- NA

 # Resample psri.s2 to psri.wv
 psri.s2.samp = resample(psri.s2, psri.wv,"bilinear")

 # Points
 s2=getValues(psri.s2.samp)
 wv=getValues(psri.wv)
 plot(s2 ~ wv, pch=".")

 # Linear regression
 abline(lm(getValues(psri.s2.samp) ~ getValues(psri.wv)))

 # (Pearson) Correlation
 legend("topleft", legend=paste("Correlation =", round(cor(getValues(psri.wv),
 getValues(psri.s2.samp), use="complete.obs"), 2)))

#### Escala de Sentinel-2
 wv.samp = resample( psri.wv,psri.s2,"bilinear")
 dat1=wv.samp@data@values
 dat2=psri.s2@data@values

 ## range01 <- function(x){(x-min(x))/(max(x)-min(x))}  
 norm.wv <- function(x){(x+2.667)/(0.369+2.667)}
 d1=norm.wv(dat1)

 norm.s2 <- function(x){(x+0.75)/(-.35+0.75)}
 d2=norm.s2(dat2)

 plot(d2 ~ d1, pch=".")
 fit=lm(d2 ~ d1)
 summary(fit)

 ## Plot contour

 c=cbind(d2,d1)
 c=data.frame(c)
 c=c[complete.cases(c),]
 plot(c, pch=".")
 z <- kde2d(c[,1],c[,2], n=30)
 
 library(RColorBrewer)
 k <- 11
 my.cols <- rev(brewer.pal(k, "RdYlBu"))
 contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
 fit=lm(c[,2]~c[,1])
 abline(fit, col="red")

fit=lm(c[,2]~poly(c[,1],2))
predicted.intervals <- predict(fit,data.frame(x=c[,1]),interval='confidence',level=0.99)
lines(c[,1],predicted.intervals[,2],col='red',lwd=1)
points(c[,1], predict(fit), col="red", pch=".")
 
##======================
## Recycle
wv=rasterToPoints(psri.wv, fun=NULL, spatial=TRUE)
