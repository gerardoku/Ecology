#------------------------------------------------------------------------------
# Disponibilidad Synoptic	
# Gerardo E. Soto - Mar 2015

# AKDE estimation for GPS data
# Gerardo E. Soto - Mar 2017

#------------------------------------------------------------------------------
# Paquetes
#------------------------------------------------------------------------------
 library(rgeos)
 library(raster) 
 library(tiff)
 library(SDMTools)
 library(maptools)
 library(rgdal)
 library(adehabitat)
 library(adehabitatHR)
 library(ctmm)
 library(move)

#------------------------------------------------------------------------------
# Input Locs
#------------------------------------------------------------------------------

 workdir = "C:/Users/Gerardo/Documents/Gerardo/Trabajos/Manuscritos/Carpintero/Nest clustering/R"
 setwd(workdir)

 #()# VHF
	d.VHF.p<-read.table(file="DatosVHF.csv", header=TRUE, sep=",")
	d.VHF.p$DAY<-gsub("/","-",d.VHF.p$DAY) 
      d.VHF.p[d.VHF.p[,"Ind"]==6,"Ind"]=3 					### MERGE TITA & LENG
	d.VHF.p[d.VHF.p[,"Ind"]==7,"Ind"]=6 					### CHANGE WACH TO #6
	time=as.POSIXct(paste(d.VHF.p$DAY,d.VHF.p$H, sep=" "), format="%d-%m-%Y %H:%M", tz="UTC")
	d.VHF<-data.frame(d.VHF.p$Ind,d.VHF.p$X,d.VHF.p$Y,time)
	d.VHF<-data.frame(d.VHF, rep("VHF", NROW(d.VHF)))
	colnames(d.VHF)<-c("id", "X", "Y", "date", "type")

 #()# ATS
	d.ATS<-read.table(file="DatosATS.csv", header=TRUE)
	d.ATS$X<-round(d.ATS$X, digits = 0);d.ATS$Y<-round(d.ATS$Y, digits = 0)
	d.ATS<-cbind(d.ATS, rep("ATS", NROW(d.ATS)))
	colnames(d.ATS)<-c("id", "X", "Y", "date", "type")
 	d.ATS[,1]=d.ATS[,1]+6

 #()# Lotek
	d.Lotek<-read.table(file="DatosLotek.csv", header=TRUE)
	d.Lotek$X<-round(d.Lotek$X, digits = 0);d.Lotek$Y<-round(d.Lotek$Y, digits = 0)
	d.Lotek[,"date"]=as.POSIXct(d.Lotek[,"date"], format="%d-%m-%Y %H:%M:%S", tz="UTC")
	d.Lotek<-data.frame(d.Lotek, rep("Lotek", NROW(d.Lotek)))
	colnames(d.Lotek)<-c("id", "X", "Y", "date", "type")
	d.Lotek[,1]=d.Lotek[,1]+6


#------------------------------------------------------------------------------
# AKDE para todos
#------------------------------------------------------------------------------
 
 ## Input function
 as.telem=function(x){
  		sputm <- SpatialPoints(x[,2:3], proj4string=CRS("+proj=utm +zone=19 +datum=WGS84"))
 		spgeo <- spTransform(sputm, CRS("+proj=longlat +datum=WGS84"))
 		x.latlong=coordinates(spgeo)
		l=length(x[,1])
		x=data.frame(eventid=x[,1], visible=rep(TRUE,l), timestamp=as.POSIXct(x[,4], format="%Y-%m-%d %H:%M:%S", tz="UTC")+3*60*60, ## diferencia de 3 horas con UTC
			latlong=x.latlong, manuallymarkedoutlier=rep(NA,l), sensortype=rep("gps",l), individualtaxoncanonicalname=rep("Campephilus magellanicus",l),
			taglocalidentifier=rep(1,l), individuallocalidentifier=rep("x",l), studyname=rep("Navarino",l), 
			utm=x[,2:3], utmzone=rep("19S",l), studytimezone=rep("Argentina Time",l), 
			studylocaltimestamp=as.POSIXct(x[,4], format="%Y-%m-%d %H:%M:%S", tz="UTC")
 			)
  		colnames(x)=c("event.id","visible","timestamp","location.long","location.lat","manually.marked.outlier",
 			"sensor.type","individual.taxon.canonical.name","tag.local.identifier","individual.local.identifier",
 			"study.name","utm.easting","utm.northing","utm.zone","study.timezone","study.local.timestamp")
 		#head(x)
		x=as.telemetry(x)
		return(x)
 } 
 
 telemetry=list() ## Input lists
 for(i in 1:6) telemetry[[i]]=as.telem(d.VHF[d.VHF[,1]==i,])
 for(i in 7:24) telemetry[[i]]=as.telem(d.ATS[d.ATS[,1]==i,])
 for(i in 25:30) telemetry[[i]]=as.telem(d.Lotek[d.Lotek[,1]==i,])


 ## AKDE VHFs
 UD=list() ## list of UDs
 for(i in 1:6){
  m2 <- ctmm(tau=c(200*24*60,60*24)*60) # ~ 365 days and 24*60 minutes autocorrelation timescales
  M2 <- ctmm.fit(telemetry[[i]])#,m2)
  UD[[i]] <- akde(telemetry[[i]],M2)
 } 

 ## AKDE Zaña id=9
 m2 <- ctmm(tau=c(7*24*60,5)*60) # ~ 7 days and 5 minutes autocorrelation timescales
 M2 <- ctmm.fit(telemetry[[9]],m2)
 UD.ATS9 <- akde(telemetry[[9]],M2)
  
 ## PLOT results
 par(mfrow=c(3,3))
 for(i in 1:6) plot(telemetry[[i]],UD=UD[[i]], main=paste("VHF",i,sep=""))  
 plot(telemetry[[9]],UD=UD.ATS9, main="ATS9")


 ## Exportar Rasters ??



#------------------------------------------------------------------------------
# AKDE para Zañartu
#------------------------------------------------------------------------------
 
 points(d.ATS[d.ATS[,1]==3+6,2:3]) 		#zañartu
 points(d.ATS[d.ATS[,1]==9,2:3], pch=2) 	#titan
 points(d.ATS[d.ATS[,1]==12,2:3]) 		#oscar

 zaña=d.ATS[d.ATS[,1]==3+6,2:4]
 zaña=data.frame(zaña)

 #Coordenadas lat/lon
 sputm <- SpatialPoints(d.ATS[d.ATS[,1]==3,2:3], proj4string=CRS("+proj=utm +zone=19 +datum=WGS84")  )
 spgeo <- spTransform(sputm, CRS("+proj=longlat +datum=WGS84"))
 zaña.latlong=coordinates(spgeo)

 #Formato MoveBank
 zaña=data.frame(eventid=1, visible=TRUE, timestamp=as.POSIXct(zaña[,3], format="%Y-%m-%d %H:%M:%S", tz="UTC")+3*60*60,
		latlong=zaña.latlong, manuallymarkedoutlier=NA, sensortype="gps", individualtaxoncanonicalname="Campephilus magellanicus",
		taglocalidentifier=1, individuallocalidentifier="zaña", studyname="Navarino", 
		utm=d.ATS[d.ATS[,1]==3,2:3], utmzone="19S", studytimezone="Argentina Time", 
		studylocaltimestamp=as.POSIXct(zaña[,3], format="%Y-%m-%d %H:%M:%S", tz="UTC")
 		)
 
 colnames(zaña)=c("event.id","visible","timestamp","location.long","location.lat","manually.marked.outlier",
 	"sensor.type","individual.taxon.canonical.name","tag.local.identifier","individual.local.identifier",
 	"study.name","utm.easting","utm.northing","utm.zone","study.timezone","study.local.timestamp")
 head(zaña)
 
 ## Input AKDE 
 zaña=as.telemetry(zaña)

 ## AKDE
 m2 <- ctmm(tau=c(7*24*60,5)*60) # ~ 7 days and 5 minutes autocorrelation timescales
 M2 <- ctmm.fit(zaña,m2)
 UD2 <- akde(zaña,M2)
 UD2w <- akde(zaña,M2,weights=TRUE)
 
 ## PLOT results
 plot(zaña,UD=UD2)



#------------------------------------------------------------------------------
# Covariates for SMSU
#------------------------------------------------------------------------------

 workdir = "C:/Users/Gerardo/Documents/Gerardo/Trabajos/Manuscritos/Carpintero/Nest clustering/shps"
 setwd(workdir)


 #------------------------------------------------------------------------------
 # Avail files
 
 nidos=raster(readGDAL("C:/Users/Gerardo/Documents/Gerardo/Trabajos/Manuscritos/Carpintero/Nest clustering/kernel_bw150_nidos.tif"))
 psri=crop(raster(readGDAL("psri.tif")), extent(nidos))
 ptos=readOGR(dsn = "C:/Users/Gerardo/Documents/Gerardo/Trabajos/Manuscritos/Carpintero/Nest clustering/shps", layer = "ptos_base")

 p.psri=extract(psri, ptos, buffer=10, fun=mean)
 p.nido=extract(nidos, ptos, buffer=10, fun=mean)
 p.nido=p.nido*10000

 #respaldo
 # p.psri.r=p.psri
 # p.nido.r=p.nido

 #Normalize whole maps
 norm.maps=function(x){
  y=(x-min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T))
  return(y)
 }
 p.nido=norm.maps(p.nido)
 p.psri=norm.maps(p.psri)


 avail=data.frame(ptos@coords,p.psri,p.nido)
 colnames(avail)=c("x", "y", "psri", "nido")
 head(avail)
 
 #plot(avail[avail[,"psri"]!="NA",1:2], pch=".") 
 #plot(avail[avail[,"nido"]>0,1:2], pch=".") 
 
 # Main Avail file
 write.table(avail,file="Avail.txt",sep=",",dec=".", row.names=F)

 #------------------------------------------------------------------------------
 # Locs files

 psri.locs=list()
 nido.locs=list()
 avail.locs=list()
 locs.covs=list()
 for(i in 1:6){
  #extract values form rasters
  coords=SpatialPoints(d.VHF[d.VHF[,1]==i,2:3], proj4string=CRS("+proj=utm +zone=19 +datum=WGS84"))
  psri.locs[[i]]=extract(psri, coords, buffer=10, fun=mean)
  nido.locs[[i]]=extract(nidos, coords, buffer=10, fun=mean)

  locs.covs[[i]]=data.frame(d.VHF[d.VHF[,1]==i,2:3],psri.locs[[i]],nido.locs[[i]])
  colnames(locs.covs[[i]])=c("x", "y", "psri", "nido")
  write.table(locs.covs[[i]],file=paste("locs.covs_",i,".txt", sep=""),sep=",",dec=".", row.names=F)

 #clip Avail files for each individual
  xmin=extent(coords)[][1]
  xmax=extent(coords)[][2]
  ymin=extent(coords)[][3]
  ymax=extent(coords)[][4]

  avail.locs[[i]]=avail[avail[,1]>(xmin-100)&avail[,1]<(xmax+100)&avail[,2]<(ymax+100)&avail[,2]>(ymin-100),]
  write.table(avail.locs[[i]],file=paste("Avail_locs_",i,".txt", sep=""),sep=",",dec=".", row.names=F)
 } 

  # GPS zaña
  coords=SpatialPoints(d.ATS[d.ATS==9,2:3], proj4string=CRS("+proj=utm +zone=19 +datum=WGS84"))
  ext=as(extent(nidos), "SpatialPolygons")
  proj4string(ext)=CRS("+proj=utm +zone=19 +datum=WGS84")
  coords=gIntersection(coords, ext, byid = F)
  psri.locs[[7]]=extract(psri, coords, buffer=10, fun=mean)
  nido.locs[[7]]=extract(nidos, coords, buffer=10, fun=mean)


  locs.covs[[7]]=data.frame(coords@coords,psri.locs[[7]],nido.locs[[7]])
  colnames(locs.covs[[7]])=c("x", "y", "psri", "nido")
  write.table(locs.covs[[7]],file=paste("locs.covs_",7,".txt", sep=""),sep=",",dec=".", row.names=F)

  #clip Avail files for zaña
  xmin=extent(coords)[][1]
  xmax=extent(coords)[][2]
  ymin=extent(coords)[][3]
  ymax=extent(coords)[][4]

  avail.locs[[7]]=avail[avail[,1]>(xmin-100)&avail[,1]<(xmax+100)&avail[,2]<(ymax+100)&avail[,2]>(ymin-100),]
  write.table(avail.locs[[7]],file="Avail_locs_7.txt",sep=",",dec=".", row.names=F)


 for(i in 1:7){
  print(summary(avail.locs[[i]]))
 }
 

 #------------------------------------------------------------------------------
 # Arreglar la omisión de normalización para locs

 for(i in 1:7){
 c.avail=avail.locs[[i]] 
  for(j in 1:nrow(locs.covs[[i]])){
   coord=locs.covs[[i]][j,1:2] 
   c.cov=c.avail[c.avail[,1]>(coord[,1]-6)&c.avail[,1]<(coord[,1]+6)&c.avail[,2]>(coord[,2]-6)&c.avail[,2]<(coord[,2]+6),][,3:4]
   locs.covs[[i]][j,3:4]=c.cov[1,] #Solo el primer registro si es que hay duplicados
  }
  write.table(locs.covs[[i]],file=paste("locs.covs_",i,".txt", sep=""),sep=",",dec=".", row.names=F)
 } 
 
 
 











