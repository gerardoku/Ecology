#setwd("C:/Users/LEC-i7/Documents/GERARDO/Trabajos/Manuscritos/Carpintero/Uso de habitat/R/17mayo")
setwd("C:/Users/LEC/Documents/GERARDO/Trabajos/Manuscritos/Carpintero/Uso de habitat/R/AvailConTodasLasVariables")
setwd("C:/Users/Gerardo/Documents/GERARDO/Trabajos/Manuscritos/Carpintero/Uso de habitat/R/AvailConTodasLasVariables")

# Require libraries
library(rgeos)
library(raster) 
library(tiff)
library(SDMTools)
library(rgdal)
library(adehabitat)
library(adehabitatHR)
library(corrgram)

#Files Avail
AvailFileList<-c("Avail_Burst1.txt", "Avail_Burst2.txt", "Avail_Burst3.txt", "Avail_Burst4.txt", "Avail_Burst5.txt","Avail_Burst6.txt")

#Files locations
LocationsFileList<-c("Locs1.txt", "Locs2.txt", "Locs3.txt", "Locs4.txt", "Locs5.txt","Locs6.txt")

# =======================================================================================
# Leer Tablas

#Open files
a=Sys.time()
Habmat<-list()
Track<-list()
for (IndividNum in 1:length(LocationsFileList)){
	file <- AvailFileList[IndividNum]
	Habmat[[IndividNum]] <- data.frame(read.table(file=file, head=TRUE, sep="	", row.names=NULL)) 
	head(Habmat[[IndividNum]])
	locationsfile <- LocationsFileList[IndividNum]
	Track[[IndividNum]] <- data.frame(read.table(file=locationsfile, head=TRUE, sep="	", row.names=NULL)) 
	head(Track[[IndividNum]])
}
Sys.time()-a

# =======================================================================================
# Exportar Tablas

for(IndividNum in 1:length(LocationsFileList)){
fileOUT=paste("Locs",IndividNum,".txt", sep="")
write.table(Track[[IndividNum]], file=fileOUT, row.names = FALSE, col.names = TRUE)
fileO=paste("Avail_Burst",IndividNum,".txt", sep="")
write.table(Habmat[[IndividNum]], file=fileO, row.names = FALSE, col.names = TRUE)
}

# =======================================================================================
# Extraer DH

NAvail <- data.frame(read.table(file="NAvail.txt", head=TRUE, sep="	", row.names=NULL)) 
head(NAvail)

wd="C:/Users/LEC/Documents/GERARDO/Trabajos/Manuscritos/Carpintero/Uso de habitat/R/AvailConTodasLasVariables"
distance=raster(readGDAL(paste(wd,"/rasternulo/distancias_P11.tif", sep="")))

a<-Sys.time()
for (IndividNum in 1:length(LocationsFileList)){
  xy<-Habmat[[IndividNum]][,1:2]
  DH=extract(distance,xy)
  Habmat[[IndividNum]]=cbind(Habmat[[IndividNum]],DH)
}
for (IndividNum in 1:length(LocationsFileList)){
  xy<-Track[[IndividNum]][,1:2]
  DH=extract(distance,xy)
  Track[[IndividNum]]=cbind(Track[[IndividNum]][,1:(NCOL(Track[[IndividNum]])-2)],DH,Track[[IndividNum]][,(NCOL(Track[[IndividNum]])-1):NCOL(Track[[IndividNum]])])
}
Sys.time()-a


# =======================================================================================
# Calcular bordes de acuerdo a distancias de árboles al borde

p=list()
for(IndividNum in 1:length(LocationsFileList)){
a=Habmat[[IndividNum]]["DH"]
p[[IndividNum]]=nrow(subset(a, DH<5))/nrow(Habmat[[IndividNum]])
}

# =======================================================================================
# Calcular bordes de acuerdo a raster de clasificación

clas<- raster(readGDAL("Supervisada_4Tierr_B_Np23_TIFF.tif"))
for(IndividNum in 1:length(LocationsFileList)){

}

# =======================================================================================
# 

for(IndividNum in 1:length(LocationsFileList)){

}


# =======================================================================================
# Interacciónes de variables medidas sobre PSRI

for (IndividNum in 1:length(LocationsFileList)){

### Cálculo de interacciones en locs
locationsfile <- LocationsFileList[IndividNum]
track <- data.frame(read.table(file=locationsfile, head=TRUE, sep="	", row.names=NULL)) 
head(track)

#Separar variables
TY1<-track$TY1
TY2<-track$TY2
SP1<-track$SP1
SP2<-track$SP2
SP3<-track$SP3
PSRI<-track$PSRI_NOR

#TY*SP=TS 
TS1<-TY1*SP1
TS2<-TY1*SP2
TS3<-TY1*SP3

#SP*PSRI=SPS
SPS1<-PSRI*SP1
SPS2<-PSRI*SP2
SPS3<-PSRI*SP3

#TS*PSRI=TSP
TSP1<-TS1*PSRI
TSP2<-TS2*PSRI
TSP3<-TS3*PSRI

#TY*PSRI=TP
TP1<-TY1*PSRI
TP2<-TY2*PSRI

track<-cbind(track[1:(NCOL(track)-2)],TS1,TS2,TS3,SPS1,SPS2,SPS3,TSP1,TSP2,TSP3,TP1,TP2,track[(NCOL(track)-1):NCOL(track)])
head(track)

fileOUT=paste("Locs",IndividNum,".txt", sep="")
write.table(track, file=fileOUT, row.names = FALSE, col.names = TRUE)

### Cálculo de interacciones en disponibilidad
habmat <- AvailFileList[IndividNum]
habmat <- data.frame(read.table(file=habmat, head=TRUE, sep="	", row.names=NULL)) 
head(habmat)

#Separar variables
TY1<-habmat$TY1
TY2<-habmat$TY2
SP1<-habmat$SP1
SP2<-habmat$SP2
SP3<-habmat$SP3
PSRI<-habmat$PSRI_NOR

#TY*SP=TS 
TS1<-TY1*SP1
TS2<-TY1*SP2
TS3<-TY1*SP3

#SP*PSRI=SPS
SPS1<-PSRI*SP1
SPS2<-PSRI*SP2
SPS3<-PSRI*SP3

#TS*PSRI=TSP
TSP1<-TS1*PSRI
TSP2<-TS2*PSRI
TSP3<-TS3*PSRI

#TY*PSRI=TP
TP1<-TY1*PSRI
TP2<-TY2*PSRI

habmat<-cbind(habmat,TS1,TS2,TS3,SPS1,SPS2,SPS3,TSP1,TSP2,TSP3,TP1,TP2)
head(habmat)

fileOUT=paste("Avail_Burst",IndividNum,".txt", sep="")
write.table(habmat, file="NAvail.txt", row.names = FALSE, col.names = TRUE)
}

# =======================================================================================
# Interacciónes de variables medidas a escala de HR

#------------------------------------------------------------------------------
# Calculo de Áreas de HRs
#------------------------------------------------------------------------------
area<-rep(NA,length(LocationsFileList))
for (IndividNum in 1:length(LocationsFileList)){
	track<-Track[[IndividNum]]
	polygons<-list()
	xy<-SpatialPoints(track[,1:2], proj4string=CRS("+proj=utm +zone=19 +datum=WGS84"))
	ker<-kernelUD(xy, h = "LSCV", grid = 400, same4all = FALSE, hlim = c(0.1, 100), kern = c("bivnorm", "epa"), extent = 5)				
	area[IndividNum]<-kernel.area(ker, percent = 95, unin = "m", unout = "ha", standardize = FALSE)	
	polygons[[IndividNum]]<-getverticeshr(ker,95)
}

#------------------------------------------------------------------------------
# Calculo de suma y promedio de PSRI
#------------------------------------------------------------------------------
sum.mean<-matrix(NA, ncol=(2*(NCOL(Habmat[[1]])-2)), nrow=length(LocationsFileList), dimnames=NULL)
for (IndividNum in 1:length(LocationsFileList)){
	habmat<-Habmat[[IndividNum]]
	mean.cov<-apply(habmat[,-(1:2)], 2, mean)
	sum.cov<-apply(habmat[,-(1:2)], 2, sum)
	sum.mean[IndividNum,]<-c(mean.cov,sum.cov)
	names<-colnames(habmat[3:NCOL(habmat)])
	names.mean<-rep(0,length(names))
	names.sum<-rep(0,length(names))
	for(i in 1:length(names)){
		names.mean[i]<-paste(names[i], "Mean", sep="")
		names.sum[i]<-paste(names[i], "Sum", sep="")
	}
}
colnames(sum.mean)<-c(names.mean,names.sum)
sum.mean

#------------------------------------------------------------------------------
# Multiplicación de cov(HR) con cov(tree)
#------------------------------------------------------------------------------

for (IndividNum in 1:length(LocationsFileList)){
	track<-Track[[IndividNum]][,-c((ncol(Track[[IndividNum]])-1):ncol(Track[[IndividNum]]))]
	habmat<-Habmat[[IndividNum]]
	inter<-cbind(habmat[-(1:2)],habmat[-(1:2)])
	intert<-cbind(track[-(1:2)],track[-(1:2)])
	iarea<-habmat[-(1:2)]*area[IndividNum]
	iareat<-track[-(1:3)]*area[IndividNum]
	imean.psri<-habmat[-(1:2)]*sum.mean[IndividNum,"PSRI_NORMean"]
	imean.slop<-habmat[-(1:2)]*sum.mean[IndividNum,"SLOPE_NORMean"]
	imeant.psri<-track[-(1:3)]*sum.mean[IndividNum,"PSRI_NORMean"]
	imeant.slop<-track[-(1:3)]*sum.mean[IndividNum,"SLOPE_NORMean"]	
	prop.ty<-cbind(habmat[,3]*sum.mean[IndividNum,1],habmat[,4]*sum.mean[IndividNum,2])
	prop.ty.names<-c("Prop.Ty1","Prop.Ty2")
	propt.ty<-cbind(track[,4]*sum.mean[IndividNum,1],track[,5]*sum.mean[IndividNum,2])
	Habmat[[IndividNum]]<-cbind(habmat,iarea,imean.psri,imean.slop,prop.ty)
	Track[[IndividNum]]<-cbind(track,iareat,imeant.psri,imeant.slop,propt.ty,Track[[IndividNum]][,c((ncol(Track[[IndividNum]])-1):ncol(Track[[IndividNum]]))])
}

covs<-colnames(track[-(1:3)])
ar<-vector();ps<-vector();sl<-vector()
for(h in 1:length(covs)){
	ar<-c(ar,paste("area.",covs[h],sep=""))
	ps<-c(ps,paste("psri.",covs[h],sep=""))
	sl<-c(sl,paste("slop.",covs[h],sep=""))
}

for (IndividNum in 1:length(LocationsFileList)){
colnames(Habmat[[IndividNum]])<-c(colnames(habmat),ar,ps,sl,prop.ty.names)
colnames(Track[[IndividNum]])<-c(colnames(track),ar,ps,sl,prop.ty.names,colnames(Track[[IndividNum]][,c((ncol(Track[[IndividNum]])-1):ncol(Track[[IndividNum]]))]))
}


#------------------------------------------------------------------------------
# Crear nuevo locs y disponibilidad
#------------------------------------------------------------------------------

for (IndividNum in 1:length(LocationsFileList)){
	fileOUT=paste("Avail_Burst",IndividNum,".txt", sep="")
	filetOUT=paste("Locs",IndividNum,".txt", sep="")
	write.table(Habmat[[IndividNum]], file=fileOUT, row.names = FALSE, col.names = TRUE)
	#write.table(Track[[IndividNum]], file=filetOUT, row.names = FALSE, col.names = TRUE)
}

# =======================================================================================
# Agregar variable familia a Locs y Disp

a<-matrix(0,1,length(LocationsFileList))
for (IndividNum in 1:length(LocationsFileList)){
	a<-matrix(0,1,length(LocationsFileList))
	a[,IndividNum]<-1
	Habmat[[IndividNum]]<-cbind(Habmat[[IndividNum]],a)
	Track[[IndividNum]]<-cbind(Track[[IndividNum]][,1:(ncol(Track[[IndividNum]])-2)],a,Track[[IndividNum]][,(ncol(Track[[IndividNum]])-1):ncol(Track[[IndividNum]])])
}

#### Alternativa con puros 0s en disponibilidad
a<-matrix(0,1,length(LocationsFileList))
for (IndividNum in 1:length(LocationsFileList)){
	Habmat[[IndividNum]][,(ncol(Habmat[[IndividNum]])-5):ncol(Habmat[[IndividNum]])]<-rep(0,length(LocationsFileList))
}


