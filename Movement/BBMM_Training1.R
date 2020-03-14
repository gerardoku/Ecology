##install.packages("adehabitatHR");install.packages("adehabitatLT");install.packages("raster");install.packages("maptools");install.packages("rgdal")
##install.packages("chron");install.packages("adehabitat")
library(adehabitatHR)
library(chron)
library(adehabitatHR)
library(raster)
library(maptools)
library(rgdal)
library(adehabitatLT)
library(adehabitat)

LISTA<-read.table(file="C:/Users/Pablo/Documents/Manuscritos/woodpecker/BBRW/LISTA.csv", sep=",", header=TRUE)
fix(LISTA)

RUTAS<-read.table(file="C:/Users/Pablo/Documents/Manuscritos/woodpecker/BBRW/Rutas.BBSM1.csv", sep=",", header=TRUE)
RUTAS<-read.table(file="C:/Users/Pablo/Documents/Manuscritos/woodpecker/BBRW/Rutas.BBSM1.txt", sep="\t", header=TRUE)
fix(RUTAS)
tapply(rep(1,NROW(RUTAS)),RUTAS$No.ruta,sum) 

MIN.X<-min(RUTAS$UTME)-100 ;MAX.X<- max(RUTAS$UTME)+100; MIN.Y<-min(RUTAS$UTMN)-100; MAX.Y<-max(RUTAS$UTMN)+100
LISTA.RES<-LISTA[LISTA$X>MIN.X & LISTA$X<MAX.X & LISTA$Y>MIN.Y & LISTA$Y<MAX.Y,]
   sum(ifelse(LISTA.RES$X<MIN.X,1,0));sum(ifelse(LISTA.RES$X>MAX.X,1,0));sum(ifelse(LISTA.RES$Y<MIN.Y,1,0));sum(ifelse(LISTA.RES$Y>MAX.Y,1,0))
   ## LISTA$X	LISTA$Y	Tree.ID	SPP	TIPO	PSRI

dts <- as.character(RUTAS$Dia); tms <- as.character(times(RUTAS$Hora))

x <- chron(dates = dts, times = tms, format = c(dates = "m-d-y", times = "h:m:s"))
times.hour <- hours(x);times.min <- minutes(x);times.sec <- seconds(x); times.dat<-round(times.hour + (times.min/60) + (times.sec/3600),4)
da <- paste(dts, tms); RUTAS<-cbind(RUTAS,times.dat,da); No.RUTAS<-max(RUTAS$No.ruta)



##################################################

FIRST.R<-11##first route
LAST.R<-15##last route
Arrival.time<-5 ##Time spent when arrival
MB<-2##Extra time when moving between trees
Max.BNorm<-100 ## 0.5 x dimension of matrix used for fit Bivariate Normal kernel, lower values make faster programming
Sig2<-4  #standard deviation observation error
veloc<-10
UMBRAL<-0.90

BBSM.OUTPUT<- function(RUTAS,LISTA.RES, FIRST.R, LAST.R, Sig2,Arrival.time,MB,Max.BNorm,veloc,UMBRAL) {

USED.TREE.INFO.TOT<-c()
BBSM.TOT<-c()
SHARED.TREES<-c()
for (h in FIRST.R:LAST.R){
##h<-14
       RUTAS1<-subset(RUTAS,RUTAS$No.ruta==h)
        XOB<-RUTAS1$UTME ; YOB<-RUTAS1$UTMN;  no.fixes <-NROW(XOB); fixes<- 1:no.fixes; tree<-RUTAS1$No.Arbol
        temp<-RUTAS1$times.dat; temp<-60*(temp-temp[1])##tiempo acumulado en min

##establish the landscape area surrounding the route
land<-matrix(0,ncol=200+max(XOB)-min(XOB),nrow=200+max(YOB)-min(YOB),byrow=F)
land<-as.asc(land)  
land.X0<-min(XOB)-100; land.Y0<-min(YOB)-100
attr(land,"xll")<-land.X0; attr(land,"yll")<-land.Y0     
attr(land,"cellsize")<-1
ha.tot<-dim(land)[1]*dim(land)[2]/10000; Num.tree<-500*ha.tot

##load trees in the study area 
MIN.X<-min(RUTAS1$UTME)-100 ;MAX.X<- max(RUTAS1$UTME)+100; MIN.Y<-min(RUTAS1$UTMN)-100; MAX.Y<-max(RUTAS1$UTMN)+100
OB.TREES<-LISTA.RES[LISTA.RES$X>MIN.X & LISTA.RES$X<MAX.X & LISTA.RES$Y>MIN.Y & LISTA.RES$Y<MAX.Y,]
sum(ifelse(OB.TREES$X<MIN.X,1,0));sum(ifelse(OB.TREES$X>MAX.X,1,0));sum(ifelse(OB.TREES$Y<MIN.Y,1,0));sum(ifelse(OB.TREES$Y>MAX.Y,1,0))
image(land, xlim=c(max(XOB)+100,min(XOB-100)), col="white",ylim=c(max(YOB)+100,min(YOB-100))); 
points(OB.TREES[,1],OB.TREES[,2], pch=21, col="gray");points(XOB,YOB, cex=1.4, col="red", pch=19); lines(XOB,YOB,lwd=2.5, col="red")


##VARIANCE

da <- as.POSIXct(strptime(paste(RUTAS1$Dia, RUTAS1$Hora), "%m-%d-%y %H:%M:%S"));
xy<-cbind(XOB,YOB); x <- as.ltraj(xy , date = da, id=rep(1,NROW(RUTAS1$No.ruta)),typeII=TRUE)
##Sig2<-4 ##Parameter required  
lik2 <- liker(x, sig2 = Sig2, rangesig1 = c(0, 2.5),plotit =FALSE); Sig1<-lik2[]$`1`$sig1
##Arrival.time<-5 ##Parameter requiredS
##veloc<-5  ##Parameter required (m / seg )
##MB<-2##Parameter required time extra for moving
Dist<-numeric(no.fixes-1); for (i in 1:(no.fixes-1)){ Dist[i] <-dist(rbind(xy[i,],xy[i+1,])) }; tt<-Arrival.time + MB*(Dist/veloc)
temp1<-cumsum(c(0,tt))##tiempo en sec
T <- alfa <- VarBB<-dist.tree <- Step.UT <- numeric(no.fixes-2);MeanBB<-cbind(T,T); 

for (i in 2:(no.fixes-1)){
T[i-1]<-temp1[i+1]-temp1[i-1]
alfa[i-1]<-(temp1[i]-temp1[i-1])/T[i-1]
MeanBB[i-1,1] <- XOB[i-1] + alfa[i-1]*(XOB[i+1]-XOB[i-1])
MeanBB[i-1,2] <- YOB[i-1] + alfa[i-1]*(YOB[i+1]-YOB[i-1])
VarBB[i-1]    <- (T[i-1] * alfa[i-1]* (1 - alfa[i-1]) * (Sig1^2)) + ((1 - alfa[i-1]) * (Sig2^2) ) + (alfa[i-1] * (Sig2^2))
dist.tree[i-1]<- dist(rbind(MeanBB[i-1,],c(XOB[i],YOB[i]))) 
Step.UT[i-1]<-i-1
}

text(xy[,1]+15, xy[,2]+15, as.character(1:no.fixes));lines(XOB,YOB);points(MeanBB[,1],MeanBB[,2], col="green",pch=19)

Dat<-cbind(MeanX=MeanBB[,1],MeanY=MeanBB[,2],VarBB,Dist.centroide=dist.tree,Step.UT=Step.UT);Dat<-rbind(rep(0,NCOL(Dat)),Dat,rep(0,NCOL(Dat)))

USED.TREE.INFO<-cbind(xy, 1:NROW(xy),Dist.trees=c(0,Dist), long.dist=c(0,ifelse(Dist>100,1,0)),Dat, 
Distance.95=rep(0,NROW(xy)),UD.USED.TREES=rep(0,NROW(xy)),NO.USED.TREES=rep(0,NROW(xy)) )
col.used<-NCOL(USED.TREE.INFO)

##########################################
#####################

BBSM<-c()

for (i in 1:NROW(MeanBB)){

##Max.BNorm<-100
DIX<-ifelse(dist.tree[i]>Max.BNorm,dist.tree[i], Max.BNorm); mu <- c(DIX,DIX)
var<-VarBB[i]
Sigma <- matrix(c(var,0.5, 0.5,  var), 2)

xmult<-cbind(rep(seq(1,DIX*2,1),(DIX*2)),rep(seq(1,DIX*2,1),each=(DIX*2)))
XDEN<-PDdMvn(xmult,mu,Sigma)
dmat<-matrix(XDEN[,3],ncol=DIX*2,nrow=DIX*2,byrow=F)
##95%
bb<-XDEN[sort.list(XDEN[,3],decreasing=TRUE), ]
cum.bb<-cumsum(bb[,3])
##plot(1:NROW(bb),cumsum(bb[,3]))
Dist.95<-(sum(ifelse(cum.bb>UMBRAL,0,1))^0.5)/2
Dist.95

##Filtering trees < Dist.95
OB.TREES1<-OB.TREES
OB.TREES1[,1]<-round(OB.TREES[,1]-MeanBB[i,1]+DIX,0)
OB.TREES1[,2]<-round(OB.TREES[,2]-MeanBB[i,2]+DIX,0)
OB.TREES1<-OB.TREES1[OB.TREES1[,1]<max(XDEN[,1])& OB.TREES1[,1]>min(XDEN[,1])&
           OB.TREES1[,2]<max(XDEN[,2])& OB.TREES1[,2]>min(XDEN[,2]),]
Dist.max<-numeric(NROW(OB.TREES1))
for (j in 1:NROW(OB.TREES1)){Dist.max[j] <-dist(rbind(mu,OB.TREES1[j,1:2]))  }
OB.TREES1<-subset(OB.TREES1,Dist.max<Dist.95); 

##Extract UD from available trees (< Dist.95)
NULL.UD0<-numeric(NROW(OB.TREES1))
for (j in 1:NROW(OB.TREES1)){NULL.UD0[j] <-dmat[OB.TREES1[j,1],OB.TREES1[j,2]]  }
USO<-rep(0,NROW(OB.TREES1)); OB.TREES1<-as.data.frame(cbind(OB.TREES1,NULL.UD0,USO)); 
OB.TREES1[,1]<-OB.TREES1[,1]+MeanBB[i,1]-DIX
OB.TREES1[,2]<-OB.TREES1[,2]+MeanBB[i,2]-DIX

##Extract UD from the used tree (< Dist.95)
USED.TREE<-numeric(dim(OB.TREES1)[2])
NULL.UD1 <-dmat[round(XOB[i+1]-MeanBB[i,1]+DIX,0),round(YOB[i+1]-MeanBB[i,2]+DIX,0)]  
if(OB.TREES[,1]==round(XOB[i+1],0)&OB.TREES[,2]==round(YOB[i+1],0)){
USED.TREE[1:2]<-OB.TREES[OB.TREES[,1]==round(XOB[i+1],0)&OB.TREES[,2]==round(YOB[i+1],0),]} else {
USED.TREE[1:2]<-c(XOB[i+1],YOB[i+1])}
USED.TREE[3]<-c(Tree.ID=0) ##add used tree ID
USED.TREE[4]<-c(SPP=0) ##add SPP
USED.TREE[5]<-c(TIPO=0) ##add TIPO
USED.TREE[6]<-c(PSRI=0) ##add PSRI
USED.TREE[(NROW(USED.TREE)-1):NROW(USED.TREE)]<-c(NULL.UD1,1)
options("scipen"=100, "digits"=4)

TREES.UA<-rbind(USED.TREE,OB.TREES1)
No.TREES<-NROW(TREES.UA)
##Standarize for Null UD
TREES.UA$NULL.UD0<- TREES.UA$NULL.UD0 / sum(TREES.UA$NULL.UD0 )

STEP<-rep(i,NROW(TREES.UA))
ROUTE<-rep(h,NROW(TREES.UA))
TREES.UA<-	cbind(TREES.UA,STEwP,ROUTE)
BBSM<-rbind(BBSM,TREES.UA)
USED.TREE.INFO[i+1,c(col.used-2,col.used-1,col.used)]<-c(Dist.95,NULL.UD1,No.TREES) 

}

BBSM.TOT<-rbind(BBSM.TOT,BBSM)
ROUTE<-rep(h,NROW(USED.TREE.INFO))
USED.TREE.INFO<-	cbind(USED.TREE.INFO,ROUTE)
USED.TREE.INFO.TOT<-rbind(USED.TREE.INFO.TOT,USED.TREE.INFO)
}
##how many trees are shaded between steps
SHARED.TREES<-c()
rut<-c();stp<-c()
for (g in min(BBSM.TOT$ROUTE):max(BBSM.TOT$ROUTE)){
AA<-subset(BBSM.TOT, BBSM.TOT$ROUTE==g)
if(max(AA$STEP)>1){
rut<-c(rut,rep(g,max(AA$STEP)-1))
for(j in 2:max(AA$STEP)){
AA1<-c(-1,subset(AA, AA$STEP==j)$Tree.ID)
AA2<-c(-1,subset(AA, AA$STEP==(j-1))$Tree.ID)
SHARED.TREES<-c(SHARED.TREES,NROW(intersect(AA1,AA2))-1)
stp<-c(stp,j)}}
else{
rut<-c(rut,g)
SHARED.TREES<-c(SHARED.TREES,0)
stp<-c(stp,1)
}}
SHARED.TREES<-cbind(SHARED.TREES,rut,stp)

return(list(BBSM.TOT=BBSM.TOT,USED.TREE.INFO.TOT=USED.TREE.INFO.TOT,SHARED.TREES=SHARED.TREES))
}

BBSm11<-BBSM.OUTPUT(RUTAS,LISTA.RES, FIRST.R=1, LAST.R=5, Sig2=3,Arrival.time=3,MB=2,Max.BNorm=100,veloc=20,UMBRAL=0.9)
BBSm12<-BBSM.OUTPUT(RUTAS,LISTA.RES, FIRST.R=6, LAST.R=10, Sig2=2,Arrival.time=3,MB=2,Max.BNorm=100,veloc=20,UMBRAL=0.9)
BBSm13<-BBSM.OUTPUT(RUTAS,LISTA.RES, FIRST.R=11, LAST.R=15, Sig2=2,Arrival.time=3,MB=2,Max.BNorm=100,veloc=20,UMBRAL=0.9)
BBSm14<-BBSM.OUTPUT(RUTAS,LISTA.RES, FIRST.R=16, LAST.R=20, Sig2=2,Arrival.time=3,MB=2,Max.BNorm=100,veloc=20,UMBRAL=0.9)
BBSm15<-BBSM.OUTPUT(RUTAS,LISTA.RES, FIRST.R=21, LAST.R=25, Sig2=2,Arrival.time=3,MB=2,Max.BNorm=100,veloc=20,UMBRAL=0.9)
BBSm16<-BBSM.OUTPUT(RUTAS,LISTA.RES, FIRST.R=26, LAST.R=30, Sig2=2,Arrival.time=3,MB=2,Max.BNorm=100,veloc=20,UMBRAL=0.9)
BBSm17<-BBSM.OUTPUT(RUTAS,LISTA.RES, FIRST.R=31, LAST.R=35, Sig2=2,Arrival.time=3,MB=2,Max.BNorm=300,veloc=20,UMBRAL=0.9)
BBSm18<-BBSM.OUTPUT(RUTAS,LISTA.RES, FIRST.R=36, LAST.R=40, Sig2=2,Arrival.time=3,MB=2,Max.BNorm=100,veloc=20,UMBRAL=0.9)
BBSm19<-BBSM.OUTPUT(RUTAS,LISTA.RES, FIRST.R=41, LAST.R=44, Sig2=2,Arrival.time=3,MB=2,Max.BNorm=100,veloc=20,UMBRAL=0.9)

#########################################################################
BBSm13<-BBSM.OUTPUT(RUTAS,LISTA.RES, FIRST.R=1, LAST.R=5, Sig2=3,Arrival.time=3,MB=2,Max.BNorm=100,veloc=20,UMBRAL=0.9)
BBSm12$BBSM.TOT
BBSm12$USED.TREE.INFO.TOT
BBSm12$SHARED.TREES

BBSm2<-BBSM.OUTPUT(RUTAS,LISTA.RES, FIRST.R=6, LAST.R=10, Sig2=2,Arrival.time=5,MB=2,Max.BNorm=100,veloc=10)
BBSm2$USED.TREE.INFO.TOT
BBSm2$BBSM.TOT
BBSm2$SHARED.TREES

BBSm3<-BBSM.OUTPUT(RUTAS,LISTA.RES, FIRST.R=11, LAST.R=15, Sig2=2,Arrival.time=5,MB=2,Max.BNorm=100,veloc=10)
BBSm3$USED.TREE.INFO.TOT
BBSm3$BBSM.TOT
BBSm3$SHARED.TREES

BBSm4<-BBSM.OUTPUT(RUTAS,LISTA.RES, FIRST.R=16, LAST.R=20, Sig2=2,Arrival.time=5,MB=2,Max.BNorm=100,veloc=10)
BBSm4$USED.TREE.INFO.TOT
BBSm4$BBSM.TOT
BBSm4$SHARED.TREES

BBSm5<-BBSM.OUTPUT(RUTAS,LISTA.RES, FIRST.R=21, LAST.R=25, Sig2=2,Arrival.time=5,MB=2,Max.BNorm=100,veloc=10)
BBSm5$USED.TREE.INFO.TOT
BBSm5$BBSM.TOT
BBSm5$SHARED.TREES

BBSm6<-BBSM.OUTPUT(RUTAS,LISTA.RES, FIRST.R=26, LAST.R=30, Sig2=2,Arrival.time=5,MB=2,Max.BNorm=100,veloc=10)
BBSm6$USED.TREE.INFO.TOT
BBSm6$BBSM.TOT
BBSm6$SHARED.TREES


BBSm7<-BBSM.OUTPUT(RUTAS,LISTA.RES, FIRST.R=33, LAST.R=33, Sig2=2,Arrival.time=5,MB=2,Max.BNorm=270,veloc=10)
BBSm7$USED.TREE.INFO.TOT
BBSm7$BBSM.TOT
BBSm7$SHARED.TREES

BBSm8<-BBSM.OUTPUT(RUTAS,LISTA.RES, FIRST.R=36, LAST.R=40, Sig2=2,Arrival.time=5,MB=2,Max.BNorm=100,veloc=10)
BBSm8$USED.TREE.INFO.TOT
BBSm8$BBSM.TOT
BBSm8$SHARED.TREES

BBSm9<-BBSM.OUTPUT(RUTAS,LISTA.RES, FIRST.R=41, LAST.R=44, Sig2=2,Arrival.time=5,MB=2,Max.BNorm=100,veloc=10)
BBSm9$USED.TREE.INFO.TOT
BBSm9$BBSM.TOT
BBSm9$SHARED.TREES

BBSm19$USED.TREE.INFO.TOT

Type.forest<-readGDAL("C:/Users/Pablo/Documents/Manuscritos/woodpecker/BBRW/RST_TIPO.TIF")
summary(Type.forest)
image(Type.forest, col = grey(1:99/100), axes = TRUE)
points(MeanBB[,1],MeanBB[,2], col="red", cex=5)
Type.forest1<-as.matrix(Type.forest)
max(Type.forest1,na.rm = TRUE);min(Type.forest1,na.rm = TRUE)

tt<-cbind(BBSm1$USED.TREE.INFO.TOT[,"NO.USED.TREES"],BBSm12$USED.TREE.INFO.TOT[,"NO.USED.TREES"],BBSm13$USED.TREE.INFO.TOT[,"NO.USED.TREES"])

apply(tt,2,sum)

cbind(1:NROW(BBSm1$BBSM.TOT),BBSm1$BBSM.TOT)
BBSm13$USED.TREE.INFO.TOT
