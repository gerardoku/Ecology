# Mod. Pablo

file="C:/Users/Gerardo/Documents/Gerardo/Trabajos/Manuscritos/Carpintero/Seleccion Puntual/Rutas_Pablo_Agosto.csv"
write.table(tabla, file = file, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
C:/Users/Gerardo/Documents/Gerardo/Trabajos/Manuscritos/Carpintero/Uso de habitat/R/SMSU GPS

RUTAS=read.table(file="C:/Users/Gerardo/Documents/Gerardo/Trabajos/Manuscritos/Carpintero/Uso de habitat/R/SMSU GPS/Rutas_GPS.txt", sep="	", header=TRUE)
fix(RUTAS)

# [1] "BURST"      "x"         
# [3] "y"          "Time_Days" 
# [5] "Hour"       "TYPE"      
# [7] "SPECIES"    "PSRI"      
# [9] "ALT"        "SLOPE_NOR" 
#[11] "PSRI_NOR"   "ALT_NOR"   
#[13] "MODEL"      "ExtentFile"

#date = as.POSIXct(strptime(paste(RUTAS$Time_Days, RUTAS$Hour), "%d-%m-%y %H:%M"))
time.lag = rep(1,NROW(RUTAS))

library(BBMM)
BBMM = brownian.bridge(RUTAS$x, RUTAS$y, time.lag=time.lag, location.error=1, area.grid = NULL, cell.size = NULL, time.step = 10, max.lag = NULL)

da <- as.POSIXct(strptime(paste(RUTAS1$Dia, RUTAS1$Hora), "%m-%d-%y %H:%M:%S"));
xy <- cbind(XOB, YOB); 
x <- as.ltraj(xy , date = da, id = rep(1,NROW(RUTAS1$No.ruta)), typeII = TRUE)
##Sig2<-4 ##Parameter required  

lik2 <- liker(x, sig2 = Sig2, rangesig1 = c(0, 2.5), plotit =FALSE)
Sig1<-lik2[]$`1`$sig1

##Arrival.time<-5 ##Parameter requiredS
##veloc<-5  ##Parameter required (m / seg )
##MB<-2##Parameter required time extra for moving

Dist <- numeric(no.fixes-1); 
for (i in 1:(no.fixes-1)){
  Dist[i] <-dist(rbind(xy[i, ], xy[i+1, ])) 
}
tt <- Arrival.time + MB*(Dist/veloc)
temp1 <- cumsum(c(0, tt))##tiempo en sec
T <- alfa <- VarBB <- dist.tree <- Step.UT <- numeric(no.fixes-2)
MeanBB <- cbind(T,T); 

for (i in 2:(no.fixes-1)){
  T[i-1]<-temp1[i+1]-temp1[i-1]
  alfa[i-1]<-(temp1[i]-temp1[i-1])/T[i-1]
  MeanBB[i-1,1] <- XOB[i-1] + alfa[i-1]*(XOB[i+1]-XOB[i-1])
  MeanBB[i-1,2] <- YOB[i-1] + alfa[i-1]*(YOB[i+1]-YOB[i-1])
  VarBB[i-1]    <- (T[i-1] * alfa[i-1]* (1 - alfa[i-1]) * (Sig1^2)) + ((1 - alfa[i-1]) * (Sig2^2) ) + (alfa[i-1] * (Sig2^2))
  dist.tree[i-1]<- dist(rbind(MeanBB[i-1,],c(XOB[i],YOB[i]))) 
  Step.UT[i-1]<-i-1
}

text(xy[,1]+15, xy[,2]+15, as.character(1:no.fixes))
lines(XOB,YOB);points(MeanBB[,1],MeanBB[,2], col="green",pch=19)

Dat <- cbind(MeanX=MeanBB[,1],MeanY=MeanBB[,2],VarBB,Dist.centroide=dist.tree,Step.UT=Step.UT)
Dat<-rbind(rep(0,NCOL(Dat)),Dat,rep(0,NCOL(Dat)))

USED.TREE.INFO<-cbind(xy, 1:NROW(xy), Dist.trees = c(0,Dist), long.dist = c(0, ifelse(Dist>100,1,0)), Dat, 
Distance.95 = rep(0,NROW(xy)), UD.USED.TREES = rep(0,NROW(xy)), NO.USED.TREES = rep(0,NROW(xy)))
col.used <- NCOL(USED.TREE.INFO)
