## Source script of algorithms for "Soto et al. Magellanic woodpecker cavities"
## Gerardo E. Soto - Feb 2019
## CC licence: BY-NC-SA

## Gaussian Process Movement Model (Original code written by Devin Johnson <devin.johnson@noaa.gov>) ----
## Jon Horne added code for Synoptic model 9-04-2009; jhorne@uidaho.edu ----

# Create data frame from raster with distance to points ----
pointsToGrid <- function(points, extent, fun, h = 0.005){
  require(raster)
  if(fun == "distance"){
    r <- raster(ncol = 1000, nrow = 1000, ext = extent(extent))
    xy <- points %>% dplyr::select(lon, lat)
    coordinates(xy) <- xy %>% dplyr::select(lon, lat)
    proj4string(xy) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    grid <- distanceFromPoints(r, xy, plot = FALSE)
    # grid <- data.frame(lon = coordinates(dist)[, 1], lat = coordinates(dist)[, 2], dist = values(dist))
  } else {
    if(fun == "kernel"){
      require(MASS)
      k <- kde2d(points$lon, points$lat, h = h, n = 1000, lims = extent)
      grid <- raster(k)
    } else {
      print("no specified function!")
    }
  }
  return(grid)
}


# Weighted distributions Model ----
bmmle = function(track, 
                 habmat, 
                 paramCovs, 
                 Hessian = TRUE, 
                 maxit = 100)
{
  require(MASS)
  require(optimParallel)
  nms = names(paramCovs)
  stime = Sys.time()
  cl <- makeCluster(3)  # Select number of clusters 
  setDefaultCluster(cl=cl)  # set clusters 
  mle.bm = optimParallel(paramCovs, 
                 sbvnLogLik, 
                 track = track, 
                 habmat = habmat,
                 method = "L-BFGS-B",
                 lower = c(-10),
                 upper = c(10),
                 control = list(maxit = maxit),
                 hessian = Hessian,
                 parallel=list(loginfo=TRUE))
  if(Hessian){
    covmat = 2*ginv(mle.bm$hessian)
    se = sqrt(diag(covmat))
  }else{
    covmat = NULL
    se = rep(NA, 
             length(mle.bm$par))
  }       
  z = mle.bm$par/se
  p.val = 2*(1-pnorm(abs(z)))
  parTable = cbind(Est = mle.bm$par, 
                   SE = se, 
                   Lower = mle.bm$par-1.96*se,
                   Upper = mle.bm$par+1.96*se, 
                   Z = z, 
                   P = p.val)
  rownames(parTable) = c(colnames(habmat)[-c(1:2)])
  list(parTable = zapsmall(parTable), 
       covmat = covmat,
       AIC = mle.bm$value + 2*length(paramCovs),
       evalTime = difftime(Sys.time(),stime),
       convergence = (mle.bm$convergence==0),
       print(mle.bm))
}

# habmat <- data.frame(x = coordinates(psriNava_sp)[,1], y = coordinates(psriNava_sp)[,2], psri = values(psriNava_sp))
# habmat <- habmat %>% filter(!is.na(psri))
# track <- useNava %>% filter(ind == 1) %>% select(lon, lat)
# psri <- extract(psriNava_sp, track)
# track <- track %>% mutate(psri = psri) %>% filter(!is.na(psri))

### el track tiene que tener las covariables también!!!
# paramCovs <- 0.5

# Likelihood function for bivariate normal ----
sbvnLogLik = function(paramCovs, track, habmat)
{
  # needed to apply pipes
  if (!require(dplyr)) install.packages('dplyr')
  library(dplyr)
  
  # . Bivariate normal joint probability density function ----
  bvnjPDF <- function(x, meanx, meany, sdx, sdy, corr){
    if(x %>% is.vector){
      bvnjPDF_firTerm = ((2 * pi) * sdx * sdy * sqrt(1 - corr ^ 2))
      bvnjPDF_secTerm = -1 / (2 * (1 - corr ^ 2))
      bvnjPDF_thiTerm = (((x[1] - meanx) / sdx) ^ 2) + (((x[2] - meany) / sdy) ^ 2) - 
        (2 * corr * ((x[1] - meanx) / sdx) * ((x[2] - meany) / sdy))
      bvnjPDF = exp(bvnjPDF_secTerm * bvnjPDF_thiTerm) / bvnjPDF_firTerm  
    }else{
      bvnjPDF_firTerm = ((2 * pi) * sdx * sdy * sqrt(1 - corr ^ 2))
      bvnjPDF_secTerm = -1 / (2 * (1 - corr ^ 2))
      bvnjPDF_thiTerm = (((x[, 1] - meanx) / sdx) ^ 2) + (((x[, 2] - meany) / sdy) ^ 2) - 
                         (2 * corr * ((x[, 1] - meanx) / sdx) * ((x[, 2] - meany) / sdy))
      bvnjPDF = exp(bvnjPDF_secTerm * bvnjPDF_thiTerm) / bvnjPDF_firTerm  
    }
    return(bvnjPDF)
  }
  
  # . Function parameters from track ----
  meanx = mean(track[, 1])
  meany = mean(track[, 2])
  sdx = sd(track[, 1])
  sdy = sd(track[, 2])
  corr = cor(track[, 1:2])[1, 2]
  # meanx = paramSBVN[1]
  # meany = paramSBVN[2]
  # sdx = exp(paramSBVN[3])
  # sdy = exp(paramSBVN[4])
  # corr = paramSBVN[5] 
  
  # . Calculate volume under non-normalized use function ----
  bvnjPDFmap <- bvnjPDF(habmat, meanx, meany, sdx, sdy, corr)
  
  # . . Exponential Selection Function ----
  if (length(paramCovs) == 0){
    wMap = 1
  } else {
    if(length(paramCovs) == 1){
      wMap = exp(habmat[, 3] * paramCovs[1])
    } else {
      wMap = exp(habmat[, 3:ncol(habmat)] %*% paramCovs[1:length(paramCovs)])
    }
  }
  MapNonNorm.g.u = bvnjPDFmap * wMap
  
  dimX <- ((habmat[, 1] %>% max) - (habmat[, 1] %>% min)) / habmat[, 1] %>% unique %>% length
  dimY <- ((habmat[, 2] %>% max) - (habmat[, 2] %>% min)) / habmat[, 2] %>% unique %>% length
  cellsize = as.numeric(dimX * dimY)  
  MapVolume = sum(MapNonNorm.g.u * cellsize)
  
  # . Loop over locations ----
  # . Likelihood elements to populate  
  SumLogLik = 0
  LogBvnjPDFi = array(0, nrow(track))
  
  for (i in 1:nrow(track)){  
    # Exponential selection function
    if(length(paramCovs) == 0){
      wi = 1
    } else {
      if(length(paramCovs) == 1){
        wi = exp(track[i, -(1:2)] * paramCovs[1])
      } else {
        wi = exp(track[i, -(1:2)] %*% paramCovs[1:length(paramCovs)])
      }
    }
    bvnjPDFi <- bvnjPDF(track[i, 1:2], meanx, meany, sdx, sdy, corr)  # what is this doing?
    density = array(bvnjPDFi * wi / MapVolume)
    if(is.na(density)) density = 10 ^ -320 
    if(density == 0) density = 10 ^ -320 
    LogBvnjPDFi[i] = log(density)
    SumLogLik = SumLogLik + LogBvnjPDFi[i]
  } #end loop through locations
  
  -2 * SumLogLik
}


# Likelihood function for brownian motion (bridges) ----
bmLogLik = function(paramBM, 
                    track, 
                    tDiff, 
                    habmat, 
                    cc.fac = 6){ # cc.fac: bm window factor
  # paramBM = (log(L), alpha1, ..., alphaK)
  Nu = nrow(track)
  wMap = exp(habmat[, -c(1:2)] * paramBM[2:length(paramBM)]) 
  w = exp(track[,-(1:3)] * paramBM[2:length(paramBM)])
  maxy = max(habmat[,'x'])
  maxx = max(habmat[,'y'])
  Cinv = matrix(c(exp(-0.5 * paramBM[1]), 0, 0, exp(-0.5 * paramBM[1])), 2, 2) ## Choice
  loglik = 0
  for(i in 2:Nu){ # Estimation of the likelihood values for each observation i
    mui = track[i-1,1:2]
    Cinvi = Cinv / sqrt(tDiff[i - 1]) # 
    win = ceiling(cc.fac * exp(0.5 * paramBM[1]) * sqrt(tDiff[i-1])) # window for the gaussian process
    cc = cbind(rep(seq(-win, win, 1), each = 2*win+1), rep(seq(-win, win, 1), 2*win+1))
    cci = cbind(cc[,1]+track[i-1,1], cc[,2]+track[i-1,2]) 
    cci = cci[(cci[,1] >= 1 & cci[,1] <= maxx) & (cci[,2] >= 1 & cci[,2] <= maxy),] 
    indi = cci[,2] + (cci[,1]-1) * maxy 
    zi = cbind(cci[,1]-as.numeric(mui[1]), cci[,2]-as.numeric(mui[2])) %*% Cinvi
    gai = exp(-0.5 * ((zi*zi) %*% c(1,1))) ## ga at present location 
    gaitr = exp(-0.5 * crossprod(Cinvi %*% t(track[i,1:2]-mui)))  ## ga at trajectory
    loglik = loglik + log(gaitr) + log(w[i]) - log(crossprod(gai, wMap[indi])) 
  }
  -2*loglik
}
