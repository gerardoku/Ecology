## STODE ----
## Example 1. A simple sediment biogeochemical model ----

model<-function(t, y, pars)
{
  
  with (as.list(c(y, pars)),{
    
    Min       = r*OM
    oxicmin   = Min*(O2/(O2+ks))
    anoxicmin = Min*(1-O2/(O2+ks))* SO4/(SO4+ks2)
    
    dOM  = Flux - oxicmin - anoxicmin
    dO2  = -oxicmin      -2*rox*HS*(O2/(O2+ks)) + D*(BO2-O2)
    dSO4 = -0.5*anoxicmin  +rox*HS*(O2/(O2+ks)) + D*(BSO4-SO4)
    dHS  = 0.5*anoxicmin   -rox*HS*(O2/(O2+ks)) + D*(BHS-HS)
    
    list(c(dOM, dO2, dSO4, dHS), SumS = SO4+HS)
  })
}

# parameter values
pars <- c(D = 1, Flux = 100, r = 0.1, rox = 1,
          ks = 1, ks2 = 1, BO2 = 100, BSO4 = 10000, BHS = 0)
# initial conditions
y<-c(OM = 1, O2 = 1, SO4 = 1, HS = 1)

# direct iteration  - enforces  positivitiy..
ST <- stode(y = y, func = model, parms = pars, pos = TRUE)
ST



## Example 2. 1000 simultaneous equations ----

model <- function (time, OC, parms, decay, ing) {
  # model describing organic Carbon (C) in a sediment, 
  # Upper boundary = imposed flux, lower boundary = zero-gradient
  Flux  <- v * c(OC[1] ,OC) +              # advection
    -Kz*diff(c(OC[1],OC,OC[N]))/dx  # diffusion;
  Flux[1]<- flux     # imposed flux
  
  # Rate of change= Flux gradient and first-order consumption
  dOC   <- -diff(Flux)/dx - decay*OC
  
  # Fraction of OC in first 5 layers is translocated to mean depth
  dOC[1:5]  <- dOC[1:5] - ing*OC[1:5]
  dOC[N/2]  <- dOC[N/2] + ing*sum(OC[1:5])
  list(dOC)
}

v    <- 0.1    # cm/yr
flux <- 10
dx   <- 0.01
N    <- 1000 
dist <- seq(dx/2,by=dx,len=N)
Kz   <- 1  #bioturbation (diffusion), cm2/yr
print( system.time(
  ss   <- stode(runif(N), func = model, parms = NULL, positive = TRUE, 
                decay = 5, ing = 20)))

plot(ss$y[1:N], dist, ylim = rev(range(dist)), type = "l", lwd = 2,
     xlab = "Nonlocal exchange", ylab = "sediment depth",
     main = "stode, full jacobian")


## Example 3. Solving a system of linear equations ----

# this example is included to demonstrate how to use the "jactype" option.
# (and that stode is quite efficient).

A <- matrix(nrow = 500, ncol = 500, runif(500*500))
B <- 1:500

# this is how one would solve this in R
print(system.time(X1 <- solve(A, B)))

# to use stode:
# 1. create a function that receives the current estimate of x
# and that returns the difference A%*%x-b, as a list:

fun <- function (t, x, p)  # t and p are dummies here...
  list(A%*%x-B)

# 2. jfun returns the Jacobian: here this equals "A"
jfun <- function (t, x, p) # all input parameters are dummies
  A

# 3. solve with jactype="fullusr" (a full Jacobian, specified by user)
print (system.time(
  X <- stode(y = 1:500, func = fun, jactype = "fullusr", jacfunc = jfun)
))

# the results are the same (within precision)
sum((X1-X$y)^2)



steady.2D

## =======================================================================
## Diffusion in 2-D; imposed boundary conditions
## =======================================================================
diffusion2D <- function(t, Y, par)   {
  
  y    <- matrix(nr=n,nc=n,data=Y)  # vector to 2-D matrix
  dY   <- -r*y        # consumption
  BND   <- rep(1,n)   # boundary concentration 
  
  #diffusion in X-direction; boundaries=imposed concentration
  Flux <- -Dx * rbind(y[1,]-BND,(y[2:n,]-y[1:(n-1),]),BND-y[n,])/dx
  dY   <- dY - (Flux[2:(n+1),]-Flux[1:n,])/dx
  
  #diffusion in Y-direction
  Flux <- -Dy * cbind(y[,1]-BND,(y[,2:n]-y[,1:(n-1)]),BND-y[,n])/dy
  dY    <- dY - (Flux[,2:(n+1)]-Flux[,1:n])/dy
  
  return(list(as.vector(dY)))
}

# parameters
dy    <- dx <- 1   # grid size
Dy    <- Dx <- 1   # diffusion coeff, X- and Y-direction
r     <- 0.025     # consumption rate

n  <- 100
y  <- matrix(nrow = n, ncol = n, 10.)

# stodes is used, so we should specify the size of the work array (lrw)
# We take a rather large value

system.time(
  ST2 <- steady.2D(y, func = diffusion2D, parms = NULL, pos = TRUE,
                   dimens = c(n, n), lrw = 1000000, 
                   atol = 1e-10, rtol = 1e-10, ctol = 1e-10)
)

## Not run:   # this takes a long time...
system.time(
  ST3 <- steady.2D(y, func = diffusion2D, parms = NULL, 
                   dimens = c(n, n), lrw = 1000000, method = "runsteady", 
                   time = c(0, 1e6), atol = 1e-10, rtol = 1e-10)
)

## End(Not run)

# the actual size of lrw is in the attributes()$dims vector.     
# it is best to set lrw as small as possible 
attributes(ST2)     

image(ST2, legend = TRUE)

# The hard way of plotting:    
#y <- matrix(nr = n, nc = n, data = ST2$y)
#     filled.contour(y, color.palette = terrain.colors)
