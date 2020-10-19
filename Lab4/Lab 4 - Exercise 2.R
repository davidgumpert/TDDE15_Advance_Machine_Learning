library(kernlab)
data <- read.csv("https://github.com/STIMALiU/AdvMLCourse/raw/master/GaussianProcess/Code/TempTullinge.csv", header=TRUE, sep=";")

time <- seq(1, nrow(data))
temp <- data$temp
count <- 1
day <- c()
for(i in 1:nrow(data)) {
  if(count > 365) count <- 1
  day <- rbind(day, count)
  count <- count + 1
}

timeSparse <- time[seq(1, length(time), 5)]
daySparse <- day[seq(1, length(day), 5)]
tempSparse <- data$temp[seq(1, length(data$temp), 5)]

# Exercise 2 - GP Regression with kernlab
# (1)
# First define your own square exponential function with parameters ell
# and sigmaf. Evaluate in point x=1, x'=2. Use the kernalMatrix function
# to compute the covariance matrixK(X,Xstar) for the input vectors
# X = t(1, 3, 4) and Xstar= t(2,3,4).

# Setting initial parameters
x <- 1
xPrime <- 2
X <- c(1,3,4)
Xstar <- c(2,3,4)

# The covariance function
SquaredExpKernel <- function(x1, x2, sigmaF, l) {
  n1 <- length(x1)
  n2 <- length(x2)
  K <- matrix(NA,n1,n2)
  for (i in 1:n2){
    K[,i] <- sigmaF^2*exp(-0.5*( (x1-x2[i])/l)^2 )
  }
  return(K)
}

# Defining the SE Kernal for the kernelMatrix function
SEKernel <- function(sigmaF=1, l=3) {
  InternalKernel <- function(X, Xstar) {
    n1 <- length(X)
    n2 <- length(Xstar)
    K <- matrix(NA,n1,n2)
    for (i in 1:n2){
      K[,i] <- sigmaF^2*exp(-0.5*( (X-Xstar[i])/l)^2 )
    }
    return(K)
  }
  class(InternalKernel) <- 'kernel'
  return(InternalKernel)
}

# Following demo in KernLabDemo
# kernel <- rbfdot(sigma = 1/(2*ell^2))
seKernel <- SEKernel()
covMat1 <- seKernel(1,2)
covMat1

# Now we use it in the kernelMatrix function, creating the comvariance matrix
# for the two vectors X and Xstar
covMat2 <- kernelMatrix(kernel=seKernel, x=X, y=Xstar )
covMat2

# (2)
# Setting initial parameters
sigmaF <- 20
ell  <- 0.2
seKernel <- SEKernel(sigmaF, ell)
quadraticRegressionFit <- lm(tempSparse ~ timeSparse + timeSparse^2)
sigmaN <- sd(quadraticRegressionFit$residuals)

# Estimating the given Gaussian process regression model using the 
# squared exponential function from (1)
fittedModel <- gausspr(timeSparse, tempSparse, kernel=seKernel, var=sigmaN^2)

# Computing posterior mean at every datapoint in the sparse data set with
# the predict function
postMean <- predict(fittedModel, timeSparse)

# Plotting the combination plot (scatter for the datapoints and posterior mean
# of f as a curve)
plot(timeSparse, tempSparse, ylab="Temp", xlab="Time",
     main="Scatter of data and curve of posterior mean of f")
lines(timeSparse, postMean, col="red", lwd=2)  

# (3)
# When plotting the 96% confidential intervals we need to use our own
# computations of the posterior variance of f due to a bug in the 
# Kernlab library.

# Reusing the PosteriorGP function from exercise 1, removed explanatory comments
PosteriorGP <- function(X, y, XStar, sigmaNoise, k, ...) {
  N <- length(X)
  K <- k(X,X,...)
  KStar <- k(X,XStar,...)
  L <- t(chol(K + sigmaNoise^2*diag(N)))
  alpha <- solve(t(L), solve(L,y))
  fBarStar <- t(KStar)%*%alpha 
  v <- solve(L,KStar)
  fStarVar <- k(XStar, XStar, ...) - t(v)%*%v
  postMeanPred <- fBarStar
  postVarPred <- fStarVar
  return(list(mean=postMeanPred, var=postVarPred))
} 

# Reusing the plotting function from exercise 1
PlotSim <- function(observations, grid, mean, var, confInt = TRUE) {
  lowerY <- mean-1.96*sqrt(diag(var))
  upperY <- mean+1.96*sqrt(diag(var))
  plot(grid, mean, type = "l",
       ylim=c(min(lowerY), max(upperY)), col="darkblue", 
       main="Posterior mean with 95% probability bands",
       xlab="Time step", ylab="Temperature"
  )
  if(confInt) {
    lines(grid, lowerY, col="red", lwd=2, lty=21)
    lines(grid, upperY, col="red", lwd=2, lty=21)
  }
  points(observations[,1], observations[,2])
}

# Initiate kernel
seKernel <- SEKernel(sigmaF, ell)

# Computing my own variance
posterior <- PosteriorGP(scale(timeSparse), scale(tempSparse),
                         scale(timeSparse), sigmaN, seKernel)

# Plotting the sparse data with own variance and posterior mean from previous
# exercise
PlotSim(data.frame(x=timeSparse, y=tempSparse), 
        timeSparse, postMean, posterior$var)

# (4)
# Here we fit the gaussian process to the day varibale instead 
# of the time variable. 
sigmaF <- 20
ell <- 0.2
seKernel <- SEKernel(sigmaF, ell)
# Computing the sigmaN^2 as in exercise 2.2
quadraticRegressionFit <- lm(tempSparse ~ daySparse + daySparse^2)
sigmaN <- sd(quadraticRegressionFit$residuals)

fittedModelDay <- gausspr(x = daySparse, y = tempSparse, kernel = seKernel,
                          var = sigmaN^2)
postMeanDay <- predict(fittedModelDay, daySparse)

# Add line for the day prediction
PlotSim(data.frame(x=timeSparse, y=tempSparse), 
        timeSparse, postMean, posterior$var, confInt = FALSE)
lines(x = timeSparse, y = postMeanDay, col = "green")


# (5)
# Here we will use a general periodic kernel
# See the GeneralPeriodicKernel function

# Setting parameters
sigmaF <- 20
l1 <- 1
l2 <- 10
d <- 365/sd(timeSparse)

# Creating the general periodic kernel
GeneralPeriodicKernel <- function(sigmaF=1, l1=1, l2=10, d) {
  InternalKernel <- function(X, Xstar) {
    diff <- X-Xstar
    sigmaF^2*exp(-(2*(sin(pi*diff/d)^2)/(l1^2)))*exp(-(0.5*diff^2)/(l2^2))
  }
  class(InternalKernel) <- 'kernel'
  return(InternalKernel)
}
periodicKernel <- GeneralPeriodicKernel(sigmaF, l1, l2, d)
fittedModelGeneral <- gausspr(x = timeSparse, y = tempSparse,
                              kernel = periodicKernel,
                              var = sigmaN^2)
postMeanPeriodic <- predict(fittedModelGeneral, timeSparse)
lines(timeSparse, postMeanPeriodic, type = "l", col = "red")

