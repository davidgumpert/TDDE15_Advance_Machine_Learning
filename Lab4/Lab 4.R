library(mvtnorm)

# Covariance function
# Takes x1 and x2 as input and computes the exponential kernal, giving the
# resulting covariance between two values
SquaredExpKernel <- function(x1,x2,sigmaF=1,l=3){
  n1 <- length(x1)
  n2 <- length(x2)
  K <- matrix(NA,n1,n2)
  for (i in 1:n2){
    K[,i] <- sigmaF^2*exp(-0.5*( (x1-x2[i])/l)^2 )
  }
  return(K)
}

# X=training inputs
# y=training labels
# XStar=inputs where the posterior distribution is evaluated
# sigmaNoise=noise standard deviation
# k=covariance function of kernel
PosteriorGP <- function(X, y, XStar, sigmaNoise, k, ...) {
  # The number of functions we need to generate, same as the number of inputs
  N <- length(X)
  
  # Computing the covariance matrix [K = K(X,X)]
  K <- k(X,X,...)
  KStar <- k(X,XStar,...)
  
  # Computing L with cholesky decomposition
  # Taking the transpose of chol because it returns the upper 
  # triangulated matrix and we want the lower
  L <- t(chol(K + sigmaNoise^2*diag(N)))
  
  # Computing the predictive mean fBar*
  # Gotten by multiplying transpose of KStar with alpha
  # alpha = t(L)\(L\y)
  alpha <- solve(t(L), solve(L,y))
  fBarStar <- t(KStar)%*%alpha 
  
  # Calculating the predictive variance
  v <- solve(L,KStar)
  fStarVar <- k(XStar, XStar, ...) - t(v)%*%v
  
  # Saving result in appropriat named variables
  postMeanPred <- fBarStar
  postVarPred <- fStarVar
  
  return(list(mean=postMeanPred, var=postVarPred))
} 

# Function for plotting the posterior simulation
PlotSim <- function(observations, grid, mean, var) {
  lowerY <- mean-1.96*sqrt(diag(var))
  upperY <- mean+1.96*sqrt(diag(var))
  plot(grid, mean, type = "l",
       ylim=c(min(lowerY), max(upperY)), col="darkblue", 
       main="Posterior mean with 95% probability bands",
       xlab="x", ylab="Posterior mean"
       )
  lines(grid, lowerY, col="red", lwd=2, lty=21)
  lines(grid, upperY, col="red", lwd=2, lty=21)
  points(observations[,1], observations[,2])
}

# Exercise 1
# (2) - Update prior with single observation

# Setting priors for hyperparameters and other variables
xGrid <- seq(-1,1,length=100)
sigmaF <- 1
sigmaN <- 0.1
l <- 0.3
observations <- data.frame(x=0.4, y=0.719)

# Drawing one simulation
postSim <- PosteriorGP(observations$x,
                       observations$y,
                       xGrid,
                       sigmaN,
                       SquaredExpKernel,
                       sigmaF,
                       l)

# Plotting one simulation
PlotSim(observations, xGrid, postSim$mean, postSim$var)


# (3) - Update the posterior with two sample observations
x <- c(0.4, -0.6)
y <- c(0.719, -0.044)
observations <- data.frame(x, y)

# Drawing two simulations
postSim <- PosteriorGP(observations$x,
                       observations$y,
                       xGrid,
                       sigmaN,
                       SquaredExpKernel,
                       sigmaF,
                       l)

# Plotting two simulations
PlotSim(observations, xGrid, postSim$mean, postSim$var)

# (4) - Updating the posterior using all 5 observations
x <- c(-1, -0.6, -0.2, 0.4, 0.8)
y <- c(0.768, -0.044, -0.940, 0.719, -0.664)
observations <- data.frame(x,y)

# Drawing from the posterior
postSim <- PosteriorGP(observations$x,
                       observations$y,
                       xGrid,
                       sigmaN,
                       SquaredExpKernel,
                       sigmaF,
                       l)
PlotSim(observations, xGrid, postSim$mean, postSim$var)

# (5) - Changing hyperparameter l
l <- 1
postSim <- PosteriorGP(observations$x,
                       observations$y,
                       xGrid,
                       sigmaN,
                       SquaredExpKernel,
                       sigmaF,
                       l)
PlotSim(observations, xGrid, postSim$mean, postSim$var)

