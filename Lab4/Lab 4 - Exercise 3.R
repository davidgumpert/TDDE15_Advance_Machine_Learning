# GP Classification with kernlab
# (1)
library(AtmRay)
data <- read.csv("https://github.com/STIMALiU/AdvMLCourse/raw/master/GaussianProcess/Code/banknoteFraud.csv", header=FALSE, sep=",")
names(data) <- c("varWave","skewWave","kurtWave","entropyWave","fraud")
data[,5] <- as.factor(data[,5])
set.seed(111);
SelectTraining <- sample(1:dim(data)[1], size = 1000, replace = FALSE)

# Assigning training and test data
trainData <- data[SelectTraining,]
testData <- data[-SelectTraining,]

# Fit gaussian Process Classifier
gpFit <- gausspr(fraud ~ varWave + skewWave, data = trainData)

# Create grid over the varWave and skewWave variables
xGrid1 <- seq(from=min(trainData$varWave), to=max(trainData$varWave), length=100)
xGrid2 <- seq(from=min(trainData$skewWave), to=max(trainData$skewWave), length=100)
gridPoints <- meshgrid(xGrid1, xGrid2)
gridPoints <- cbind(c(gridPoints$x), c(gridPoints$y))
gridPoints <- data.frame(gridPoints)
names(gridPoints) <- c("varWave", "skewWave")

# predict the fitted gaussian process classification on the grid
predGrid <- predict(gpFit, gridPoints, type="probabilities")

# Get the indencies classified as fraud
fraud <- which(trainData$fraud == 1)

# Set up a contour of the varWave and skewWave
contour(x=xGrid1,
        y=xGrid2,
        z=matrix(predGrid[,2],100,byrow=TRUE),
        20,
        xlab="varWave",
        ylab="skewWave")

# Add the points to the grid, cases of fraud = blue, else = red
points(x=trainData$varWave[fraud],
       y=trainData$skewWave[fraud],
       col="blue")
points(x=trainData$varWave[-fraud],
       y=trainData$skewWave[-fraud],
       col="red")

# Predict through the fitted GP classification on the training data
predTrain <- predict(gpFit, trainData)

# Confusion matix
CM <- table(predictions=predTrain, true=trainData$fraud )

# Computing the accuracy
acc <- sum(diag(CM))/sum(CM)


# (2)
# Making predictions on the test set
predTest <- predict(gpFit, testData)

# Confusion matix
CM2 <- table(predictions=predTest, true=testData$fraud )

# Computing the accuracy
acc2 <- sum(diag(CM2))/sum(CM2)


# (3)
# Training a model using all 4 covariates
gpFit2 <- gausspr(fraud ~., data = trainData)
predTest2 <- predict(gpFit2, testData)
CM3 <- table(predictions=predTest2, true=testData$fraud)
acc3 <- sum(diag(CM3))/sum(CM3)
