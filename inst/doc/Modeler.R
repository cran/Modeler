### R code from vignette source 'Modeler.Rnw'

###################################################
### code chunk number 1: lib
###################################################
library(Modeler)


###################################################
### code chunk number 2: seed
###################################################
set.seed(234843)


###################################################
### code chunk number 3: params
###################################################
nFeatures <- 10000
nSignif <- 100
pB <- 0.4
delta <- 1
sigma <- 0.3
nTrain <- 100
nTest <- 100


###################################################
### code chunk number 4: paramlist
###################################################
paramlist <- c("nFeatures", "nSignif", "pB",
               "delta", "sigma", "nTrain", "nTest")


###################################################
### code chunk number 5: alpha.beta
###################################################
alpha <- 0.75
beta <- 0.95
round(100*pbeta(seq(0.1, 0.9, 0.1), alpha, beta), 1)
xx <- seq(0, 1, length=300)
yy <- dbeta(xx, alpha, beta)


###################################################
### code chunk number 6: Modeler.Rnw:81-82
###################################################
plot(xx, yy, type='l')


###################################################
### code chunk number 7: signed
###################################################
signed <- -1 + 2*rbinom(nSignif, 1, 0.5)


###################################################
### code chunk number 8: offsets
###################################################
offsets <- c(signed*rnorm(nSignif, delta, sigma), # can change in either direction
             rep(0, nFeatures - nSignif))         # but most don't change at all


###################################################
### code chunk number 9: out
###################################################
lp <- function(p) log(p/(1-p))
ea <- function(a) 1/(1+exp(-a))
pOut <- rbeta(nTrain, alpha, beta)
trainOutcome <- lp(pOut)


###################################################
### code chunk number 10: class
###################################################
# TODO: Fix this so it looks at correlation with the continuous outcome
# instead of just diferential expression between classes
trainClass <- factor(c("cyan", "magenta")[1 + 1*(pOut > 0.5)])
summary(trainClass)
isB <- trainClass=="magenta"
summary(isB)


###################################################
### code chunk number 11: train
###################################################
trainData <- matrix(rnorm(nFeatures*nTrain), ncol=nTrain) # pure noise
trainData[,isB] <- sweep(trainData[,isB], 1, offsets, "+")
trainData <- t(scale(t(trainData)))
dimnames(trainData) <- list(paste("gene", 1:nFeatures, sep=''),
                            paste("trainsamp", 1:nTrain, sep=''))



###################################################
### code chunk number 12: Modeler.Rnw:133-141
###################################################
K <- 3
truncData <- trainData[1:500,]
truncData[truncData < -K] <- -K
truncData[truncData > K] <- K
heatmap(truncData, col=redgreen(64),
        ColSideColors=as.character(trainClass), 
        scale='none', zlim=c(-K, K), main="Training Data")



###################################################
### code chunk number 13: test.out
###################################################
pOut <- rbeta(nTest, alpha, beta)
testOutcome <- lp(pOut)


###################################################
### code chunk number 14: test.class
###################################################
testClass <- factor(c("cyan", "magenta")[1 + 1*(pOut > 0.5)])
summary(testClass)
isB <- testClass=="magenta"
summary(isB)


###################################################
### code chunk number 15: testData
###################################################
testData <- matrix(rnorm(nFeatures*nTest), ncol=nTest) # pure noise
testData[,isB] <- sweep(testData[,isB], 1, offsets, "+")
testData <- t(scale(t(testData)))
dimnames(testData) <- list(paste("gene", 1:nFeatures, sep=''),
                            paste("testsamp", 1:nTest, sep=''))


###################################################
### code chunk number 16: Modeler.Rnw:167-174
###################################################
K <- 3
truncData <- testData[1:500,]
truncData[truncData < -K] <- -K
truncData[truncData > K] <- K
heatmap(truncData, col=redgreen(64),
        ColSideColors=as.character(testClass), 
        scale='none', zlim=c(-K, K), main="Test Data")


###################################################
### code chunk number 17: clean
###################################################
rm(list=paramlist)
rm(pOut, isB, signed, offsets)
rm(xx, yy, alpha, beta)
rm(paramlist)


###################################################
### code chunk number 18: t.tests
###################################################
library(ClassComparison)
mtt <- MultiTtest(trainData, trainClass)


###################################################
### code chunk number 19: Modeler.Rnw:194-198
###################################################
opar <- par(mfrow=c(2,1))
plot(mtt)
hist(mtt)
par(opar)


###################################################
### code chunk number 20: bum
###################################################
bum <- Bum(mtt@p.values)
countSignificant(bum, alpha=0.01, by="FDR")
countSignificant(bum, alpha=0.05, by="FDR")



###################################################
### code chunk number 21: Modeler.Rnw:209-210
###################################################
hist(bum)


###################################################
### code chunk number 22: geneset
###################################################
geneset <- rownames(trainData)[selectSignificant(bum, alpha=0.05, by="FDR")]
length(geneset)
trainSubset <- trainData[geneset,]
testSubset  <- testData[geneset,]


###################################################
### code chunk number 23: 3nn
###################################################
knnFitted <- learn(modeler3NN, trainSubset, trainClass)
knnPredictions <- predict(knnFitted, testSubset)
table(knnPredictions, testClass)



###################################################
### code chunk number 24: 5nn
###################################################
knnFitted <- learn(modeler5NN, trainSubset, trainClass)
knnPredictions <- predict(knnFitted, testSubset)
table(knnPredictions, testClass)


###################################################
### code chunk number 25: cart
###################################################
rpartFitted <- learn(modelerRPART, trainSubset, trainClass)
rpartPredictions <- predict(rpartFitted, testSubset, type='class')
table(rpartPredictions, testClass)


###################################################
### code chunk number 26: cart.reg
###################################################
rpartFitted <- learn(modelerRPART, trainSubset, trainOutcome)
rpartPredictions <- predict(rpartFitted, testSubset)
table(rpartPredictions > 0, testClass)
cor(rpartPredictions, testOutcome)
temp <- lm(testOutcome ~ rpartPredictions)


###################################################
### code chunk number 27: Modeler.Rnw:256-258
###################################################
plot(rpartPredictions, testOutcome)
abline(coef(temp))


###################################################
### code chunk number 28: lr (eval = FALSE)
###################################################
## # takes too long for the vignette, because of the "step"
## # across glm fits.
## lrFitted <- learn(modelerLR, trainSubset, trainClass)
## lrPredictions <- predict(lrFitted, testSubset)
## table(lrPredictions, testClass)


###################################################
### code chunk number 29: lr.reg
###################################################
lrFitted <- learn(modelerLR, trainSubset, trainOutcome)
lrPredictions <- predict(lrFitted, testSubset)
table(lrPredictions > 0, testClass)
cor(lrPredictions, testOutcome)
temp <- lm(testOutcome ~ lrPredictions)


###################################################
### code chunk number 30: Modeler.Rnw:281-283
###################################################
plot(lrPredictions, testOutcome)
abline(coef(temp))


###################################################
### code chunk number 31: ccp
###################################################
ccpFitted <- learn(modelerCCP, trainSubset, trainClass)
ccpPredictions <- predict(ccpFitted, testSubset)
table(ccpPredictions, testClass)


###################################################
### code chunk number 32: svm
###################################################
# takes too long for the vignette, because of the "step"
# across glm fits.
svmFitted <- learn(modelerSVM, trainSubset, trainClass)
svmPredictions <- predict(svmFitted, testSubset)
table(svmPredictions, testClass)


###################################################
### code chunk number 33: svm.reg
###################################################
svmFitted <- learn(modelerSVM, trainSubset, trainOutcome)
svmPredictions <- predict(svmFitted, testSubset)
table(svmPredictions > 0, testClass)
cor(svmPredictions, testOutcome)
temp <- lm(testOutcome ~ svmPredictions)


###################################################
### code chunk number 34: Modeler.Rnw:314-316
###################################################
plot(svmPredictions, testOutcome)
abline(coef(temp))


###################################################
### code chunk number 35: nnet
###################################################
nnetFitted <- learn(modelerNNET, trainSubset, trainClass)
nnetPredictions <- predict(nnetFitted, testSubset)
table(nnetPredictions, testClass)


###################################################
### code chunk number 36: nnet.reg
###################################################
nnetFitted <- learn(modelerNNET, trainSubset, trainOutcome)
nnetPredictions <- predict(nnetFitted, testSubset)
table(nnetPredictions > 0, testClass)
cor(nnetPredictions, testOutcome)
temp <- lm(testOutcome ~ nnetPredictions)


###################################################
### code chunk number 37: Modeler.Rnw:337-339
###################################################
plot(nnetPredictions, testOutcome)
#abline(coef(temp))


###################################################
### code chunk number 38: rf
###################################################
rfFitted <- learn(modelerRF, trainSubset, trainClass)
rfPredictions <- predict(rfFitted, testSubset)
table(rfPredictions, testClass)


###################################################
### code chunk number 39: rf.reg
###################################################
rfFitted <- learn(modelerRF, trainSubset, trainOutcome)
rfPredictions <- predict(rfFitted, testSubset)
table(rfPredictions > 0, testClass)
cor(rfPredictions, testOutcome)
temp <- lm(testOutcome ~ rfPredictions)


###################################################
### code chunk number 40: Modeler.Rnw:360-362
###################################################
plot(rfPredictions, testOutcome)
abline(coef(temp))


