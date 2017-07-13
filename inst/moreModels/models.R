#########################################################
# LDA
library(MASS)

learnLDA <- function(data, status, params, pfun) {
  tdata <- data.frame(Stat=status, t(data))
  model <- lda(Stat ~ ., data=tdata)
  FittedModel(pfun, data, status, details=list(model=model))
}

predictLDA <- function(newdata, details, status, ...) {
  preds <- predict(details$model, newdata=data.frame(t(newdata)))
  preds$class
}

ldaModel <- Modeler(learnLDA, predictLDA)
ldaFitted <- learn(ldaModel, trainSubset, trainClass)
ldaPredictions <- predict(ldaFitted, testSubset)
table(ldaPredictions, testClass)

###################
# alternative LDA

require(MASS)
learnLDA <- function(data, status, params, pfun) {
  if (is.null(params$alpha)) {
    params$alpha <- 0.10
  }
  if (is.null(params$minNGenes)) {
    params$minNGenes <- 10
  }
  if (is.null(params$maxNGenes)) {
    params$maxNGenes <- 100
  }
  if (is.null(params$verbose)) {
    params$verbose=TRUE
  }
  if (is.null(params$prior)) {
    params$prior <- c(0.5, 0.5)
  }
  # perform two-sample t-tests w.r.t status
  mtt1 <- MultiTtest(data, status)
  bum1 <- Bum(mtt1@p.values)
  # select features that have small p-values
  sel <- selectSignificant(bum1, alpha=params$alpha, by="FDR")
  if(sum(sel) < params$minNGenes) {
    tgt <- sort(mtt1@p.values)[1+params$minNGenes]
    sel <- mtt1@p.values <= tgt
  }
  if(sum(sel) > params$maxNGenes) {
    tgt <- sort(mtt1@p.values)[1+params$maxNGenes]
    sel <- mtt1@p.values <= tgt
  }
  # remember how many genes were used
  nGenesSelected <- sum(sel)
  trdata <- data.frame(Stat=status, t(data[sel,]))
  model <- lda(Stat ~ ., trdata)
  fm <- FittedModel(pfun, data, status,
                    details=list(prior=params$prior, sel=sel,
                      model=model),
                    nGenesSelected=nGenesSelected)
  fm
}

predictLDA <- function(newdata, details, status, ...) {
  # project the test set into the principal component space so we know
  # the values of the predictors in the test set
  temp <- data.frame(t(newdata[details$sel,]))
  preds <- predict(details$model, newdata=temp, prior=details$prior)
  values <- data.frame(Class=preds$class, Posterior=preds$posterior, Coef=preds$x)
  values
} 

modelerLDA <- Modeler(learnLDA, predictLDA)


###################
# QDA
library(MASS)
learnQDA <- function(data, status, params, pfun) {
  tdata <- data.frame(Stat=status, t(data))
  model <- qda(Stat ~ ., data=tdata)
  FittedModel(pfun, data, status, details=list(model=model))
}
predictQDA <- predictLDA

qgeneset <- rownames(trainData)[mtt@p.values < sort(mtt@p.values)[10]]
dataForQDA <- trainData[qgeneset,]

qdaModel <- Modeler(learnQDA, predictQDA)
qdaFitted <- learn(qdaModel, dataForQDA, trainClass)
qdaPredictions <- predict(qdaFitted, testSubset)
table(qdaPredictions, testClass)

###################
# alternative QDA

require(MASS)
learnQDA <- function(data, status, params, pfun) {
  if (is.null(params$alpha)) {
    params$alpha <- 0.10
  }
  if (is.null(params$minNGenes)) {
    params$minNGenes <- 10
  }
  if (is.null(params$maxNGenes)) {
    params$maxNGenes <- 100
  }
  if (is.null(params$verbose)) {
    params$verbose=TRUE
  }
  if (is.null(params$prior)) {
    params$prior <- c(0.5, 0.5)
  }
  # perform two-sample t-tests w.r.t status
  mtt1 <- MultiTtest(data, status)
  bum1 <- Bum(mtt1@p.values)
  # select features that have small p-values
  sel <- selectSignificant(bum1, alpha=params$alpha, by="FDR")
  if(sum(sel) < params$minNGenes) {
    tgt <- sort(mtt1@p.values)[1+params$minNGenes]
    sel <- mtt1@p.values <= tgt
  }
  if(sum(sel) > params$maxNGenes) {
    tgt <- sort(mtt1@p.values)[1+params$maxNGenes]
    sel <- mtt1@p.values <= tgt
  }
  # remember how many genes were used
  nGenesSelected <- sum(sel)
  trdata <- data.frame(Stat=status, t(data[sel,]))
  model <- qda(Stat ~ ., trdata)
  fm <- FittedModel(pfun, data, status,
                    details=list(prior=params$prior, sel=sel,
                      model=model),
                    nGenesSelected=nGenesSelected)
  fm
}

predictQDA <- function(newdata, details, status, ...) {
  # project the test set into the principal component space so we know
  # the values of the predictors in the test set
  temp <- data.frame(t(newdata[details$sel,]))
  preds <- predict(details$model, newdata=temp, prior=details$prior)
  values <- data.frame(Class=preds$class, Posterior=preds$posterior, Coef=preds$x)
  values
} 

modelerQDA <- Modeler(learnQDA, predictQDA)


###################
# PAMR (classification only)
library(pamr)
learnPAMR <- function(data, status, params, pfun) {
  inputs <- list(y=status, x = data)
  model <- pamr.train(inputs)
  res <- data.frame(threshold=model$threshold,
                    errors=model$errors,
                    nonzero=model$nonzero)
  n <-  1 + sum(res >= 50)
  th <- res$threshold[n]
  FittedModel(pfun, data, status,
              details=list(model=model, 
                threshold=th))
}
predictPAMR <- function(newdata, details, status, ...) {
  pamr.predict(details$model,
               newdata,
               threshold=details$threshold, ...)
}

pamrModel <- Modeler(learnPAMR, predictPAMR)
pamrFitted <- learn(pamrModel, trainData, trainClass)
pamrPredictions <- predict(pamrFitted, testData, type='class')
table(pamrPredictions, testClass)

###################
# FirstSign
library(FirstSign)
fsModel <- FirstSign(trainData, geneset)

fs.selfpred <- predict(fsModel, dichot=TRUE)
table(pred=fs.selfpred, trainClass)

fs.testpred <- predict(fsModel, newdata=testData, dichot=TRUE)
table(pred=fs.testpred, testClass)


#########################################################
# RDA
library(rda)
rdaModel <- rda(trainData,
                as.numeric(trainClass))
dang <- rda.cv(rdaModel, trainData,
               as.numeric(trainClass))

# this is an entire bloody array of predictions for various tunable parameters.
rda.selfpred <- predict(rdaModel, trainData, as.numeric(trainClass),
                        xnew=trainData, type='class')
table(pred=rda.selfpred)

###################
# superpc -- FAIL COMPLETELY
library(superpc)
forSUPER <- list(x = trainData, y = trainClass=="cyan")
superpcModel <- superpc.train(forSUPER, type="regression")
cv <- superpc.cv(superpcModel, forSUPER)
data.frame(cv$thresh, cv$nonzero, t(cv$scor))
w <- which(cv$scor[3,]==max(cv$scor[3,]))
th <- cv$threshold[w]
super.pred <- superpc.predict(superpcModel, forSUPER,
                              list(x=testData),
                              threshold=th)
v <- super.pred$v.pred[,3]
