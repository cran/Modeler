# Copyright (C) Kevin R. Coombes, 2007-2013

# LR, linear or logistic regression
learnLR <- function(data, status, params, pfun) {
  if (is.null(params$prior)) {
    params$prior <- 0.5
  }
  tdata <- data.frame(Stat=status, t(data))
  # switch on the class of the status vector
  if (is.numeric(status)) {       # linear regression
#    cat("learning from numeric\n", file=stderr())
    model <- step(lm(Stat ~ ., data=tdata), trace=0)
  } else if (is.factor(status)) { # logistic regression
#    cat("learning from factor\n", file=stderr())
    model <- step(glm(Stat ~ ., data=tdata, family=binomial), trace=0)
  } else {
    warning("'status' must be a factor or a numeric vector, not a", class(status))
    model <- NULL
  }
#  cat("passing", class(status), "to FittedModel\n", file=stderr())
  FittedModel(pfun, data, status, details=list(model=model, prior=params$prior))
}

predictLR <- function(newdata, details, status, type="response", ...) {
#  print(summary(status))
  preds <- predict(details$model, newdata=data.frame(t(newdata)), type=type)
  if (is.factor(status)) {
#    cat("binarizing\n", file=stderr())
    values <- rep(levels(status)[2], length(preds))
    values[preds < details$prior] <- levels(status)[1]
    preds <- values
  }
  preds
}

modelerLR <- Modeler(learnLR, predictLR)

# PCALR, do principal components analysis first,
# followed by linear or logistic regression
#
# currently only works for binary classification, because
# of the way it uses t-tests to select features.
learnPCALR <- function(data, status, params, pfun) {
  # manually set default values for the parameters
  if (is.null(params$alpha)) {
    params$alpha <- 0.10
  }
  if (is.null(params$minNGenes)) {
    params$minNGenes <- 10
  }
  if (is.null(params$perVar)) {
    params$perVar <- 0.8
  }
  if (is.null(params$verbose)) {
    params$verbose=FALSE
  }
  if (is.null(params$prior)) {
    params$prior <- 0.5
  }
  # perform two-sample t-tests w.r.t status
  mtt1 <- MultiTtest(data, status)
  bum1 <- Bum(mtt1@p.values)
  # select features that have small p-values
  sel <- selectSignificant(bum1, alpha=params$alpha, by="FDR")
  if(sum(sel) < params$minNGenes) {
    tgt <- sort(mtt1@p.values)[1+params$minNGenes]
    sel <- mtt1@p.values < tgt
  }
  # remember how many genes were used
  nGenesSelected <- sum(sel)
  # compute principal components on training set with selected features
  spca <- SamplePCA(data[sel,])
  # decide how many PCs to use
  NC <- sum(cumsum(spca@variances)/sum(spca@variances) < params$perVar)
  if(NC < 2) NC <- 2
  if(params$verbose) cat(paste("n comps:", NC, "\n"))
  # remember how many components were used
  nCompAvail <- NC
  # assemble training status and PC features into a data frame
  trdata <- spca@scores[, 1:NC]
  # rely on the existing LR code to fit the regression model
  fmBase <- learnLR(t(trdata), status, params, predictLR)
  mmod <- fmBase@details$model
  nCompUsed <- length(mmod$coefficients)-1 #$
  FittedModel(pfun, data, status,
              details=list(prior=params$prior,
                sel=sel,
                spca=spca,
                baseModel=fmBase,
                nCompAvail=nCompAvail))
}

predictPCALR <- function(newdata, details, status, ...) {
  # project the test set into the principal component space so we know
  # the values of the predictors in the test set
  proj <- predict(details$spca, newdata=newdata[details$sel,])
  temp <- data.frame(proj)[, 1:details$nCompAvail]
  base <- details$baseModel
  base@predictFunction(t(temp), base@details, base@trainStatus, ...)
} 

modelerPCALR <- Modeler(learnPCALR, predictPCALR)

learnSelectedLR <- function(data, status, params, pfun) {
  if (is.null(params$alpha)) {
    params$alpha <- 0.10
  }
  if (is.null(params$minNGenes)) {
    params$minNGenes <- 10
  }
  if (is.null(params$perVar)) {
    params$perVar <- 0.8
  }
  if (is.null(params$verbose)) {
    params$verbose=TRUE
  }
  if (is.null(params$prior)) {
    params$prior <- 0.5
  }
  # perform two-sample t-tests w.r.t status
  mtt1 <- MultiTtest(data, status)
  bum1 <- Bum(mtt1@p.values)
  # select features that have small p-values
  sel <- selectSignificant(bum1, alpha=params$alpha, by="FDR")
  if(sum(sel) < params$minNGenes) {
    tgt <- sort(mtt1@p.values)[1+params$minNGenes]
    sel <- mtt1@p.values < tgt
  }
  # remember how many genes were used
  nGenesSelected <- sum(sel)
  trdata <- data.frame(Stat=status, t(data[sel,]))
  # fit a logistic prediction model
  # use step-wise AIC to get the "optimal" model
  mmod <- step(glm(Stat ~ ., data=trdata, family=binomial), trace=0)
  nFeaturesUsed <- length(mmod$coefficients)-1 #$
  fm <- FittedModel(pfun, data, status,
                    details=list(prior=params$prior, sel=sel,
                      mmod=mmod),
                    nGenesSelected=nGenesSelected,
                    nFeaturesUsed=nFeaturesUsed)
  fm
}

predictSelectedLR <- function(newdata, details, status, ...) {
  # project the test set into the principal component space so we know
  # the values of the predictors in the test set
  temp <- data.frame(t(newdata[details$sel,]))
  preds <- predict(details$mmod, newdata=temp, type='response')
  values <- rep(levels(status)[2], length(preds))
  values[preds < details$prior] <- levels(status)[1]
  factor(values)
} 

modelerSelectedLR <- Modeler(learnSelectedLR, predictSelectedLR)
