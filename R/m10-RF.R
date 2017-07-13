# Copyright (C) Kevin R. Coombes, 2007-2013

# random forest

learnRF <- function(data, status, params, pfun) {
  tdata <- data.frame(Stat=status, t(data))
  model <- randomForest(Stat ~ ., data=tdata)
  FittedModel(pfun, data, status,
              details=list(model=model))
}
predictRF <- function(newdata, details, status, ...) {
  predict(details$model, t(newdata), ...)
}

modelerRF <- Modeler(learnRF, predictRF)
