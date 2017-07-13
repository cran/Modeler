# Copyright (C) Kevin R. Coombes, 2007-2013

learnRPART <- function(data, status, params, pfun) {
  tda <- data.frame(status, t(data))
  rp <- rpart(status ~ ., data=tda, control=params)
  fm <- FittedModel(pfun, data, status,
                    details=list(rp=rp))
}

predictRPART <- function(newdata, details, status, ...) {
  nd <- data.frame(t(newdata))
  predict(details$rp, newdata=nd, ...)
}

modelerRPART <- Modeler(learnRPART, predictRPART)
