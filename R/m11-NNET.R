# Copyright (C) Kevin R. Coombes, 2007-2013

learnNNET <- function(data, status, params, pfun) {
  tdata <- data.frame(Stat=status, t(data))
  arglist <- c(list(formula = Stat ~ ., data=tdata),  params)
  model <- do.call(nnet, arglist)
  FittedModel(pfun, data, status,
              details=list(model=model))
}
predictNNET <- function(newdata, details, status, ...) {
  predict(details$model, t(newdata), ...)
}

modelerNNET <- Modeler(learnNNET, predictNNET, size=5)
