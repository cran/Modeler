# Copyright (C) Kevin R. Coombes, 2007-2013

# Support Vector Machines

learnSVM <- function(data, status, params, pfun) {
  tdata <- data.frame(Stat=status, t(data))
  arglist <- c(list(formula = Stat ~ ., data=tdata),  params)
  model <- do.call(svm, arglist)
  FittedModel(pfun, data, status,
              details=list(model=model))
}

predictSVM <- function(newdata, details, status, ...) {
  predict(details$model, t(newdata), ...)
}

modelerSVM <- Modeler(learnSVM, predictSVM)
