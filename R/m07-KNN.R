# Copyright (C) Kevin R. Coombes, 2007-2013

learnKNN <- function(data, status, params, pfun) {
  K <- ifelse(is.null(params$K), 3, params$K)
  fm <- FittedModel(pfun, data, status,
                    details=list(train=data,
                      status=status,
                      K=K))
}

predictKNN <- function(newdata, details, status, ...) {
  data <- details$train
  K <- details$K
  status <- details$status
  knn(t(data), t(newdata), status, k=K)
}

modeler3NN <- Modeler(learnKNN, predictKNN, params=list(K=3))
modeler5NN <- Modeler(learnKNN, predictKNN, params=list(K=5))
