# Copyright (C) Kevin R. Coombes, 2016

learnNNET2 <- function(data, status, params, pfun) {
  f <- as.formula( paste("Stat ~ ",
                         paste(rownames(data), collapse="+")))
  status <- as.numeric(status)-1
  tdata <- data.frame(Stat=status, t(data))
  arglist <- c(list(formula = f, data=tdata),  params)
  model <- do.call(neuralnet, arglist)
  FittedModel(pfun, data, status,
              details=list(model=model))
}
predictNNET2 <- function(newdata, details, status, ...) {
  compute(details$model, t(newdata), ...)$net.result[,1]
}

modelerNNET2 <- Modeler(learnNNET2, predictNNET2, hidden=c(8,5))
