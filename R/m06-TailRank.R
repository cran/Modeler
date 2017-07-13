# Copyright (C) Kevin R. Coombes, 2007-2013

learnTailRank <- function(data, status, params, pfun) {
  if (is.null(params$spec)) {
    params$spec <- 0.8
  }
  if (is.null(params$conf)) {
    params$conf <- 0.9
  }  
  tr <- TailRankTest(data, status, specificity=params$spec,
                     confidence=params$conf, direction="two")
  tdata <- data.frame(Stat=status, t(data))
  model <- step(glm(Stat ~ ., data=tdata, family=binomial), trace=0)
  FittedModel(pfun, data, status, details=list(model=model, prior=params$prior))
}

predictTailRank <- function(newdata, details, status, ...) {
  preds <- predict(details$model, newdata=data.frame(t(newdata)), type='response')
  values <- rep(levels(status)[2], length(preds))
  values[preds < details$prior] <- levels(status)[1]
  values
}

modelerTailRank <- Modeler(learnTailRank, predictTailRank)

