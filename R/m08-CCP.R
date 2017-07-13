# Copyright (C) Kevin R. Coombes, 2007-2013

learnCCP <- function(data, status, params, pfun) {
  norm.t <- MultiTtest(data, status) #$
  ccp <- matrix(norm.t@t.statistics, nrow=1) %*% data
  temp <- sapply(levels(status), function(f) {mean(ccp[status==f])})
  cutpoint <- mean(temp)
  FittedModel(pfun, data, status,
              details=list(norm.t=norm.t, ccp=ccp,
                cutpoint=cutpoint, big=which(temp==max(temp))))
}

predictCCP <- function(newdata, details, status, ...) {
  new.ccp <- matrix(details$norm.t@t.statistics, nrow=1) %*% newdata
  big <- details$big # must be 1 or 2
  # implies 3-big = 2 or 1
  pred <- rep(levels(status)[3-big], length(new.ccp))
  pred[new.ccp > details$cutpoint] <- levels(status)[big]
  factor(pred)
}

modelerCCP <- Modeler(learnCCP, predictCCP)
