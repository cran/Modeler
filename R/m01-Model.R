# Copyright (C) Kevin R. Coombes, 2007-2016

###
### m01-model.R
###

setClass("Modeler",
         slots = c(learnFunction="function",
                   predictFunction="function",
                   paramList="list"))

## Generates a Modeler object
Modeler <- function(learn, predict, ...) {
    new("Modeler",
        learnFunction=learn,
        predictFunction=predict,
        paramList=list(...))
}

learn <- function(model, data, status, prune=keepAll) {
  keep <- prune(data, status)
  data <- data[keep, , drop = FALSE]
  fitted <- model@learnFunction(data,
                                status,
                                model@paramList,
                                model@predictFunction)
  fitted@fsVector <- keep
  fitted
}

