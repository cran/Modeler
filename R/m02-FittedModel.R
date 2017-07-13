# Copyright (C) Kevin R. Coombes, 2007-2013

###
### m02-fittedModel.R
###

setClassUnion("numericOrFactor", c("numeric", "factor"))

setClass("FittedModel",
         representation(predictFunction="function",
                        trainData="matrix",
                        trainStatus="numericOrFactor",
                        details="list",
                        extras="list",
                        fsVector = "logical"))

## Generates a FittedModel object
FittedModel <- function(predict, data, status, details, ...) {
#  cat("FittedModel received", class(status), "\n", file=stderr())
    new("FittedModel",
        predictFunction=predict,
        trainData=data,
        trainStatus=status,
        details=details,
        extras=list(...),
        fsVector=logical())
}

setMethod("predict", signature(object="FittedModel"),
          function(object,
                   newdata=object@trainData,
                   ...) {
    if(nrow(newdata) > sum(object@fsVector)) {
      newdata <- newdata[object@fsVector,]
    }
    object@predictFunction(newdata=newdata,
                           object@details,
                           object@trainStatus,
                           ...)
})

