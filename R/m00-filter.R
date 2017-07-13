# Copyright (C) Kevin R. Coombes, 2007-2016

###
### m00-filter.R
###

# Filter interface
# INPUT:
#	(training) data
#	optional parameter list
# OUTPUT:
#	logical vector equal in length to number of data rows

filterMean <- function(cutoff) {
  function(data) {
    m <- as.vector(oompaBase::matrixMean(data))
    m > cutoff
  }
}

filterSD <- function(cutoff) {
  function(data) {
    m <- oompaBase::matrixMean(data)
    s <- as.vector(sqrt(oompaBase::matrixVar(data, m)))
    s > cutoff
  }
}

filterRange <- function(cutoff) {
  function(data) {
    r <- apply(data, 1, function(x) diff(range(x)))
    r > cutoff
  }
}

filterIQR <- function(cutoff) {
  function(data) {
    r <- apply(data, 1, function(x) diff(quantile(x,c(0.25, 0.75))))
    r > cutoff
  }
}

filterMin <- function(cutoff) {
  function(data) {
    r <- apply(data, 1, min, na.rm=TRUE)
    r > cutoff
  }
}

filterMax <- function(cutoff) {
  function(data) {
    r <- apply(data, 1, max, na.rm=TRUE)
    r > cutoff
  }
}

filterMedian <- function(cutoff) {
  function(data) {
    r <- apply(data, 1, median, na.rm=TRUE)
    r > cutoff
  }
}
