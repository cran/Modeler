\name{Modeler-package}
\alias{Modeler-package}
\docType{package}
\title{Model for classification and regression}
\description{
  The \code{Modeler} package provides generic tools for learning models
  and making predictions in the context of high-throughput data sets
  such as those produced by gene expression microarrays. In order to use
  this package , you must first prepare a pair of functions (one for
  learning models from training data, and one for making predictions on
  test data). These functions, along with any required meta-parameters,
  are used to create an object of the \code{\link{Modeler-class}}. That
  object is then passed to the \code{\link{Modeler}} function along with
  the full training data set.
}
\details{
  \tabular{ll}{
    Package: \tab Modeler\cr
    Type: \tab Package\cr
    Version: \tab 2.0.0\cr
    Date: \tab 2013-07-05\cr
    License: \tab Artistic-2.0\cr
    LazyLoad: \tab yes\cr
  }
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com}
}
\keyword{package}
\keyword{classif}
\keyword{multivariate}

