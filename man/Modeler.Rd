\name{Modeler}
\alias{Modeler}
\title{
  Constructor for "Modeler" objects
}
\description{
  The \code{\link{Modeler-class}} represents (parametrized but not yet
  fit) statistical models that can predict binary outcomes. The
  \code{Modeler} function is used to construct objects of this class.
}
\usage{
Modeler(learn, predict, ...)
}
\arguments{
  \item{learn}{
    Object of class \code{"function"} that will be used to fit the model
    to a data set.  See \code{\link{learn}} for details.
}
\item{predict}{
  Object of class \code{"function"} that will be used to make
  predictions on new data from a fitted model.  See
  \code{\link{predict}} for details.
}
  \item{\dots}{
    Additional parameters required for the specific kind of
    classificaiton model that will be constructed.  See Details. 
}
}
\details{
  Objects of the \code{\link{Modeler-class}} provide a general
  abstraction for classification models that can be learned from one
  data set and then applied to a new data set.  Each type of classifier
  is likely to have its own specific parameters.  For instance, a
  K-nearest neighbors classifier requires you to specify \code{k}.  The
  more complex classifier, PCA-LR has many more parameters, including
  the false discovery rate (\code{alpha}) used to select features and
  the percentage of variance (\code{perVar}) that should be explained by
  the number of principal components created from those features.  All
  additional parameters should be suplied as named arguments to the
  \code{Modeler} constructor; these additional parameters will be
  bundled into a list and inserted into the \code{params} slot of the
  resulting object of the \code{\link{Modeler-class}}.
}
\value{
  Returns an object of the \code{\link{Modeler-class}}.
}
\author{
  Kevin R. Coombes <krc@silicovore.com>
}
\seealso{
  See the descriptions of the \code{\link{learn}} function and
  the \code{\link{predict}} method for details on how to fit models on
  training data and make predictions on new test data. 
  
  See the description of the \code{\link{FittedModel-class}} for details
  about the kinds of objects produced by \code{\link{learn}}.
}
\examples{
learnNNET
predictNNET
modelerNNET <- Modeler(learnNNET, predictNNET, size=5)
modelerNNET
}
\keyword{ classif }
\keyword{ multivariate }
