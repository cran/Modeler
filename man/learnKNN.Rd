\name{learnKNN}
\alias{learnKNN}
\alias{predictKNN}
\alias{modeler3NN}
\alias{modeler5NN}
\title{
  Fit models and make predictions with a KNN classifier
}
\description{
  These functions are used to apply the generic train-and-test
  mechanism to a K-nearest neighbors (KNN) classifier. 
}
\usage{
learnKNN(data, status, params, pfun)
predictKNN(newdata, details, status, ...)
}
\arguments{
  \item{data}{
    The data matrix, with rows as features and columns as the samples to
    be classified.
}
  \item{status}{
    A factor, with two levels, classifying the samples. The length must
    equal the number of \code{data} columns.
}
  \item{params}{
    A list of additional parameters used by the classifier; see Details.
}
  \item{pfun}{
    The function used to make predictions on new data, using the
    trained classifier.
}
  \item{newdata}{
    Another data matrix, with the same number of rows as \code{data}.
}
  \item{details}{
    A list of additional parameters describing details about the
    particular classifier; see Details.
}
  \item{\dots}{
    Optional extra parameters required by the generic "predict" method.
}
}
\details{
  The input arguments to both \code{learnKNN} and \code{predictKNN}
  are dictated by the requirements of the general train-and-test
  mechanism provided by the \code{\link{Modeler-class}}.
  
  The implementation uses the \code{\link[class]{knn}} method from the
  \code{\link{class}} package.  The \code{params} argument to
  \code{learnKNN} must be alist that at least includes the component
  \code{k} that specifies the number of neighbors used.
}
\value{
  The \code{learnKNN} function returns an object of the
  \code{\link{FittedModel-class}}, logically representing a KNN
  classifier that has been fitted on a training \code{data} set.

  The \code{predictKNN} function returns a factor containing the
  predictions of the model when applied to the new data set.

}
\references{
  Ripley, B. D. (1996) \emph{Pattern Recognition and Neural Networks}.
  Cambridge.

  Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied
  Statistics with S}. Fourth edition.  Springer.
}
\author{
  Kevin R. Coombes <krc@silicovore.com>
}
\seealso{
  See \code{\link{Modeler-class}} and \code{\link{Modeler}} for details
  about how to train and test models.  See
  \code{\link{FittedModel-class}} and \code{\link{FittedModel}} for
  details about the structure of the object returned by \code{learnPCALR}.
}
\examples{
# simulate some data
data <- matrix(rnorm(100*20), ncol=20)
status <- factor(rep(c("A", "B"), each=10))

# set up the parameter list
knn.params <- list(k=5)

# learn the model
fm <- learnKNN(data, status, knn.params, predictKNN)

# Make predictions on some new simulated data
newdata <- matrix(rnorm(100*30), ncol=30)
predictKNN(newdata, fm@details, status)
}
\keyword{ classif }
\keyword{ multivariate }
