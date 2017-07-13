\name{FittedModel}
\alias{FittedModel}
\title{
Creating FittedModel objects
}
\description{
  Construct an object of the \code{\link{FittedModel-class}}.
}
\usage{
FittedModel(predict, data, status, details, ...)
}
\arguments{
  \item{predict}{
A function that applies the model to predict outcomes on new test data.
}
  \item{data}{
A matrix containing the training data.
}
  \item{status}{
A vector containing the training outcomes, which should either be a
binary-valued factor or a numeric vector of contiuous outcomes.
}
  \item{details}{
A list of the fitted parameters for the specified model.
}
  \item{\dots}{
Any extra information that is produced while learning the model; these
wil be saved in the \code{extras} slot of the \code{FittedModel} object.
}
}
\details{
Most users will never need to use this function; instead, they will
first use an existing object of the \code{\link{Modeler-class}}, 
call the \code{\link{learn}} method of that object with the training data
to obtain a \code{FittedModel} object, and then apply its
\code{\link{predict}} method to test data.  Only people who want to
implement the learn-predict interface for a new classification algorithm
are likely to need to call this function directly.
}
\value{
  Returns an object of the \code{\link{FittedModel-class}}.
}
\author{
  Kevin R. Coombes <krc@silicovore.com.
}
\seealso{
  See the descriptions of the \code{\link{learn}} function and
  the \code{\link{predict}} method for details on how to fit models on
  training data and make predictions on new test data. 
  
  See the description of the \code{\link{Modeler-class}} for details
  about the kinds of objects produced by \code{\link{learn}}.
}
\examples{
# see the examples for learn and predict and for specific
# implementations of classifiers.
}
\keyword{ classif }
\keyword{ multivariate }
