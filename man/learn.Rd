\name{learn}
\alias{learn}
\title{
Learning models from data
}
\description{
  The \code{learn} function provides an abstraction that can be used to
  fit a binary classification model to a training data set.
}
\usage{
learn(model, data, status, prune=keepAll)
}
\arguments{
  \item{model}{
    An object of the \code{\link{Modeler-class}}
}
  \item{data}{
    A matrix containing the training data, with rows as features and
    columns as samples to be classified.
}
  \item{status}{
    A factor, with two levels, containing the known classification of
    the training data.
  }
  \item{prune}{A "pruning" funciton; that is, a funciton that takes two
    arguments (a data matrix and a class factor) and returns a logical
    vector indicating which features to retain.
  }
}
\details{
  Objects of the \code{\link{Modeler-class}} contain functions to learn
  models from training data to make predictions on new test data.  These
  functions have to be prepared as pairs, since they have a shared
  opinion about how to record and use specific details about the
  parameters of the model.  As a result, the learn function is
  implemented by:
\preformatted{
  learn <- function(model, data, status) {
    model@learn(data, status, model@params, model@predict)
  }
}
}
\value{
  An object of the \code{\link{FittedModel-class}}.
}
\author{
  Kevin R. Coombes <krc@silicovore.com>
}
\seealso{
  See \code{\link[stats]{predict}} for how to make predictions on new test data
  from an object of the \code{\link{FittedModel-class}}.
}
\examples{
# set up a generic RPART model
rpart.mod <- Modeler(learnRPART, predictRPART, minsplit=2, minbucket=1)

# simulate fake data
data <- matrix(rnorm(100*20), ncol=20)
status <- factor(rep(c("A", "B"), each=10))

# learn the specific RPART model
fm <- learn(rpart.mod, data, status)

# show the predicted results from the model on the trianing data
predict(fm)

# set up a nearest neighbor model
knn.mod <- Modeler(learnKNN, predictKNN, k=3)

# fit the 3NN model on the same data
fm3 <- learn(knn.mod, data, status)
# show its performance
predict(fm3)
}
\keyword{ classif }
\keyword{ multivariate }
