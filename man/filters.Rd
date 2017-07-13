\name{feature.filters}
\alias{filterMean}
\alias{filterSD}
\alias{filterRange}
\alias{filterIQR}
\alias{filterMin}
\alias{filterMax}
\alias{filterMedian}

\title{Feature Filtering}
\description{
  Functions to create functions that filter potential predictive
  features using statistics that do not access class labels.
}
\usage{
filterMean(cutoff)
filterMedian(cutoff)
filterSD(cutoff)
filterMin(cutoff)
filterMax(cutoff)
filterRange(cutoff)
filterIQR(cutoff)
}
\arguments{
  \item{cutoff}{
    A real number, the level above which features with this statistic
    should be retained and below which should be discarded.
  }
}
\details{
  Following the usual conventions introduced from the world of
  gene expression microarrays, a typical data matrix is constructed from
  columns representing samples on which we want to make predictions
  amd rows representing the features used to construct the predictive
  model. In this context, we define a \emph{filter} to be a function
  that accepts a data matrix as its only argument and returns a logical
  vector, whose length equals the number of rows in the matrix, where
  'TRUE' indicates features that should be retrained. Most filtering
  functions belong to parametrized families, with one of the most common
  examples being
  "retain all features whose mean is above some pre-specified cutoff".
  We implement this idea using a set of function-generating functions,
  whose arguments are the parameters that pick out the desired member
  of the family.  The return value is an instantiation of a particular
  filtering function.  The decison to define things this way is to be
  able to apply the methods in cross-validation (or other) loops where
  we want to ensure that we use the same filtering rule each time.
}
\value{
  Each of the seven functions described here return a filter function,
  \code{f}, that can be used by code that basically looks like
  \code{logicalVector <- filter(data)}.
}
\author{
  Kevin R. Coombes <krc@silicovore.com>
}
\seealso{
  See \code{\link{Modeler-class}} and \code{\link{Modeler}} for details
  about how to train and test models.
}
\examples{
set.seed(246391)
data <- matrix(rnorm(1000*30), nrow=1000, ncol=30)
fm <- filterMean(1)
summary(fm(data))

summary(filterMedian(1)(data))
summary(filterSD(1)(data))
}
\keyword{ multivariate }
