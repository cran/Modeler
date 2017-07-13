\name{learnRF}
\alias{learnRF}
\alias{predictRF}
\alias{modelerRF}
\title{
  Fit models and make predictions with a PCA-LR classifier
}
\description{
  These functions are used to apply the generic train-and-test
  mechanism to a classifier that combines principal component analysis
  (PCA) with logistic regression (LR). 
}
\usage{
learnRF(data, status, params, pfun)
predictRF(newdata, details, status, ...)
}
\arguments{
  \item{data}{
    The data matrix, with rows as features ("genes") and columns as the
    samples to be classified.
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
    trained classifier.  Should always be set to
    \code{predictRF}.
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
  The input arguments to both \code{learnRF} and \code{predictRF}
  are dictated by the requirements of the general train-and-test
  mechanism provided by the \code{\link{Modeler-class}}.
  
  The RF classifier is similar in spirit to the "supervised principal
  components" method implemented in the \code{superpc} package.  We
  start by performing univariate two-sample t-tests to identify features
  that are differentially expressed between two groups of training
  samples.  We then set a cutoff to select features using a bound
  (\code{alpha}) on the false discovery rate (FDR).  If the number of
  selected features is smaller than a prespecified goal
  (\code{minNgenes}), then we increase the FDR until we get the desired
  number of features.  Next, we perform PCA on the selected features
  from the trqining data.  we retain enough principal components (PCs)
  to explain a prespecified fraction of the variance (\code{perVar}).
  We then fit a logistic regression model using these PCs to predict the
  binary class of the training data.  In order to use this model to make
  binary predictions, you must specify a \code{prior} probability that a
  sample belongs to the first of the two groups (where the ordering is
  determined by the levels of the classification factor, \code{status}).

  In order to fit the model to data, the \code{params} argument to the
  \code{learnRF} function should be a list containing components
  named \code{alpha}, \code{minNgenes}, \code{perVar}, and \code{prior}.
  It may also contain a logical value called \code{verbose}, which
  controls the amount of information that is output as the algorithm runs.

  The result of fitting the model using \code{learnRF} is a member of
  the \code{\link{FittedModel-class}}.  In additon to storing the
  prediction function (\code{pfun}) and the training data and status,
  the FittedModel stores those details about the model that are required
  in order to make predictions of the outcome on new data.  In this
  acse, the details are: the \code{prior} probability, the set of 
  selected features (\code{sel}, a logical vector), the principal
  component decomposition (\code{spca}, an object of the
  \code{\link[ClassDiscovery]{SamplePCA}} class), the logistic
  regression model (\code{mmod}, of class \code{\link{glm}}), the number
  of PCs used (\code{nCompUsed}) as well as the number of components
  available (\code{nCompAvail}) and the number of gene-features selected
  (\code{nGenesSelecets}).  The \code{details} object is appropriate for
  sending as the second argument to the \code{predictRF} function in
  order to make predictions with the model on new data.  Note that the
  status vector here is the one used for the \emph{training} data, since
  the prediction function only uses the \emph{levels} of this factor to
  make sure that the direction of the predicitons is interpreted
  correctly.
}
\value{
  The \code{learnRF} function returns an object of the
  \code{\link{FittedModel-class}}, representing a RF classifier
  that has been fitted on a training \code{data} set.

  The \code{predictRF} function returns a factor containing the
  predictions of the model when applied to the new data set.
}
\author{
  Kevin R. Coombes <krc@silicovore.com>
}
\seealso{
  See \code{\link{Modeler-class}} and \code{\link{Modeler}} for details
  about how to train and test models.  See
  \code{\link{FittedModel-class}} and \code{\link{FittedModel}} for
  details about the structure of the object returned by \code{learnRF}.
}
\examples{
# simulate some data
data <- matrix(rnorm(100*20), ncol=20)
status <- factor(rep(c("A", "B"), each=10))

# set up the parameter list
svm.params <- list(minNgenes=10, alpha=0.10, perVar=0.80, prior=0.5)

# learn the model
#fm <- learnRF(data, status, svm.params, predictRF)

# Make predictions on some new simulated data
#newdata <- matrix(rnorm(100*30), ncol=30)
#predictRF(newdata, fm@details, status)
}
\keyword{ classif }
\keyword{ multivariate }
