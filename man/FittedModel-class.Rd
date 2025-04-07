\name{FittedModel-class}
\Rdversion{1.1}
\docType{class}
\alias{FittedModel-class}
\alias{predict,FittedModel-method}

\title{Class "FittedModel" }
\description{	Objects of this class represent parametrized statistical
  models (of the \code{\link{Modeler-class}}) after they have been fit
  to a training data set.  These objects can be used to
  \code{\link[stats]{predict}} binary outcomes on new test data sets.}
\section{Objects from the Class}{
Objects can be created by calls to the constructor function,
\code{\link{FittedModel}}.  In practice, however, most
\code{FittedModel} objects are created as the result of applying the
\code{\link{learn}} function to an object of the
\code{\link{Modeler-class}}. 
}
\section{Slots}{
  \describe{
    \item{\code{predictFunction}:}{Object of class \code{"function"} that
      implemnts the ability to make predictions using the fitted model. }
    \item{\code{trainData}:}{Object of class \code{"matrix"} containing
      the trainng data set.  Rowes are features and columns are samples. }
    \item{\code{trainStatus}:}{Object of class \code{"vector"}. Should
      either be a numeric vector representing outcome or a factor with two
      levels, containing the classes of the training data set.}
    \item{\code{details}:}{Object of class \code{"list"} containing the
      fitted parameters for the specific model. }
    \item{\code{extras}:}{Object of class \code{"list"} containing any
      extra information (such as diagnostics) produced a a result of
      learning the model from the training data set.}
    \item{\code{fsVector}:}{Logical vector indicating which features
      should be retained (TRUE) of discared (FALSE) after performing
      featgure selection on the training data.}
  }
}
\section{Methods}{
  \describe{
    \item{predict}{\code{signature(object = "FittedModel")}: Predict the
    binary outcome on a new data set.}
	 }
}
\author{
  Kevin R. Coombes <krcoombes@mdanderson.org>
}
\seealso{
  See \code{\link{Modeler-class}} and \code{\link{learn}} for details on
  how to fit a model to data.
}
\examples{
showClass("FittedModel")
}
\keyword{classes}
