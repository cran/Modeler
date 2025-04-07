\name{Modeler-class}
\Rdversion{1.1}
\docType{class}
\alias{Modeler-class}

\title{Class "Modeler" }
\description{The \code{Modeler} class represents (parametrized but not
  yet fit) statistical models that can predict binary outcomes.}
\section{Objects from the Class}{
Objects can be created by calls to the constructor fuinction, \code{\link{Modeler}}.
}
\section{Slots}{
  \describe{
    \item{\code{learnFunction}:}{Object of class \code{"function"} that is used
      to fit the model to a data set.  See \code{\link{learn}} for details.}
    \item{\code{predictFunction}:}{Object of class \code{"function"} that is
      used to make predictions on new data from a fitted model.  See
      \code{\link{predict,FittedModel-method}} for details. }
    \item{\code{paramList}:}{Object of class \code{"list"} that contains
      parameters that are specific for one type of classifier.}
  }
}
\section{Methods}{
No methods are defined with class "Modeler" in the signature.  The only
function that can be applied to a \code{Modeler} object is
\code{\link{learn}}, which has not been made into a generic funtion.
}
\author{
  Kevin R. Coombes <krc@silicovore.com>
}
\seealso{
  See the description of the \code{\link{FittedModel-class}} for details
  about the kinds of objects produced by \code{\link{learn}}.
}
\examples{
showClass("Modeler")
}
\keyword{classes}
