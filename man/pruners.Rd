\name{feature.selection}
\alias{keepAll}
\alias{fsTtest}
\alias{fsModifiedFisher}
\alias{fsPearson}
\alias{fsSpearman}
\alias{fsMedSplitOddsRatio}
\alias{fsChisquared}
\alias{fsEntropy}
\alias{fsFisherRandomForest}
\alias{fsTailRank}

\title{Feature Selection}
\description{
  Functions to create functions that perform feature selection (or at
  least feature reduction)  using statistics that access class labels.
}
\usage{
keepAll(data, group)
fsTtest(fdr, ming=500)
fsModifiedFisher(q)
fsPearson(q = NULL, rho)
fsSpearman(q = NULL, rho)
fsMedSplitOddsRatio(q = NULL, OR)
fsChisquared(q = NULL, cutoff)
fsEntropy(q = 0.9, kind=c("information.gain", "gain.ratio", "symmetric.uncertainty"))
fsFisherRandomForest(q)
fsTailRank(specificity=0.9, tolerance=0.5, confidence=0.5)
}
\arguments{
  \item{data}{
    A matrix containng the data; columns are samples and rows are features.
  }
  \item{group}{
    A factor with two levels defining the sample classes.
  }
  \item{fdr}{
    A real number between 0 and 1 specifying the target false discovery rate (FDR).
  }
  \item{ming}{
    An integer specifing the minimum number of features to return;
    overrides the FDR.
  }
  \item{q}{
    A real number between 0.5 and 1 specifiying the fraction of features to discard.
  }
  \item{rho}{
    A real number between 0 and 1 specifying the absolute value of the
    correlation coefficient used to filter features.
  }
  \item{OR}{
    A real number specifying the desired odds ratio for filtering features.
  }
  \item{cutoff}{
    A real number specifiyng the targeted cutoff rate when using the
    statistic to filter features.
  }
  \item{kind}{
    The kind of information metric to use for filtering features.
  }
  \item{specificity}{
    See \code{\link[TailRank]{TailRankTest}}.
  }
  \item{tolerance}{
    See \code{\link[TailRank]{TailRankTest}}.
  }
  \item{confidence}{
    See \code{\link[TailRank]{TailRankTest}}.
  }
}
\details{
  Following the usual conventions introduced from the world of
  gene expression microarrays, a typical data matrix is constructed from
  columns representing samples on which we want to make predictions
  amd rows representing the features used to construct the predictive
  model. In this context, we define a \emph{feature selector} or
  \emph{pruner} to be a function that accepts a data matrix and a
  two-level factor as its only arguments and returns a logical 
  vector, whose length equals the number of rows in the matrix, where
  'TRUE' indicates features that should be retrained. Most pruning
  functions belong to parametrized families.  We implement this idea
  using a set of function-generating functions, whose arguments are the
  parameters that pick out the desired member of the family.  The return
  value is an instantiation of a particular filtering function.  The
  decison to define things this way is to be able to apply the methods
  in cross-validaiton (or other) loops where we want to ensure that we
  use the same feature selection rule each time.

  We have implemented the following algorithms:
  \itemize{
    \item \code{keepAll}: retain all features; do nothing.
    \item \code{fsTtest}: Keep features based on the false discovery rate
    from a two-goup t-test, but always retain a specified minimum number
    of genes.
    \item \code{fsModifiedFisher} Retain the top quantile of features
    for the statistic \deqn{\frac{(m_A - m)^2 + (m_B - m)^2}{v_A + v_B}}
    where m is the mean and v is the variance.
    \item \code{fsPearson}: Retain the top quantile of features based on
    the absolute value of the Pearson correlation with the binary outcome.
    \item \code{fsSpearman}: Retain the top quantile of features based on
    the absolute value of the Spearman correlation with the binary outcome.
    \item \code{fsMedSplitOddsRatio}: Retain the top quantile of
    features based on the odds ratio to predict the binary outcome,
    after first dichotomizing the continuous predictor using a split at
    the median value. 
    \item \code{fsChisquared}: retain the top quantile of features based
    on a chi-squared test comparing the binary outcome to continous
    predictors discretized into ten bins.
    \item \code{fsEntropy}: retain the top quantile of features based on
    one of three information-theoretic measures of entropy. 
    \item \code{fsFisherRandomForest}: retain the top features based on
    their importance in a random forest analysis, after first filtering
    using the modified Fisher statistic. 
    \item \code{fsTailRank}: Retain features that are significant based
    on the TailRank test, which is a measure of whether the tails of the
    distributions are different.
  }
}
\value{
  The \code{keepAll} function is a "pruner"; it takes the data matrix and
  grouping factor as arguments, and returns a logical vector indicating
  which features to retain.
  
  Each of the other nine functions described here return uses its
  arguments to contruct and return a pruning function,
  \code{f}, that has the same interface as \code{keepAll}.
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
data <- matrix(rnorm(1000*36), nrow=1000, ncol=36)
data[1:50, 1:18] <- data[1:50, 1:18] + 1
status <- factor(rep(c("A", "B"), each=18))

fsel <- fsPearson(q = 0.9)
summary(fsel(data, status))
fsel <- fsPearson(rho=0.3)
summary(fsel(data, status))

fsel <- fsEntropy(kind="gain.ratio")
summary(fsel(data, status))

}
\keyword{ multivariate }

