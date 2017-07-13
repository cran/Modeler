# Copyright (C) Kevin R. Coombes, 2007-2016

###
### m00-pruner.R
###

# Pruner interface
# INPUT:
#	(training) data
#	(training) class-status
#	optional parameter list
# OUTPUT:
#	logical vector equal in length to number of data rows

# the keep-everything feature selection method
keepAll <- function(data, group) {
  rep(TRUE, length=nrow(data))
}

# pruning by a t-statistic
fsTtest <- function(fdr, ming=500) {
  function(data, group) {
    mtt <- MultiTtest(data, group)
    bum <- Bum(mtt@p.values)
    x <- selectSignificant(bum, alpha=fdr, by="FDR")
    if (sum(x) < ming) {
      p <- sort(mtt@p.values)[ming]
      x <- mtt@p.values <= p
    }
    x
  }
}

# feature selection or "pruning" by a modified Fisher criterion F-score
fsModifiedFisher <- function(q) {
  function(data, group) {
    isA <- group == levels(group)[1]
    isB <- !isA
    m  <- oompaBase::matrixMean(data)
    mA <- oompaBase::matrixMean(data[,isA])
    mB <- oompaBase::matrixMean(data[,isB])
    vA <- oompaBase::matrixVar(data[,isA], mA)
    vB <- oompaBase::matrixVar(data[,isB], mB)
    fval <- as.vector( ((mA - m)^2 + (mB - m)^2)/(vA + vB) )
    fval > quantile(fval, q)
  }
}

## hidden helper function
.corstat <- function(data, group, q, rho) {
  target <- matrix(scale(1 * (group == levels(group)[1])), ncol=1)
  scaled <- t(scale(t(data)))
  x <- abs(as.vector(scaled %*% target) / length(group))
  if (!is.null(q)) {
    rho <- quantile(x, q)
  }
  x > rho
}

# pruning by Pearson correlation
fsPearson <- function(q = NULL, rho) {
  function(data, group) {
    .corstat(data, group, q, rho)
  }
}

# pruning by Spearman correlation
fsSpearman <- function(q = NULL, rho) {
  function(data, group) {
    data <- t(apply(data, 1, rank))
    .corstat(data, group, q, rho)
  }
}


# pruning by median-split odds-ratio
fsMedSplitOddsRatio <- function(q = NULL, OR) {
  function(data, group) {
    meds <- apply(data, 1, median)
    x <- 1*(sweep(data, 1, meds, '-') > 0)
    A <- levels(group)[1]
    pHA <- as.vector( x %*% matrix(group == A, ncol=1) ) / sum(group == A)
    pHB <- 1 - pHA
    pLA <- as.vector((1 - x) %*% matrix(group == A, ncol=1) ) / sum(group == A)
    pLB <- 1 - pLA
    logodds <- log( (pHA/pLA) / (pHB/pLB) )
    lor <- ifelse(is.null(q), log(OR), quantile(abs(logodds), q))
    abs(logodds) > lor
  }
}

# discretize continuous data
# the 'simple' part of binning inspires the name...
disco <- function(x, k=10, type=c("width", "freq")) {
  type <- match.arg(type)
  cutpts <- switch(type,
                   width = seq(min(x), max(x), length=k+1),
                   freq = quantile(x, seq(0, 1, by=1/k)))
  cutpts <- sort(unique(cutpts)) # may be needed for frequency-based cuts
  cut(x, cutpts, include.lowest = TRUE, labels=FALSE)
}

discodata <- function(data, type=c("width", "freq")) {
  k <- ceiling(min(10, sqrt(ncol(data))))
  t(apply(data, 1, disco, k=k, type=type))
}

fsChisquared <- function(q = NULL, cutoff) {
  function(data, group) {
    dd <- discodata(data, type='freq')
    stats <- apply(dd, 1, function(x) {
      tab <- table(x, group)
        rowsums = apply(tab, 1, sum)
        colsums = apply(tab, 2, sum)
        allsum = sum(colsums)
        expectedmatrix = t(as.matrix(colsums) %*% t(as.matrix(rowsums)))/allsum
        chis = sum((tab - expectedmatrix)^2/expectedmatrix)
        if (chis == 0 || length(colsums) < 2 || length(rowsums) < 2) {
            return(0)
        } else {
          return(sqrt(chis/(allsum * min(length(colsums) - 1,
                                         length(rowsums) - 1))))
        }
    })
    if (!is.null(q)) cutoff <- quantile(stats, q)
    stats > cutoff
  }
}

entropic <- function(x) {
  y <- table(x, useNA='always')
  freqs <- y / sum(y)
  -sum( ifelse(freqs > 0, freqs*log(freqs), 0) )
}

fsEntropy <- function(q = 0.9, kind=c("information.gain",
                                         "gain.ratio",
                                         "symmetric.uncertainty")) {
  kind <- match.arg(kind)
  function(data, group) {
    dd <- discodata(data, 'width')
    itemEnt <- apply(dd, 1, entropic)
    groupEnt <- entropic(group)
    jointEnt <- apply(dd, 1, function(x) {
      entropic(data.frame(group, x))
    })
    IG <- groupEnt + itemEnt - jointEnt
    results <- switch(kind,
                      information.gain = IG,
                      gain.ratio = IG / itemEnt,
                      symmetric.uncertainty = 2 * IG / (itemEnt + groupEnt))
    results > quantile(results, q)
  }
}

fsFisherRandomForest <- function(q) {
  function(data, group) {
    temp <- fsModifiedFisher(q=2*q-1)(data, group)
    daft <- data.frame(Y=group, t(data[temp,]))
    orchard <- randomForest(Y ~ ., daft, ntree = 1000, keep.forest = FALSE, 
                            importance = TRUE)
    imp <- importance(orchard, type=1)
    rating <- rep(min(imp)-1, length(temp))
    rating[temp] <- imp
    rating > quantile(rating, q)
  }
}

fsTailRank <- function(specificity=0.9, tolerance=0.5, confidence=0.5) {
  function(data, group) {
    trt <- TailRankTest(data, group, direction="two",
                        specificity=specificity,
                        tolerance=tolerance,
                        confidence=confidence)
    as.logical(trt)
  }
}
