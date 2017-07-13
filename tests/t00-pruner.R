library(Modeler)
set.seed(246391)
data <- matrix(rnorm(1000*36), nrow=1000, ncol=36)
data[1:50, 1:18] <- data[1:50, 1:18] + 1
status <- factor(rep(c("A", "B"), each=18))

fsel <- fsTtest(fdr=0.10, ming=125)
summary(fsel(data, status))

fsel <- fsModifiedFisher(q = 0.75)
summary(fsel(data, status))

fsel <- fsPearson(q = 0.9)
summary(fsel(data, status))
fsel <- fsPearson(rho=0.3)
summary(fsel(data, status))

fsel <- fsSpearman(q=0.9)
summary(fsel(data, status))
fsel <- fsSpearman(rho=0.3)
summary(fsel(data, status))

fsel <- fsMedSplitOddsRatio(OR = 5)
summary(fsel(data, status))
fsel <- fsMedSplitOddsRatio(q = 0.9) # note effect opf ties
summary(fsel(data, status))

fsel <- fsChisquared(cutoff=0.5)
summary(fsel(data, status))
fsel <- fsChisquared(q = 0.9)
summary(fsel(data, status))

fsel <- fsEntropy()
summary(fsel(data, status))
fsel <- fsEntropy(kind="gain.ratio")
summary(fsel(data, status))
fsel <- fsEntropy(kind="symm")
summary(fsel(data, status))

fsel <- fsFisherRandomForest(q = 0.9)
summary(fsel(data, status))

fsel <- fsTailRank(spec=0.9, conf=0.5)
summary(fsel(data, status))

