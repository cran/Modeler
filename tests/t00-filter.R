library(Modeler)
set.seed(246391)
data <- matrix(rnorm(1000*30), nrow=1000, ncol=30)

fm <- filterMean(0.2)
summary(fm(data))
summary(fm(as.data.frame(data)))

fd <- filterMedian(0.2)
summary(fd(data))

fs <- filterSD(1)
summary(fs(data))

fr <- filterRange(3)
fx <- filterMax(2)
fn <- filterMin(-1)
lap <- lapply(list(fm, fd, fs, fr, fx, fn), function(f) summary(f(data)))

try(fm(1:3)) # unclear whether this should work
try(fs(1:3)) # unclear whether this should work
try(fd(1:3)) # fails
try(fr(1:3)) # fails
try(fx(1:3)) # fails
try(fn(1:3)) # fails

try(fx("A")) # fails

