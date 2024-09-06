Rfiles <- dir(".", pattern="\\.[R|r]$", full.names=TRUE)
library(parallel)
clus <- makeCluster(length(Rfiles))
parLapply(clus, Rfiles, source)
