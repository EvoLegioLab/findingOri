
# check for zero crossings ------------------------------------------------

require("changepoint")

#Uses gccount

candidates <- function(gcdata){
  
  data <- c(gcdata, gcdata[1:10000])
  # datalength <- length(data)
  # chunklength <- 100
  # i <- 1
  # j <- 1
  # gcvec <- numeric(datalength/chunklength)
  # while (i < datalength - chunklength){
  #   gcvec[j]=mean(data[i:i+chunklength])
  #   i <- i + chunklength
  #   j <- j + 1
  # }
  gcvec <- data
  ts.plot(gcvec, xlab = "Index")
  gcvec.pelt <- cpt.mean(gcvec, method = "BinSeg", Q = 2)
  plot(gcvec.pelt, type = "l", cpt.col = "blue", xlab = "Index", cpt.width = 4)
  print(cpts(gcvec.pelt))
}
