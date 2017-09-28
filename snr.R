
# snr by window -----------------------------------------------------------

snr <- function(data){
  
  datalength = length(data)
  windowlength = 100
  i = 1
  j = 1
  windowend = 0
  veclength = floor(datalength/windowlength)
  skewvec.index <- numeric(veclength)
  skewvec.mean <- numeric(veclength)
  skewvec.std <- numeric(veclength)
  while (i < datalength - windowlength){
    windowend <- i + windowlength
    skewvec.index[j] <- j
    skewvec.mean[j] <- mean(data[i:windowend])
    skewvec.std[j] <- sd(data[i:windowend])
    i <- i + windowlength
    j <- j + 1
  }
  
  # print(skewvec.mean)
  # print(skewvec.std)
  # min_y = min(c(min(skewvec.mean), min(skewvec.std)))
  # max_y = max(c(max(skewvec.mean), max(skewvec.std)))
  # 
  # plot(skewvec.index, skewvec.mean, type="l", main = "skew statistics", xlab="Index",ylab="Black=Mean, Green=Std.dev", ylim=(c(min_y,max_y)))
  # lines(skewvec.std, type = "l", col = "green")

  hist(skewvec.std)

  snr <- data.frame(mean=mean(skewvec.std), std=sd(skewvec.std), snr=mean(skewvec.std)/sd(skewvec.std))
  return(snr)

}
