
# snr by window -----------------------------------------------------------

snr_windows<- function(seqdata, seqpos){

  data <- c(tail(seqdata, 5000), seqdata, head(seqdata, 5000))
  xpos <- c(tail(seqpos, 5000), seqpos, head(seqpos, 5000))
  
  # # kmers:
  # wrapdata <- seqdata
  # wrappos <- seqpos
  
  datalength = length(data)
  
  windowlength = 2000
  # # kmers:
  # windowlength = 10
  
  i = 1
  j = 1
  windowend = 0
  veclength = floor(datalength/windowlength)
  skewvec.mid_window <- numeric(veclength)
  skewvec.mean <- numeric(veclength)
  skewvec.std <- numeric(veclength)
  skewvec.snr <- numeric(veclength)
  while (i < datalength - windowlength){
    windowend <- i + windowlength
    skewvec.mid_window[j] <- xpos[i+floor(windowlength/2)]
    skewvec.mean[j] <- mean(data[i:windowend])
    skewvec.std[j] <- sd(data[i:windowend])
    #skewvec.snr[j] <- skewvec.mean[j]^2/skewvec.std[j]^2
    skewvec.snr[j] <- skewvec.mean[j]/skewvec.std[j]
    i <- i + windowlength
    j <- j + 1
  }
  
  # print(skewvec.mid_window)
  # print(skewvec.mean)
  # print(skewvec.std)
  # min_y = min(c(min(skewvec.mean), min(skewvec.std)))
  # max_y = max(c(max(skewvec.mean), max(skewvec.std)))
  # 
  # plot(skewvec.index, skewvec.mean, type="l", main = "skew statistics", xlab="Index",ylab="Black=Mean, Green=Std.dev", ylim=(c(min_y,max_y)))
  # lines(skewvec.std, type = "l", col = "green")

  # hist(skewvec.snr)

  snr_windows <- list(xpos=skewvec.mid_window, mean=skewvec.mean, std=skewvec.std, snr=skewvec.snr)
  return(snr_windows)

}
