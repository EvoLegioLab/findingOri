
# snr pairwise by window -----------------------------------------------------------

snr_pairwise<- function(seqdata, seqpos){

  #data <- c(seqdata, seqdata[1:10000])
  data <- seqdata
  #xpos <- c(seqpos, seqpos[1:10000])
  xpos <- seqpos
  datalength = length(data)

  windowlength = 200
  # # kmers:
  # windowlength = 10
  
  i = 1
  j = 1
  windowend = 0
  veclength = floor(datalength/windowlength)
  skewvec.pos <- numeric(veclength)
  skewvec.mean <- numeric(veclength)
  skewvec.std <- numeric(veclength)
  skewvec.snr <- numeric(veclength)
  
  while (i < datalength - 2*windowlength){

    window1start <- i
    window1stop <- window1start + windowlength
    window2start <- window1stop + 1
    window2stop <- window2start + windowlength
    
    skewvec.pos[j] <- xpos[window1stop]
    mean1 <- mean(data[window1start:window1stop])
    mean2 <- mean(data[window2start:window2stop])
    skewvec.mean[j] <- mean2-mean1
    skewvec.std[j] <- sd((data[window1start:window1stop]-mean1) + (data[window2start:window2stop]-mean2))
    
    skewvec.snr[j] <- skewvec.mean[j]^2/skewvec.std[j]^2
    # skewvec.snr[j] <- skewvec.mean[j]/skewvec.std[j]
    
    i <- i + windowlength
    j <- j + 1
  }
  
  #hist(skewvec.snr)

  outlen <- veclength-1
  snr_windows <- list(xpos=skewvec.pos[1:outlen], mean=skewvec.mean[1:outlen], std=skewvec.std[1:outlen], snr=skewvec.snr[1:outlen])
  return(snr_windows)

}
