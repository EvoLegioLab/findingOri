#Summary:
#Data is summarized by snr in pairwise sliding windows
#Input arguments:
#seqdata = vector with sequence data
#seqpos = vector with x-positions
#label = histogram label
#verb = T/F for progress reports and plotting

# snr pairwise by sliding windows ------------------------------------------

snr_sliding <- function(seqdata, seqpos, label, verb){

  data <- seqdata
  xpos <- seqpos
  datalength = length(data)

  windowlength = 200
  steplength = 10

  i = 1
  j = 1
  windowend = 0
  veclength = floor(windowlength/steplength)*floor(datalength/windowlength)
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

    i <- i + steplength
    j <- j + 1
  }
  
  if(verb){hist(skewvec.snr, main="Histogram for snr", xlab=label)}

  outlen <- veclength-1
  snr_windows <- list(xpos=skewvec.pos[1:outlen], mean=skewvec.mean[1:outlen], std=skewvec.std[1:outlen], snr=skewvec.snr[1:outlen])
  return(snr_windows)

}
