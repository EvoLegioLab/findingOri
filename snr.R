#Summary:
#Data is summarized by 
#Input arguments:
#data = vector with data
#label = histogram label
#verb = T/F for progress reports and plotting

# snr by window -----------------------------------------------------------

snr <- function(data, label, verb){
  
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
  
  if (verb){hist(skewvec.std, main = label)}
  snr <- data.frame(mean=mean(skewvec.std), std=sd(skewvec.std), snr=mean(skewvec.std)/sd(skewvec.std))
  
  return(snr)

}
