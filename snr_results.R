
# plot raw & cum skew & best pairwise snr windows ----------------------------------------------

snr_results <- function(){
  
  ## 1 raw gc
  # source("GC_skew.R")
  # seqdata <- gccount
  # skewdata <- seqdata$gc
  # skewpos <- seqdata$xpos
  # cumgc <- seqdata$cumgc
  
  ## 2 gc3s gc
  # source("GC_skew.R")
  seqdata <- gc3s
  skewdata <- seqdata$gc
  skewpos <- seqdata$xpos
  cumgc <- seqdata$cumgc
  
  ## 3 raw ta
  # source("GC_skew.R")
  # seqdata <- gccount
  # skewdata <- seqdata$ta
  # skewpos <- seqdata$xpos
  # cumgc <- seqdata$cumgc
  
  # 4 gc3s ta
  # source("GC_skew.R")
  # seqdata <- gc3s
  # skewdata <- seqdata$ta
  # skewpos <- seqdata$xpos
  # cumgc <- seqdata$cumgc
  
  # 5 kmers
  # seqdata <- source("kmer_both_strands.R")
  # skewdata <- seq_skews[[3]]
  # skewpos <- seq_skews[[1]]
  # cumgc <- seq_skews[[2]]

  source("snr_rank.R")
  ranked_windows <- snr_rank(skewdata, skewpos)
  pos_change <- ranked_windows$pos
  neg_change <- ranked_windows$neg
  
  source("snr_windows.R")
  raw_windows <- snr_windows(skewdata, skewpos)
  pos_snr <- which.max(raw_windows$snr)
  neg_snr <- which.min(raw_windows$snr)
  
  min_y = min(c(min(cumgc), min(skewdata)))
  max_y = max(c(max(cumgc), max(skewdata)))

  plot(skewpos, cumgc, type="l", main = "skew statistics", xlab="Position",ylab="Black=cumGC, Blue=skew", ylim=(c(min_y,max_y)))
  lines(skewpos, 200*skewdata, type = "l", col = "blue")
  points(pos_change[1:1,1], pos_change[1:1,2], type = "p", col = "red", lwd=4)
  points(neg_change[1:1,1], neg_change[1:1,2], type = "p", col = "green", lwd=4)
  
  #kmers:
  # plot(skewpos, cumgc, type="l", main = "skew statistics", xlab="Position",ylab="Black=cumGC, Blue=skew", ylim=(c(min_y,max_y)))
  # lines(skewpos, skewdata, type = "l", col = "blue")
  # points(pos_change[1], pos_change[2], type = "p", col = "red", lwd=4)
  # points(neg_change[1], neg_change[2], type = "p", col = "green", lwd=4)
  
  #both:
  points(raw_windows$xpos[pos_snr], raw_windows$snr[pos_snr], type = "p", col = "black", lwd=4)
  points(raw_windows$xpos[neg_snr], raw_windows$snr[neg_snr], type = "p", col = "cadetblue", lwd=4)
}