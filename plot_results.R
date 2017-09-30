
# plot raw & cum skew & window rank ----------------------------------------------

plotresults <- function(seqdata){
  
  source("snr_rank.R")
  ranked_windows <- snr_rank(seqdata$gc,seqdata$xpos)
  pos_change <- ranked_windows$pos
  neg_change <- ranked_windows$neg
  
  source("snr_windows.R")
  raw_windows <- snr_windows(seqdata$gc,seqdata$xpos)
  pos_snr <- which.max(raw_windows$snr)
  neg_snr <- which.min(raw_windows$snr)
  
  min_y = min(c(min(seqdata$cumgc), min(seqdata$gc)))
  max_y = max(c(max(seqdata$cumgc), max(seqdata$gc)))

  plot(seqdata$xpos, seqdata$cumgc, type="l", main = "skew statistics", xlab="Position",ylab="Black=cumGC, Blue=GC", ylim=(c(min_y,max_y)))
  lines(seqdata$xpos, 500*seqdata$gc, type = "l", col = "blue")
  points(pos_change[1:1,1], pos_change[1:1,2], type = "p", col = "red", lwd=4)
  points(neg_change[1:1,1], neg_change[1:1,2], type = "p", col = "green", lwd=4)
  points(raw_windows$xpos[pos_snr], raw_windows$snr[pos_snr], type = "p", col = "black", lwd=4)
  points(raw_windows$xpos[neg_snr], raw_windows$snr[neg_snr], type = "p", col = "cadetblue", lwd=4)
}