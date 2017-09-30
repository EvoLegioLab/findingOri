
# rank window change by snr -----------------------------------------------------------

snr_rank <- function(seqdata, seqpos){

  wrapdata <- c(tail(seqdata, 500), seqdata, head(seqdata, 500))
  wrappos <- c(tail(seqpos, 500), seqpos, head(seqpos, 500))
  
  # # kmers:
  # wrapdata <- seqdata
  # wrappos <- seqpos
  
  # source("snr_pairwise.R")
  # snrdata <- snr_pairwise(wrapdata, wrappos)
  source("snr_sliding.R")
  snrdata <- snr_sliding(wrapdata, wrappos)
  veclength = length(snrdata$xpos)
  pos_change <- list(xpos=NULL, snr=NULL)
  neg_change <- list(xpos=NULL, snr=NULL)
  
  # list positive signals by position
  for(i in 1:veclength) {
       if(snrdata$mean[i] < 0){
         pos_change$xpos <- c(pos_change$xpos, snrdata$xpos[i])
         pos_change$snr <- c(pos_change$snr, snrdata$snr[i])
       }
    }

  # list negative signals by position
  for(i in 1:veclength) {
    if(snrdata$mean[i] > 0){
      neg_change$xpos <- c(neg_change$xpos, snrdata$xpos[i])
      neg_change$snr <- c(neg_change$snr, snrdata$snr[i])
    }
  }
  
  # sort positive signals by snr
  m <- matrix(unlist(pos_change), ncol=2)
  pos_sorted <- m[sort.list(-m[,2]),]
  #print(pos_sorted)
  
  # sort negative signals by snr
  m <- matrix(unlist(neg_change), ncol=2)
  neg_sorted <- m[sort.list(-m[,2]),]
  #print(neg_sorted)
  
  snr_rank <- list(pos=pos_sorted, neg=neg_sorted)
  return(snr_rank)

}
