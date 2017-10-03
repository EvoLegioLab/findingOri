#Summary:
#Data is summarized by snr in pairwise sliding windows
#Input arguments:
#seqdata = vector with sequence data
#seqpos = vector with x-positions
#label = histogram label
#verb = T/F for progress reports and plotting

# rank window change by snr -----------------------------------------------------------

snr_ori_ter <- function(seqdata, seqpos, label, verb){

  # Wrap at least two windows over sequence ends which can contain Ori or Ter
  # Winow size is 200 xpos
  wrapdata <- c(tail(seqdata, 500), seqdata, head(seqdata, 500))
  wrappos <- c(tail(seqpos, 500), seqpos, head(seqpos, 500))
  
  # Collect snr by sliding windows
  source("snr_sliding.R")
  snrdata <- snr_sliding(wrapdata, wrappos, label, verb)
  
  # List ori and ter candidates
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
  
  # Return ori and ter candidates
  # For gc ori will be positive skew change, for ta negative
  snr_candidates <- list(pos=max(pos_change), neg=max(neg_change)
  return(snr_candidates)

}
