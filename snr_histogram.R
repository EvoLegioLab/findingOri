
# snr by window -----------------------------------------------------------

snr <- function(data){

  require("mixtools")
  
  skewdata <- data
  hist(skewdata)
  
  set.seed(100)
  system.time(out2<-normalmixEM(skewdata, arbvar = FALSE, epsilon = 1e-03, fast=TRUE))
  plot(out2,density = TRUE, w = 1.1)
  
  ##  Attach mclust library.
  library(mclust)
  
  ##  Fit bimodal mixture model.
  yBIC = mclustBIC(skewdata, modelNames="V")
  yModel = mclustModel(skewdata, yBIC)
  
  ##  Print model parameters.
  print(yModel$parameters$mean)
  
  print(yModel$parameters$variance$sigmasq)

  #plot(yBIC)
  dens = densityMclust(skewdata)
  plot(dens, skewdata, what = "density")
}
