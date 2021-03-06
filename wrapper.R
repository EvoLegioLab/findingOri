#Summary:
#Import a FASTA file to be analyzed for statistics of Ori and Ter
#based on GC, TA, kmer skews and strand bias.
#Input arguments:
#1. Either of:
#access = "NCBI accession number"
#fi = "local file path" (UNIX standard)
#2. Optional:
#verbose = "v" for progress reports and plotting

findingOri<-function(access, fi, verbose){
  require(seqinr)
  require(ape)
  require(stringr)
  require(pracma)
  
  cat("-------------------------------------\n")
  cat("findingOri v1 [September, 2017]\n")
  cat("Uppsala University\n")
  cat("Lionel Guy, et al.\n")
  cat("-------------------------------------\n")
  
  #------------------Calling
  
  # Parse arguments
  if (!missing(access) & missing(fi)){
    readfile = FALSE
    cat ("Reading in fasta from accession number", access, "\n")
  } else if (missing(access) & !missing(fi)) {
    readfile = TRUE
    cat("Reading in fasta from file", fi, "\n")
  } else {
    cat("Unspecified arguments!\n")
    cat("Call findingOri in at least one of the following ways:\n")
    cat("findingOri(accessnumber)\n")
    cat("or findingOri(,fastafile)\n")
    cat("Add","'v(erbose)'"," for slower calculation with figures\n")
    cat("on the third parameter position in the function call.")
    stop()
  }
  
  if (missing(verbose)){
    verb = FALSE
    cat("No figures will be produced!\n")
  } else {
    verb = TRUE
    cat("Figures will be produced!\n")
  }
  
  # Read fasta file
  if (readfile){
    
    # # Test files
    # source("readFasta.R")
    # myDNAbin <- readFasta(verb)
    
    # User files
    source("readseqfile.R")
    myDNAbin <- readseqfile(fi, verb)
    #print(myDNAbin)
    if (verb){cat("Processing sequence data...")}
    myFasta<-write.fasta(sequences=myDNAbin, names=names(myDNAbin),file.out ="myFasta.fasta")
    mySeq<-read.fasta("myFasta.fasta")
    if (verb){cat("done!\n")}
  }
  
  # Fetch fasta by accession number
  if (!readfile) {mySeq<-getfastafiles(access)}
  
  # Assemble nucleotide and gene statistics 
  source('GC_skew.R')
  rawoutput<-GCall(mySeq,verb)
  
  #Calculate position of ori and ter from snr
  if (verb){cat("Processing snr data...")}
  source("snr_ori_ter.R")
  snrgc <- snr_ori_ter(rawoutput$gc, rawoutput$gcxpos, 'gc sliding snr', verb)
  snrta <- snr_ori_ter(rawoutput$ta, rawoutput$gcxpos, 'ta sliding snr', verb)
  snrgc3 <- snr_ori_ter(rawoutput$gc3, rawoutput$gc3xpos, 'gc3 sliding snr', verb)
  snrta3 <- snr_ori_ter(rawoutput$ta3, rawoutput$gc3xpos, 'ta3 sliding snr', verb)
  if (verb){cat("done!\n")}
  
  #snr output
  cat('Based on gc skew snr, the origin is located at', snrgc$neg, 'and the terminus at',snrgc$pos,'\n')
  cat('Based on ta skew snr, the origin is located at', snrta$pos,' and the terminus at',snrta$neg,'\n')
  cat('Based on gc3 skew snr, the origin is located at', snrgc3$neg, 'and the terminus at',snrgc3$pos,'\n')
  cat('Based on ta3 skew snr, the origin is located at', snrta3$pos, 'and the terminus at',snrta3$neg,'\n')
  
  # Assemble kmer statistics
  source('kmer_both_strands_updated.R')
  kmerout<-kmers(mySeq, verb)
  #rawoutput<-list(gc=gccount$gc,gc3=gc3s$gc,gcxpos=gccount$xpos,gc3xpos=gc3s$xpos,genebias=cumgenebias$genes,genexpos=cumgenebias$xpos)
  oriout<-list(gcori=rawoutput$gcori, gc3ori=rawoutput$gc3ori, gcter=rawoutput$gcter, gc3ter=rawoutput$gc3ter, taori=rawoutput$taori, ta3ori=rawoutput$ta3ori,tater=rawoutput$tater,ta3ter=rawoutput$ta3ter,geneori=rawoutput$geneori, geneter=rawoutput$geneter)
  #kmerout<-seq_skews
    
  #calculating the position of ori and ter from kmers
  kmer1<-kmerout$pos_interp[match(max(kmerout$tot_change_inter),kmerout$tot_change_inter)]
  #print(kmer1)
  seqlen<-length(getSequence(mySeq[[1]]))
  #print(seqlen)
  if(kmer1+(seqlen/2)>seqlen){
    #peak_diff<-(kmer1-seqlen)-round((seqlen/2))
    #kmer2<-1+peak_diff
    kmer2 <- kmer1 + round(seqlen/2) - seqlen
  }else{
    kmer2<-kmer1+round((seqlen/2))
  }
  margin<-round((0.1*seqlen))
  #print(margin)
  if(kmer2+margin>seqlen){
    #kmer2<-max(kmerout$tot_change_inter[which(kmerout$pos_interp%in%c((kmer2-margin):seqlen))])
    kmer2<-max(kmerout$tot_change_inter[which(kmerout$pos_interp%in%c((kmer2-margin):seqlen,1:(kmer2+margin-seqlen)))])
    kmer2pos<-kmerout$pos_interp[which(kmerout$tot_change_inter%in%kmer2)]
  }else if(kmer2-margin<1){
    #kmer2<-max(kmerout$tot_change_inter[which(kmerout$pos_interp %in% c(1:(kmer2+margin)))])
    kmer2<-max(kmerout$tot_change_inter[which(kmerout$pos_interp %in% c((seqlen+kmer2-margin):seqlen,1:(kmer2+margin)))])
    kmer2pos<-kmerout$pos_interp[which(kmerout$tot_change_inter %in% kmer2)]
  }else{
    kmer2<-max(kmerout$tot_change_inter[which(kmerout$pos_interp %in% c((kmer2-margin):(kmer2+margin)))])
    kmer2pos<-kmerout$pos_interp[which(kmerout$tot_change_inter %in% kmer2)]
  }
  #print(kmer2)
  #print(kmer2pos)
  if (oriout$gc3ori>oriout$gc3ter){
    if(kmer1>kmer2pos){
      kmerter<-kmer2pos
      kmerori<-kmer1
    }else{
      kmerter<-kmer1
      kmerori<-kmer2pos
    }
  }else{
    if(kmer1>kmer2pos){
      kmerter<-kmer1
      kmerori<-kmer2pos
    }else{
      kmerter<-kmer2pos
      kmerori<-kmer1
    }
  }
  #kmer output
  cat('The position of the origin based on k-mer skew is at',kmerori,'\n')
  cat('The position of the terminus based on k-mer skew is at',kmerter,'\n')

  #Calling snr for weigthed averages
  source('snr.R')
  gcsnr<-snr(rawoutput$gc, 'gc snr', verb)[,'snr']
  gc3snr<-snr(rawoutput$gc3, 'gc3 snr', verb)[,'snr']
  #genesnr<-snr(rawoutput$genebias, 'gene snr', verb)[,'snr'] cannot be calculated for genebias in a comparable way
  #kmersnr<-snr(kmerout$tot_change_inter, 'kmer snr', verb)[,'snr']
  tasnr<-snr(rawoutput$ta, 'ta snr', verb)[,'snr']
  ta3snr<-snr(rawoutput$ta3, 'ta3 snr', verb)[,'snr']
  #print(gcsnr)
  #print(gc3snr)
  #print(genesnr)
  #print(kmersnr)
  #gcweight<-gcsnr/sum(c(gcsnr,gc3snr, kmersnr, tasnr, ta3snr))
  #gc3weight<-gc3snr/sum(c(gcsnr,gc3snr, kmersnr, tasnr, ta3snr))
  gcweight<-gcsnr/sum(c(gcsnr,gc3snr, tasnr, ta3snr))
  gc3weight<-gc3snr/sum(c(gcsnr,gc3snr, tasnr, ta3snr))
  #kmerweight<-kmersnr/sum(c(gcsnr,gc3snr, kmersnr, tasnr, ta3snr))
  #taweight<-tasnr/sum(c(gcsnr,gc3snr, kmersnr, tasnr, ta3snr))
  #ta3weight<-ta3snr/sum(c(gcsnr,gc3snr, kmersnr, tasnr, ta3snr))
  taweight<-tasnr/sum(c(gcsnr,gc3snr, tasnr, ta3snr))
  ta3weight<-ta3snr/sum(c(gcsnr,gc3snr, tasnr, ta3snr))
  #geneweight<-1
  #Simple statistics with all metrics 
  combOri<-c(oriout$gcori,oriout$gc3ori,kmerori, oriout$geneori, oriout$taori, oriout$ta3ori)
  SEori<-sd(combOri)/sqrt(length(combOri))*1.96
  #Ori_upper<-mean(combOri)+SEori
  #Ori_lower<-mean(combOri)-SEori
  combTer<-c(oriout$gcter,oriout$gc3ter,kmerter, oriout$geneter,oriout$tater, oriout$ta3ter)
  SEter<-sd(combTer)/sqrt(length(combTer))*1.96
  #Ter_upper<-mean(combTer)+SEter
  #Ter_lower<-mean(combTer)-SEter
  cat('The combined measure of origin is at position',mean(combOri),'+-',SEori,'\n')
  cat('The combined measure of terminus is at position',mean(combTer),'+-',SEter,'\n')
  #for plotting
  orilimits<-SEori/(1/4*seqlen)
  terlimits<-SEter/(1/4*seqlen)
  if(verb==TRUE){
    plot(c(-1,1),c(1,-1), type='n', asp=1)
    title(main='Graphical representation of calculated ori and ter positions (black dots) with 95% confidence intervals (red line)')
    text(0.5,0.96,paste('Origin of replication at',round(mean(combOri))), pos=4,cex=1.2)
    text(0.5,-0.96,paste('Terminus of replication at',round(mean(combTer))),pos=4, cex=1.2)
    radius<-1
    theta<-seq(0,2*pi, length=100)
    x<-radius*cos(theta)
    y<-radius*sin(theta)
    lines(x,y)
    points(mean(x),max(y), pch=19, cex=2 )
    points(mean(x),min(y), pch=19, cex=2)
    lines(c(-orilimits,orilimits),c(1,1), type='l',col='red', lwd=3)
    lines(c(-terlimits,terlimits),c(-1,-1), type='l',col='red', lwd=3)
    
  }
  #Stats with weights, kmers and gene bias excluded 
#  weightedcombOri<-c(gcweight*oriout$gcori,gc3weight*oriout$gc3ori,kmerweight*kmerori, taweight*oriout$taori, ta3weight*oriout$ta3ori)
#  weightedcombOri<-c(gcweight*oriout$gcori, gc3weight*oriout$gc3ori, taweight*oriout$taori, ta3weight*oriout$ta3ori)
#  wSEori<-sd(weightedcombOri)/sqrt(length(weightedcombOri))*1.96
  weightedcombOri <- sum(gcweight*oriout$gcori, gc3weight*oriout$gc3ori, taweight*oriout$taori, ta3weight*oriout$ta3ori)
# Assume variance is the same in all factors to derive SE:
  weightfactorsOri <- c(oriout$gcori, oriout$gc3ori, oriout$taori, oriout$ta3ori)
  wSEori <- sd(weightfactorsOri)/sqrt(length(weightfactorsOri))*1.96
  #wOri_upper<-mean(weightedcombOri)+wSEori
  #wOri_lower<-mean(weightedcombOri)-wSEori
#  weightedcombTer<-c(gcweight*oriout$gcter,gc3weight*oriout$gc3ter,kmerweight*kmerter, taweight*oriout$tater, ta3weight*oriout$ta3ter)
#  weightedcombTer<-c(gcweight*oriout$gcter, gc3weight*oriout$gc3ter, taweight*oriout$tater, ta3weight*oriout$ta3ter)
#  wSEter<-sd(weightedcombTer)/sqrt(length(weightedcombTer))*1.96
  weightedcombTer <- sum(gcweight*oriout$gcter, gc3weight*oriout$gc3ter, taweight*oriout$tater, ta3weight*oriout$ta3ter)
  # Assume variance is the same in all factors to derive SE:
  weightfactorsTer <- c(oriout$gcter, oriout$gc3ter, oriout$tater, oriout$ta3ter)
  wSEter <- sd(weightfactorsTer)/sqrt(length(weightfactorsTer))*1.96
  #wTer_upper<-mean(weightedcombTer)+wSEter
  #wTer_lower<-mean(weightedcombTer)-wSEter
  cat('The combined measure of origin from weighted means is at position', weightedcombOri,'+-',wSEori,'\n')
  cat('The combined measure of terminus from weighted means is at position', weightedcombTer,'+-',wSEter,'\n')
  #for plotting
  worilimits<-wSEori/(1/4*seqlen)
  wterlimits<-wSEter/(1/4*seqlen)
  if(verb==TRUE){
    plot(c(-1,1),c(1,-1), type='n', asp=1)
    title(main='Graphical representation of calculated ori and ter positions from \n weighted means (black dots) with 95% confidence intervals (red line)')
    text(0.5,0.96,paste('Weighted origin of replication at',round(weightedcombOri)), pos=4,cex=1.2)
    text(0.5,-0.96,paste('Weighted terminus of replication at',round(weightedcombTer)),pos=4,cex=1.2)
    radius<-1
    theta<-seq(0,2*pi, length=100)
    x<-radius*cos(theta)
    y<-radius*sin(theta)
    lines(x,y)
    points(mean(x),max(y), pch=19, cex=2 )
    points(mean(x),min(y), pch=19, cex=2)
    lines(c(-worilimits,worilimits),c(1,1), type='l',col='red', lwd=3)
    lines(c(-wterlimits,wterlimits),c(-1,-1), type='l',col='red', lwd=3)
    
  }
}