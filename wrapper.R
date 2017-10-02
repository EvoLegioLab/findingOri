#Summary:
#Import a FASTA file to be analyzed for statistics of Ori and Ter
#based on GC, TA, kmer skews and strand bias.
#Input arguments:
#1. Either of:
#access = "NCBI acession number"
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
    if (verb){cat("Processing sequence data...")}
    myFasta<-write.dna(myDNAbin,file ="myFasta.fasta", format = "fasta")
    mySeq<-read.fasta("myFasta.fasta")
    if (verb){cat("done!\n")}
  }
  
    source('GC_skew.R')
    source('kmer_both_strands_updated.R')
    if (!readfile) {mySeq<-getfastafiles(access)}

    rawoutput<-GCall(mySeq,verb)
    kmerout<-kmers(mySeq, verb)
    #rawoutput<-list(gc=gccount$gc,gc3=gc3s$gc,gcxpos=gccount$xpos,gc3xpos=gc3s$xpos,genebias=cumgenebias$genes,genexpos=cumgenebias$xpos)
    oriout<-list(gcori=rawoutput$gcori, gc3ori=rawoutput$gc3ori, gcter=rawoutput$gcter, gc3ter=rawoutput$gc3ter, geneori=rawoutput$geneori, geneter=rawoutput$geneter)
    #kmerout<-seq_skews

  #calculating the position of ori and ter from kmers
  kmer1<-kmerout$pos_interp[match(max(kmerout$tot_change_inter),kmerout$tot_change_inter)]
  #print(kmer1)
  seqlen<-length(getSequence(mySeq[[1]]))
  #print(seqlen)
  if(kmer1+seqlen/2>seqlen){
    peak_diff<-(kmer1-seqlen)-(seqlen/2)
    kmer2<-1+peak_diff
  }else{
    kmer2<-kmer1+(seqlen/2)
  }
  margin<-500000
  if(kmer2+margin>seqlen){
    kmer2<-max(kmerout$tot_change_inter[which(kmerout$pos_interp%in%c((kmer2-margin):seqlen))])
    kmer2pos<-kmerout$pos_interp[which(kmerout$tot_change_inter%in%kmer2)]
  }else if(kmer2-margin<1){
    kmer2<-max(kmerout$tot_change_inter[which(kmerout$pos_interp %in% c(1:(kmer2+margin)))])
    kmer2pos<-kmerout$pos_interp[which(kmerout$tot_change_inter %in% kmer2)]
  }else{
    kmer2<-max(kmerout$tot_change_inter[which( kmerout$pos_interp %in% c((kmer2-margin):(kmer2+margin)))])
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
  #Simple statistics with all metrics 
  combOri<-c(oriout$gcori,oriout$gc3ori,kmerori, oriout$geneori)
  SEori<-sd(combOri)/sqrt(length(combOri))*1.96
  Ori_upper<-mean(combOri)+SEori
  Ori_lower<-mean(combOri)+SEori
  combTer<-c(oriout$gcter,oriout$gc3ter,kmerter, oriout$geneter)
  SEter<-sd(combTer)/sqrt(length(combTer))*1.96
  Ter_upper<-mean(combTer)+SEter
  Ter_lower<-mean(combTer)+SEter
  cat('The combined measure of origin is at position',mean(combOri),'+-',SEori,'\n')
  cat('The combined measure of terminus is at position',mean(combTer),'+-',SEter,'\n')
  #for plotting
  orilimits<-SEori/(1/4*seqlen)
  terlimits<-SEter/(1/4*seqlen)
  if(verb==TRUE){
    plot(c(-1,1),c(1,-1), type='n', asp=1)
    title(main='Gaphical representation of calculated ori and ter positions (black dots) with 95% confidence intervals (red line)')
    text(0.5,1,paste('Origin of replication at',round(mean(combOri))), pos=4)
    text(0.5,-1,paste('Terminus of replication at',round(mean(combTer))),pos=4)
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
}