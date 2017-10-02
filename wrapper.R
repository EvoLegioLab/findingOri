#Just starting to create the wrapper. Need to change plotting in both GC skew and k-mer so that it is optional.
#Also need to link file fetching here, and then just input that to sources
findingOri<-function(access, fi, verbose){
  require(seqinr)
  require(ape)
  require(stringr)
  require(pracma)
  #------------------Calling with file selection
  if (missing(fi)==!TRUE){
    #cat ('calling', fi)
    #fi<-fi
    if (missing(verbose)){
      print ('figures will be produced!')
      verb<-'n'
    } else {
      print('no figures will be produced!')
      verb<-'y'
    }
    source("readFasta.R")
    myDNAbin<-readFasta()
    myFasta<-write.dna(myDNAbin,file ="myFasta.fasta", format = "fasta")
    mySeq<-read.fasta("myFasta.fasta")
    #mySeq<-readFasta()
    #print(mySeq[1])
    source('GC_skew.R')
    source('kmer_both_strands_updated.R')
    #mySeq<-getDNAbin(myDNAbin)
    rawoutput<-GCall(mySeq,verb)
    kmerout<-kmers(mySeq, verb)
    #rawoutput<-list(gc=gccount$gc,gc3=gc3s$gc,gcxpos=gccount$xpos,gc3xpos=gc3s$xpos,genebias=cumgenebias$genes,genexpos=cumgenebias$xpos)
    oriout<-list(gcori=rawoutput$gcori, gc3ori=rawoutput$gc3ori, gcter=rawoutput$gcter, gc3ter=rawoutput$gc3ter, geneori=rawoutput$geneori, geneter=rawoutput$geneter)
    #kmerout<-seq_skews
    #----------------Calling with accecssion number
    } else if (missing(access)==!TRUE){
    cat ('reading in fasta from accession number', access)
    access<-access
      if (missing(verbose)){
      print ('figures will be produced!')
      verb<-'n'
      } else {
      print('no figures will be produced!')
      verb<-'y'
      }
    source('GC_skew.R')
    source('kmer_both_strands_updated.R')
    mySeq<-getfastafiles(access)
    #print(mySeq[1])
    rawoutput<-GCall(mySeq,verb)
    kmerout<-kmers(mySeq,verb)
    #rawoutput<-list(gc=gccount$gc,gc3=gc3s$gc,gcxpos=gccount$xpos,gc3xpos=gc3s$xpos,genebias=cumgenebias$genes,genexpos=cumgenebias$xpos)
    oriout<-list(gcori=rawoutput$gcori, gc3ori=rawoutput$gc3ori, gcter=rawoutput$gcter, gc3ter=rawoutput$gc3ter, geneori=rawoutput$geneori, geneter=rawoutput$geneter)
    #kmerout<-seq_skews
    } else {
    cat('Call findingOri with at least one of the following ways:\n findingOri(accessnumber)\n or findingOri(,fastafile)\n add call to verbose for faster calculations without figures \n on the third parameter position in the function call')
    stop()
    }
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
  if(missing(verb)){
    plot(c(-1,1),c(1,-1), type='n', asp=1)
    radius<-1
    theta<-seq(0,2*pi, length=200)
    x<-radius*cos(theta)
    y<-radius*sin(theta)
    lines(x,y)
    points(mean(x),testori)
    points(mean(x),min(y))
    lines(c(Oriup,orilow),c(1,1), type='l',col='red')
  }
}