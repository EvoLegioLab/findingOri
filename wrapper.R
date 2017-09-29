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
  print(kmer1)
  seqlen<-length(getSequence(mySeq[[1]]))
  if(kmer1+seqlen/2>seqlen){
    peak_diff<-(kmer1-seqlen)-(seqlen/2)
    kmer2<-1+peak_diff
  }else{
    kmer2<-kmer1+(seqlen/2)
  }
  margin<-500000
  if(kmer2+margin>seqlen){
    kmer2<-max(kmerout$tot_change_inter[match(c((kmer2-margin):seqlen),kmerout$pos_interp)])
    kmer2pos<-kmerout$pos_interp[match(kmer2,kmerout$tot_change)]
  }else if(kmer2-margin<1){
    kmer2<-max(kmerout$tot_change_inter[match(c(1:(kmer2+margin)),kmerout$pos_interp)])
    kmer2pos<-kmerout$pos_interp[kmer2%in%kmerout$tot_change]
  }else{
    kmer2<-max(kmerout$tot_change_inter[c((kmer2-margin):(kmer2+margin)) %in% kmerout$pos_interp])
    kmer2pos<-kmerout$pos_interp[match(kmer2,kmerout$tot_change)]
  }
  print(kmer2)
  print(kmer2pos)
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
  combOri<-c(gcori,gc3ori,kmerori, geneori)
  SEori<-sd(combOri)/sqrt(length(combOri))*1.96
  Ori_upper<-mean(combOri)+SEori
  Ori_lower<-mean(combOri)+SEori
  combTer<-c(gcter,gc3ter,kmerter, geneter)
  SEter<-sd(combTer)/sqrt(length(combTer))*1.96
  Ter_upper<-mean(combTer)+SEter
  Ter_lower<-mean(combTer)+SEter
}
#-----------------Output from scripts and modelling, TRIALS WITH MANOVA! NOT IN USE IN THE LAST VERSION----------------
  #output<-list(gc=gccount$cumgc,gc3=gc3s$cumgc,gcxpos=gccount$xpos,gc3xpos=gc3s$xpos,genebias=cumgenebias$cumgene,genexpos=cumgenebias$xpos)#,kmer=...kmerxpos...)
  output<-list(gc=gccount$gc,gc3=gc3s$gc,gcxpos=gccount$xpos,gc3xpos=gc3s$xpos,genebias=cumgenebias$genes,genexpos=cumgenebias$xpos)
  genelen<-list()
  smallwinlen<-list()
  for (i in 1:length(annotinfo)){
    genelen[i]<-list((annotinfo[[i]][2]-annotinfo[[i]][1]))
  }
  genelen<-unlist(genelen, recursive=TRUE)
  if(2*max(genelen)>(0.02*length(getSequence(mySeq[[1]])))){
    smallwinlen<-ceiling(2*max(genelen))
  }else{
    smallwinlen<-ceiling(0.02*length(getSequence(mySeq[[1]])))
  }
  #alternative
  #smallwinlen<-ceiling(0.002*length(getSequence(mySeq[[1]])))
  bigwin<-4*smallwinlen
  bigstart<-1
  bigend<-bigstart+bigwin
  bigstep<-0.2*bigwin
  #for (k in 1:length(getSequence(mySeq[1]))){
    #for (i in bigstart:bigend){
      window<-list()
      #start<-bigstart
      start<-1
      count<-0
      for (j in 1:4){
        #print(j)
        wingc<-list()
        wingc3<-list()
        wingene<-list()
        start<-start
        stop<-start+smallwinlen
        #print(start)
        #print(stop)
        midpoint<-median(c(start:stop))
        #Calculating parameter values for each window
        wingc<-output$gc[which(output$gcxpos %in% c(start:stop))]
        wingc3<-output$gc3[which(output$gc3xpos %in% c(start:stop))]
        wingene<-output$genebias[which(output$genexpos %in% c(start:stop))]
        #dividing each window to 10 to get n=10
        window$gcav[(count+1):(count+10)]<-partaverages(wingc)
        window$gc3av[(count+1):(count+10)]<-partaverages(wingc3)
        window$geneav[(count+1):(count+10)]<-partaverages(wingene)
        window$point[(count+1):(count+10)]<-paste('w',j)
        #print(window$point)
        count<-count+10
        #setting start for next small window
        start<-stop+1
      }
        #window$position<-midpoint
    #}
      #model<-manova(cbind(gcav,gc3av, geneav)~point, data=window)
      #model2<-manova(cbind(gc,gc3, gene)~point, data=window, subset=point %in% c('w 1', 'w 2'))
      #model3<-manova(cbind(gc,gc3, gene)~point, data=window, subset=point %in% c('w 2', 'w 3'))
      #model4<-manova(cbind(gc,gc3, gene)~point, data=window, subset=point %in% c('w 3', 'w 4'))
      #model5<-manova(cbind(gc,gc3, gene)~point, data=window, subset=point %in% c('w 1', 'w 4'))
      #summary<-summary(model, test='Wilks')
      #summary2<-summary(model2, test='Wilks')
      #summary3<-summary(model3, test='Wilks')
      #summary4<-summary(model4, test='Wilks')
      #summary5<-summary(model5, test='Wilks')
      #bigstart<-bigstart+bigstep
    #}
   # 
  #}
  #
  #
  #summary$stats...
  #print(summary)
#}
filltable<-function(mydata){
  myfill<-list()
  #filling shorter variable lists to same length in a list
  if(length(mydata)<repeat_len){
    fill<-vector('list',repeat_len-length(mydata))
    fill<-sapply(fill, function(x) sample(mydata,1))
    #print(fill)
    myfill<-unlist(c(mydata,fill), recursive=TRUE)
    #print(myfill)
    mydata<-myfill
    }
    return(mydata)
}
partaverages<-function(mydata){
  partav<-list()
  len<-round(length(mydata)/10)
  part_start<-1
  part_stop<-len
  #part<-0
  for (i in 1:10){
    #print(part_start)
    #print(len)
    partav[i]<-mean(mydata[part_start:part_stop])
    part_start<-part_start+len
    part_stop<-part_start+len-1
  }
  partav<-unlist(partav,recursive = TRUE)
  return(partav)
}