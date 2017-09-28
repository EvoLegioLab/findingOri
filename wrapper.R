#Just starting to create the wrapper. Need to change plotting in both GC skew and k-mer so that it is optional.
#Also need to link file fetching here, and then just input that to sources
findingOri<-function(access, fi, verbose){
  if (missing(access)){
    cat ('calling', fi)
  }else if(missing(fi)){
    cat ('reading in fasta from accession number', access)
    access<-access
  }
  if (missing(verbose)){
    print ('figures will be produced!')
    verb<-'n'
  } else {
    print('no figures will be produced!')
    verb<-'y'
  }
  #if(!exists(c("getfastafiles",'countgc', 'threes', 'simplegenecount', 'cumgenecount'), mode="function")){
  source('GC_skew.R')
  #mySeq<-getfastafiles(access)
  #print(verb)
  GCall(access, ,verb)
 # }else{
  #  print('R code not available')
  #}
#-----------------Output from scripts and modelling----------------
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
    for (i in bigstart:bigend){
      window<-list()
      #start<-bigstart
      start<-1900000
      counter<-0
      for (j in 1:4){
        #print(j)
        wingc<-list()
        wingc3<-list()
        wingene<-list()
        filledwingc<-list()
        filledwingc3<-list()
        filledwingene<-list()
        start<-start
        stop<-start+smallwinlen
        #print(start)
        #print(stop)
        midpoint<-median(c(start:stop))
        #Calculating parameter values for each window
        #print(which(output$gcxpos %in% c(start:stop)))
        wingc<-output$gc[which(output$gcxpos %in% c(start:stop))]
        lengc<-length(wingc)
        #print(lengc)
        wingc3<-output$gc3[which(output$gc3xpos %in% c(start:stop))]
        lengc3<-length(wingc3)
        #print(lengc3)
        wingene<-output$genebias[which(output$genexpos %in% c(start:stop))]
        lengene<-length(wingene)
        repeat_len<-max(c(lengc, lengc3, lengene))
        #filling list to same lengths with filltable function (defined below)
        filledwingc<-filltable(wingc)
        filledwingc3<-filltable(wingc3)
        filledwingene<-filltable(wingene)
        #within each window
        for (n in 1:repeat_len){
          #print(counter+n)
          #print(repeat_len)
          window$gc[(counter+n)]<-filledwingc[n]
          window$gc3[(counter+n)]<-filledwingc3[n]
          #window$kmer[j]<-
          window$gene[(counter+n)]<-filledwingene[n]
          window$point[(counter+n)]<-paste('w',j)
          window$position[(counter+n)]<-midpoint
        }
        counter<-counter+repeat_len
        #print(counter)
        #setting start for next small window
        start<-stop+1
        #print(start)
      }
      model<-manova(cbind(gc,gc3, gene)~point, data=window)
      model2<-manova(cbind(gc,gc3, gene)~point, data=window, subset=point %in% c('w 1', 'w 2'))
      model3<-manova(cbind(gc,gc3, gene)~point, data=window, subset=point %in% c('w 2', 'w 3'))
      model4<-manova(cbind(gc,gc3, gene)~point, data=window, subset=point %in% c('w 3', 'w 4'))
      model5<-manova(cbind(gc,gc3, gene)~point, data=window, subset=point %in% c('w 1', 'w 4'))
      summary<-summary(model, test='Wilks')
      summary2<-summary(model2, test='Wilks')
      summary3<-summary(model3, test='Wilks')
      summary4<-summary(model4, test='Wilks')
      summary5<-summary(model5, test='Wilks')
      #bigstart<-bigstart+bigstep
    }
   # 
  #}
  #
  #
  #summary$stats...
  print(summary)
}
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
