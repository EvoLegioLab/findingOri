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
  print(verb)
  GCall(access,,verb)
 # }else{
  #  print('R code not available')
  #}
  
  #output<-list(gc=...gc3=...gcxpos=...gc3xpos=...kmer=...kmerxpos...genebias=...genexpos=...)
  #for (each big window){
  #output$window<- set a bunch of sequence neighbourhoods to test, e.g. len=big window/10
  #for (each small window){
  #calculate variable averages per small window, minimum 1 observation/variable
  #}
  #}
  #model<-manova(cbind(gc, gc3,kmer,genebias)~window, data=output)
  #summary<-summary(model, test=c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"))
  #summary$stats...
}