require(seqinr)
require(ape)
#Calculating GC skew from whole genome and genomic 3rd codon positions. 
#NOTE! Prodigal MUST be installed on the computer and working through bash!
#-----------------------------Import data-------------------------------------
#Importing data from binary DNA file or with an accession number
getfastafiles<-function(accession){
  myDNAbin<-read.GenBank(accession)
  myFasta<-write.dna(myDNAbin,file ="myFasta.fasta", format = "fasta")
  mySeq<-read.fasta("myFasta.fasta")
  return(mySeq)
}
getDNAbin<-function(mybinfile){
  myFasta<-write.dna(mybinfile,file ="myFasta.fasta", format = "fasta")
  mySeq<-read.fasta("myFasta.fasta")
  return(mySeq)
}
#-----------------------GC calculations (functions)----------------------------------------
#Counting numbers of g and c per window
#Input to countgc either a fasta file (for whole genome gc), pure sequence list or fasta list and sequence position list (for gc3) 
countgc<-function(seqfile, seqpos){
  if (length(seqfile)==1){
    mydna<-getSequence(seqfile[[1]])
  }else{
    mydna<-getSequence(seqfile)
  }
  cC<-list()
  cG<-list()
  x<-list()
  G_C<-list()
  windowlen<-1000
  step<-100
  start<-1
  end<-start+windowlen
  if(missing(seqpos)){
    xpos<-c(1:length(mydna))
  } else {
    xpos<-seqpos
  }
  for (i in 1:(length(mydna)/step)){
    if (end <= length(mydna)){
      tabseq<-table(mydna[start:end])
      cC[i]<-tabseq[['c']]
      cG[i]<-tabseq[['g']]
      G_C[i]<-tabseq[['g']]-tabseq[['c']]
      x[i]<-xpos[median(c(start:end))]
      start<-start+step-1
      end<-start+windowlen
    } else {
      tabseq<-table(c(mydna[start:length(mydna)],seqfile[1:(end-length(mydna))]))
      cC[i]<-tabseq[['c']]
      cG[i]<-tabseq[['g']]
      G_C[i]<-tabseq[['g']]-tabseq[['c']]
      x[i]<-xpos[median(c(start:end))]
      start<-start+step-1
      end<-start+windowlen
    }
  }
  gccount<-list(cC=unlist(cC,recursive=TRUE), cG=unlist(cG, recursive=TRUE), gc=unlist(G_C,recursive=TRUE),xpos=unlist(x, recursive = TRUE))
}
#--------------------------------------GC3 positions and figure------------------------------------
threes<-function(mydata, myannot){
  tempposit<-list()
  tempcod<-list()
  positcod<-list()
  tempnegcod<-list()
  tempcodonseqpos_neg<-list()
  negcod<-list()
  seqpos<-list()
  codonseqpos<-list()
  codonseqpos_neg<-list()
  #This function gets the third codons and their positions for the strand marked (+1)(or -1) in prodigal
  #it returns the 3rd nucleotides and their positions for both strands
  for (i in 1:length(mydata)){
    if (names(mydata[i])==names(myannot[i])){
      mydna<-getSequence(mydata[[i]])
      seqpos<-c(myannot[[i]][1]:myannot[[i]][2])
      #print(length(seqpos))
      #getting third nucleotides from sequence on the +1 strand
      if (myannot[[i]][3]==1){
        #with the foor loop code takes forever and has not looked right
        #for (j in 1:(length(mydna)/3)){
         # positcod<-c(positcod,mydna[3*j]) 
        #Setting the positions of the nucleotides to another list
          #codonseqpos<-c(codonseqpos,seqpos[3*j])
          #print(length(positcod));print(length(codonseqpos))
        #with this code has worked before
        positcod<-c(positcod,mydna[seq(3,length(mydna),3)])
        codonseqpos<-c(codonseqpos,seqpos[seq(3,length(seqpos),3)])
        #}
      } else {
        #Getting third nucleotides from sequence on the -1 strand
        #for (j in 1:(length(mydna)/3)){
        #  negcod<-c(negcod,mydna[3*j])
        #  codonseqpos_neg<-c(codonseqpos_neg,seqpos[3*j])
        negcod<-c(negcod,mydna[seq(3,length(mydna),3)])
        codonseqpos_neg<-c(codonseqpos_neg,seqpos[seq(3,length(seqpos),3)])
        #}
      }
      #If names do not match
      } else {  
      print('Datasets do not match!')
    }
  }
  #print(head(positcod))
  positcod<-unlist(positcod,recursive=TRUE)
  codonseqpos<-unlist(codonseqpos, recursive=TRUE)
  negcod<-unlist(negcod, recursive = TRUE)
  codonseqpos_neg<-unlist(codonseqpos_neg, recursive=TRUE)
  thirds<-list(plusncl=positcod,pluspos=codonseqpos,minusncl=negcod)
  #thirds<-list(plusncl=positcod,pluspos=codonseqpos)
  return(thirds)
}
##-----------------RUnning with input--------------
#GCskews<-function(accession){
#change your accession code here
accession<-"CP011051.1"
mySeq<-getfastafiles(accession)
#------------GC-----------------
gccount<-countgc(mySeq)
gccount$cumgc<-cumsum(gccount$gc)
plot(gccount$xpos, gccount$cumgc, type='l')
#-------------gc3-------------
system('prodigal -i myFasta.fasta -o myGenCoord.fasta -d myGenSeqs.fasta')
mygenseq<-read.fasta('myGenSeqs.fasta')
myAnnot<-sapply(mygenseq,function(x) strsplit(getAnnot(x),'#'))
annotinfo<- lapply(lapply(myAnnot,'[', 1:4),function(x) as.integer(x[2:4]))
thirds<-threes(mygenseq,annotinfo)
gc3s<-countgc(thirds$plusncl, thirds$pluspos)
gc3s$cumgc<-cumsum(gc3s$gc)
plot(gc3s$xpos, gc3s$cumgc, type='l')
#}