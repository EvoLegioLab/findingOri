library(seqinr)
library(ape)
#Calculating GC skew from whole genome and genomic 3rd codon positions. 
#NOTE! Prodigal MUST be installed on the computer and working through bash!
#-----------------------------Import data-------------------------------------
#Importing data from binary DNA file or with an accession number
getfastafiles<-function(accession){
  myDNAbin<-read.GenBank("accession")
  myFasta<-write.dna(myDNAbin,file ="myFasta.fasta", format = "fasta")
  mySeq<-read.fasta("myFasta.fasta")
}
getDNAbin<-function(mybinfile){
  myFasta<-write.dna(mybinfile,file ="myFasta.fasta", format = "fasta")
  mySeq<-read.fasta("myFasta.fasta")
}
#-----------------------GC calculations----------------------------------------
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
gccount<-countgc(mySeq)
gccount$cumgc<-cumsum(gccount$gc)
plot(gccount$xpos, gccount$cumgc, type='l')
#--------------------------------------GC3 positions and figure------------------------------------
system('prodigal -i myFasta.fasta -o myGenCoord.fasta -d myGenSeqs.fasta')
mygenseq<-read.fasta('myGenSeqs.fasta')
myannot<-sapply(mygenseq,function(x) strsplit(getAnnot(x),'#'))
annotinfo<- lapply(lapply(myannot,'[', 1:4),function(x) as.integer(x[2:4]))
positcod<-list()
negcod<-list()
seqpos<-list()
codonseqpos<-list()
seqpos_neg<-list()
codonseqpos_neg<-list()
threes<-function(mydata, myannot){
  #This function gets the third codons and their positions for the strand marked (+1)(or -1) in prodigal
  #it returns the strand vectors
  for (i in 1:length(mydata)){
    if (names(mydata[i])==names(myannot[i])){
      mydna<-getSequence(mydata[[i]])
      #getting third nucleotides from sequence on the +1 strand
      if (myannot[[i]][3]==1){
        #print('plus strand')
        #positcod list creation and concatenation works, tried with two elements
        positcod<-c(positcod,mydna[seq(3,length(mydna),3)])
        #Setting the positions of the nucleotides to another list
        seqpos<-c(myannot[[i]][1]:myannot[[i]][2])
        codonseqpos<-c(codonseqpos,seqpos[seq(3,length(seqpos),3)])
      } else {
        #print('minus strand')
        #Getting third nucleotides from sequence on the -1 strand
        negcod<-c(negcod,mydna[seq(3,length(mydna),3)])
        seqpos_neg<-c(myannot[[i]][1]:myannot[[i]][2])
        codonseqpos_neg<-c(codonseqpos_neg,seqpos_neg[seq(3,length(seqpos_neg),3)])
        #print(unlist(codonseqpos_neg, recursive=TRUE)) #works!
      }
    } else {  #If names do not match
      print('Datasets do not match!')
    }
  }
  positcod<-unlist(positcod,recursive=TRUE)
  codonseqpos<-unlist(codonseqpos, recursive=TRUE)
  negcod<-unlist(negcod, recursive = TRUE)
  codonseqpos_neg<-unlist(codonseqpos_neg, recursive=TRUE)
  thirds<-list(plusncl=positcod,pluspos=codonseqpos,minusncl=negcod,minuspos=codonseqpos_neg)
  #thirds<-list(plusncl=positcod,pluspos=codonseqpos)
  return(thirds)
}
thirds<-threes(mygenseq,annotinfo)
gc3s<-countgc(thirds$plusncl, thirds$pluspos)
gc3s$cumgc<-cumsum(gc3s$gc)
plot(gc3s$xpos, gc3s$cumgc, type='l')