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
  cA<-list()
  cT<-list()
  x<-list()
  G_C<-list()
  T_A<-list()
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
    #If window length does not go over the last nucleotides
    if (end <= length(mydna)){
      tabseq<-table(mydna[start:end])
      cC[i]<-tabseq[['c']]
      cG[i]<-tabseq[['g']]
      cA[i]<-tabseq[['a']]
      cT[i]<-tabseq[['t']]
      G_C[i]<-tabseq[['g']]-tabseq[['c']]
      T_A[i]<-tabseq[['t']]-tabseq[['a']]
      x[i]<-xpos[median(c(start:end))]
      start<-start+step-1
      end<-start+windowlen
    } else {
      #When we are at the last window (and it goes over the last nucleotide) we allow it to go to the beginning again
      tabseq<-table(c(mydna[start:length(mydna)],seqfile[1:(end-length(mydna))]))
      cC[i]<-tabseq[['c']]
      cG[i]<-tabseq[['g']]
      cA[i]<-tabseq[['a']]
      cT[i]<-tabseq[['t']]
      G_C[i]<-tabseq[['g']]-tabseq[['c']]
      T_A[i]<-tabseq[['t']]-tabseq[['a']]
      x[i]<-xpos[median(c(start:end))]
      start<-start+step-1
      end<-start+windowlen
    }
  }
  gccount<-list(cC=unlist(cC,recursive=TRUE), cG=unlist(cG, recursive=TRUE), cT=unlist(cT, recursive=TRUE), cA=unlist(cA, recursive=TRUE),gc=unlist(G_C,recursive=TRUE),ta=unlist(T_A, recursive=TRUE), xpos=unlist(x, recursive = TRUE))
}
#--------------------------------------GC3 positions and figure------------------------------------
#This function first checks that the sequence file and annotation file annotations are compatible. Then for both strands, it goes through the codons
#takes every third nucleotide and respective positions from the created sequence position list
threes<-function(mydata, myannot){
  positcod<-list()
  seqpos<-list()
  codonseqpos<-list()
  #This function gets the third codons and their positions for the strand marked (+1)(or -1) in prodigal
  #it returns the 3rd nucleotides and their positions for both strands
  for (i in 1:length(mydata)){
    if (names(mydata[i])==names(myannot[i])){
      mydna<-getSequence(mydata[[i]])
      seqpos<-c(myannot[[i]][1]:myannot[[i]][2])
      #getting third nucleotides from sequence on the +1 strand
      if (myannot[[i]][3]==1){
        positcod<-c(positcod,mydna[seq(3,length(mydna),3)])
        #Setting the positions of the nucleotides to another list
        codonseqpos<-c(codonseqpos,seqpos[seq(3,length(seqpos),3)])
        #}
      } else {
        #Getting third nucleotides from sequence on the -1 strand
        #taking the complement of the nucleotides to work with one strand only
        mydna=comp(mydna)
        positcod<-c(positcod,mydna[seq(3,length(mydna),3)])
        codonseqpos<-c(codonseqpos,seqpos[seq(3,length(seqpos),3)])
      }
      #If names do not match
      } else {  
      print('Datasets do not match!')
    }
  }
  positcod<-unlist(positcod,recursive=TRUE)
  codonseqpos<-unlist(codonseqpos, recursive=TRUE)
  thirds<-list(plusncl=positcod,pluspos=codonseqpos)
  return(thirds)
}
##-----------------RUnning with input--------------
#GCskews<-function(accession){
#change your accession code here
# accession<-"NZ_CM000488.1"
# mySeq<-getfastafiles(accession)

#fetch FASTA file
source("readFasta.R")
myDNAbin<-readFasta()
myFasta<-write.dna(myDNAbin,file ="myFasta.fasta", format = "fasta")
mySeq<-read.fasta("myFasta.fasta")

#NOTE! Fasta used in the following calls MUST be in the variable mySeq!
#------------GC-----------------
#Running gc counts on whole sequence
gccount<-countgc(mySeq)
gccount$cumgc<-cumsum(gccount$gc)
gccount$cumta<-cumsum(gccount$ta)
plot(gccount$xpos, gccount$cumgc, type='l')
plot(gccount$xpos, gccount$cumta, type='l')
#-------------gc3-------------
#Running gene prediction with prodigal and gc3 counts for 3rd codon positions
system('prodigal -i myFasta.fasta -o myGenCoord.fasta -d myGenSeqs.fasta')
mygenseq<-read.fasta('myGenSeqs.fasta')
myAnnot<-sapply(mygenseq,function(x) strsplit(getAnnot(x),'#'))
annotinfo<- lapply(lapply(myAnnot,'[', 1:4),function(x) as.integer(x[2:4]))
thirds<-threes(mygenseq,annotinfo)
gc3s<-countgc(thirds$plusncl, thirds$pluspos)
gc3s$cumgc<-cumsum(gc3s$gc)
gc3s$cumta<-cumsum(gc3s$ta)
plot(gc3s$xpos, gc3s$cumgc, type='l')
plot(gc3s$xpos,gc3s$cumta, type='l')
ter<-gc3s$xpos[match(max(gc3s$cumgc),gc3s$cumgc)]
ori<-gc3s$xpos[match(min(gc3s$cumgc),gc3s$cumgc)]
cat('The origin is located at', ori ,' and the terminus at', ter)

#---Counting genes------------------
#Functions
simplegenecount<-function(annotations){
  lespos<-list()
  laspos<-list()
  lescount<-0
  lascount<-0
  if(ori<ter){
    lespos<-c(ori:ter)
    laspos<-c(ter:length(mySeq), 1:ori)
  }else{
    lespos<-c(ori:length(mySeq), 1:ter)
    laspos<-c(ter:ori)
  }
  for (i in 1:length(annotations)){
    if(annotations[[i]][2] %in% lespos){
      lescount<-lescount+1
    } else{
      lascount<-lascount+1
    }
  }
  genecounts<-list(lesc=lescount, lasc=lascount)
}
#------------Running gene counts-----------------
genecount<-simplegenecount(annotinfo)
genecount
genebias<-(genecount$lesc/(genecount$lesc+genecount$lasc))*100
cat('\n', genebias, '% of genes are on the leading strand')
#more elaborated try with windows and plotting

return(gc3s)
