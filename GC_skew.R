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
GCall<-function(mysequence,verb){
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
  for (i in 1:ceiling((length(mydna)/step))){
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
      tabseq<-table(c(mydna[start:length(mydna)],mydna[1:(end-length(mydna))]))
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
mySeq<-mysequence
#print(mySeq)
#NOTE! Fasta used in the following calls MUST be in the variable mySeq!
#------------GC-----------------
#Running gc counts on whole sequence
if (verb){cat("Running raw nucleotide calculations...")}
gccount<-countgc(mySeq)
gccount$cumgc<-cumsum(gccount$gc)
gccount$cumta<-cumsum(gccount$ta)
if (verb){
  plot(gccount$xpos, gccount$cumgc, type='l', main='Cumulative GC-skew based whole sequence', xlab='Position on chromosome', ylab='Cumulative GC')
  plot(gccount$xpos, gccount$cumta, type='l', main='Cumulative TA-skew based on whole sequence',xlab='Position on chromosome',ylab='Cumulative TA')
}
ter<-gccount$xpos[match(max(gccount$cumgc),gccount$cumgc)]
ori<-gccount$xpos[match(min(gccount$cumgc),gccount$cumgc)]
tamin<-gccount$xpos[match(min(gccount$cumta),gccount$cumta)]
tamax<-gccount$xpos[match(max(gccount$cumta),gccount$cumta)]
if (ter>ori){
  if(tamin>tamax){
    tater<-tamin
    taori<-tamax
  }else{
    ta3ter<-ta3max
    ta3ori<-ta3min
  }
}else{
  if(tamin>tamax){
    tater<-tamax
    taori<-tamin
  }else{
    tater<-tamim
    taori<-tamax
  }
}
if (verb){cat("done!\n")}
cat('Based on gc, the origin is located at', ori ,' and the terminus at', ter,'\n')
#-------------gc3-------------
#Running gene prediction with prodigal and gc3 counts for 3rd codon positions
if (verb){cat("Running gene prediction...\n")
system('prodigal -i myFasta.fasta -o myGenCoord.fasta -d myGenSeqs.fasta')}
system('prodigal -q -i myFasta.fasta -o myGenCoord.fasta -d myGenSeqs.fasta')
if (verb){cat("...gene prediction done!\n")
cat("Running 3rd codon position calculations...\n")}
mygenseq<-read.fasta('myGenSeqs.fasta')
myAnnot<-sapply(mygenseq,function(x) strsplit(getAnnot(x),'#'))
annotinfo<- lapply(lapply(myAnnot,'[', 1:4),function(x) as.integer(x[2:4]))
thirds<-threes(mygenseq,annotinfo)
gc3s<-countgc(thirds$plusncl, thirds$pluspos)
gc3s$cumgc<-cumsum(gc3s$gc)
gc3s$cumta<-cumsum(gc3s$ta)
if (verb){
  plot(gc3s$xpos, gc3s$cumgc, type='l', main='Cumulative GC-skew based on 3rd codon position',xlab='Position on chromosome', ylab='Cumulative GC')
  plot(gc3s$xpos,gc3s$cumta, type='l', main='Cumulative TA-skew based on 3rd codon position',xlab='Position on chromosome',ylab='Cumulative TA')
}
gc3ter<-gc3s$xpos[match(max(gc3s$cumgc),gc3s$cumgc)]
gc3ori<-gc3s$xpos[match(min(gc3s$cumgc),gc3s$cumgc)]
ta3min<-gc3s$xpos[match(min(gc3s$cumta),gc3s$cumta)]
ta3max<-gc3s$xpos[match(max(gc3s$cumta),gc3s$cumta)]
if (gc3ter>gc3ori){
  if(ta3min>ta3max){
    ta3ter<-ta3min
    ta3ori<-ta3max
  }else{
    ta3ter<-ta3max
    ta3ori<-ta3min
  }
}else{
  if(ta3min>ta3max){
    ta3ter<-ta3max
    ta3ori<-ta3min
  }else{
    ta3ter<-ta3mim
    ta3ori<-ta3max
  }
}
if (verb){cat("...3rd codon position calculations done!\n")}
cat('Based on gc3, the origin is located at', gc3ori ,' and the terminus at', gc3ter, "\n")

#---Counting genes------------------
#Functions
simplegenecount<-function(annotations){
  lespos<-list()
  laspos<-list()
  lescount<-0
  lascount<-0
  if(ori<ter){
    lespos<-c(gc3ori:gc3ter)
    laspos<-c(gc3ter:length(mySeq), 1:gc3ori)
    #assuming prodigal +1 is the lagging strand when ori comes before ter
    for (i in 1:length(annotations)){
      if((annotations[[i]][2] %in% lespos)&&(annotations[[i]][3]==-1)|(annotations[[i]][2] %in% laspos)&&(annotations[[i]][3]==1)){
        lescount<-lescount+1
      } else{
        lascount<-lascount+1
      }
    }
  }else{
    lespos<-c(gc3ori:length(mySeq), 1:gc3ter)
    laspos<-c(gc3ter:gc3ori)
    #assuming prodigal +1 is the leading strand when ter comes before ori
    for (i in 1:length(annotations)){
      if((annotations[[i]][2] %in% lespos)&&(annotations[[i]][3]==1)|(annotations[[i]][2] %in% laspos)&&(annotations[[i]][3]==-1)){
        lescount<-lescount+1
      } else{
        lascount<-lascount+1
      }
    }
  }
  genecounts<-list(lesc=lescount, lasc=lascount)
}

#------------Running gene counts-----------------
if (verb){cat("Running gene bias calculations...\n")}
genecount<-simplegenecount(annotinfo)
#genecount
genebias<-(genecount$lesc/(genecount$lesc+genecount$lasc))*100
cat(genebias, '% of genes are on the leading strand\n')

#more elaborated gene bias cumulative calculation with windows and plotting
cumgenecount<-function(myannot){
  xpos<-list()
  ypos<-list()
  cumgene<-list()
  for (i in 1:length(myannot)){
    xpos[i]<-median(c(myannot[[i]][1]:myannot[[i]][2]))
    ypos[i]<-myannot[[i]][3]
  }
  xpos<-unlist(xpos,recursive=TRUE)
  ypos<-unlist(ypos,recursive=TRUE)
  #print(head(xpos))
  #print(head(ypos))
  cumgene<-cumsum(ypos)
  cumgenebias<-list(cumgene=cumgene, genes=ypos,xpos=xpos)
}

cumgenebias<-cumgenecount(annotinfo)
geneter<-cumgenebias$xpos[match(max(cumgenebias$cumgene),cumgenebias$cumgene)]
geneori<-cumgenebias$xpos[match(min(cumgenebias$cumgene), cumgenebias$cumgene)]
if (verb){cat("...gene bias calculations done!\n")}
cat('Based on gene bias, the origin is located at', geneori ,' and the terminus at', geneter, '\n')
if (verb){
  plot(cumgenebias$xpos, cumgenebias$cumgene, type='l',main='Cumulative gene bias over the whole chromosome',xlab='Position on chromosome',ylab='Cumulative gene count')
}
rawoutput<-list(gc=gccount$gc,gc3=gc3s$gc,ta=gccount$ta,ta3=gc3s$ta,cumgc=gccount$cumgc,cumgc3=gc3s$cumgc, cumta=gccount$cumta, cumta3=gc3s$cumta,gcxpos=gccount$xpos,gc3xpos=gc3s$xpos,genebias=cumgenebias$genes,genexpos=cumgenebias$xpos,gcori=ori, gc3ori=gc3ori, gcter=ter, gc3ter=gc3ter, taori=taori, tater=tater, geneori=geneori, geneter=geneter)
#print(head(rawoutput))
return(rawoutput)
}
