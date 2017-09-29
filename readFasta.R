
# Load Packages -----------------------------------------------------------

#source("http://bioconductor.org/biocLite.R")
#biocLite()
#install.packages("seqinr")
#require(seqinr)


# read seqinr test file ---------------------------------------------------

# dnafile <- system.file("sequences/malM.fasta", package = "seqinr")
# testseq <- read.fasta(file = dnafile, as.string = TRUE)
# print(testseq)


# read test file ----------------------------------------------------------

# dnafile <- file("test.fna")
# testseq <- read.fasta(file = dnafile, as.string = TRUE)
# print(testseq)


# read fasta file ---------------------------------------------------------

# dnafile <- file("FASTA/Bacillus subtilis txid1423/GCF_000009045.1_ASM904v1_genomic.fna")
# seq <- read.fasta(file = dnafile, as.string = TRUE)
# print(seq)



# list dirs and get file name ---------------------------------------------

readFasta <- function() {
  print(list.files("FASTA", all.files = TRUE, recursive = TRUE, pattern="*.faa"))
  print(list.files("FASTA", all.files = TRUE, recursive = TRUE, pattern="*.fna"))
  name <- readline("Give FASTA file name: ")
  dna_file <- paste("FASTA/", name, sep="")
  print(dna_file)
  dna_temp <- read.fasta(dna_file)
  
#  dna_seq <- getSequence(dna_temp[[1]])
#  print(dna_seq)
#  return(dna_seq)
  
  # Return the frame:
  return(dna_temp)
}
