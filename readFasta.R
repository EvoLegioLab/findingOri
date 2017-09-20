
# Load Packages -----------------------------------------------------------

#source("http://bioconductor.org/biocLite.R")
#biocLite()
#install.packages("seqinr")
require(seqinr)


# read seqinr test file ---------------------------------------------------

# dnafile <- system.file("sequences/malM.fasta", package = "seqinr")
# testseq <- read.fasta(file = dnafile, as.string = TRUE, forceDNAtolower = FALSE)
# print(testseq)


# read test file ----------------------------------------------------------

# dnafile <- file("test.fna")
# testseq <- read.fasta(file = dnafile, as.string = TRUE, forceDNAtolower = FALSE)
# print(testseq)


# read fasta file ---------------------------------------------------------

# dnafile <- file("FASTA/Bacillus subtilis txid1423/GCF_000009045.1_ASM904v1_genomic.fna")
# seq <- read.fasta(file = dnafile, as.string = TRUE, forceDNAtolower = FALSE)
# print(seq)



# list dirs and get file name ---------------------------------------------

readFasta <- function() {
  print(list.files("FASTA", all.files = TRUE, recursive = TRUE, pattern="*.faa"))
  print(list.files("FASTA", all.files = TRUE, recursive = TRUE, pattern="*.fna"))
  name <- readline("Give FASTA file name: ")
  dnafile <- paste("FASTA/", name, sep="")
  print(dnafile)
  seqstring <- read.fasta(file = dnafile,  forceDNAtolower = FALSE, seqonly = TRUE)
  seq <- strsplit(noquote(seqstring[[1]]), NULL)
  # print(seq)
  return(seq)
}
