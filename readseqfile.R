#Summary:
#Import a FASTA file to be analyzed
#Input arguments:
#dna_file = "local file path" (UNIX standard)
#verb = T/F for progress reports and plotting

# Load Packages -----------------------------------------------------------

#source("http://bioconductor.org/biocLite.R")
#biocLite()
#install.packages("seqinr")
#require(seqinr)

# Read by file name ---------------------------------------------
readseqfile <- function(dna_file, verb) {
  if (verb){cat("Reading", dna_file, "sequence...")}
  dna_seq <- read.fasta(dna_file)
  if (verb){cat("done!\n")}
  cat("Sequence is", getName(dna_seq), "\n")

  # Return the frame:
  return(dna_seq)
}
