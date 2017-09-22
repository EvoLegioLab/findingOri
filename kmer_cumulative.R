#Summary:
#Import a FASTA file to be analyzed for kmer bias and GC skew
#First go through the sequence once using smaller, non-overlapping windows to find interesting kmers
#From each window, select the top x kmers with highest difference between observed counts and expected counts based on nucleotide composition
#Go through the chromosome again, using larger, overlapping windows. Count the occurrences of the selected kmers and the complementary kmers.
#The counts of complementary kmers is equal to the counts of the selected kmers on the complementary strand
#Plot the (kmer count) - (complementary kmer count) as a cumulative graph. Also calculate and plot the cumulative GC skew.


#Required packages: seqinr, stringr, ape
require("seqinr")
require("stringr")
require("ape")

#Download sequence from database. Function from Little book of R for bioinf
getncbiseq <- function(accession) 
{
  # first find which ACNUC database the accession is stored in:
  dbs <- c("genbank","refseq","refseqViruses","bacterial")
  numdbs <- length(dbs)
  for (i in 1:numdbs)
  {
    db <- dbs[i]
    choosebank(db)
    # check if the sequence is in ACNUC database 'db':
    resquery <- try(query(".tmpquery", paste("AC=", accession)), silent = TRUE)
    if (!(inherits(resquery, "try-error")))
    {
      query2 <- query("query2", paste("AC=",accession,sep=""))
      # see if a sequence was retrieved:
      seq <- getSequence(query2$req[[1]])
      closebank()
      return(seq)
    }
    closebank()
  }
  print(paste("ERROR: accession",accession,"was not found"))
}

#Get the expected number of kmers based on nucleotide composition in a certain chunk
get_expected <- function(kmer, chunk_comp1, window_factor) {
  #kmer = kmer to be analyzed (string)
  #chunk_comp1 = count of each base in the chunk to be analyzed (integer vector with names)
  #window_factor = total number of kmers in the chunk (integer)
  #Returns: exp_count = expected count of the input kmer in the chunk
  #Count of each nucleotide:
  A_w = chunk_comp1[1]
  C_w = chunk_comp1[2]
  G_w = chunk_comp1[3]
  T_w = chunk_comp1[4]
  all_w = A_w + C_w + G_w + T_w
  #Frequency of each nucleotide:
  ACGT_frac = c(A_w/all_w, C_w/all_w, G_w/all_w, T_w/all_w)
  #Translate kmer to vector of indeces for the vector ACGT_frac:
  kmer_split <- unlist(strsplit(kmer,""))
  kmer_transl <- as.integer(chartr("acgt","1234",kmer_split)) 
  exp_count <- window_factor
  #For each base in the kmer, multiply exp_count with frequency of that base:
  for (i in 1:length(kmer_transl)) { 
    j = kmer_transl[i]
    exp_count <- exp_count*ACGT_frac[j]
  }
  return(exp_count)
}

#Convert a vector of kmers to a vector of the complementary, inverted kmers
get_compl_rev <- function(kmer_vec) {
  num_kmers <- length(kmer_vec)
  for (i in 1:num_kmers) {
    kmer_vec[i] <- chartr("acgt","tgca",kmer_vec[i])
    kmer_temp <- rev(unlist(strsplit(kmer_vec[i],"")))
    kmer_vec[i] <- paste(kmer_temp, collapse = "", sep = "")
  }
  return(kmer_vec)
}

#Calculate and plot metrics in sliding (overlapping) windows:
slidingwindowplot <- function(windowsize, find_kmers_wind, step_len, inputseq, accession, oligo, top_kmers)
{
  #windowsize = length of window to be used in the main analysis
  #find_kmers_wind = length of window when identifying most frequent kmers
  #step_len = step length between consecutive overlapping windows in main analysis
  #inputseq = chromosomal sequence to be analyzed
  #accession = name of the sequence, to be printed in the title of the plot
  #oligo = length of kmer
  #top_kmers = How many of the most frequent kmers in each window should be included?
  seq_length <- length(inputseq)
  oct_per_wind <- find_kmers_wind-oligo+1 #number of kmers in each window in the first step
  starts <- seq(1, seq_length, by = find_kmers_wind) #For stepping through the chromosome the first time (no overlaps between windows)
  n <- length(starts) 
  kmer_vec <- numeric(n*top_kmers) #Store selected kmers in this vector
  next_row = 1 #Counting the number of unique kmers selected below
  #Go through the chromosome to select kmers (windows not overlapping):
  for (i in 1:n) {
    #Define each chunk (window):
    start <- starts[i]
    end <- start + find_kmers_wind-1
    if (end <= seq_length){
      chunk <- inputseq[start:end]
    }
    else { #Last chunk includes the beginning of the sequence also to become as long as the other chunks
      chunk <- inputseq[start:seq_length]
      chunk <- c(chunk,inputseq[1:(end-seq_length)])
    }
    chunk_comp1 <- count(chunk,1) #Count number of each nucleotide in the window
    chunk_comp8 <- count(chunk,oligo) #Count all kmers in the window
    chunk_comp8_len <- length(chunk_comp8)
    #Calculate the expected number of each kmer based on nucleotide composition:
    expected_comp8 <- numeric(length(chunk_comp8))
    for (ii in 1:chunk_comp8_len) {
      expected_comp8[ii] <- get_expected(names(chunk_comp8[ii]), chunk_comp1, oct_per_wind)
    }
    #Take the difference between observed and expected counts of each kmer:
    diff_comp8 <- chunk_comp8 - expected_comp8 
    #Sort kmers to get the kmers with largest difference between observed and expected frequency first:
    sorted_kmers <- sort(diff_comp8, decreasing=TRUE, na.last=TRUE) 
    #Save the first top_kmers of the sorted kmers to kmer_vec
    for (j in 1:top_kmers) { 
      kmer <- names(sorted_kmers[j]) 
      #Only add kmer if it and its complement don't already exist in the vector
      if (! (kmer %in% kmer_vec)) { #((! (kmer %in% kmer_vec)) & (! (chartr("acgt","tgca",kmer) %in% kmer_vec))) { 
        kmer_vec[next_row] <- kmer
        next_row = next_row + 1 
      }
    }
  }
  tot <- next_row-1 #Total number of unique kmers found above
  kmers_found <- sort(kmer_vec[1:tot]) #Put kmers in same order as they come from the count() function
  kmers_compl_rev <- get_compl_rev(kmers_found) #Get the complementary, inverted kmers
  overlapping = tot/(n*top_kmers) #To what extent did the top kmers from different windows overlap each other?
  cat("(Unique selected kmers)/(All selected kmers): ", overlapping, "\n")
  cat("Selected kmers: \n")
  print(kmers_found)
  #Sort the complementary kmers so that they come in the same orders as in the output from the count() function
  #But before sorting, save information about how it will be sorted in compl_index, so that the complementary kmers
  #can be rearranged to the initial order again later to match the order in the kmers_found vector
  compl_index <- sort.list(kmers_compl_rev) #Get indeces for sorting the complementary kmers
  meta_index <- sort.list(compl_index) #Get indeces for sorting it back to the initial order
  sorted_compl_rev <- sort(kmers_compl_rev) #Actually sort the complementary kmers
  
  #####Now use larger windows, overlapping each other######
  oct_per_wind <- windowsize-oligo+1 #Update number of kmers per window
  starts <- seq(1, seq_length, by = step_len) #Update starts vector for overlapping windows
  n <- length(starts) 
  kmer_matrix <- matrix(0L,nrow=tot,ncol=n) #Matrix for cumulative kmer count
  cum_gc_vec <- numeric(n) #Vector for cumulative GC skew
  #Start by looking only at the first window, 
  #to get the indices of the selected kmers in the output from the count() function
  chunk <- inputseq[1:windowsize-1] 
  chunk_comp8 <- count(chunk,oligo)
  chunk_comp1 <- count(chunk,1)
  C_w = chunk_comp1[2]
  G_w = chunk_comp1[3]
  cum_gc_vec[1] <- (G_w - C_w)/(G_w + C_w) #Cumulative GC skew
  k = 1
  h = 1
  index_vec = numeric(tot)
  index_vec_compl = numeric(tot)
  for (j in 1:chunk_comp8_len){ 
    #Go through all kmers from count() to search for the selected kmers and complementary kmers
    #Since the kmers come in the same order in all three cases, the search becomes easier
    #(The kmers_found and sorted_compl_rev are sorted alphabetically and the output from count() comes in alphabetical order)
    if (k <= tot & (names(chunk_comp8[j])==kmers_found[k])) {
      #For each kmer, save the difference between the count of this kmer and the count of the complementary kmer
      kmer_matrix[k,1] <- kmer_matrix[k,1] + chunk_comp8[j] 
      index_vec[k] = j #Save the chunk_comp8-index of each kmer
      k = k+1
    }
    if (h <= tot & (names(chunk_comp8[j])==sorted_compl_rev[h])) {
      #For each kmer, save the difference between the count of this kmer and the count of the complementary kmer
      kmer_matrix[h,1] <- kmer_matrix[h,1] - chunk_comp8[j] 
      index_vec_compl[h] = j #Save the chunk_comp8-index of each complementary kmer
      h = h+1
    }
    if (k>tot & h>tot){ #When all selected kmers and their complements are found
      break
    }
  }
  #Rearrange the complementary kmers so that their order again matches the order of the selected kmers
  index_vec_compl <- index_vec_compl[meta_index] 
  
  #### Now look at all the remaining windows #####
  for (i in 2:n) { 
    start <- starts[i]
    end <- start + windowsize-1
    if (end <= seq_length){
      chunk <- inputseq[start:end]
    }
    else {
      chunk <- inputseq[start:seq_length]
      chunk <- c(chunk,inputseq[1:(end-seq_length)])
    }
    chunk_comp8 <- count(chunk,oligo)
    chunk_comp1 <- count(chunk,1)
    C_w = chunk_comp1[2]
    G_w = chunk_comp1[3]
    cum_gc_vec[i] <- (G_w - C_w)/(G_w + C_w) + cum_gc_vec[i-1] #Cumulative GC skew
    for (j in 1:tot) { #Go through all selected kmers
      jj <- index_vec[j] #Get index of the kmer in chunk_comp8
      jj_compl <- index_vec_compl[j] #Get index of the complementary kmer in chunk_comp8
      #For each kmer, add the difference between the count of this kmer and the count of the complementary kmer to the previous value (cumulatively)
      kmer_matrix[j,i] <- kmer_matrix[j,i-1] + chunk_comp8[jj] - chunk_comp8[jj_compl]
    }
  }
  
  ##### For the plot: #####
  window_mid <- starts + (windowsize/2)
  title <- paste("Kmer:", oligo, "/File:", accession, "/Window:", windowsize, "/Tot length:", seq_length, sep = " ")
  scaled_gc <- 1000*cum_gc_vec
  min_y = min(c(min(kmer_matrix), min(scaled_gc)))
  max_y = max(c(max(kmer_matrix), max(scaled_gc)))
  plot(window_mid,scaled_gc,type="l", col="blue", main = title,  xlab="Middle pos of each window",ylab="Red=Cum kmer count, Blue=Cum GC-skew", ylim=(c(min_y,max_y)))
  for (m in 1:tot) {
    lines(window_mid, kmer_matrix[m,], type = "l", col = "red")
  }
}

#***************************** Test different values of the parameters below ****************************

##### Import DNA sequence in one of three possible ways: ###############

#accession_nb <-("NC_004307")
#dna_seq <- getncbiseq(accession_nb)

#dna_bin <- read.GenBank(accession_nb)
#write.dna(dna_bin, file="dna_temp.fasta", format = "fasta", append = FALSE)
#fasta_name <- "dna_temp.fasta"

#fasta_name <- "fasta_files/BacillusSubtilisT30.fasta" #Change to your file name and path
fasta_name <- "fasta_files/BifidobacteriumLongum105A.fasta"
#fasta_name <- "fasta_files/EColiBstrREL606.fasta"
#fasta_name <- "fasta_files/ParvibaculumLavamentivoransDS1.fasta"
#fasta_name <- "fasta_files/RickettsiaSlovaca13B.fasta"
dna_temp <- read.fasta(fasta_name)
dna_seq <- getSequence(dna_temp[[1]])
oligo = 8 #Length of kmer
window <- ceiling(length(dna_seq)*0.2) #Windowsize = X % of sequence length
step <- ceiling(window*0.1) #Steplength = 10 % of windowsize
find_kmers_wind <- ceiling(length(dna_seq)*0.2) #Window to use when selecting kmers
top_kmers = 10 #How many of the most frequent kmers in each window should be included?
seq_name <- strsplit(fasta_name, "/", fixed=TRUE)
seq_name <- strsplit(unlist(seq_name)[2], ".", fixed = TRUE)
seq_name <- unlist(seq_name)[1]
slidingwindowplot(window, find_kmers_wind, step, dna_seq, seq_name, oligo, top_kmers)




