#Required packages: seqinr, stringr, ape

getncbiseq <- function(accession) #function from Little book of R for bioinf
{
  require("seqinr") # this function requires the SeqinR R package
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

get_expected <- function(kmer, chunk_comp1, window_factor) {
  kmer_split <- unlist(strsplit(kmer,""))
  kmer_transl <- as.integer(chartr("acgt","1234",kmer_split)) 
  A_w = chunk_comp1[1]
  C_w = chunk_comp1[2]
  G_w = chunk_comp1[3]
  T_w = chunk_comp1[4]
  all_w = A_w + C_w + G_w + T_w
  ACGT_frac = c(A_w/all_w, C_w/all_w, G_w/all_w, T_w/all_w)
  exp_count <- window_factor
  for (i in 1:length(kmer_transl)) {
    j = kmer_transl[i]
    exp_count <- exp_count*ACGT_frac[j]
  }
  return(exp_count)
}

slidingwindowplot <- function(windowsize, find_kmers_wind, step_len, inputseq, accession, oligo)
{
  seq_len <- length(inputseq)
  oct_per_wind <- find_kmers_wind-oligo+1
  top_kmers = 50
  starts <- seq(1, seq_len, by = find_kmers_wind)
  n <- length(starts)    # Find the length of the vector "starts"
  kmer_vec <- numeric(n*top_kmers)
  next_row = 1
  #cum_gc_vec <- numeric(n)
  #cum_gcskew = 0.0
  for (i in 1:n) {
    start <- starts[i]
    end <- start + find_kmers_wind-1
    if (end <= seq_len){
      chunk <- inputseq[start:end]
    }
    else { #Last chunk includes the beginning of the sequence also to become as long as the other chunks
      chunk <- inputseq[start:seq_len]
      chunk <- c(chunk,inputseq[1:(end-seq_len)])
    }
    chunk_comp1 <- count(chunk,1) #Count number of each nucleotide in window
    #C_w = chunk_comp1[2]
    #G_w = chunk_comp1[3]
    #chunk_gcskew <- (G_w - C_w)/(G_w + C_w) #GC-skew = (number of G)/(number of C)
    #cum_gcskew <- cum_gcskew + chunk_gcskew #Cumulative GC-skew, pos when G>C, neg when C>G
    #cum_gc_vec[i] <- cum_gcskew
    chunk_comp8 <- count(chunk,oligo) #Count all octamers in the window
    expected_comp8 <- numeric(length(chunk_comp8))
    for (ii in 1:length(chunk_comp8)) {
      expected_comp8[ii] <- get_expected(names(chunk_comp8[ii]), chunk_comp1, oct_per_wind)
    }
    
    diff_comp8 <- chunk_comp8 - expected_comp8
    sorted_kmers <- sort(diff_comp8, decreasing=TRUE, na.last=TRUE) #Most frequent octamers first
    for (j in 1:top_kmers) { #Save top_kmers most frequent octamers in a vector
      kmer <- names(sorted_kmers[j]) #don't count kmers with >= 7 of same nucleotide
      if (! (kmer %in% kmer_vec)) {
        kmer_vec[next_row] <- kmer
        next_row = next_row + 1
      }
    }
  }
  tot <- next_row-1 #Total number of unique octamers found among 5 most frequent in each window
  kmers_found <- sort(kmer_vec[1:tot]) #Put kmers in same order as they come from the count() function
  overlapping = tot/(n*top_kmers)
  print(overlapping)
  print(kmers_found)
  num_windows <- ceiling(seq_len/step_len) #Number of windows
  kmer_matrix <- matrix(0L,nrow=tot,ncol=num_windows)
  kmer_not_cum <- matrix(0L,nrow=tot,ncol=num_windows)
  kmer_max_diff <- numeric(num_windows)
  starts <- seq(1, seq_len, by = step_len)
  n <- length(starts) 
  cum_gc_vec <- numeric(n)
  diff_gc_vec <- numeric(n)
  oct_per_wind <- windowsize-oligo+1
  window_half <- windowsize/2
  chunk <- inputseq[1:window_half-1] #Look at only the first chunk first
  chunk_comp8 <- count(chunk,oligo)
  chunk_comp1 <- count(chunk,1)
  C_w = chunk_comp1[2]
  G_w = chunk_comp1[3]
  k = 1
  index_vec = numeric(tot)
  for (j in 1:length(chunk_comp8)){ 
    #Go through all kmers from count() to search for the kmers identified above
    #Since the kmers come in the same order in both cases, the search becomes easier
    if (names(chunk_comp8[j])==kmers_found[k]) {
      expected <- get_expected(names(chunk_comp8[j]), chunk_comp1, oct_per_wind)
      kmer_not_cum[k,1] <- chunk_comp8[j] #- expected
      #kmer_not_cum[k,1] <- expected
      index_vec[k] = j #Save the index of each kmer in the output from count()
      k = k+1
    }
    if (k>tot){
      break
    }
  }
  chunk <- inputseq[window_half:windowsize-1] #Look at only the first chunk first
  chunk_comp8 <- count(chunk,oligo)
  chunk_comp1 <- count(chunk,1)
  C_w2 = chunk_comp1[2]
  G_w2 = chunk_comp1[3]
  cum_gc_vec[1] <- (G_w + G_w2 - C_w - C_w2)/(G_w + G_w2 + C_w + C_w2) 
  diff_gc_vec[1] <- abs(((G_w - C_w)/(G_w + C_w)) - ((G_w2 - C_w2)/(G_w2 + C_w2))) 
  for (j in 1:tot) { 
    jj <- index_vec[j]
    expected <- get_expected(names(chunk_comp8[jj]), chunk_comp1, oct_per_wind)
    kmer_not_cum[j,1] <- abs(kmer_not_cum[j,1] - chunk_comp8[jj])
  }
  for (i in 2:n) { #Look at all the remaining chunks
    start <- starts[i]
    end <- start + window_half-1
    if (end <= seq_len){
      chunk <- inputseq[start:end]
    }
    else {
      chunk <- inputseq[start:seq_len]
      end <- end - seq_len #redefine end to avoid problems in next part
      chunk <- c(chunk,inputseq[1:end])
    }
    chunk_comp8 <- count(chunk,oligo)
    chunk_comp1 <- count(chunk,1)
    C_w = chunk_comp1[2]
    G_w = chunk_comp1[3]
    for (j in 1:tot) { #For all selected kmers, get the frequency in each window
      jj <- index_vec[j]
      expected <- get_expected(names(chunk_comp8[jj]), chunk_comp1, oct_per_wind)
      kmer_not_cum[j,i] <- chunk_comp8[jj] #- expected
      #kmer_not_cum[j,i] <- expected
    }
    start2 <- end+1
    end2 <- start2 + window_half-1
    if (end2 <= seq_len){
      chunk <- inputseq[start2:end2]
    }
    else {
      chunk <- inputseq[start2:seq_len]
      chunk <- c(chunk,inputseq[1:(end2-seq_len)])
    }
    chunk_comp8 <- count(chunk,oligo)
    chunk_comp1 <- count(chunk,1)
    C_w2 = chunk_comp1[2]
    G_w2 = chunk_comp1[3]
    cum_gc_vec[i] <- (G_w + G_w2 - C_w - C_w2)/(G_w + G_w2 + C_w + C_w2) + cum_gc_vec[i-1]
    diff_gc_vec[i] <- abs(((G_w - C_w)/(G_w + C_w)) - ((G_w2 - C_w2)/(G_w2 + C_w2))) 
    for (j in 1:tot) { 
      jj <- index_vec[j]
      expected <- get_expected(names(chunk_comp8[jj]), chunk_comp1, oct_per_wind)
      kmer_not_cum[j,i] <- abs(kmer_not_cum[j,i] - chunk_comp8[jj])
    }
  }
  window_mid <- starts + (windowsize/2)
  title <- paste(oligo, accession, windowsize, step_len, sep = ", ")
  scaled_gc <- 10000*cum_gc_vec + 1000
  scaled_diff_gc_vec <- 10000*diff_gc_vec
  tot_change <- colSums(kmer_not_cum)
  min_y = min(c(min(kmer_not_cum), min(scaled_gc), min(tot_change)))
  max_y = max(c(max(kmer_not_cum), max(scaled_gc), max(tot_change)))
  plot(window_mid,scaled_gc,type="l", col="blue", main = title, xlab="Middle pos of each window",ylab="Red=Count per 8mer, Blue=Cum GC-skew", ylim=(c(min_y,max_y)))
  for (m in 1:tot) {
    lines(window_mid, kmer_not_cum[m,], type = "l", col = "red")
  }
  lines(window_mid, tot_change, type = "l", col = "green")
  #lines(window_mid, scaled_diff_gc_vec, type = "l", col = "red")
}

oligo = 8

#accession_nb <-("NC_004307")
#dna_seq <- getncbiseq(accession_nb)
#dna_bin <- read.GenBank(accession_nb)
#write.dna(dna_bin, file="dna_temp.fasta", format = "fasta", append = FALSE)

fasta_name <- "fasta_files/BacillusSubtilisT30.fasta"
#fasta_name <- "fasta_files/BifidobacteriumLongum105A.fasta"
#fasta_name <- "fasta_files/EColiBstrREL606.fasta"
#fasta_name <- "fasta_files/ParvibaculumLavamentivoransDS1.fasta"
#fasta_name <- "fasta_files/RickettsiaSlovaca13B.fasta"
dna_temp <- read.fasta(fasta_name)
dna_seq <- getSequence(dna_temp[[1]])
window <- length(dna_seq)*0.2
step <- window*0.1
find_kmers_wind <- length(dna_seq)*0.1
slidingwindowplot(window, find_kmers_wind, step, dna_seq, fasta_name, oligo)




