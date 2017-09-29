#Summary:
#Import a FASTA file to be analyzed for kmer bias and GC skew
#First go through the sequence once using smaller, non-overlapping windows to find interesting kmers
#From each window, select the kmers with the most unlikely distribution between the strands. Then only save the kmers that were found in more than 1 window.
#Go through the chromosome again, using larger, overlapping windows, where each window is divided into two halves
#Count the occurences of the kmers and the complementary kmers (=kmers on the other strand) in both window halves 
#Calculate the difference between (the kmer counts of one window half + the complementary kmer counts on the other window half) and the opposite.
#This corresponds to (count_leading - count_lagging) below, even though it is not known which segment is on the leading or lagging strand, or if both seqments include both strands.
#Use the formula from Worning et al. 2006 for each kmer: (count_leading - count_lagging)*log2((count_leading+r)/(count_lagging+r)), where r is a constant and r=5 seems to be a good choice.
#Sum this value for all kmers in each window and plot the sum. Also interpolate the values to get as many data points as one gets from the GC3 function.
#Also calculate and plot the cumulative GC skew. 


#Required packages: seqinr, stringr, ape, pracma
#require("seqinr")
#require("stringr")
#require("ape")
#require("pracma")
kmers<-function(mydna, verb){
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
slidingwindowplot <- function(windowsize, find_kmers_wind, step_len, inputseq, accession, oligo, max_kmers, gc_step)
{
  #windowsize = length of window to be used in the main analysis
  #find_kmers_wind = length of window when identifying most frequent kmers
  #step_len = step length between consecutive overlapping windows in main analysis
  #inputseq = chromosomal sequence to be analyzed
  #accession = name of the sequence, to be printed in the title of the plot
  #oligo = length of kmer
  #max_kmers = Maximum number of kmers to include in the analysis
  #gc_step = Step length in the GC3 function
  seq_length <- length(inputseq)
  oct_per_wind <- find_kmers_wind-oligo+1 #number of kmers in each window in the first step
  starts <- seq(1, seq_length, by = find_kmers_wind) #For stepping through the chromosome the first time (no overlaps between windows)
  n <- length(starts) 
  chunk_comp8_len <- 4^oligo
  kmer_vec <- numeric(chunk_comp8_len) #Store selected kmers in this vector
  kmer_p <- numeric(chunk_comp8_len) +1
  kmer_g <- numeric(chunk_comp8_len)
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
    chunk_comp1 <- count(chunk,1) #Count number of each nucleotide in window
    chunk_comp8 <- count(chunk,oligo) #Count all octamers in the window
    compl_comp8_names <- get_compl_rev(names(chunk_comp8)) #Get all the complementary kmers
    index_comp8 <- sort.list(compl_comp8_names) #Get indeces for sorting complementary kmers alphabetically
    meta_index_comp8 <- sort.list(index_comp8) #Get indeces for sorting the alphabetical list back again
    compl_comp8 <- chunk_comp8[meta_index_comp8] #Get the counts of the complementary kmers by reordering the chunk_comp8
    
    #Calculate the probability of finding the observed distribution of each kmer by chance
    for (p in 1:chunk_comp8_len) { #Go through all kmers
      kmer <- chunk_comp8[p]
      compl <- compl_comp8[p]
      tot_count <- kmer + compl
      p1 <- 2*nchoosek(tot_count, compl)*(0.5^tot_count) #Probability of finding the observed distribution of the kmer on both strands by chance
      if (p1 < 0.01) { #Keep the kmers with low probability of getting their distribution by chance
        found <- 0
        for (g in 1:next_row) {
          #If kmer or compl has already been added to the vector, increment the counter in kmer_g and update kmer_p to minimum p1 value
          if (kmer_vec[g] == names(kmer) | kmer_vec[g] == names(compl)) { 
            kmer_g[g] <- kmer_g[g] + 1 #How many times has this kmer been found?
            kmer_p[g] <- min(kmer_p[g], p1) #What is the minimum p1 value of the kmer?
            found <- 1
            break
          }
        }
        if (found == 0) { #If kmer or compl was not found in the vector, add it and increment the row counter
          kmer_vec[next_row] <- names(kmer)
          kmer_g[next_row] <- kmer_g[next_row] + 1
          kmer_p[next_row] <- p1
          next_row <- next_row+1
        }
      }
    }
  }
  if (next_row == 1) {
    cat("No significantly skewed kmers found! \n")
    stop()
  }
  kmers_found <- numeric(next_row-1)
  kmer_p_found <- numeric(next_row-1)
  kmer_g_found <- numeric(next_row-1)
  tot <- 0
  for (t in 1:(next_row-1)) {
    #If a kmer was found in more than g_min of the windows, save it in the kmers_found vector
    #Also save its kmer_g count and its p1 value 
    g_min <- 1
    if (kmer_g[t] > g_min) {
      tot <- tot + 1
      kmers_found[tot] <- kmer_vec[t]
      kmer_p_found[tot] <- kmer_p[t]
      kmer_g_found[tot] <- kmer_g[t]
    }
  }
  if (tot == 0) {
    cat("No skewed kmers existed in more than ",g_min," windows! \n")
    stop()
  }
  cat("Number of kmers found in more than ", g_min, " windows: ", tot, "\n" )
  #Order the kmers by the p1 values divided by the kmer_g count values
  pg_vec <- kmer_p_found / kmer_g_found
  order_by_pg <- sort.list(pg_vec, decreasing = FALSE) 
  kmers_found <- kmers_found[order_by_pg]
  tot <- min(max_kmers, tot) #If more than max_kmers were found, keep only the top max_kmers based on the pg values
  kmers_found <- sort(kmers_found[1:tot]) #Put kmers in same order as they come from the count() function
  
  kmers_compl_rev <- get_compl_rev(kmers_found) #Get complementary, inverted kmers
  #Sort the complementary kmers so that they come in the same orders as in the output from the count() function
  #But before sorting, save information about how it will be sorted in compl_index, so that the complementary kmers
  #can be rearranged to the initial order again later to match the order in the kmers_found vector
  compl_index <- sort.list(kmers_compl_rev) #Get indeces for sorting the complementary kmers
  meta_index <- sort.list(compl_index) #Get indeces for sorting it back to the initial order
  sorted_compl_rev <- sort(kmers_compl_rev) #Actually sort the complementary kmers
  
  ##### Now use larger windows, divided into two halves, and overlapping each other ########
  oct_per_wind <- windowsize-oligo+1 #Update number of kmers per window
  window_half <- windowsize/2
  mid_points <- seq(1, seq_length, by = step_len)
  n <- length(mid_points) 
  #Matrix for storing differences in counts between pairs of window halves for each selected kmer:
  kmer_not_cum <- matrix(0L,nrow=tot,ncol=n)
  cum_gc_vec <- numeric(n) #For cumulative GC skew
  #Start by looking only at the first window half, to get the indices of the selected kmers in the output from the count() function
  chunk <- inputseq[(seq_length-window_half):seq_length]
  chunk_comp8 <- count(chunk,oligo)
  chunk_comp1 <- count(chunk,1)
  C_w = chunk_comp1[2]
  G_w = chunk_comp1[3]
  k = 1
  h=1
  index_vec = numeric(tot)
  index_vec_compl = numeric(tot)
  #Vectors for saving counts of kmers and complementary kmers temporarily
  temp_kmers = numeric(tot)
  temp_compl = numeric(tot)
  for (j in 1:chunk_comp8_len){ 
    #Go through all kmers from count() to search for the selected kmers and complementary kmers
    #Since the kmers come in the same order in all three cases, the search becomes easier
    #(The kmers_found and sorted_compl_rev are sorted alphabetically and the output from count() comes in alphabetical order)
    if (k <= tot & (names(chunk_comp8[j])==kmers_found[k])) {
      #Save the counts of the selected kmers
      temp_kmers[k] <- chunk_comp8[j] 
      index_vec[k] = j #Save the chunk_comp8-index of each kmer
      k = k+1
    }
    if ( h <= tot & (names(chunk_comp8[j])==sorted_compl_rev[h])) {
      #Save the counts of the complementary kmers
      temp_compl[h] <- chunk_comp8[j] 
      index_vec_compl[h] = j #Save the chunk_comp8-index of each complementary kmer
      h = h+1
    }
    if (k>tot & h>tot){ #When all selected kmers and complementary kmers are found
      break
    }
  }
  #Rearrange the complementary kmers and indeces so that their order again matches the order of the selected kmers
  index_vec_compl <- index_vec_compl[meta_index]
  temp_compl <- temp_compl[meta_index]
  
  #Look at the second window half
  chunk <- inputseq[1:window_half]
  chunk_comp8 <- count(chunk,oligo)
  chunk_comp1 <- count(chunk,1)
  C_w2 = chunk_comp1[2]
  G_w2 = chunk_comp1[3]
  cum_gc_vec[1] <- (G_w + G_w2 - C_w - C_w2)/(G_w + G_w2 + C_w + C_w2) #Cumulative GC skew
  for (j in 1:tot) { #Go through all selected kmers
    jj <- index_vec[j] #Get index of the kmer in chunk_comp8
    jj_compl <- index_vec_compl[j] #Get the index of the compementary kmer in chunk_comp8
    #Calculate the difference between kmer counts and complementary kmer counts in both window halves:
    str1_temp <- temp_kmers[j] + chunk_comp8[jj_compl]
    str2_temp <- temp_compl[j]+chunk_comp8[jj]
    kmer_not_cum[j,1] <- (str1_temp - str2_temp)*log2((str1_temp+5)/(str2_temp+5)) #According to formula used by Worning et al. 2006
  }
  
  #### Now look at all the remaining windows #####
  for (i in 2:n) {
    #First half of the window
    start <- mid_points[i] - window_half
    end <- mid_points[i]
    if (start < 1) { #A negative start position means start from the end of the sequence
      start <- start + seq_length
    }
    if (end > seq_length) { #If the end is after the end of the sequence, start from the beginning of the sequence again
      end <- end - seq_length
    }
    if (end < start) { #Combine the last part of the sequence and the first part of the sequence
      chunk <- c(inputseq[start:seq_length],inputseq[1:end])
    }
    else {
      chunk <- inputseq[start:end]
    }
    chunk_comp8 <- count(chunk,oligo)
    chunk_comp1 <- count(chunk,1)
    C_w = chunk_comp1[2]
    G_w = chunk_comp1[3]
    for (j in 1:tot) { #Get the counts of all selected kmers and their complements
      jj <- index_vec[j]
      jj_compl <- index_vec_compl[j]
      temp_kmers[j] <- chunk_comp8[jj] 
      temp_compl[j] <- chunk_comp8[jj_compl]
    }
    #Second half of the window:
    start2 <- end
    end2 <- start2 + window_half 
    if (end2 > seq_length){
      end2 <- end2 - seq_length
    }
    if (end2 < start2) {
      chunk <- c(inputseq[start2:seq_length],inputseq[1:end2])
    }
    else {
      chunk <- inputseq[start2:end2]
    }
    chunk_comp8 <- count(chunk,oligo)
    chunk_comp1 <- count(chunk,1)
    C_w2 = chunk_comp1[2]
    G_w2 = chunk_comp1[3]
    #Cumulative GC-skew:
    cum_gc_vec[i] <- (G_w + G_w2 - C_w - C_w2)/(G_w + G_w2 + C_w + C_w2) + cum_gc_vec[i-1]
    for (j in 1:tot) { 
      jj <- index_vec[j]
      jj_compl <- index_vec_compl[j]
      #Calculate the difference between kmer counts and complementary kmer counts in both window halves:
      str1_temp <- temp_kmers[j] + chunk_comp8[jj_compl]
      str2_temp <- temp_compl[j]+chunk_comp8[jj]
      kmer_not_cum[j,i] <- (str1_temp - str2_temp)*log2((str1_temp+5)/(str2_temp+5)) #According to formula used by Worning et al. 2006
    }
  }
  
  tot_change <- colSums(kmer_not_cum) #Sum the difference between window halves over all kmers in each window
  #Interpolate to get as many data points as from the GC3 calculation
  interpol_len <- ceiling(seq_length/gc_step)
  tot_change_interpolated <- numeric(interpol_len)
  pos_interpolated <- seq(1,seq_length,by=gc_step) #All positions included in the inerpolated results
  temp_distances <- numeric(n)
  temp_index <- numeric(n)
  closest <- numeric(2)
  vec_pos <- 1
  for (seq_pos in pos_interpolated) { #Go through all positions of the interpolated data
    temp_distances <- abs(mid_points - seq_pos) #Distances between the position and all midpoints
    #Get the two closest midpoints:
    temp_index <- sort.list(temp_distances,decreasing = FALSE)
    first <- temp_index[1]
    #If the position is after the last midpoint, the first midpoint of the sequence is one of the two closest:
    if (seq_pos > mid_points[length(mid_points)]) {
      second <- 1
      temp_distances[second] <- seq_length - seq_pos + mid_points[1] #Distance betw. position and the first midpoint of the sequence
    }
    else {
      second <- temp_index[2]
    }
    closest <- tot_change[c(first,second)] #Use the obtained indeces to find the elements of tot_change that correspond to the two closest midpoints
    #Multiply the values of tot_change by a factor corresponding to the distance between the position and each of the midpoints
    second_factor <- temp_distances[first]/(temp_distances[first]+temp_distances[second])
    first_factor <- temp_distances[second]/(temp_distances[first]+temp_distances[second])
    tot_change_interpolated[vec_pos] <- (closest[1]*first_factor) + (closest[2]*second_factor)
    vec_pos <- vec_pos +1
  }
  #### For plotting #####
  title <- paste(oligo, "mers ", "/Seq: ", accession, " /Wind: ", windowsize, " /Tot len: ", seq_length, sep = "")
  scaling_factor <- max(tot_change) - min(tot_change)
  shifting_const <- mean(tot_change)
  scaled_gc <- scaling_factor*cum_gc_vec + shifting_const #So that it's visible in the plot
  min_y = min(c(min(scaled_gc), min(tot_change)))
  max_y = max(c(max(scaled_gc), max(tot_change)))
  if(verb=='n'){
    plot(mid_points,scaled_gc,type="l", col="blue", main = title, xlab="Middle pos of each window",ylab="Blue=Scaled cum GCskew, Green=Sum(kmer change), Red=Interpol green", ylim=(c(min_y,max_y)))
    lines(mid_points, tot_change, type = "p", col = "green")
    lines(pos_interpolated,tot_change_interpolated, type = "l", col = "red")
  }
  return(list(mid_points, cum_gc_vec, tot_change, pos_interpolated, tot_change_interpolated))
}

#***************************** Test different values of the parameters below ****************************

##### Import DNA sequence in one of three possible ways: ###############

#accession_nb <-("NC_004307")
#dna_seq <- getncbiseq(accession_nb)

#dna_bin <- read.GenBank(accession_nb)
#write.dna(dna_bin, file="dna_temp.fasta", format = "fasta", append = FALSE)
#fasta_name <- "dna_temp.fasta"


#fasta_name <- "fasta_files/BacillusSubtilisT30.fasta" #Change to your file name and path
#fasta_name <- "fasta_files/BifidobacteriumLongum105A.fasta"
#fasta_name <- "fasta_files/EColiBstrREL606.fasta"
#fasta_name <- "fasta_files/ParvibaculumLavamentivoransDS1.fasta"
#fasta_name <- "fasta_files/RickettsiaSlovaca13B.fasta"
#fasta_name <- "fasta_files/BacteroidesThetaiotaomicron7330.fasta"

#dna_temp <- read.fasta(fasta_name)
#dna_seq <- getSequence(dna_temp[[1]])

#4:
#source("readFasta.R")
#dna_temp <- readFasta()
dna_temp<-mydna
seq_name <- getName(dna_temp[[1]])
dna_seq <- getSequence(dna_temp[[1]])

oligo = 8 #Length of kmer
window <- ceiling(length(dna_seq)*0.2) #Windowsize = X % of sequence length
step <- ceiling(window*0.1) #Steplength = Y % of windowsize
find_kmers_wind <- ceiling(length(dna_seq)*0.2) #Window to use when selecting kmers
max_kmers <- 2000L #Maximum number of kmers to include in the analysis
gc_step <- 100L #Step length in the GC3 function
#seq_name <- strsplit(fasta_name, "/", fixed=TRUE)
#seq_name <- strsplit(unlist(seq_name)[2], ".", fixed = TRUE)
#seq_name <- unlist(seq_name)[1]
seq_skews <- slidingwindowplot(window, find_kmers_wind, step, dna_seq, seq_name, oligo, max_kmers, gc_step)

#return(seq_skews)
}
