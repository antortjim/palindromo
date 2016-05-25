consensus.displayer <- function(seqs, n.mismatch,  mismatch.threshold)
{
  ## Select only those seqs that have less mismatches than stated in the treshold
  ## if threshold is 2, only 0 and 1 mismatches are selected
  seqs <- seqs[n.mismatch < mismatch.threshold,]
  
  ## Initialize a vector that will store the consensus score of every position i from
  ## the sequence in the i position of the vector
  conservados <- numeric(length = 10)
  
  # for every position in the sequence
  for(i in 1:10)
  {
    #compute how many rows have the same letter as the one stated in the first row
    conserva <- sum(seqs[1,i] == seqs[-1,i])
    #add this sum to the corresponding position on conservados
    conservados[i] <- conserva
  }
  #normalize the vector so that scores range between 0 and 1
  conservados <- conservados / nrow(seqs)
  names(conservados) <- seqs[1,] %>% as.list() %>% unlist() %>% as.character
  return(conservados)
} 

  


