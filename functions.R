library("purrr") #convierte factores en characters
library("Biostrings")
library("ggplot2")
library("dplyr")
library(devtools)
library("seqinr")
library("seqLogo")
source_gist(4676064)

consensus.displayer <- function(seqs, mismatch.threshold)
{
  library("dplyr")
  ## Select only those seqs that have less mismatches than stated in the treshold
  ## if threshold is 2, only 0 and 1 mismatches are selected
  seqs <- seqs %>% dplyr::filter(n.mismatch < mismatch.threshold) #seqs[n.mismatch < mismatch.threshold,]
  seqs <- dplyr::select(seqs, -n.mismatch)
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


#Function that produces a matrix that seqlogo understands
alignment.analyzer <- function(alignment.df)
{
  number.of.columns <- ncol(alignment.df)
  conservation <- matrix(0, ncol = number.of.columns, nrow = 4)
  base.code <- 1:4
  names(base.code) <- c("A","C","G","T")
  for (column in 1:number.of.columns)
  {
    #print(column)
    current.col <- alignment.df[,column] %>% as.character() %>% table() / nrow(alignment.df)
    present.bases <- names(current.col)
    
    for (i in 1:length(present.bases))
    {
      current.base <- present.bases[i]
      #print(current.base)
      conservation[base.code[current.base], column] <- current.col[[current.base]]
      
    }
  }
  return(conservation)
}



select.by.mismatches <- function(df, mypattern, mysubject, mismatch)
{
  vmatchPattern.result <- vmatchPattern(pattern  = mypattern,
                                        #subject is an object that stores the sequences returned by FIMO
                                        subject = mysubject,
                                        max.mismatch = mismatch, min.mismatch = mismatch)
  
  
  df.rows <- vmatchPattern.result %>% elementNROWS() %>% as.logical()
  if(sum(df.rows) == 0)
  {
    print(paste("No sequences with", mismatch, "mismatches", sep = " "))
    return(F)
  }
  else
  {
  mydf <- df[df.rows, ]

  mydf$n.mismatch <- mismatch
  return(mydf)
  }
}

break.by.break.fields <- function(big.set, underscore.fields.column, seqId.fields, numeric.fields, break.char)
{
  ## Break strings by break.char and generate a new dataframe where each row is the result of breaking each string
  small.set <- data.frame(do.call('rbind', strsplit(as.character(big.set[[underscore.fields.column]]), break.char, fixed=TRUE)))
  

  ## Add the colnames given by seqId.fields
  colnames(small.set) <- seqId.fields
  ## Set TSS as numeric, instead of as.character (cant be helped by command above)
  small.set[,numeric.fields]<- as.numeric(small.set[,numeric.fields])
  return(small.set)
}

count.reps <- function(df, feature.containing.repeats)
{
  repeated.positions <- df[[feature.containing.repeats]] %>% as.list %>% duplicated
  if(repeated.positions %>% sum != 0)
  {
   mydf <- data.frame((df[[feature.containing.repeats]][repeated.positions] %>% table) +1)
   colnames(mydf) <- c(feature.containing.repeats, "total.count")
   mydf <- full_join(df, mydf, by = feature.containing.repeats)
   mydf[["total.count"]][is.na(mydf$total.count)] <- 1
  }
  else
  {
   df[["total.count"]] <- 1
   mydf <- df
  }
    return(mydf)
}


# subset.by.underscore.fields <- function(big.set, small.set)
# {
# }



complementary.reversal <- function(df, sequence.var = "matched.sequence", strand.var = "Strand")
{
  sequences <- df[[sequence.var]]
  strand <- df[[strand.var]]
  names(sequences) <- strand
  for(i in 1:length(sequences))
  {
    if (names(sequences[i]) == "rev")
    {
      sequences[i] <- sequences[i] %>% s2c %>% rev %>% comp %>% toupper %>% c2s
    }
  }
  df[[sequence.var]] <- sequences
  return(df)
}



fimo.subject.sequences <- function(fimo)
{
   #Make match.sequence a character vector where every element i is the ith sequence detected
   match.sequence        <- dplyr::select(fimo, matched.sequence)[,1] %>% as.character()
   names(match.sequence) <- as.character(fimo[["sequence.name"]])
   
   #Transform match.sequence into Biostrings data
   DNAS.match.sequence <- DNAStringSet(match.sequence)
   DNA.subject <- DNAS.match.sequence
   return(DNA.subject)
}


make.numeric.those.numeric <- function(df)
{
  right.order <- colnames(df)
  df[] <- lapply(df, as.character)
  digits.positions <- grepl("^[[:digit:]]+$", df[1,])
  digits.df <- df[, digits.positions]
  digits.columns <- colnames(digits.df)
  digits.df[] <- lapply(digits.df, as.numeric)
  df <- df[, !digits.positions]
  df <- cbind(df, digits.df)
  df <- df[,right.order]
  return(df)
}

##IMPORTANTE: DA ERROR SI FIMO CONTIENE INSTANCIAS DE TSS EN PLASMIDOS Y NO EN EL CROMOSOMA PRINCIPAL
##El error consiste en que quiere añadr columnas al data.frame relevant.tss.expression con mas
## filas que las que tiene la propia df

## This function organizes the fimo output by number of mismatches, and adds a significant expression tag based on threshold value
fimo.output.integrator <- function(fimo.df, seqId.fields, numeric.fields, DNA.pattern, mismatch, transcriptomics = NULL, threshold = 3, pattern.id, break.char, join.fields)
{
  fimo.df$pattern.id <- pattern.id
  fimo.df <- fimo.df %>% dplyr::select(-strand)
  DNA.subject <- fimo.subject.sequences(fimo.df)
  
  #La columna de ids contiene varios datos en la misma cadena. Separa esos datos por _ y
  ## Genera una data.frame en la que cada fila representa un TSS y hay 4 columnas:
  ## source, coordenada, gen y tipo de tss, y ponles el nombre adecuado
  fimo.tss <- break.by.break.fields(fimo.df, "sequence.name", seqId.fields, numeric.fields, break.char)

  # Añade esta data.frame a fimo (la principal) 
  fimo.df <- cbind(fimo.tss, fimo.df %>% dplyr::select(-sequence.name))
  rm(fimo.tss)
 #fimo.df <- make.numeric.those.numeric(fimo.df)
  fimo.df$Annotation <- as.character(fimo.df$Annotation)
  if(sum(colnames(fimo.df) %in% "Source"))
    {
    fimo.df <- fimo.df %>% filter(Source == "chr")
   }
  ## Genera una data.frame analoga con los TSS que tienen un palindromo con el numero de mismatches predefinido
  ## en su entorno
  fimo.subset <- select.by.mismatches(fimo.df, DNA.pattern, DNA.subject, mismatch)
  #Return false if fimo.subset is not a data.frame, which happens when there are no instances
  #with the supplied number of mismatches
  if(fimo.subset %>% class != "data.frame")
  {
    return(F)
  }
  else
  {
  ## elimina aquellas filas en las que el nombre del gen y la secuencia detectada sean exactamente las mismas
  ## esto ocurre cuando varios TSS se encuentran asociados al mismo gen (desde el sd01.xslx) y ambos resultan
  ## estar en el entorno de un palindromo. Les sale que es el mismo
  ## Este codigo considera duplicado si la secuencia del palindromo detectado coincide porque asume que se
  ## trata de la misma instancia, pero podría ser que se dé la casualidad de que coincida y sin embargo no
  ## sean la misma instancia (son dos secuencias iguales pero en sitios distintos, aunque sean cercanos)
  pasted.seq.ann <- paste(fimo.subset$matched.sequence, fimo.subset$Annotation, sep = "_")
  fimo.subset$pasted.seq.ann <- pasted.seq.ann
  fimo.subset <- count.reps(df = fimo.subset, "pasted.seq.ann") 
  fimo.subset <- fimo.subset %>% dplyr::select(-pasted.seq.ann)
  

  #expression.df <- subset(expression.df, expression.df$TSS %in% fimo.subset$TSS)
  # fimo.subset$TSS <- as.numeric(fimo.subset$TSS)
  # fimo.subset$start <- as.numeric(fimo.subset$start)
  # fimo.subset$stop <- as.numeric(fimo.subset$stop)
  if(class(transcriptomics)  == "data.frame")
  { 
  #haz un subset de expression.df que se quede con tss que no estan en los plasmidos (mejorar codigo para soportar plasmidos)
  #y haz un subset que solo muestre informacion sobre los tss que hemos seleccionado de fimo
  #con eso, filtramos la informacion transcriptomica de los tss que ha detectado fimo
  transcriptomics <- transcriptomics %>% filter(Source == "chr")
  df <- dplyr::inner_join(transcriptomics, fimo.subset, by = join.fields)
  }
  else
  {
  df <- fimo.subset
  }
  
 
  
  ##Correct the sequences found in the - strand to show their complementary (defaults to the + strand)
  df <- complementary.reversal(df, sequence.var = "matched.sequence", strand.var = "Strand")
  ##Correct the relative position of these sequences with respect to the TSS
  df[df[["Strand"]] == "rev", "start"] <- df[df[["Strand"]] == "rev", "start"] + 10
  df[df[["Strand"]] == "rev", "stop"]  <- df[df[["Strand"]] == "rev", "stop"]  + 10
  
  return(df)
  }

}



addDistanceAsCharacter <- function(df, var = "start", distance = upstream.distance.prom, length = 9, new.var = "Distance.to.TSS", formatted.var = "pretty.distance")
{
  ## Create a new variale that stores the distance to the TSS. Negative numbers mean that the sequence is upstream of the TSS,
  ## whereas a positive number means downstream of the TSS
  df[[new.var]] <- df[[var]] - distance
  
  #Create a new vector called symbol where + means that the sequence is beyond the TSS and - that the palindrom is before the TSS
  symbol <- ifelse(df[[new.var]] > 0, "+", "-")
  
  #Get the absolute value
  new.var.abs <- abs(df[[new.var]])
  #Paste the symbol and the absolute values into a character vector (thus 3 is displayed as +3, making it prettier in the table)
  df[[formatted.var]] <- paste(symbol, new.var.abs)
  return(df)
}

transcriptomics.data <- expression.df
read.fimo.data <- function(fimo.file, seqId.fields, numeric.fields, num.ran, break.char, transcriptomics = transcriptomics.data, join.fields = join.fields)
{
  ## Initialize the data.frame that will store extracted data
  data <- data.frame(stringsAsFactors = FALSE)
  ## For every possible mismatch, extract data and add it to that of all other mismatches
  ## fimo.output.integrator works row-based. fimo.output.integrator must be executed once per row. Each iteration adds a new row corresponding
  ## to a differnet hit in fimo.txt
  cmd <- paste("sed -i.backup 's/^#//' ", fimo.file, sep ="")
  system(cmd)
  #Lee la tabla y guardala en la variable fimo
  fimo <- read.table(fimo.file, header = T, sep = "\t", quote = "")
  my.motifs.names <- c("palindromo", paste("random_", 1:num.ran, sep = ""))
  
  for(i in 1:length(my.motifs.names))
  {
    print(i)
    pattern.id <- i -1
    print(my.motifs.names[i])
    mysubset <- filter(fimo, pattern.name == my.motifs.names[i])
    DNA.pattern <- my.motifs.seqs[i] %>% DNAString()
    for(mismatch in min.mismatch:max.mismatch)
    {
      print(mismatch)
      current.mismatch <- fimo.output.integrator(fimo.df = mysubset, seqId.fields = seqId.fields, numeric.fields = numeric.fields,
                                                 DNA.pattern, mismatch, transcriptomics = transcriptomics, threshold,
                                                 pattern.id = pattern.id, break.char = break.char, join.fields)
      if(current.mismatch %>% class == "data.frame")
      {
        data <- rbind(data, current.mismatch)
      }
      rm(current.mismatch)
    }
  }
  data <- addDistanceAsCharacter(data)
  sequences <- data$matched.sequence %>% as.list() %>% lapply(s2c) %>% as.data.frame.list
  data <- cbind(data, sequences); rm(sequences)
  return(data)
}


