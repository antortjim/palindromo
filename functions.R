library("purrr") #convierte factores en characters
library("Biostrings")
library("ggplot2")
library("dplyr")
library(devtools)
library("seqinr")
library("seqLogo")
source_gist(4676064)

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



select.by.mismatches <- function(df, mypattern, mysubject, mismatch, colnames)
{
  vmatchPattern.result <- vmatchPattern(pattern  = mypattern,
                                        #subject is an object that stores the sequences returned by FIMO
                                        subject = mysubject,
                                        max.mismatch = mismatch, min.mismatch = mismatch)
  
  df.rows <- vmatchPattern.result %>% elementNROWS() %>% as.logical()
  mydf <- df[df.rows, ]
  mydf$n.mismatch <- mismatch
  return(mydf)
}

break.by.underscore.fields <- function(big.set, underscore.fields.column, seqId.fields)
{
  small.set <- big.set[[underscore.fields.column]] %>% as.character() %>% strsplit(split = "_") %>% as.data.frame.list()
  colnames(small.set) <- seqId.fields
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



complementary.reversal <- function(data, sequence.var = "matched.sequence", strand.var = "Strand")
{
  sequences <- data[[sequence.var]]
  strand <- data[[strand.var]]
  names(sequences) <- strand
  for(i in 1:length(sequences))
  {
    if (names(sequences[i]) == "rev")
    {
      sequences[i] <- sequences[i] %>% s2c %>% rev %>% comp %>% toupper %>% c2s
    }
  }
  data[[sequence.var]] <- sequences
  return(data)
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
fimo.output.integrator <- function(fimo.file, seqId.fields, DNA.pattern, mismatch, expression, threshold = 3)
{
  #Lee la tabla y guardala en la variable fimo
  fimo <- read.table(fimo.file, header = T, sep = "\t", quote = "")
  fimo <- fimo %>% dplyr::select(-pattern.name, -strand)
  DNA.subject <- fimo.subject.sequences(fimo)
  
  #La columna de ids contiene varios datos en la misma cadena. Separa esos datos por _ y
  ## Genera una data.frame en la que cada fila representa un TSS y hay 4 columnas:
  ## source, coordenada, gen y tipo de tss, y ponles el nombre adecuado
  fimo.tss <- break.by.underscore.fields(fimo, "sequence.name", seqId.fields)

  # Añade esta data.frame a fimo (la principal) 
  fimo <- cbind(fimo.tss, fimo %>% dplyr::select(-sequence.name))
  rm(fimo.tss)
  fimo <- make.numeric.those.numeric(fimo)
  fimo <- fimo %>% filter(Source == "chr")
  ## Genera una data.frame analoga con los TSS que tienen un palindromo con el numero de mismatches predefinido
  ## en su entorno
  fimo.subset <- select.by.mismatches(fimo, DNA.pattern, DNA.subject, mismatch, seqId.fields)
  
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
  
  
  
  
  #haz un subset de expression.df que se quede con tss que no estan en los plasmidos (mejorar codigo para soportar plasmidos)
  #y haz un subset que solo muestre informacion sobre los tss que hemos seleccionado de fimo
  #con eso, filtramos la informacion transcriptomica de los tss que ha detectado fimo
  expression.df <- expression.df %>% filter(Source == "chr")
  #expression.df <- subset(expression.df, expression.df$TSS %in% fimo.subset$TSS)


  df <- dplyr::inner_join(expression.df, fimo.subset, by = c("TSS", "Annotation", "TSS.class","Source"))
  ##Correct the sequences found in the - strand to show their complementary (defaults to the + strand)
  df <- complementary.reversal(df, sequence.var = "matched.sequence", strand.var = "Strand")
  ##Correct the relative position of these sequences with respect to the TSS
  df[df[["Strand"]] == "rev", "start"] <- df[df[["Strand"]] == "rev", "start"] + 10
  df[df[["Strand"]] == "rev", "stop"]  <- df[df[["Strand"]] == "rev", "stop"]  + 10
  
  return(df)
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
