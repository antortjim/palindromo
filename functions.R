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
  small.set <- big.set[[underscore.fields.column]] %>% as.character() %>% strsplit(split = "_") %>% as.data.frame()
  colnames(small.set) <- seqId.fields
  return(small.set)
}


count.reps <- function(df, feature.containing.repeats)
{
  repeated.positions <- df[[feature.containing.repeats]] %>% as.list %>% duplicated
  mydf <- data.frame((df[[feature.containing.repeats]][repeated.positions] %>% table) +1)
  colnames(mydf) <- c(feature.containing.repeats, "total.count")
  mydf <- full_join(df, mydf, by = feature.containing.repeats)
  mydf[["total.count"]][is.na(mydf$total.count)] <- 1
  return(mydf)
}


# subset.by.underscore.fields <- function(big.set, small.set)
# {
# }

#Define the palindromic sequence
palindromo <- "TACCCGGGTA"   # palindromo presente en el promotor de all1873
# muy conservado en todas las cianos
DNAS.palindromo <-  DNAString(palindromo)  # pasamos la secuencia del palindromo
# a DNAString (Biostrings package)
DNA.pattern <- DNAS.palindromo

#Parse tss metadata structured by _ (Source_TSS_Annotation_TSS.class)
seqId.fields <- c("Source","TSS","Annotation","TSS.class")


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



##IMPORTANTE: DA ERROR SI FIMO CONTIENE INSTANCIAS DE TSS EN PLASMIDOS Y NO EN EL CROMOSOMA PRINCIPAL
##El error consiste en que quiere añadr columnas al data.frame relevant.tss.expression con mas
## filas que las que tiene la propia df

## This function organizes the fimo output by number of mismatches, and adds a significant expression tag based on threshold value
fimo.output.integrator <- function(fimo.file, seqId.fields, DNA.pattern, DNA.subject, mismatch, expression, threshold = 3)
{
  #Lee la tabla y guardala en la variable fimo
  fimo <- read.table(fimo.file, header = T, sep = "\t", quote = "")
  fimo <- fimo %>% dplyr::select(-pattern.name, -strand)
  DNA.subject <- fimo.subject.sequences(fimo)
  
  #La columna de ids contiene varios datos en la misma cadena. Separa esos datos por _ y
  ## Genera una data.frame en la que cada fila representa un TSS y hay 4 columnas:
  ## source, coordenada, gen y tipo de tss, y ponles el nombre adecuado
  fimo.tss <- break.by.underscore.fields(fimo, "sequence.name", seqId.fields)
  ## Añade esta data.frame a fimo (la principal) 
  fimo <- cbind(fimo.tss, fimo %>% dplyr::select(-sequence.name))
  fimo$TSS <- fimo$TSS %>% as.numeric
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
  
  
  ## Extrae datos de los TSS
  ##Coordenadas. El resto esta guardado en fimo.subset porque son datos que no van dentro del nombre del TSS
  coordinates <- fimo.subset$TSS
  ## Nombre del gen
  geneNames <- fimo.subset$Annotation
  
  ## Start del palindromo
  startPositions <- fimo.subset$start
  ## Stop del palindromo
  stopPositions  <- fimo.subset$stop

  #haz un subset de expression.df que se quede con tss que no estan en los plasmidos (mejorar codigo para soportar plasmidos)
  #y haz un subset que solo muestre informacion sobre los tss que hemos seleccionado de fimo
  #con eso, filtramos la informacion transcriptomica de los tss que ha detectado fimo
  expression.df <- expression.df %>% filter(Source == "chr")
  #expression.df <- subset(expression.df, expression.df$TSS %in% fimo.subset$TSS)
  
  df <- dplyr::full_join(expression.df, fimo.subset, by = c("TSS", "Annotation"))
  expression.df[with(expression.df, order(TSS)), ]
  fimo.subset[with(fimo.subset, order(TSS)), ]
  
  fimo.subset$TSS %>% class
  expression.df$TSS %>% class
  
  vector1 <- relevant.tss.expression$TSS %>% as.character() %>% as.numeric()
  vector2 <- hits.df.tss.subset$TSS %>% as.character() %>% as.numeric()
  mynumbers <- match(vector1, vector2)
  
  
  relevant.tss.expression$start <- relevant.tss.start[mynumbers]
  relevant.tss.expression$stop  <- relevant.tss.stop[mynumbers]
  relevant.tss.expression$sequence <- hits.df[relevant.tss.index,]$matched.sequence[mynumbers]
  relevant.tss.expression$TSS.class <- hits.df.tss [relevant.tss.index,]$TSS.class[mynumbers]
  relevant.tss.expression$Annotation <- relevant.tss.gene[mynumbers]
  relevant.tss.expression$n.mismatch <- mismatch
  
  return(relevant.tss.expression)
}


## Initialize the data.frame that will store extracted data
data <- data.frame()

## For every possible mismatch, extract data and add it to that of all other mismatches
## fimo.output.integrator works row-based. fimo.output.integrator must be executed once per row. Each iteration adds a new row corresponding
## to a differnet hit in fimo.txt
for(mismatch in min.mismatch:max.mismatch)
{
  print(mismatch)
  current.mismatch <- fimo.output.integrator(DNAS.palindromo, DNAS.match.sequence, mismatch, hits.df = fimo, expression = expression.df[,-(3:4)], threshold)
  data <- rbind(data, current.mismatch)
}

##Correct the sequences found in the - strand to show their complementary (defaults to the + strand)
sequences <- data$sequence %>% as.character()
for (row in 1:nrow(data))
{
  current.strand <- data$Strand[row]
  if(current.strand == "rev")
  {
    current.sequence <- sequences[row] %>% as.character() %>% strsplit(split = "") %>% unlist()
    current.sequence <- current.sequence %>% seqinr::comp() %>% rev() %>% toupper()
    current.sequence <- paste(current.sequence, collapse = "") 
    sequences[row] <- current.sequence
  }
}

## Add them to data passing each position as a different column. Sequences must all be of equal length
sequences <- lapply(sequences %>% as.list, function(x) unlist(strsplit(x, split = ""))) %>% as.data.frame
data <- select(data, -sequence) %>% cbind(sequences)


## Make the number of mismatches a factor (n.mismatch) so that ggplot can handle it
data$n.mismatch <- factor(data$n.mismatch, levels = min.mismatch:max.mismatch)

## Create a new variale that stores the distance to the TSS. Negative numbers mean that the sequence is upstream of the TSS,
## whereas a positive number means downstream of the TSS
data$Distance.to.TSS <- data$start - upstream.distance.prom
## remove rownames because they are non sense
rownames(data) <- NULL





## REMOVING DUPLICATED PALINDROMES
# This block of code removes sequences that are duplicated because there is more than one palindrom close. Only the first occurence in data is kept,
##Get the coordinates of the palindrom (that of the TSS + the distance to the TSS)
palindromes.coordinates <- data$TSS + data$Distance.to.TSS
#Determine which of them are duplicated
pal.coord <- !(palindromes.coordinates %>% duplicated.Vector())


#Create a matrix where the first column is the vector palindromes coordinates, and the second column is the same +9 (encoding the end of the palindromes)
#Call it the same (machaca el vector)
palindromes.coordinates <- matrix(c(palindromes.coordinates, palindromes.coordinates + 9), ncol = 2)
#Create a new vector called symbol where + means that the sequence is beyond the TSS and - that the palindrom is before the TSS
symbol <- ifelse(data$Distance.to.TSS > 0, "+", "-")
#Get the absolute value
Distance.to.TSS.abs <- abs(data$Distance.to.TSS)
#Paste the symbol and the absolute values into a character vector (thus 3 is displayed as +3, making it prettier in the table)
data$Distance.to.TSS.pretty <- paste(symbol, Distance.to.TSS.abs)


#Remove duplicates and order data. data.no.duplicates has no duplicates and is the one what will be used 
data.no.duplicates <- data[pal.coord,]
data <- data[with(data, order(n.mismatch, start)), ]
data.no.duplicates <- data.no.duplicates[with(data.no.duplicates, order(n.mismatch, start)), ]

