#Load packages and set working space
setwd("~/MEGA/CuartoCurso/TFG/Bioinformatica/")
library("Biostrings")
library("dplyr")
library(devtools)
library("seqinr")
source_gist(4676064)

if(interactive())
{
  max.mismatch <- 2
  min.mismatch <- 1
  threshold <- 3
  upstream.distance <- 100
}
else
{
 arguments <- as.integer(commandArgs(trailingOnly = T))
 max.mismatch <- arguments[1]
 min.mismatch <- arguments[2]
 threshold <- arguments[3]
 upstream.distance <- arguments[4]
}

commandArgs <- function() threshold
source("rna-seq-accession.R")


 
##IMPORTANTE: DA ERROR SI FIMO CONTIENE INSTANCIAS DE TSS EN PLASMIDOS Y NO EN EL CROMOSOMA PRINCIPAL
##El error consiste en que quiere añadr columnas al data.frame relevant.tss.expression con mas
## filas que las que tiene la propia df
comparator <- function(DNA.pattern, DNA.subject, mismatch = 1, hits.df, expression, threshold = 3)
{
  ## Genera una data.frame en la que cada fila representa un TSS y hay 4 columnas:
  ## source, coordenada, gen y tipo de tss
  hits.df.tss <- hits.df[["sequence.name"]] %>% as.character() %>% strsplit(split = "_") %>% as.data.frame()
  colnames(hits.df.tss) <- c("source","coordinate","gene","TSS.type")
  
  ## Genera una data.frame analoga con los TSS que tienen un palindromo con el numero de mismatches predefinido
  ## en su entorno
  mismatches <- unlist(vmatchPattern(pattern  = DNA.pattern,
                                       subject = DNA.subject,
                                       max.mismatch = mismatch, min.mismatch = mismatch)) %>% 
  names() %>% 
  strsplit(split = "_") %>% 
  as.data.frame.list()
  colnames(mismatches) <- c("source","coordinate","gene","TSS.type")
  
  ##Obten el vector de caracteres relevant.tss, que no es mas que:
  #para los tss que encuentre que tienen un palindromo con los mismatches que sean, recompon su nombre (source_coord_gen_tipoTSS)
  relevant.tss <- with(mismatches, paste0(source, "_", coordinate, "_", gene, "_", TSS.type))
  
  #determina qué posiciones de hits.df hay que extraer: aquellas que correspondan a la primera ocurrencia
  #de un relevant.tss.
  relevant.tss.index <- hits.df$sequence.name %in% relevant.tss & ! duplicated(hits.df$sequence.name)
  #haz un subset de fimo que contenga solo esos tss
  hits.df.subset <- subset(hits.df, relevant.tss.index)
  
  #genera el subset equivalente (las mismas filas) de la data.frame de los tss y extrae las coordenadas
  hits.df.tss.subset <- subset(hits.df.tss, relevant.tss.index)
  relevant.tss.coordinates <- hits.df.tss.subset$coordinate
  
  relevant.tss.start <- hits.df.subset$start
  relevant.tss.stop  <- hits.df.subset$stop
  relevant.tss.gene <- hits.df.tss.subset$gene %>% as.character()
  
  #extrae las coordenadas de todos los tss en expression.df y quedate con las posiciones de los que estan
  #en relevant.tss.coordinates para saber qué posiciones de expression.df nos interesan (ok.positions)
  all.coordinates <- expression.df$coordinate
  relevant.tss.positions <- which(all.coordinates %in% relevant.tss.coordinates)
  

  relevant.tss.expression <- expression.df[relevant.tss.positions,]
  
  vector1 <- relevant.tss.expression$coordinate %>% as.character() %>% as.numeric()
  vector2 <- hits.df.tss.subset$coordinate %>% as.character() %>% as.numeric()
  mynumbers <- match(vector1, vector2)
  
  
  relevant.tss.expression$start <- relevant.tss.start[mynumbers]
  relevant.tss.expression$stop  <- relevant.tss.stop[mynumbers]
  relevant.tss.expression$sequence <- hits.df[relevant.tss.index,]$matched.sequence[mynumbers]
  relevant.tss.expression$TSS.type <- hits.df.tss [relevant.tss.index,]$TSS.type[mynumbers]
  relevant.tss.expression$gene <- relevant.tss.gene[mynumbers]
  relevant.tss.expression$n.mismatch <- mismatch
  
  return(relevant.tss.expression)
}




#Define the palindromic sequence
palindromo <- "TACCCGGGTA"   # palindromo presente en el promotor de all1873
# muy conservado en todas las cianos
DNAS.palindromo <-  DNAString(palindromo)  # pasamos la secuencia del palindromo
# a DNAString (Biostrings package)


#Incorporate FIMO output as data.frame
fimo <- read.table("fimo.txt", header = T, sep = "\t", quote = "")
hits.df <- fimo
DNA.pattern <- DNAS.palindromo

#Parse tss metadata structured by _
fimo.tss <- fimo[["sequence.name"]] %>% as.character() %>% strsplit(split = "_") %>% as.data.frame()
colnames(fimo.tss) <- c("source","coordinate","gene","TSS.type")

#Make match.sequence a character vector where every element i is the ith sequence detected
match.sequence <-dplyr::select(fimo, matched.sequence)[,1] %>% as.character()

names(match.sequence) <- as.character(fimo[["sequence.name"]])

#Transform it into Biostrings data
DNAS.match.sequence <- DNAStringSet(match.sequence)

DNA.subject <- DNAS.match.sequence
## Initialize the data.frame that will store extracted data
data <- data.frame()

## For every possible mismatch, extract data and add it to that of all other mismatches
for(mismatch in min.mismatch:max.mismatch)
{
  print(mismatch)
  current.mismatch <- comparator(DNAS.palindromo, DNAS.match.sequence, mismatch, hits.df = fimo, expression = expression.df[,-(3:4)], threshold)
  data <- rbind(data, current.mismatch)
}

##Correct the sequences found in the complementary strand
sequences <- data$sequence %>% as.character()
for (row in 1:nrow(data))
{
  current.strand <- data$strand.character[row]
  if(current.strand == "-")
  {
    print(row)
    current.sequence <- sequences[row] %>% as.character() %>% strsplit(split = "") %>% unlist()
    current.sequence <- current.sequence %>% seqinr::comp() %>% rev() %>% toupper()
    current.sequence <- paste(current.sequence, collapse = "") 
    print(current.sequence)
    sequences[row] <- current.sequence
  }
}

data$sequence <- sequences


## Make the number of mismatches a factor (n.mismatch) so that ggplot can handle it
data$n.mismatch <- factor(data$n.mismatch, levels = min.mismatch:max.mismatch)

data <- data[with(data, order(n.mismatch, significant, start)), ]

data$Distance.to.TSS <- data$start - upstream.distance

Alicia <- data %>% dplyr::select(coordinate, strand.character, gene, TSS.type, Distance.to.TSS, fold.change, significant, sequence, n.mismatch)


data$symbol <- ifelse(data$Distance.to.TSS > 0, "+", "-")
data$Distance.to.TSS <- abs(data$Distance.to.TSS)



write.table(data, file = "fimo_data.txt", row.names = F)

##Porcentaje de TSS con un palindromo de 1 mismatch en torno a el
100*nrow(data)/nrow(expression.df)

##Visualiza la distribucion de los tipos de TSS en el fondo y en los tss de interes
plot(table(fimo.tss[["TSS.type"]]))
plot(table(data[["TSS.type"]]))


## Graphics production
ggplot(data = data, aes(x = log10(WT0), y = log10(WT8), col = significant)) + geom_point()
ggsave("plot1_1.pdf")

ggplot(data = data, aes(x = log10(WT0), y = log10(WT8), col = significant, shape = n.mismatch)) + 
  geom_point() +
  scale_shape_manual(values = c(17, 19))
ggsave("plot1_2.pdf")
##Perfil de la distribucion del comienzo de los palindromos en el entorno de los TSS de interes
##Cuanto más alta sea la señal, es que la frecuencia con la que el TSS se posiciona ahí es mayor
ggplot(data = data, aes(x = start - upstream.distance)) + geom_histogram(binwidth = 5)
ggsave("plot2.pdf")
ggplot(data = data, aes(x = start - upstream.distance, fill = significant)) + geom_histogram()
ggsave("plot3.pdf")
ggplot(data = data, aes(x = start - upstream.distance)) + geom_density(size = 1.5)
ggsave("plot4.pdf")
ggplot(data = data, aes(x = start - upstream.distance, col = significant)) + geom_density(size = 1.5)
ggsave("plot5.pdf")
ggplot(data = data, aes(x = significant, fill = significant)) + geom_bar()
ggsave("plot6.pdf")
