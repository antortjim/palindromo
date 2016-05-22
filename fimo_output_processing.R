#Load packages and set working space
setwd("~/MEGA/CuartoCurso/TFG/Bioinformatica/")
library("Biostrings")
library("ggplot2")
library("dplyr")
library(devtools)
library("seqinr")
library("seqLogo")
source_gist(4676064)

 arguments <- as.integer(commandArgs(trailingOnly = T))
 max.mismatch <- arguments[1]
 min.mismatch <- arguments[2]
 threshold <- arguments[3]
 upstream.distance <- arguments[4]


if(interactive())
{
  max.mismatch <- 2
  min.mismatch <- 1
  threshold <- 3
  upstream.distance <- 100
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
  colnames(hits.df.tss) <- c("Source","TSS","Annotation","TSS.class")
  
  ## Genera una data.frame analoga con los TSS que tienen un palindromo con el numero de mismatches predefinido
  ## en su entorno
  mismatches <- unlist(vmatchPattern(pattern  = DNA.pattern,
                                       subject = DNA.subject,
                                       max.mismatch = mismatch, min.mismatch = mismatch)) %>% 
  names() %>% 
  strsplit(split = "_") %>% 
  as.data.frame.list()
  colnames(mismatches) <- c("Source","TSS","Annotation","TSS.class")
  
  ##Obten el vector de caracteres relevant.tss, que no es mas que:
  #para los tss que encuentre que tienen un palindromo con los mismatches que sean, recompon su nombre (source_coord_gen_tipoTSS)
  relevant.tss <- with(mismatches, paste0(Source, "_", TSS, "_", Annotation, "_", TSS.class))
  
  #determina qué posiciones de hits.df hay que extraer: aquellas que correspondan a la primera ocurrencia
  #de un relevant.tss.
  relevant.tss.index <- hits.df$sequence.name %in% relevant.tss & ! duplicated(hits.df$sequence.name)
  #haz un subset de fimo que contenga solo esos tss
  hits.df.subset <- subset(hits.df, relevant.tss.index)
  
  #genera el subset equivalente (las mismas filas) de la data.frame de los tss y extrae las coordenadas
  hits.df.tss.subset <- subset(hits.df.tss, relevant.tss.index)
  relevant.tss.coordinates <- hits.df.tss.subset$TSS
  
  relevant.tss.start <- hits.df.subset$start
  relevant.tss.stop  <- hits.df.subset$stop
  relevant.tss.gene <- hits.df.tss.subset$Annotation %>% as.character()
  
  #extrae las coordenadas de todos los tss en expression.df y quedate con las posiciones de los que estan
  #en relevant.tss.coordinates para saber qué posiciones de expression.df nos interesan (ok.positions)
  all.coordinates <- expression.df$TSS
  relevant.tss.positions <- which(all.coordinates %in% relevant.tss.coordinates)
  

  relevant.tss.expression <- expression.df[relevant.tss.positions,]
  
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
colnames(fimo.tss) <- c("Source","TSS","Annotation","TSS.class")

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
  current.strand <- data$Strand[row]
  if(current.strand == "rev")
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

data <- rbind(data.frame(filter(expression.df, Annotation == "all1873"), start = 101, stop = 111,
           sequence = palindromo, n.mismatch = 0, Distance.to.TSS = 1), data)

color.sequences <- data$sequence %>%
                   as.character() %>%
                   strsplit(split = "") %>% 
                   as.data.frame 
color.sequences

# Problema: resolver problema de longitud de current.col
# cuando la columna analizada no contiene en ninguna 
# secuncia a una de las letras (resulta ser de longitud
# menor al numero de filas de la matrix y da error

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

pdf(file = "~/MEGA/CuartoCurso/TFG/Documento/TFG/Figures/seqlogo.pdf", height = 5)
pwm <- alignment.analyzer(color.sequences) %>% makePWM() %>% seqLogo(xaxis = F, yaxis = F)
dev.off()
system("inkscape -f ~/MEGA/CuartoCurso/TFG/Documento/TFG/Figures/seqlogo.pdf -A ~/MEGA/CuartoCurso/TFG/Documento/TFG/Figures/seqlogo.pdf")

#write.table(file = "palindromos_color.txt", row.names = F)

#Alicia <- data %>% dplyr::select(TSS, Strand, Annotation, TSS.class, Distance.to.TSS, fold.change.WT, significant, sequence, n.mismatch)

#write.table(Alicia, file = "palindromos_hasta_2_mismatches.txt", row.names = F)

data$symbol <- ifelse(data$Distance.to.TSS > 0, "+", "-")


data$Distance.to.TSS <- abs(data$Distance.to.TSS)
Distance.to.TSS <- paste(data$symbol, data$Distance.to.TSS)

tabular.data <- data %>% select(TSS.class:stop, n.mismatch) %>% cbind(Distance.to.TSS, color.sequences)


#write.table(dplyr::select(data, TSS, Strand, Annotation, TSS.class, symbol, Distance.to.TSS, fold.change.WT, sequence), file = "fimo_data.txt", row.names = F)
tabular.data <- select_(tabular.data, "TSS", "Strand", "Annotation", "TSS.class", "Distance.to.TSS", "fold.change.WT", "X1:X10", "n.mismatch")

commandArgs <- function() list(tabular.data, 20)
source("palindromo/tableR.R")


blast <- readline('Would you like to blast interesting genes? Yes/No:\n')
if(blast == "Yes")
{
  query.genes <- readline("Enter interesting genes: ") %>% strsplit(split = " ") %>% unlist()
  commandArgs <- function() list(query.genes)
  source("palindromo/blast.R")
}




sequences <- data$sequence %>% as.character() %>% as.list()
seqs.id <- paste(data$TSS, data$Annotation, data$TSS.class, sep ="_")
write.fasta(file = "sequences.fa", sequences = sequences, names = seqs.id)  

##Porcentaje de TSS con un palindromo de 1 mismatch en torno a el
100*nrow(data)/nrow(expression.df)

##Visualiza la distribucion de los tipos de TSS en el fondo y en los tss de interes
plot(table(fimo.tss[["TSS.class"]]))
plot(table(data[["TSS.class"]]))


## Graphics production
ggplot(data = data, aes(x = log10(WT0), y = log10(WT8), col = significant)) + geom_point()
ggsave("plot1_1.pdf")

ggplot(data = data, aes(x = log10(WT0), y = log10(WT8), col = significant, shape = n.mismatch)) + 
  geom_point() +
  scale_shape_manual(values = c(17, 19,1))
ggsave("plot1_2.pdf")
##Perfil de la distribucion del comienzo de los palindromos en el entorno de los TSS de interes
##Cuanto más alta sea la señal, es que la frecuencia con la que el TSS se posiciona ahí es mayor
ggplot(data = data, aes(x = start - upstream.distance)) + geom_histogram(binwidth = 5)
ggsave("plot2.pdf")
ggplot(data = data, aes(x = start - upstream.distance, fill = significant)) + geom_histogram()
ggsave("plot3.pdf")
ggplot(data = data, aes(x = start - upstream.distance)) + geom_density(size = 1.5)
ggsave("plot4.pdf")

ggplot(data = data, aes(x = start - upstream.distance, col = significant)) + geom_density(size = 1.5) +
  xlab("Intervalo (-100, +50) respecto a los TSS") + ylab("Densidad") + theme_bw() +
  scale_color_manual(name = "Actividad transcripcional\nen el estrés por carencia\nde nitrógeno",
                     labels = c("Disminuye","Aumenta","No significativo"),
                     breaks = c("decrease","increase","no"),
                     values = c("red","green","blue"))

ggsave("../Documento/TFG/Figures/plot5.pdf",  width = 10, height = 5, dpi = 120)

ggplot(data = data, aes(x = significant, fill = significant)) + geom_bar()
ggsave("plot6.pdf")


## Significant genes
significant.genes <- rep(1, times = length(data$Annotation))
names(significant.genes) <- data$Annotation



## Background
background.genes <- rep(0, times = length(expression.df$Annotation))
names(background.genes) <-  expression.df$Annotation


## Enriquecimiento
source("GSEA.R")
GeneSetEnrichmentAnalysis()
hist(as.numeric(cast(nGO[,1:2], . ~ GID)))
