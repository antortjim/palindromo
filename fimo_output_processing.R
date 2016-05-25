#Load packages and set working space
setwd("~/MEGA/CuartoCurso/TFG/Bioinformatica/")
library("Biostrings")
library("ggplot2")
library("dplyr")
library(devtools)
library("seqinr")
library("seqLogo")
source_gist(4676064)


# Uncomment when script is completed
 # arguments <- as.integer(commandArgs(trailingOnly = T))
 # max.mismatch <- arguments[1]
 # min.mismatch <- arguments[2]
 # threshold <- arguments[3]
 # upstream.distance.prom <- arguments[4]


# if(interactive())
#{
  max.mismatch <- 3
  min.mismatch <- 0
  threshold <- 3
  upstream.distance.prom <- 100
#}


# Load transcriptomic data
commandArgs <- function() threshold
source("rna-seq-accession.R")


 
##IMPORTANTE: DA ERROR SI FIMO CONTIENE INSTANCIAS DE TSS EN PLASMIDOS Y NO EN EL CROMOSOMA PRINCIPAL
##El error consiste en que quiere añadr columnas al data.frame relevant.tss.expression con mas
## filas que las que tiene la propia df

## This function organizes the fimo output by number of mismatches, and adds a significant expression tag based on threshold value
comparator <- function(DNA.pattern, DNA.subject, mismatch, hits.df, expression, threshold = 3)
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
##IMPORTANT. The leading hash on the first line must be removed manually before introducing it into R!!!!!!!
##If you dont the script breaks because names are used to acces data and they are ignored by R if the hash remains
answer <- readline("Have you removed the leading hash from fimo.txt? If you haven't, the script won't work. Input Y if you have: \n")
if(answer != "Y")
{
  break
}

##Prepare the data to be organized in the for loop
# Read
fimo <- read.table("fimo.txt", header = T, sep = "\t", quote = "")

# #Redefine data to include all1872. YOU SHOULD WRITE THE CODE SO THAT all1873  COMES ALREADY
# all1873 <- data.frame(pattern.name = 1, sequence.name = "chr_2236056_all1873_gTSS", start = )
# head(fimo)
# fimo <- rbind(all1873, data)


hits.df <- fimo
DNA.pattern <- DNAS.palindromo

#Parse tss metadata structured by _ (Source_TSS_Annotation_TSS.class)
fimo.tss <- fimo[["sequence.name"]] %>% as.character() %>% strsplit(split = "_") %>% as.data.frame()
colnames(fimo.tss) <- c("Source","TSS","Annotation","TSS.class")

#Make match.sequence a character vector where every element i is the ith sequence detected
match.sequence        <- dplyr::select(fimo, matched.sequence)[,1] %>% as.character()
names(match.sequence) <- as.character(fimo[["sequence.name"]])

#Transform match.sequence into Biostrings data
DNAS.match.sequence <- DNAStringSet(match.sequence)
DNA.subject <- DNAS.match.sequence

## Initialize the data.frame that will store extracted data
data <- data.frame()

## For every possible mismatch, extract data and add it to that of all other mismatches
## Comparator works row-based. Comparator must be executed once per row. Each iteration adds a new row corresponding
## to a differnet hit in fimo.txt
for(mismatch in min.mismatch:max.mismatch)
{
  print(mismatch)
  current.mismatch <- comparator(DNAS.palindromo, DNAS.match.sequence, mismatch, hits.df = fimo, expression = expression.df[,-(3:4)], threshold)
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

# Machaca el data.frame sequences antiguo para que recoja la reordenacion y la eliminacion de duplicados
sequences <- select(data.no.duplicates, X1:X10)
n.mismatch <- data.no.duplicates$n.mismatch %>% as.character() %>% as.numeric()


#Generate the seqlogo image 
pdf(file = "~/MEGA/CuartoCurso/TFG/Documento/TFG/Figures/seqlogo_corrected.pdf", height = 5)
pwm <- alignment.analyzer(sequences[n.mismatch < 3,]) %>% makePWM() %>% seqLogo(xaxis = F, yaxis = F)
dev.off()
system("inkscape -f ~/MEGA/CuartoCurso/TFG/Documento/TFG/Figures/seqlogo_corrected.pdf -A ~/MEGA/CuartoCurso/TFG/Documento/TFG/Figures/seqlogo_corrected.pdf")


# Run tableR.R, which generates xlsx and pdf files formatted as a beautiful table with the features defined in select_
tabular.data <- select_(data.no.duplicates, "TSS", "Strand", "Annotation", "TSS.class", "Distance.to.TSS.pretty", "fold.change.WT", "X1:X10", "n.mismatch")
commandArgs <- function() list(tabular.data, ntss = 25)
source("palindromo/tableR.R")

seqs <- tabular.data[,7:16]
n.mismatch <- tabular.data[,17] %>% as.character %>% as.numeric
source("palindromo/consenso.R")
consensus.displayer(seqs = seqs, n.mismatch = n.mismatch, mismatch.threshold = 3)
                               
colnames(tabular.data)[17] <- "N"

# query.genes <- tabular.data %>%
#   filter(N < 3) %>%
#   select(Gen) %>% distinct()
# query.genes <- query.genes[,1]

query.genes <- c("alr2431","asr3168","alr0296")
blast <- readline('Would you like to blast interesting genes? Yes/No:\n')
if(blast == "Yes")
{
  #query.genes <- readline("Enter interesting genes: ") %>% strsplit(split = " ") %>% unlist()
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
ggplot(data = data, aes(x = start - upstream.distance.prom)) + geom_histogram(binwidth = 5)
ggsave("plot2.pdf")
ggplot(data = data, aes(x = start - upstream.distance.prom, fill = significant)) + geom_histogram()
ggsave("plot3.pdf")
ggplot(data = data, aes(x = start - upstream.distance.prom)) + geom_density(size = 1.5)
ggsave("plot4.pdf")

ggplot(data = data, aes(x = start - upstream.distance.prom, col = significant)) + geom_density(size = 1.5) +
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

save.image("fimo.RData")
