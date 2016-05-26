#Load packages and set working space
setwd("~/MEGA/CuartoCurso/TFG/Bioinformatica/")


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
source("functions.R")

#Define the palindromic sequence
"TACCCGGGTA" %>% DNAString -> DNAS.palindromo 

#Parse tss metadata structured by _ (Source_TSS_Annotation_TSS.class)
seqId.fields <- c("Source","TSS","Annotation","TSS.class")

## Initialize the data.frame that will store extracted data
data <- data.frame()
## For every possible mismatch, extract data and add it to that of all other mismatches
## fimo.output.integrator works row-based. fimo.output.integrator must be executed once per row. Each iteration adds a new row corresponding
## to a differnet hit in fimo.txt
for(mismatch in min.mismatch:max.mismatch)
{
  print(mismatch)
  current.mismatch <- fimo.output.integrator(fimo.file = "fimo.txt", seqId.fields = seqId.fields, DNAS.palindromo, mismatch, expression = expression.df, threshold)
  data <- rbind(data, current.mismatch)
}

data <- addDistanceAsCharacter(data)

sequences <- data$matched.sequence %>% as.list() %>% lapply(s2c) %>% as.data.frame.list

#Generate the seqlogo image 
pdf(file = "~/MEGA/CuartoCurso/TFG/Documento/TFG/Figures/seqlogo_corrected.pdf", height = 5)
pwm <- alignment.analyzer(sequences[data$n.mismatch < 3,]) %>% makePWM() %>% seqLogo(xaxis = F, yaxis = F)
dev.off()
system("inkscape -f ~/MEGA/CuartoCurso/TFG/Documento/TFG/Figures/seqlogo.pdf -A ~/MEGA/CuartoCurso/TFG/Documento/TFG/Figures/seqlogo.pdf")



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



## Make the number of mismatches a factor (n.mismatch) so that ggplot can handle it
data$n.mismatch <- factor(data$n.mismatch, levels = min.mismatch:max.mismatch)

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


# ## Significant genes
# significant.genes <- rep(1, times = length(data$Annotation))
# names(significant.genes) <- data$Annotation
# 
# 
# 
# ## Background
# background.genes <- rep(0, times = length(expression.df$Annotation))
# names(background.genes) <-  expression.df$Annotation
# 
# 
# ## Enriquecimiento
# source("GSEA.R")
# GeneSetEnrichmentAnalysis()
# hist(as.numeric(cast(nGO[,1:2], . ~ GID)))

save.image("fimo.RData")
