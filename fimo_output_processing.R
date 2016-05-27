#Load packages and set working space
setwd("~/MEGA/CuartoCurso/TFG/Bioinformatica/")
library("reshape")

# Uncomment when script is completed
 # arguments <- as.integer(commandArgs(trailingOnly = T))
 # max.mismatch <- arguments[1]
 # min.mismatch <- arguments[2]
 # threshold <- arguments[3]
 # upstream.distance.prom <- arguments[4]


# if(interactive())
#{
  max.mismatch <- 2
  min.mismatch <- 0
  threshold <- 3
  upstream.distance.prom <- 100
#}


# Load transcriptomic data
commandArgs <- function() threshold
source("rna-seq-accession.R")
source("functions.R")

#Define the palindromic sequence
palindromo <- "TACCCGGGTA"
palindromo %>% DNAString -> DNAS.palindromo 


write(palindromo, file = "sites_files/palindromo.txt")
## Create a random population of motifs
palindromo %>% s2c -> palindromo.char
num.ran <- 5000


random.motifs <- character(length = num.ran)
random.motifs.names <- character(length = num.ran)
for(i in 1:num.ran)
{
  motif.name <- paste("random_", i, sep = "")
  motif.file <- paste("sites_files/random_", i, ".txt", sep = "")
  random.motif <- sample(palindromo.char, 10) %>% c2s
  random.motifs[i] <- random.motif
  write(random.motif, file = motif.file)
}

my.motifs.seqs <- c(palindromo, random.motifs)



system("~/meme/bin/sites2meme sites_files > motifs.meme")
system("~/meme/bin/fimo --oc . --verbosity 1 --thresh 0.01 -o fimo_out/ --norc motifs.meme Anabaena_150_promotor.fasta")
##Borra la hash que inicia la primera linea (el comando borra todas)
system("sed -i.backup 's/^#//' fimo_out/fimo.txt")


#Parse tss metadata structured by _ (Source_TSS_Annotation_TSS.class)
seqId.fields <- c("Source","TSS","Annotation","TSS.class")


## Initialize the data.frame that will store extracted data
data <- data.frame()
## For every possible mismatch, extract data and add it to that of all other mismatches
## fimo.output.integrator works row-based. fimo.output.integrator must be executed once per row. Each iteration adds a new row corresponding
## to a differnet hit in fimo.txt

#Lee la tabla y guardala en la variable fimo
fimo <- read.table("fimo_out/fimo.txt", header = T, sep = "\t", quote = "")
my.motifs.names <- c("palindromo", paste("random_", 1:num.ran, sep = ""))

for(i in 1:length(my.motifs.names))
{
  pattern.id <- i -1
  print(my.motifs.names[i])
  mysubset <- filter(fimo, pattern.name == my.motifs.names[i])
  DNA.pattern <- my.motifs.seqs[i] %>% DNAString()
for(mismatch in min.mismatch:max.mismatch)
{
  print(mismatch)
  current.mismatch <- fimo.output.integrator(fimo.df = mysubset, seqId.fields = seqId.fields, DNA.pattern, mismatch,
                                             expression = expression.df, threshold, pattern.id = pattern.id)
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



# aggregate(n.mismatch ~ pattern.id, FUN=table, data=data)
#  
# aggregate(n.mismatch ~ pattern.id, FUN=sum, data=data)
# 
# with(data, tapply(n.mismatch, pattern.id, FUN=table)) %>% as.data.frame()

library(doBy)
mydata <- summaryBy(n.mismatch ~ pattern.id, data = data, 
          FUN = function(x) { c(t = table(x)) } )


colnames(mydata) <- c("pal","zero","one","two")
mydata[,4] %>% which.max()

my.motifs.seqs[442]
ggplot(data = mydata) + geom_histogram(aes(x=zero))
ggplot(data = mydata) + geom_histogram(aes(x=one))
ggplot(data = mydata) + geom_histogram(aes(x=two))


filter(data, pattern.id == 3) %>% filter(n.mismatch == 2)
# produces mpg.m wt.m mpg.s wt.s for each 
# combination of the levels of cyl and vs



data$n.mismatch %>% table



write.table(data, file = "fimo_data_and_expresion.txt", col.names = T, row.names = F)


#Generate the seqlogo image 
seqlogo.file <- "~/MEGA/CuartoCurso/TFG/Documento/TFG/Figures/seqlogo_2.pdf"
pdf(file = seqlogo.file, height = 5)
pwm <- alignment.analyzer(data[data$n.mismatch < 3,] %>% dplyr::select(X1:X10)) %>% makePWM() %>% seqLogo(xaxis = F, yaxis = F)
dev.off()
system(paste("inkscape -f ", seqlogo.file, " -A ", seqlogo.file, sep = ""))
system(paste("evince ", seqlogo.file, sep = ""))


# Run tableR.R, which generates xlsx and pdf files formatted as a beautiful table with the features defined in select_
tabular.data <- select_(data, "TSS", "Strand", "Annotation", "TSS.class", "pretty.distance", "fold.change.WT", "X1:X10", "n.mismatch")
commandArgs <- function() list(tabular.data, ntss = 14)
source("palindromo/tableR.R")

Uno <- consensus.displayer(seqs =  dplyr::select(data, X1:X10, n.mismatch), mismatch.threshold = 2)
Dos <- consensus.displayer(seqs =  dplyr::select(data, X1:X10, n.mismatch), mismatch.threshold = 3)
consenso <- data.frame(Uno, Dos, letra = "TACCCGGGTA" %>% s2c %>% as.factor)

for(i in 1:nrow(consenso))
{
 if(consenso$Uno[i] - consenso$Dos[i] > 0)
 {
  consenso$Uno[i] <- consenso$Uno[i] - consenso$Dos[i]
 }
 else
 {
   print(i)
   consenso$Uno[i] <- 0
 }
}
consenso <- melt(consenso, variable_name = "consenso", id.vars = "letra")
consenso$consenso <- factor(consenso$consenso, levels = c("Dos", "Uno"))
ggplot(data = arrange(consenso, consenso, rep(c("Dos", "Uno"), each = 10)) , 
       aes(x = rep(1:10,times = 2), y = value, fill = consenso)) + 
  geom_bar(stat = "identity", position = "stack", width = 1) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic() +
  theme(axis.text.y = element_blank()) +
  scale_y_continuous(breaks=NULL, name = NULL) + 
  scale_x_continuous(breaks=NULL, name = NULL)

ggsave("../Documento/TFG/Figures/consenso.pdf")  
  

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

n = 30
cols = gg_color_hue(n)
dev.new(width = 4, height = 4)
plot(1:n, pch = 16, cex = 2, col = cols)
cols


query.genes <- c("alr2431","asr3168","alr0296")
blast <- readline('Would you like to blast interesting genes? Yes/No:\n')
if(blast == "Yes")
{
  #query.genes <- readline("Enter interesting genes: ") %>% strsplit(split = " ") %>% unlist()
  commandArgs <- function() list(query.genes)
  source("palindromo/blast.R")
}

##Porcentaje de TSS con un palindromo de 1 mismatch en torno a el
100*nrow(data)/nrow(expression.df)

##Visualiza la distribucion de los tipos de TSS en el fondo y en los tss de interes
plot(table(fimo.tss[["TSS.class"]]))
plot(table(data[["TSS.class"]]))



## Make the number of mismatches a factor (n.mismatch) so that ggplot can handle it
data$n.mismatch <- factor(data$n.mismatch, levels = min.mismatch:max.mismatch)
data$significant <- factor(data$significant, levels = "no, increase, decrease" %>% asCharacter)

mydata <- data
## Graphics production
ggplot(filter(data, n.mismatch < 3)) + geom_point(aes(x = WT0, y = WT8))


       
ggplot(data = data, aes(x = log10(WT0), y = log10(WT8), col = significant, shape = n.mismatch)) + 
  geom_point() +
  scale_shape_manual(values = c(17, 19, 1, 3))

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
