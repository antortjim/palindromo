#Load packages and set working space
setwd("~/MEGA/CuartoCurso/TFG/Bioinformatica/")
library("Biostrings")
library("dplyr")
library(devtools)
source_gist(4676064)
source("rna-seq-accession.R")


#Define the palindromic sequence
palindromo <- "TACCCGGGTA"   # palindromo presente en el promotor de all1873
# muy conservado en todas las cianos
DNAS.palindromo <-  DNAString(palindromo)  # pasamos la secuencia del palindromo
# a DNAString (Biostrings package)


#Incorporate FIMO output as data.frame
fimo <- read.table("fimo.txt", header = T, sep = "\t", quote = "")
#Parse tss metadata structured by _
fimo.tss <- fimo[["sequence.name"]] %>% as.character() %>% strsplit(split = "_") %>% as.data.frame()
colnames(fimo.tss) <- c("source","coordinate","gene","TSS.type")

#Make match.sequence a character vector where every element i is the ith sequence detected
match.sequence <-dplyr::select(fimo, matched.sequence)[,1] %>% as.character()

names(match.sequence) <- as.character(fimo[["sequence.name"]])

#Transform it into Biostrings data
DNAS.match.sequence <- DNAStringSet(match.sequence)


oneMatch <- unlist(vmatchPattern(pattern = DNAS.palindromo,
                                 subject = DNAS.match.sequence,
                                 max.mismatch = 1, min.mismatch = 1)) %>% 
   names() %>% 
   strsplit(split = "_") %>% 
   as.data.frame.list()

colnames(oneMatch) <- c("source","coordinate","gene","TSS.type")

relevant.tss <- with(oneMatch, paste0(source, "_", coordinate, "_", gene, "_", TSS.type))

oneMatch <- subset(fimo, fimo$sequence.name %in% relevant.tss)
oneMatch$start - upstream.distance

ggplot(data = oneMatch, aes(x = start - upstream.distance)) + geom_histogram(binwidth = 5)
ggplot(data = oneMatch, aes(x = start - upstream.distance)) + geom_density()


fimo.tss[["coordinate"]] %>% as.character() %>% as.numeric()

oneMatch[,1] %>% as.character() %>% as.numeric() %>% unique()
interesting.genes <- oneMatch[,2] %>% as.character() %>% unique()




##Porcentaje de TSS con un palindromo de 1 mismatch en torno a el
100*nrow(oneMatch)/40540


plot(table(fimo.tss[,3]))
plot(table(oneMatch[,3]))















positions <- expression.df[["TSS"]] %>% as.character() %in% as.character(oneMatch[,2])
one.mismatch.subset <- subset(expression.df, subset = positions)

comparison.df <- dplyr::filter(one.mismatch.subset, genotype == "WT") %>% dplyr::select(-genotype)
WT0 <- comparison.df[comparison.df[["time"]] == 0, "score"]
WT8 <- comparison.df[comparison.df[["time"]] == 8, "score"]
comparison.df <- data.frame(seqs.ids, WT0, WT8) 
comparison.df <- na.omit(comparison.df)

number.genes <- nrow(comparison.df)/2

ratiowt8    <- numeric(length=number.genes)
ratiohetr0  <- numeric(length=number.genes)
ratiohetr8  <- numeric(length=number.genes)

index <- 1
for (i in 4*(1:number.genes)-3)
{
  print(index)
  current.ratiowt8   <- expression.df[i+1,4] / expression.df[i,4]
  current.ratiohetr0 <- expression.df[i+2,4] / expression.df[i,4]
  current.ratiohetr8 <- expression.df[i+3,4] / expression.df[i,4]
  ratiowt8[index] <- current.ratiowt8
  ratiohetr0[index] <- current.ratiohetr0
  ratiohetr8[index] <- current.ratiohetr8
  index <- index +1
}
#unos segundos y listo

ratio.matrix <- matrix(c(ratiowt8, ratiohetr0, ratiohetr8), nrow = 3, byrow = T)
colnames(ratio.matrix) <- ids

ratio.df <- data.frame(rep(ids, each = 3), as.numeric(ratio.matrix))

ggplot(data = comparison.df) + geom_point(aes(x = log10(WT0), y = log10(WT8)))


WT8.fold.change <- comparison.df[["WT8"]] / comparison.df[["WT0"]]
comparison.df[["fold_change"]] <- WT8.fold.change
comparison.df[["significant"]] <- ifelse(WT8.fold.change > 3 | WT8.fold.change < 1/3, "yes", "no" )

p <- ggplot(data = comparison.df)
p <- p + geom_point(aes(x = log10(WT0), y = log10(WT8), colour = significant))
#p <- p + labs(x = NULL, y = NULL, title = "RNA-seq de la transición a ausencia de nitrógeno 8h")
p <- p + labs(y = expression(log[1*0]*"(8h)"), x = expression(log[1*0]*"(0h)"))
p <- p + geom_text(aes(x = log10(2198), y = log10(21648) + 0.15, label = "all1873"))
p <- p + guides(colour = FALSE)
p <- p + scale_color_manual(values = c("#000099", "#E69F00"))
p <- p + theme_set(theme_bw(base_size = 15))
p <- p + theme(plot.title = element_text(hjust = 0)) 
p


