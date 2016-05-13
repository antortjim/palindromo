## Fija el espacio de trabajo
setwd("~/MEGA/CuartoCurso/TFG/Bioinformatica")

##Carga de paquetes necesarios (no especificos de Nostoc)
library("ggplot2")
library("reshape")
library("AnnotationForge")
library("Biostrings")
library("openxlsx")
library("BSgenome.PCC7120.Ensembl.29")
library(org.NPCC7120.eg.db)
library("seqinr")
library("dplyr")

bsgenome <- PCC7120_Ensembl_v29
upstream.distance <- 100
downstream.distance <- 49

library(devtools)
source_gist(4676064)


annotation.file <- "Nostoc_sp_pcc_7120.GCA_000009705.1.29.gtf"
tss.df <- read.xlsx("sd01.xlsx", sheet = 1, startRow = 2, colNames = TRUE)
chr.source <- tss.df[["Source"]]
chr.source[chr.source == "chr"]

seqs.id <- tss.df %>% dplyr::select(Source, TSS, Annotation, TSS.class)
seqs.id <- with(seqs.id, paste0(Source, "_", TSS, "_", Annotation, "_", TSS.class))


unlist(strsplit(getSeq(bsgenome), split = ""))

tss.coordinates <- tss.df$TSS
tss.source <- tss.df$Source
tss.source[tss.source == "chr"]     <- 1
tss.source[tss.source == "alpha"]   <- 2
tss.source[tss.source == "beta"]    <- 3
tss.source[tss.source == "gamma"]   <- 4
tss.source[tss.source == "delta"]   <- 5
tss.source[tss.source == "epsilon"] <- 6
tss.source[tss.source == "zeta"]    <- 7
tss.source <- as.numeric(tss.source)

genome <- getSeq(bsgenome)
genome <- lapply(genome, function(x) strsplit(as.character(x), split = ""))

seqs.list <- list()

for(i in 1:nrow(tss.df))
{
  print(i)
  current.source <- genome[[tss.source[i]]][[1]]
  current.coordinate <- tss.coordinates[i]
  window.start <- current.coordinate - upstream.distance
  window.stop  <- current.coordinate + downstream.distance
  window <- window.start:window.stop
  
  if(window.start < 0)
  {
    window.start <- window.start + length(current.source)
    window <- c(window.start:length(current.source), 1:window.stop)
  }
  else if(window.stop < window.start)
  {
   window <- c(window.start:length(current.source), 1:window.stop)
  }
    seqs.list[i]  <- paste(current.source[window], collapse = "")
}

write.fasta(file.out = "Anabaena_total_150_promotor.fasta", sequences = seqs.list, names = seqs.id, nbchar = 50, as.string = T)

