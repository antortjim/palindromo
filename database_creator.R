## Fija el espacio de trabajo
setwd("~/MEGA/CuartoCurso/TFG/Bioinformatica")

##Carga de paquetes
library("AnnotationForge")
library("Biostrings")
library("openxlsx")
library("BSgenome.PCC7120.Ensembl.29")
library("org.NPCC7120.eg.db")
library("seqinr")
library("dplyr")
library("devtools")
source_gist(4676064)
upstream.distance <- 100
downstream.distance <- 49

## Introduce las secuencias genomicas
bsgenome <- PCC7120_Ensembl_v29
genome <- getSeq(bsgenome)
genome <- lapply(genome, function(x) strsplit(as.character(x), split = ""))

## Introduce las coordenadas de los TSS
tss.df <- read.xlsx("sd01.xlsx", sheet = 1, startRow = 2, colNames = TRUE)

## Construye el identificador de todas las secuencias (seqs.id)
seqs.id <- tss.df %>% dplyr::select(Source, TSS, Annotation, TSS.class)
seqs.id <- with(seqs.id, paste0(Source, "_", TSS, "_", Annotation, "_", TSS.class))



## Extrae las coordenadas dentro de la mega data.frame con todos los datos del PNAS
## para facilitar su acceso en cada iteracion del bucle for
tss.coordinates <- tss.df$TSS

## Codifica los cromosomas bacterianos con un numero
tss.source <- tss.df$Source
tss.source[tss.source == "chr"]     <- 1
tss.source[tss.source == "alpha"]   <- 2
tss.source[tss.source == "beta"]    <- 3
tss.source[tss.source == "gamma"]   <- 4
tss.source[tss.source == "delta"]   <- 5
tss.source[tss.source == "epsilon"] <- 6
tss.source[tss.source == "zeta"]    <- 7
tss.source <- as.numeric(tss.source)


## Inicializa una lista que guarda en cada elemento una secuencia de 150 pb
seqs.list <- list()

## Para cada tss
for(i in 1:nrow(tss.df))
{
  print(i)
  # Averigua en que cromosoma/plasmido esta el tss actual
  current.source <- genome[[tss.source[i]]][[1]]
  # Implementar en un futuro el soporte a los plasmidos
  if(current.source != 1)
  {
    next
  }
  ## Extrae la coordenada del tss
  current.coordinate <- tss.coordinates[i]
  
  ## Define el window.start y el window.stop como el intervalo -upstream.distance +downstream.distance
  ## respecto a la coordenada del TSS
  window.start <- current.coordinate - upstream.distance
  window.stop  <- current.coordinate + downstream.distance
  window <- window.start:window.stop
 ## Estos if y else if implementan la circularidad para solventar el problema de un window.start negativo
 ## o un window.stop más pequeño que el window.start
  if(window.start < 0)
  {
    window.start <- window.start + length(current.source)
    window <- c(window.start:length(current.source), 1:window.stop)
  }
  else if(window.stop < window.start)
  {
   window <- c(window.start:length(current.source), 1:window.stop)
  }
  
  ## Accede a la secuencia, conviertela en una cadena, y guardala en la posicion correspondiente
  ## de la lista de secuencias
  seqs.list[i]  <- paste(current.source[window], collapse = "")
}


# Exporta todo a un archivo fasta
write.fasta(file.out = "Anabaena_150_promotor.fasta", sequences = seqs.list, names = seqs.id, nbchar = 50, as.string = T)

