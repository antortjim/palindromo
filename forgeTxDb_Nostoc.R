setwd("~/MEGAsync/CuartoCurso/TFG/Bioinformatica/")
## Generamos el paquete TxDb de Nostoc PCC 7120
## a partir de su gtf (Ensembl) con la funcion
## makeTxDbFromGFF

## Genera las 3 columnas que necesita el data.frame de chrominfo al generar
## el txdb con la funcion makeTxDbFomGFF (justo abajo)
## Requiere la presencia del fichero gtf
library("GenomicFeatures")

## Este source tarda mas o menos 1 minuto.
## Modifica el GTF de Ensembl para introducir los TSS determinados por Alicia
## De momento introduce 1/3 de todos ellos (no esta mal). Solo hace falta una vez
#source("GTFtransformer.R")

nostoc.chromosomes <- c("Chromosome", paste("pCC7120", c("alpha","beta","delta","epsilon","gamma", "zeta"), sep = ""))
nostoc.lengths <- as.numeric(seqlengths(PCC7120_Ensembl_v29))
nostoc.circularity <- rep(T, 7)
nostoc.chrominfo <- data.frame(chrom = nostoc.chromosomes, length = nostoc.lengths, is_circular = nostoc.circularity)


txdb <- GenomicFeatures::makeTxDbFromGFF(file = "Nostoc_sp_pcc_7120.v29.AOJ.gtf",
                                         chrominfo = nostoc.chrominfo, dbxrefTag = T, format = "gtf")

source("loadGTF_Nostoc.R")
