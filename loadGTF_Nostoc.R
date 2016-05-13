## Carga del gtf de Nostoc PCC7120 en R con el paquete refGenome
## El gtf de Ensembl tiene un problema: considera los TSS = inicio traducción
## Es decir, no incluye utr. Todo el transcrito se traduce, y no es verdad.
## Asi, las regiones 5'UTR se quedan como la región upstream justo al TSS, es decir, el final del promotor

## Fijamos el espacio de trabajo
setwd("~/MEGAsync/CuartoCurso/TFG/Bioinformatica/")

## Lectura del archivo .gtf de PCC7120 y almacenamiento en nostoc.gtf
annotation.file <- "Nostoc_sp_pcc_7120.v29.AOJ.gtf"
library("refGenome")                              #Cargamos el paquete refGenome
nostoc.gtf <- ensemblGenome()                     #Creamos un objeto refGenome vacío con 
#la siguiente instrucción
read.gtf(nostoc.gtf, annotation.file)             #Asignamos nuestro fichero gtf a ese objeto vacío.
nostoc.gtf <- getGtf(nostoc.gtf)

gene2tx <- nostoc.gtf[c("gene_name","gene_id")]
gene2tx <- unique(gene2tx)
rownames(gene2tx) <- NULL
