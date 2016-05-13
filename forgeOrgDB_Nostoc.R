## Este script genera el paquete de anotaciones OrgDB para Anabaena sp. PCC 7120
## Para ello accede a 

## Requiere el paquete AnnotationForge
# source("http://bioconductor.org/biocLite.R")
# biocLite("AnnotationForge")
library("AnnotationForge")

## Fijamos el espacio de trabajo
setwd("~/MEGAsync/CuartoCurso/TFG/Bioinformatica/")

## Carga el fichero gtf
source("loadGTF_Nostoc.R")

## Quedate con las columnas gene_name y gene_id para obtener un diccionario
## que asocia a cada nombre de gen (all0002) el de su transcrito (BAB77526)
gene2tx <- nostoc.gtf[c("gene_name","gene_id")]
## Elimina los duplicados (generados por las distintas anotaciones presentes el el gtf
## para el mismo gen: TSS, ATG, terminador, etc)
gene2tx <- unique(gene2tx)


## Now prepare the data.frames that makeOrgPackage receives as mandatory arguments
## nSym es un marco de datos que empareja el nombre de los genes con su función, de ser conocida
## Esta información la vamos a estraer del fichero Nostoc_genes.txt, descargado de
## Cyanobase en la url ()

## Leemos el archivo y lo guardamos en la variable nostoc
nostocFile <- "Nostoc_genes.txt"
cat(head(readLines(nostocFile)), sep = "\n")
nostoc <- read.table(nostocFile, sep="\t",  comment.char = "", quote = "")

## Vamos a generar 3 data.frames nSym, nChr y nGO que makeOrgPackage demanda:

## nSym guarda la correlacion entre id de los genes (GID), su funcion (GENENAME),
## y el nombre del transcrito (TXID)
 
nSym <- nostoc[,c(2,7)]
nSym <- nSym[nSym[,1]!="-",]
nSym <- nSym[nSym[,2]!="-",]
nSym[,1] <- as.character(nSym[,1])
nSym[,2] <- as.character(nSym[,2])

# ## Algunos de los genes no estan guardados con la nomenclatura sistemática
# ## sino con un nombre que informa sobre su funcion
# OK.positions <- as.numeric(grep("all", gene2tx[,1]))
# OK.positions <-  c(OK.positions, as.numeric(grep("alr", gene2tx[,1])))
# OK.positions <-  c(OK.positions, as.numeric(grep("asl", gene2tx[,1])))
# OK.positions <-  c(OK.positions, as.numeric(grep("asr", gene2tx[,1])))
# 
# nrow(gene2tx[!1:nrow(gene2tx) %in% OK.positions,])
# sort(OK.positions)

row <- 1
## Este bucle incorpora a nSym el nombre de los transcritos asociados a cada gen
for (row in 1:nrow(nSym))
{
  print(row)
  current.gene <- nSym[row,1]
  
  ## Con esto estamos añadiendo a cada gen (en cada fila de nSym), el nombre de su transcrito
  ## Si el nombre del gen esta guardado en la notacion sistematica (all0001) o un nombre symbol (AtpC)
  if(length(unlist(strsplit(current.gene, split = " "))) == 1)
    ## Compara gene2tx con current.gene y eso devuelve un logico con la posicion
    ## de la fila de gene2tx en la que el current.gene esta guardado.
    ## Accede a la segunda columna de esa fila en el diccionario gene2tx y otorga
    ## ese valor a la fila correspondiente de nSym
    {
    nSym[row,3] <- gene2tx[which(gene2tx[["gene_name"]] == current.gene), 2]
    }
  ## Si no, hay que encontrar alguno forma de hacer en entren
}

## Renombra las columnas tal y como la funcion makeOrgPackage lo demanda
colnames(nSym) <- c("GID", "GENENAME", "TXID")

## nChr guarda la correlacion entre id de los genes (GID), y el cromosoma en el que están (CHROMOSOME)
## y el nombre del transcrito (TXID)

nChr <- nostoc[,c(2,3)]
nChr <- nChr[nChr[,2]!="-",]
nChr[,1] <- as.character(nChr[,1])
nChr[,2] <- as.character(nChr[,2])
colnames(nChr) <- c("GID","CHROMOSOME")

## nGO guarda la conexion entre GID y términos de ontología de genes asociados (GO). Ademas
## hay que añadir una columna EVIDENCE (en todos goterm, da igual lo que ponga)

nostocGOFile <- "Nostoc_GO.txt"
genes.id <- read.table(nostocGOFile, sep="\t",  comment.char = "", quote = "")[,2]
GO.ids <- readLines("GO_ids")
GO.ids <- paste("GO", GO.ids, sep = ":")
evindence <- read.table(nostocGOFile, sep="\t",  comment.char = "", quote = "")[,3]
nGO <- data.frame(genes.id, GO.ids, evindence)
nGO[,1] <- as.character(nGO[,1])
nGO[,2] <- as.character(nGO[,2])
nGO[,3] <- as.character(nGO[,3])
colnames(nGO) <- c("GID", "GO", "EVIDENCE")

## Ya tenemos las 3 data.frames necesarias, ahora llamamos a la funcion, que genera
## una carpeta org.NPCC7120.eg.db que es compilable con R CMD
makeOrgPackage(gene_info=nSym, chromosome=nChr, go=nGO,
               version="0.1",
               author = "Antonio Ortega <antortjim@alum.us.es>",
               maintainer = "Antonio Ortega <antortjim@alum.us.es>",
               outputDir = ".",
               tax_id = "103690",
               genus = "Nostoc",
               species = "PCC7120",
               goTable="go")

## Compila el paquete (R CMD build)
system("R CMD build org.NPCC7120.eg.db")
## Instala el paquete compilado
system("R CMD INSTALL org.NPCC7120.eg.db_0.1.tar.gz")