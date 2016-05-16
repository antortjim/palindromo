## Este script proporciona el codigo necesario para realizar enriquecimiento
## de ontologia de terminos GO en Nostoc PCC 7120
## Fija el espacio de trabajo
setwd("~/MEGA/CuartoCurso/TFG/Bioinformatica")
source("http://bioconductor.org/biocLite.R")
biocLite("topGO")

# Carga el paquete topGO que proporciona el codigo base necesario
library("topGO")
library("GO.db")

# ## Macro que indica el umbral de puntuacion de alineamiento
# threshold <- commandArgs()[1]
# print(threshold)

## Define la funcion de selección de genes de interes dentro del background
topDiffGene <- function(data)
{
  return(data$Annotation)
}


# printGenes <- function(geneID2GO, goID, gene.score, threshold = NULL)
# {
#  ##genera una lista en la que cada elemento es un vector logico con algun T si el gen
#  ## al que corresponde tiene la GO indicada
#  mylist <- lapply(geneID2GO, function(x) x == goID)
#  ##genera una lista en la que cada elemento es un 1 o un 0 segun si el gen al que
#  ## corresponde tiene la GO indicada
#  mylist <- lapply(mylist, function(x) sum(x))
#  ##genera una character en el que cada elemento es la posicion de un gen con goID en
#  ##la lista geneID2GO y su name es el nombre del gen
#  mygenes <- which(unlist(lapply(mylist, function(x) x == 1)))
#  ##este bucle while permite epscar genes en categorias mas especificas
#  ##si la categoria actual (mas generica) no tiene ningun gen asociado
#  while(length(mygenes) == 0)
#  {
#    
#    goID <- xx[[goID]]
#    mygenes <- character(0)
#    for (i in 1:length(goID))
#    {
#    current.goID <- goID[i]
#    mylist <- lapply(geneID2GO, function(x) x == current.goID)
#    ##genera una lista en la que cada elemento es un 1 o un 0 segun si el gen al que
#    ## corresponde tiene la GO indicada
#    mylist <- lapply(mylist, function(x) sum(x))
#    ##genera una character en el que cada elemento es la posicion de un gen con goID en
#    ##la lista geneID2GO y su name es el nombre del gen
#    mygenes <- c(mygenes, which(unlist(lapply(mylist, function(x) x == 1))))
#    }
#    goID <- current.goID
#  }
# 
#  names(mygenes)[names(mygenes) == "all1291"] <- "cynS"
#  if(is.null(threshold))
#  { 
#    return(list(goID = goID, genes = mygenes))
#  }
#  else
#  {
#  ## extrae aquellos que ademas estan por encima del threshold
#  mygenes <- which(subset(gene.score, names(gene.score) %in%  names(mygenes)) > threshold)
#  mygenes <- gene.score[names(mygenes)]
#  ##aplica restriccion por puntuacion en el alineamiento
#  return(list(goID = goID, genes = mygenes))
#  }
# }


# Carga las relaciones de ontología de terminos de genes
# Convert the object to a list
xx <- as.list(GOMFCHILDREN)
# Remove GO identifiers that do not have any children
xx <- xx[!is.na(xx)]
if(length(xx) > 0){
  # Get the children GO identifiers for the first elents of xx
  goids <- xx[[1]]
  # Find out the GO terms for the first parent goid
  GOID(GOTERM[[goids[1]]])
  Term(GOTERM[[goids[1]]])
  Synonym(GOTERM[[goids[1]]])
  Secondary(GOTERM[[goids[1]]])
  Definition(GOTERM[[goids[1]]])
  Ontology(GOTERM[[goids[1]]])
}



## Introduce el mapping gen - termino GO para que al ver un gen veamos tambien los GOs asociados
## Esta información proviene de un archivo disponible en Cyanobase
## Importantisimoa veriguar como he obtenido el GO_terms_wide.txt
## Viene de procesar el Nostoc_GO.txt pero no se como
# Terminos GO

go.mappings <- "GO_terms_wide.txt"
## Introduce los GO.mappings en R
geneID2GO   <- readMappings(go.mappings, IDsep = ";", sep = ",")

GeneSetEnrichmentAnalysis <- function(gene.score)
{
 GOdata <- new("topGOdata",
                      description = "Simple session",
                      ontology = "MF",
                      allGenes = background.genes,
                      geneSel = topDiffGene,
                      nodeSize = 10,
                      annot = annFUN.gene2GO,
                      gene2GO = geneID2GO)
 
 resultFisher  <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
 # resultKS      <- runTest(GOdata, algorithm = "classic", statistic = "ks")
 # resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
 
 return(list(GOData = GOdata, Fisher = resultFisher,
             #classic.KS = resultKS, elim.KS = resultKS.elim
             ))
}

