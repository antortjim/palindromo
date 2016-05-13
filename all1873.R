## Fija el espacio de trabajo
setwd("~/MEGA/CuartoCurso/TFG/Bioinformatica")

##Carga de paquetes necesarios (no especificos de Nostoc)
library("ggplot2")
library("reshape")
library("dplyr")
library("AnnotationForge")
library("Biostrings")



## Implementa las funciones converter, alignment.characterization, infoDisplay, seqDisplay y seqWrite
## Facilitan el manejo de las secuencias, alineamientos y puntuaciones
converter <- function(grange)
{
  tx_name <- mcols(grange)[,2]
  return(gene2tx[gene2tx[["gene_id"]] %in% tx_name,])
}


alignment.characterization <- function(pairwiseAlignmentsResults)
{
  sorted.score <- sort(score(pairwiseAlignmentsResults))
  hist.breaks <- hist(score(pairwiseAlignmentsResults))$breaks
  plot(score(pairwiseAlignmentsResults))
  return(list(score = unique(sorted.score), breaks = hist.breaks))
}

infoDisplay <- function(threshold, bsgenome.object, promoters.grange, alignment)
{
  #guarda los indices de los alineamientos cuya puntacion supera el umbral y cuya
  ##longitud es 10
  index <- which(score(alignment) >= threshold & nchar(alignment) == 10)

  #extrae la secuencia de los promotores correspondientes
  myseqs <- getSeq(bsgenome.object, promoters.grange[index])
  myseqs <- strsplit(as.character(myseqs), split = "")
  
  #subset del grange original que solo guarda los promotores que superaron el umbral
  mygrange <- promoters.grange
  
  #pasa este grange a una data.frame en el que para cada promotor aparece el nombre de su gen y de su tx
  myref <- converter(mygrange)
  
  ## añade la puntuacion de cada gen y tambien la posicion que ocupa en el grange final
  myref[, 3:4] <- 0
  colnames(myref)[4] <- "score"
  for (i in 1:nrow(myref))
  # para cada fila de myref
  {
    print(i)
    # averigua en que posicion del grange original esta el gen que toca ahora mismo
    grange.index <- which(mcols(mygrange)[,2] == myref[i,2])
    # coloca ese numero en la tercera columna
    myref[i, 3] <- grange.index
    current.score <- score(alignment)[grange.index]
    # obten la puntuacion del alineamiento de ese gen
    #current.score <- score(alignment)[index][grange.index]
    myref[i, 4] <- current.score
  }
  return(list(alignments = alignment, reference = myref))
}



seqDisplay <- function(infoDisplay.object, gene.name)
{
 #El indice 1 indica que posicion del reference guarda informacion sobre el gen correspondiente
 index.1 <- which(infoDisplay.object$reference[["gene_name"]] == gene.name)
 #El indice 2 indica que numero del grange guarda informacion sobre el gen
 ## Ese indice es tambien la posicion del alineamiento del gen en el conjunto de alineamientos
 index.2 <- infoDisplay.object$reference[index.1,3]
 return(infoDisplay.object$alignments[index.2])
}


seqWrite <- function(infoDisplay.object, gene.name, alignment.file)
{ 
  print(gene.name)
  aln <- subject(seqDisplay(universe, gene.name))
  
  ## filtra a la secuencia mas larga en caso de que haya varias
  if (length(aln) > 1)
  {
   right.pos <- which.max(nchar(aln))
   aln <- aln[right.pos]
   seq.start <- upstream.length - start(subject(seqDisplay(universe, gene.name)))[right.pos] 
   print(-seq.start)
   id <- paste(">", gene.name, sep = "")
   #id.line <- paste(name, seq.start, sep = " ") el nombre acaba en el primer espacio tras >
   #por lo que texshade no va a imprimir el start.
   id.plus.aln <- paste(id, aln, sep="\n")
   print(id.plus.aln)
   write.table(id.plus.aln,
               file = alignment.file, row.names = F,
               quote = F, col.names = F, append = T)
   return(list(sequence.start = -seq.start, gene.name = gene.name))
   
  }
  print(aln)
  if (nchar(aln) == 9)
  {
    grange.index <- subset(universe$reference[["V3"]], universe$reference[["gene_name"]] == gene.name)
    print(grange.index)
    seq.missing.pos <- start(subject(palindromo.vs.proms[grange.index])) + 9
    missing.letter <- unlist(strsplit(as.character(getSeq(BSgenome.PCC7120.Ensembl.29, nostoc.promoters)[grange.index]), split = ""))[seq.missing.pos]
    print(missing.letter)
    aln <- paste(aln, missing.letter, sep = "")
  }
  
  seq.start <- upstream.length - start(subject(seqDisplay(universe, gene.name))) 
  print(-seq.start)
  id <- paste(">", gene.name, sep = "")
  #id.line <- paste(name, seq.start, sep = " ") el nombre acaba en el primer espacio tras >
  #por lo que texshade no va a imprimir el start.
  id.plus.aln <- paste(id, aln, sep="\n")
  print(id.plus.aln)
  write.table(id.plus.aln,
              file = alignment.file, row.names = F,
              quote = F, col.names = F, append = T)
  return(list(sequence.start = -seq.start, gene.name = gene.name))
  print(aln)
}

##Carga funciones de analisis biologico
## Forja el paquete BSgenome para Nostoc sp. PCC 7120. Solo es necesario la primera vez
#source("forgeBSgenome_Nostoc.R")
## Cargalo
library("BSgenome.PCC7120.Ensembl.29")
bsgenome <- PCC7120_Ensembl_v29

## Genera la libreria TxDb para Nostoc sp. PCC 7120. Se devuelve en la variable txdb
## No hace falta generar un paquete, se emplea directamente la salida de makeTxDbfromGFF()
## Ademas se devuelve nostoc.gtf, la lectura del gtf de Nostoc
source("forgeTxDb_Nostoc.R")


## Forja el paquete OrgDB para Nostoc sp. PCC 7120
#source("forgeOrgDB_Nostoc.R")
#Cargalo
#library("org.NPCC7120.eg.db")
#orgdb <- org.NPCC7120.eg.db

## Carga datos
data(BLOSUM62)
all1873.id <- "BAB73572"     # id del transcrito de all1873
palindromo <- "TACCCGGGTA"   # palindromo presente en el promotor de all1873
                             # muy conservado en todas las cianos
upstream.length   <- 200
downstream.length <- 100
DNAS.palindromo <-  DNAString(palindromo)  # pasamos la secuencia del palindromo
                                           # a DNAString (Biostrings package)


## Lee la tabla de terminos GO no estructurada
Nostoc.GO <- read.table("Nostoc_GO.txt", sep = "\t", quote = "")

## Extrae el nombre de los genes de Nostoc (posiblemente se puedan sacar de otro de los
## archivos que se estan importando)
Nostoc.genes <- read.table("Nostoc_genes.txt", sep = "\t", quote = "")[,-1]
colnames(Nostoc.genes) <- c("gene.name", "chromosome", "start",
                            "end","strand","mol.function")
Nostoc.genes$gene.name <- as.character(Nostoc.genes$gene.name)
Nostoc.genes$gene.name[Nostoc.genes$gene.name == "all1291"] <- "cynS"


## Analisis de regiones promotoras
## Genera un GRange que guarde los intervalos asociados a los promotores
nostoc.promoters <- promoters(txdb,
                              upstream = upstream.length,
                              downstream = downstream.length)

## Extrae las secuencias de los todos los promotores
## y pasalas al formato de Biostrings (clase DNAString)
promoterSeq <- getSeq(BSgenome.PCC7120.Ensembl.29, nostoc.promoters)
promoterSeqString <- lapply(promoterSeq, DNAString)

## Alinear directamente el palindromo de interes a todos los promotores, y ver
## qué secuencias promotoras se parecen más a la nuestra
palindromo.vs.proms <- pairwiseAlignment(pattern = rep(palindromo, length(promoterSeq)), subject = promoterSeq,
                                         type = "local", gapOpening=-16, gapExtension=-16, substitutionMatrix = BLOSUM62)





## Extrae la lista de puntuaciones en forma de data.frame de forma que cada elemento
## lleve el indice del gen
alignment.score    <- score(palindromo.vs.proms)
alignment.score.df <- as.data.frame(alignment.score)
alignment.score.df[,2] <- row.names(alignment.score.df)
colnames(alignment.score.df)[2] <- "gen"

## Fija el umbral de signifiancia en el 95%. Se empleara en la generacion de los
## graficos y en el Gene Set Enrichment Analysis
threshold <- quantile(score(palindromo.vs.proms), .95)

## Generacion de graficos de interes: Nube de puntos e histograma
## Se plotean solo 10 mil puntos para simplificar la visualizacion
## el grafico casi no cambia por eso
ggplot(data = dplyr::sample_n(alignment.score.df, size = 10000)) +
  theme_bw() +
  geom_point(aes(y = alignment.score, x = gen), size = 0.5) +
  geom_hline(aes(yintercept=threshold), col = "red") + 
  theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank())

ggsave("alignment_distribution.pdf", width = 5, height = 5)
##Guardar 5 x 5

  ggplot(data = alignment.score.df) +
  theme_bw() +
  geom_histogram(aes(x = alignment.score), binwidth = 1) +
  geom_vline(aes(xintercept=threshold), col = "red") + 
  labs(x="Puntuación", y="Frecuencia")
    
ggsave("../Documento/TFG/Figures/alignment_histogram.pdf", width = 5, height = 5)

##infoDisplay es una funcion que permite correlacionar en un solo objeto los datos guardados en
##varios objetos distintos
## Pide un umbral de puntuacion en el alineamiento, un genoma, grange que guarde los promotores y 
## el resultado de un alineamiento de pairwiseAlignment
## Para cada gen que superase el umbral muestra su nombre, el de su transcrito y la posicion que ocupa en
## el Local PairwiseAlignments object. Los transcritos que vengan del mismo gen solo aparecen una vez, por lo
## que el tamaño de reference puede ser menor que el del grange.


## Solo basado en secuencia del gen
#universe <- infoDisplay(52, BSgenome.PCC7120.Ensembl.29, nostoc.promoters, palindromo.vs.proms)

## Genera el objeto universe, que es una lista de dos elementos guardando informacion
## sobre el universo de genes de Anabaena y ordenalo
## (esto ultimo deberia implementarse en la propia infoDisplay como una opcion)
universe <- infoDisplay(0, BSgenome.PCC7120.Ensembl.29, nostoc.promoters, palindromo.vs.proms)
gene.score <- universe$reference[["score"]]
names(gene.score) <- universe$reference[["gene_name"]]
universe$reference <- arrange(universe$reference, -score)

#Analisis del objeto universe y extraccion de secuencias con maxima puntuacion de alineamiento
## Generaremos una tabla en la que se emparejan los genes y su palindromo

## Determina cuales son los 30 genes con el palindromo mas parecido y extrae la seq del alineamiento
best.genes <- as.list(universe$reference[3:32,][["gene_name"]])
alignment.file <- "../Documento/TFG/bestGenes.fa"
## Exporta cada instancia y el nombre del gen a alignment.file (fasta)
write.table(paste(">all1873", palindromo, sep="\n"),
            file = alignment.file, row.names = F,
            quote = F, col.names = F)




## Guarda las coordenadas de inicio respecto al TSS de cada instancia y exportalas a un fichero.txt 
start.positions <- lapply(best.genes, function(x) seqWrite(universe, x, alignment.file))

start.positions.file <- "../Documento/TFG/Figures/start_positions.txt"
positions.vector <- unlist(start.positions)
names(positions.vector) <- NULL
positions.matrix <- matrix(positions.vector, byrow = T, ncol = 2)
write.table(positions.matrix, file = start.positions.file,
            row.names = F, # no escribe una columna adicional a la izquierda
            quote = F,     # quita las comillas
            col.names = F) # quita la x de la primera fila que siempre pone por defecto

## Esta informacion se utilizo para elaborar la figura de los 30 alineamientos



##Carga las funciones de enriquecimiento de genes. Necesita que se haya definido el threshold
## GSEA analysis. Dado el resultado del alineamiento, almacenado en gene.score, realiza
## un enriquecimiento de terminos GO asociados a los genes significativos,
## que se fijan como aquellos cuya puntuacion es mayor que el threshold. El threshold
## es el percentil 95 de la distribución de puntuaciones (definido mas arriba)
## GSEA.R carga todas las funciones y datos necesarios, entre ellas GeneSetEnrichmentAnalysis y GenTable
commandArgs <- function() threshold
source("GSEA.R")
gsea.result <- GeneSetEnrichmentAnalysis(gene.score)
allRes <- GenTable(gsea.result$GOData, classicFisher = gsea.result$Fisher,
                   classicKS = gsea.result$classic.KS, elimKS = gsea.result$elim.KS,
                   orderBy = "classicKS", ranksOf = "classicKS", topNodes = 10)



for(i in 1:10){
  print(i)
  goID <- allRes[i, "GO.ID"]
  print(printGenes(geneID2GO, goID, gene.score, threshold)$genes)
}

## Completa allRes añadiendo una columna que muestre algun gen asociado al GO correspondiente
Genes <- character()
for(i in 1:10){
  print(i)
  goID <- allRes[i, "GO.ID"]
  mygenes <- printGenes(geneID2GO, goID, gene.score, threshold)$genes
  go.genes <- names(mygenes)
  if (length(go.genes) == 0)
  {
    Genes[i] <- NA
  }
  else
  {
    go.genes <- paste(go.genes, collapse = " ")
    Genes[i] <- go.genes
  }
}

Genes
allRes <- cbind(allRes, Genes)


## Exportacion a latex
tfg.table <- allRes %>% dplyr::select(GO.ID, Term, Significant, Expected, classicKS,  Genes)
colnames(tfg.table) <- c("GO", "Termino",  "Presente", "Esperado", "p-valor", "Genes relevantes")
Hmisc::latex(tfg.table, rowname = NULL)

printGenes(geneID2GO, allRes[5, "GO.ID"], gene.score, threshold)
#sel.terms <- sample(usedGO(gsea.result$GOData), 10)
#termStat(gsea.result$GOData, sel.terms)
#graph(gsea.result$GOData)
#printGraph(gsea.result$GOData, gsea.result$classic.KS,
#firstSigNodes = 15, gsea.result$Fisher, fn.prefix = "tGO", useInfo = "all")


#png(file = "node.tree.png", width = 2000, height = 1750, res = 500)
#showSigOfNodes(gsea.result$GOData, score(resultKS.elim), firstSigNodes = 5, useInfo = 'all')
#dev.off()


#test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
#resultWeight <- getSigGroups(gsea.result$GOData, test.stat)
#head(score(resultWeight))
#
#pvalFis <- score(gsea.result$Fisher)
#hist(pvalFis, 50, xlab = "p-values")
#sort(pvalFis)[1:10]

# pval.classic.KS <- score(gsea.result$classic.KS)
# hist(pval.classic.KS, 50, xlab = "p-values")
# sort(pval.classic.KS)[1:10]
# printGraph(gsea.result$GOData, gsea.result$classic.KS, firstSigNodes = 10,
#            fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)


# dna.binding.genes <- printGenes(geneID2GO, "GO:0003677", gene.score, threshold)
# alignment.file <- "../Documento/TFG/alignment.fa"