## Analisis de regiones codificantes
## Secuencia de las cds asociadas a los transcritos de nostoc
nostoc.cds <- cds(txdb, columns = "gene_id")
gene_id <- unlist(mcols(nostoc.cds)[["gene_id"]])
nostoc.cdsSeq <- getSeq(BSgenome.PCC7120.Ensembl.29, nostoc.cds)

#nostoc.cdsSeq <- extractTranscriptSeqs(BSgenome.PCC7120.Ensembl.29, cdsBy(txdb, by="tx"))

## Alineamiento del palindromo a todas las cds de Nostoc.
## Devuelve un objeto de clase Local PairwiseAlignments que muestra el mejor alineamiento
palindromo.vs.cds <- pairwiseAlignment(pattern = rep(palindromo, length(nostoc.cdsSeq)), subject = nostoc.cdsSeq,
                                       type = "local", gapOpening=-4, gapExtension=-1, substitutionMatrix = BLOSUM62)

plot(score(palindromo.vs.cds))
hist(score(palindromo.vs.cds))$breaks
alignment.score <- score(palindromo.vs.cds)
gene2tx[]
names(alignment.score) <- gene_id
significant.cds <- subset(alignment.score, alignment.score > 55)


## Background
genes <- nostoc.gtf[,17]
names(genes) <-  nostoc.gtf[,14]

nostocFile <- "Nostoc_genes.txt"
cat(head(readLines(nostocFile)), sep = "\n")
nostoc <- read.table(nostocFile, sep="\t",  comment.char = "", quote = "")
names(alignment.score) <- mcols(nostoc.cds)[,2]


## Enriquecimiento
source("GSEA.R")
GeneSetEnrichmentAnalysis(alignment.score)
hist(as.numeric(cast(nGO[,1:2], . ~ GID)))
