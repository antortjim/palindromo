library("BSgenome.PCC7120.Ensembl.29")
upstream.distance <- 500
query.genes <- commandArgs()[[1]]
print("Will analyze following genes")
print(query.genes)
print(class(query.genes))
setwd("~/MEGA/CuartoCurso/TFG/Bioinformatica")

source("palindromo/loadGTF_Nostoc.R")

query.seqs <- list()

#system("mkdir BLASTquery")

#extrae la secuencia de cada gen desde el bsgenome de PCC7120, teniendo en cuenta la hebra
#y exportala a un archivo con el mismo nombre que el gen.fasta en la carpeta BLASTquery
for(gene in 1:length(query.genes))
{
 gene.name <- query.genes[gene]
 print(gene.name)
 positions <- nostoc.gtf %>% filter(feature == "CDS", gene_name == gene.name) %>% select(start, end, strand)
 start <- positions[1,1]
 end   <- positions[1,2]
 strand <- positions[1,3]
 if(strand == "+")
 {
   seq <- getSeq(BSgenome.PCC7120.Ensembl.29)[[1]][start:end] %>% as.character()
 }
 else if(strand == "-")
 {
   seq <- getSeq(BSgenome.PCC7120.Ensembl.29)[[1]][start:end] %>% as.character()
   seq <- seq %>% strsplit(split = "") %>% unlist() %>% rev() %>% comp() %>% toupper() %>% paste(collapse = "")
 } 
 fasta.file <- paste("BLASTquery/BLAST_", gene.name, ".fasta", sep = "")
 write.fasta(seq, names = gene.name , file = fasta.file)
}

infileNms <- list.files(path = "BLASTquery", pattern = "*.fasta")
infileNms


outnames <- paste(unlist(sapply(infileNms, strsplit, split = "*.fasta")), 
                  ".out", sep = "")


# create commands function
cmdCreate <- function(infile, outfile){
  paste("echo BLASTING ", infile, "; ", "blastn -db nr -query BLASTquery/", infile, " -remote -out BLASTout/",  outfile,
        ' -outfmt "7 qacc sacc sseqid sgi staxids sscinames evalue qstart qend sstrand sstart send qseq sseq"', sep = "")
}

# create actual commands
cmds <- mapply(FUN = cmdCreate, infile = infileNms, outfile = outnames)
##Blast all sequences stored in the fasta files
sapply(cmds, system)

## Concatenate all BLAST output into one R data.frame
system('cat BLASTout/* | grep -v "#" > "blast_output"')
blast.output <- read.table(file = "blast_output", sep = "\t")
colnames(blast.output) <- c("qacc", "sacc", "sseqid", "sgi", "staxids",
                            "sscinames", "evalue", "qstart",
                            "qend", "strand", "sstart",
                            "send", "qseq", "sseq")



## For every genome involved, extract the sequence going from -upstream.distance
#up to the start of each BLAST subject hit

sgis <- blast.output %>% select(sgi) %>% distinct()
#system("mkdir nucleotide")

reference.seqs <- list()
upstream.sequences <- list()
seq.ids <- character()

for (sgi in sgis[,1])
#sgis[,1]
{
  print(sgi)
  #Download the genome
  sgid <- sgi %>% as.character()
  fasta.file <- paste("nucleotide/", sgid, ".fasta" , sep = "")
  if(!file.exists(fasta.file))
  { 
    cmd <- paste("wget -O ", fasta.file,  ' "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=', sgid,  '&rettype=fasta"', sep = "")
    system(cmd)
  }
  
  # Read the genome and assign the genome to the list reference.seqs with id sgi
  sequence <- read.fasta(fasta.file, seqtype = "DNA") %>% getSequence() %>% unlist()
  reference.seqs[[sgid]] <- sequence
  
  # Extract start and end positions of every hit associated with genome sgi 
  positions <- blast.output %>% filter(sgi == sgid) %>% select(sstart, send)

  # For every hit, extract the sequence 500 upstream of it
  for(i in 1:nrow(positions))
  {
    seq.id <- paste(sgid, i, sep = "_")
    #print(seq.id)
    upstream.position <- positions[i, 1] - upstream.distance
    if(upstream.position < 0)
    {
      window <- c((length(sequence) - upstream.position):length(sequence), 1:positions[i, 1])
      }
    else
    {
    window <- upstream.position:positions[i, 1]
    }
    upstream.sequence <- sequence[window]
    upstream.sequences[[seq.id]] <- upstream.sequence
    seq.ids <- c(seq.ids, seq.id)
  }
}

str(upstream.sequences)

str(reference.seqs)
str(upstream.sequences)
