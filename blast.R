library("BSgenome.PCC7120.Ensembl.29")
upstream.distance <- 500
query.genes <- commandArgs()[[1]]
query.genes <- c("alr2431","asr3168","alr0296")
print("Will analyze following genes")
print(query.genes)
print(class(query.genes))
setwd("~/MEGA/CuartoCurso/TFG/Bioinformatica")

source("palindromo/loadGTF_Nostoc.R")

query.seqs <- list()

#system("mkdir BLASTquery")

#extrae la secuencia de cada uno de los query.genes desde el bsgenome de PCC7120, teniendo en cuenta la hebra
#y exportala a un archivo con el mismo nombre que el gen.fasta en la carpeta BLASTquery
for(gene in 1:length(query.genes))
{
  gene.name <- query.genes[gene]
  fasta.file <- paste("BLASTquery/BLAST_", gene.name, ".fasta", sep = "")
  if(file.exists(fasta.file))
  {
    next
  }
  print(paste("Exporting ", gene.name, '"s fasta sequence as extracted from the BSgenome object', sep = ""))
  positions <- nostoc.gtf %>% filter(feature == "CDS", gene_name == gene.name) %>% select(start, end, strand)
  if(nrow(positions) == 0)
  {
    print(paste(gene.name, " is not present in nostoc.gtf by that name", sep = ""))
    next
  }
  start <- positions[1,1]
  end   <- positions[1,2]
  dna.strand <- positions[1,3]
  if(dna.strand == "+")
  {
    seq <- getSeq(BSgenome.PCC7120.Ensembl.29)[[1]][start:end] %>% as.character()
  }
  else if(dna.strand == "-")
  {
    seq <- getSeq(BSgenome.PCC7120.Ensembl.29)[[1]][start:end] %>% as.character()
    seq <- seq %>% strsplit(split = "") %>% unlist() %>% rev() %>% comp() %>% toupper() %>% paste(collapse = "")
  } 
  write.fasta(seq, names = gene.name , file = fasta.file)
}

infileNms <- paste("BLAST_", query.genes, ".fasta", sep = "")



outnames <- paste(unlist(sapply(infileNms, strsplit, split = "*.fasta")), 
                  ".out", sep = "")


# create commands function, This function creates the command required to blast each sequence in BLASTquery
cmdCreate <- function(infile, outfile){
  if(file.exists(paste("BLASTout/", outfile, sep = "")))
  {
    paste("echo File ", infile, " already downloaded and blasted, skipping", sep = "")
  }
  else
  {
    paste("echo BLASTING ", infile, "; ", "blastn -db nr -query BLASTquery/", infile, " -remote -out BLASTout/",  outfile,
          ' -outfmt "7 qacc stitle sacc sseqid sgi staxids sscinames evalue qstart qend strand sstart send"', sep = "")
  }
}

# create actual commands
cmds <- mapply(FUN = cmdCreate, infile = infileNms, outfile = outnames)
## CALLING BLASTN
## Blast all sequences stored in the fasta files
sapply(cmds, system)

## Concatenate all BLAST output into one R data.frame
system('cat BLASTout/* | grep -v "#" > "blast_output"')

cmd <- paste("cat ", paste("BLASTout/BLAST_", query.genes, ".out", sep = "", collapse = " "), ' | grep -v "#" > "blast_output"')
system(cmd)

blast.output <- read.table(file = "blast_output", sep = "\t")
colnames(blast.output) <- c("qacc", "stitle", "sacc", "sseqid", "sgi", "staxids",
                            "sscinames", "evalue", "qstart",
                            "qend", "strand", "sstart",
                            "send")

## Remove all entries from PCC 7120
blast.output <- blast.output %>% filter(sscinames != "Nostoc sp. PCC 7120")

## Define a new column that states how many letters upstream will be explored. It will be upstream.distance + the distance between the query start and the query alignment (qstart)
blast.output$upstream.distance <- upstream.distance + blast.output$qstart


## Initialize two lists that will store upstrem sequences and their identifiers
upstream.sequences <- list()
seq.ids <- character()

# For each row in blast.output
for(i in 1:nrow(blast.output))
{
  ## Extract the sequence's NCBI identifier (sgid)
  sgi <- (blast.output %>% select(sgi))[i,]
  
  ## Extract the sequence's strand
  dna.strand <- select(blast.output, strand)[i,] %>% as.character()
  
  ## Check if it has already been downloaded
  fasta.file <- paste("nucleotide/", sgi, ".fasta" , sep = "")
  
  # If it hasn't, create the command that downloads it and execute it
  if(!file.exists(fasta.file))
  { 
    cmd <- paste("wget -O ", fasta.file,  ' "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=', sgi,  '&rettype=fasta"', sep = "")
    print(paste("Downloading ", sgi, ' "s fasta file', sep = ""))
    system(cmd)
  }
  # If it has, warn so
  else
  {
    print(paste("Sequence ", sgi, " already downloaded", sep = ""))
  }
  
  # Read the sequence just downloaded and assign it to "sequence"
  sequence <- read.fasta(fasta.file, seqtype = "DNA") %>% getSequence() %>% unlist()
  
  # Extract the elements that will make the seq.id 
  seq.id <- select(blast.output, qacc, stitle)[i,] %>% as.list() %>% unlist() %>% as.character() 
  seq.id[2] <- seq.id[2] %>% strsplit(split = " ") %>% unlist() %>% paste(collapse = "_")
  # Collapse them into one string
  seq.id <- paste(seq.id, collapse = "|")
  
  
  
  if(dna.strand == "plus")
  {
   ## Compute the subject upstream position
   upstream.position <- blast.output$sstart[i] - blast.output$upstream.distance[i]
   if(length(sequence) < upstream.distance)
   {
     ## Window cannot be bigger than sequence
     window <- 1:length(sequence)
   }
   else if(upstream.position < 0)
   {
     ## Circular
     window <- c((length(sequence) + upstream.position):length(sequence), 1:blast.output$sstart[i])
   }
   else
   {
     # Regular window
     window <- upstream.position:blast.output$sstart[i]
   }
  }
  
  
  else if(dna.strand == "minus")
  {
   ## Compute the subject upstream position
   upstream.position <- blast.output$sstart[i] + blast.output$upstream.distance[i]
   if(length(sequence) < upstream.distance)
   {
     ## Window cannot be bigger than sequence
     window <- 1:length(sequence)
   }
   else if(upstream.position > length(sequence))
   {
     ## Circular
     window <- c(blast.output$sstart[i]:length(sequence), 1:upstream.position)
   }
   else
   {
     # Regular window
     window <- blast.output$sstart[i]:upstream.position
   }
  }
   print(length(window))
   
   ## Subset the sequence and retain only the one spanning the window
   upstream.sequence <- sequence[window]
   
   ## Add the seq.id and the upstream.sequence to its corresponding list
   upstream.sequences[[seq.id]] <- upstream.sequence
   seq.ids <- c(seq.ids, seq.id)

}

  
  
write.fasta(upstream.sequences, file.out = "500pb_upstream_homologues.fa", names = seq.ids)

write.table(file = "blast_output_processed.txt", blast.output, col.names = T)
