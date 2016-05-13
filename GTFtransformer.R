##Este script introduce en el gtf los datos de la transcriptomica de Alicia
## Asi, los transcritos quedan definidos segun su TSS, y no el ATG del gen que codifican (que esta mal la mayoria de las veces)
## Falta por implementar la insercion de TSS antisense (aTSS)

library("openxlsx")
annotation.file <- "Nostoc_sp_pcc_7120.v29.gtf"
tss.df <- read.xlsx("sd01.xlsx", sheet = 1, startRow = 2, colNames = TRUE)
tss.replacer <- function(gtf.file, tss.df)
{
  #gtf.file <- annotation.file
  gtf.readLines <- readLines(gtf.file)
  gtf.Readlines.length <- length(gtf.readLines)
  tuned.lines <- 0
  added.lines <- 0
  
  for(i in 6:gtf.Readlines.length)
  {
    #i <- 3781
    current.line <- unlist(strsplit(gtf.readLines[[i]], split = "\t"))
    former.line <- current.line
    transcription.strand <- ifelse(current.line[7] == "+", 4, 5)
    if(current.line[3] == "transcript")
    {
      print(i)
      print(transcription.strand)
      current.metadata <- unlist(strsplit(current.line[9], split = ";"))
      current.gene     <- unlist(strsplit(current.metadata[5], split = '\"'))[2]
      ## No queremos incluir los TSS antisense porque no son compatibles con el formato GTF
      ## Habria que ponerlos en una nueva linea y cambiar la hebra en la que se codifica
      tss.df.row <- tss.df[["Annotation"]] == current.gene & tss.df[["TSS.class"]] != "aTSS"
      tss.number <- sum(tss.df.row)
      if(tss.number == 0)
      {
        next
      }
      else if(tss.number > 1)
      {
        ##en el caso de que tenga mas de un tss asociado
        ##vas a repetir la operacion para todos los demas y pero añadiras
        ##la nueva linea al final de readLines.
        tss.vector <- formatC(tss.df[tss.df.row, "TSS"], width = 1, format = "d", flag = "0")
        current.line[transcription.strand] <- tss.vector[1]
        current.end <- as.integer(current.line[5])
        current.start <- as.integer(current.line[4])
        if(current.end - current.start < 0)
        {
          current.line <- former.line
          next
        }
        gtf.readLines[[i]] <- paste(current.line, collapse ="\t")
        tuned.lines <- tuned.lines + 1
        for(tss in 2:tss.number)
        { 
          current.line[transcription.strand] <- tss.vector[tss]
          current.end <- as.integer(current.line[5])
          current.start <- as.integer(current.line[4])
          if(current.end - current.start < 0)
          {
            current.line <- former.line
            next
          }
          ##Modificacion del nombre del transcrito, para que cada uno de los transcritos
          ## con TSS distinto no entre en conflicto con los demas (se llamaran tx-2, tx-3 sucesivamente)
          current.tx <- unlist(strsplit(unlist(strsplit(current.line[9], split = ";"))[3], split = '\"'))[2]
          new.tx <- paste(current.tx, tss, sep ="-")
          current.line.9 <- unlist(strsplit(current.line[9], split = ";"))
          current.line.9[3] <- paste(' transcript_id \"', new.tx, '\"', sep = "")
          current.line[9] <- paste(current.line.9, collapse=";")
          
          gtf.readLines   <- c(gtf.readLines, paste(current.line, collapse ="\t"))
          added.lines <- added.lines + 1
        }
        next
      }
      
      current.tss <- formatC(tss.df[tss.df.row, "TSS"], width = 1, format = "d", flag = "0")
      current.line[transcription.strand] <- current.tss
      current.end <- as.integer(current.line[5])
      current.start <- as.integer(current.line[4])
      if(current.end - current.start < 0)
      {
        current.line <- former.line
        next
      }
      
      gtf.readLines[[i]] <- paste(current.line, collapse ="\t")
      tuned.lines <- tuned.lines + 1
    }
  }
  write(file = "Nostoc_sp_pcc_7120.v29.AOJ.gtf" , gtf.readLines)
  return(list(tuned.lines = tuned.lines,
              added.lines = added.lines ))
}

tss.replacer(annotation.file, tss.df)