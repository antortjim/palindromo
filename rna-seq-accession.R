#This script parses the .gbk file retrievable at http://www.cyanolab.de/suppdata/7120/Ana7120TSS.gbk
#into a data.frame
#The transcriptomic information is stored as a data.frame

## Fija el espacio de trabajo
setwd("~/MEGA/CuartoCurso/TFG/Bioinformatica")

##Carga de paquetes necesarios
library("ggplot2")
library("reshape")
library("dplyr")
library("grid")
library("gtable")
library("showtext")
threshold <- commandArgs()[1] %>% as.numeric()
## Accede a los datos
data <- readLines(con = "Ana7120TSS.gbk")


## Introducelos en R: parsea el archivo
## El parseador inicializa datos que guardan la informacion correspondiente
## a las 12797  entradas del fichero de Alicia (el número de TSS).
#grep -c "^     TSS" Ana7120TSS.gbk 
## En este caso son 12797
index <- 1
ids <- character(length = 12797)
expression.matrix <- matrix(nrow = 12797, ncol = 5)
colnames(expression.matrix) <- c("WT0", "WT8", "hetR0", "hetR8", "coordinate")
strand.character <- character(length = 12797)
for (i in 1:length(data))
{
  ## El bucle for itera en cada linea, y si encuentra la palabra clave TSS ejecuta
  ## el bloque siguiente. La palabra TSS marca el inicio de una nueva entrada
  ## (Se ha dejado atras la anterior y empieza una nueva de otro gen)
  if (length(grep("TSS", data[i])) == 1)
  {
    if(index %% 1000 == 0)
    {
    print(index)
    }
    ## Accede al nombre del gen (3a fila contando con la primera,
    ## la que marca la entrada). Para acceder al gen hay que
    ## separar esa linea por el caracter " y quedarse con 
    ## el segundo elemento (lo que queda despues de ")
    current.id <- unlist(strsplit(data[i+2], split = "\"", fixed = T))[2]
    
    ## Guarda las 4 lineas en las que se gurda la informacion deseada.
    ## Estas son la 4, 5, 6, y 8 lineas desde la primera linea
    ## (la marcada por la keyword TSS)
    current.wt0 <- data[i+3]
    current.hetR0 <- data[i+4]
    current.wt8 <- data[i+5]
    current.hetR8 <- data[i+8]
    current.coordinate <- data[i]
    current.coordinate <- current.coordinate %>% strsplit(split = " ") %>% unlist() %>% tail(n = 1)
    strand <- "+"
    if(length(grep("complement", current.coordinate)) == 1)
    {
      current.coordinate <- current.coordinate %>% strsplit(split = "") %>% unlist()
      ok.pos <- current.coordinate %in% 0:9
      current.coordinate <- current.coordinate[ok.pos] %>% paste(collapse = "") %>% as.numeric()
      strand <- "-"
    }
    
    ## Inicializa una lista con las 4 cadenas de caracteres y separa cada una por
    ## el simbolo de igual "=". Quedate con lo que queda despues de = en cada cadena
    ##(el segundo elemento) y finalmente pasalo a numerico
    ## Obtienes un vector numerico que guarda la expresion del mismo TSS en
    ## WT0, WT8, hetR0 y hetR8 por ese orden
    current.tss <- list(current.wt0, current.wt8, current.hetR0, current.hetR8)
    current.tss <- as.numeric(unlist(lapply(current.tss, function(x) strsplit(x, split= "=")))[seq(from = 2, to = 8, by = 2)])
    ids[index] <- current.id
    expression.matrix[index,] <- c(current.tss, current.coordinate)
    strand.character[index] <- strand
    index <- index + 1
  }
  ## Esto corta el bucle for porque. Cada vez que se encuentra una nueva entrada,
  ## una linea que contiene "TSS", index suma 1. Como sabemos que hay 12797
  ## (grep -c "TSS" fichero .gbk) sabremos que cuando se haya encontrado esa
  ## cantidad no tendra sentido que el bucle continue
  if (index > 12797) break
}


expression.df <- as.data.frame(expression.matrix)
expression.df <- cbind(expression.df, strand.character)

expression.df$WT0   <- expression.df$WT0 %>% as.character() %>% as.numeric()
expression.df$WT8   <- expression.df$WT8 %>% as.character() %>% as.numeric()
expression.df$fold.change <- expression.df$WT8 / expression.df$WT0

expression.df$significant <- rep("no", times = nrow(expression.df))
expression.df$significant[expression.df$fold.change > threshold]   <- "increase"
expression.df$significant[expression.df$fold.change < 1/threshold] <- "decrease"


rm(strand, strand.character, ids, i, current.tss, current.wt0, current.wt8, current.hetR0, current.hetR8, index, data, current.id, current.coordinate)
