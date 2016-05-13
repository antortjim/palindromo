
## Fija el espacio de trabajo
setwd("~/MEGA/CuartoCurso/TFG/Bioinformatica")
load("parser.RData")
## Parseador del fichero.gbk que genero el rna-seq de Alicia

##Carga de paquetes necesarios
library("ggplot2")
library("reshape")
library("dplyr")
library("grid")
library("gtable")
library("showtext")
ggplot_with_subtitle <- function(gg, 
                                 label="", 
                                 fontfamily=NULL,
                                 fontsize=10,
                                 hjust=0, vjust=0, 
                                 bottom_margin=5.5,
                                 newpage=is.null(vp),
                                 vp=NULL,
                                 ...) {
  
  if (is.null(fontfamily)) {
    gpr <- gpar(fontsize=fontsize, ...)
  } else {
    gpr <- gpar(fontfamily=fontfamily, fontsize=fontsize, ...)
  }
  
  subtitle <- textGrob(label, x=unit(hjust, "npc"), y=unit(hjust, "npc"), 
                       hjust=hjust, vjust=vjust,
                       gp=gpr)
  
  data <- ggplot_build(gg)
  
  gt <- ggplot_gtable(data)
  gt <- gtable_add_rows(gt, grobHeight(subtitle), 2)
  gt <- gtable_add_grob(gt, subtitle, 3, 4, 3, 4, 8, "off", "subtitle")
  gt <- gtable_add_rows(gt, grid::unit(bottom_margin, "pt"), 3)
  
  if (newpage) grid.newpage()
  
  if (is.null(vp)) {
    grid.draw(gt)
  } else {
    if (is.character(vp)) seekViewport(vp) else pushViewport(vp)
    grid.draw(gt)
    upViewport()
  }
  
  invisible(data)
  
}


## Accede a los datos
data <- readLines(con = "Ana7120TSS.gbk")


## Introducelos en R: parsea el archivo
## El parseador inicializa datos que guardan la informacion correspondiente
## a las 12797  entradas del fichero de Alicia (el número de TSS).
## En este caso son 12797
index <- 1
ids <- character(length = 12797)
expression.matrix <- matrix(nrow = 12797, ncol = 4)
colnames(expression.matrix) <- c("WT0", "WT8", "hetR0", "hetR8")
for (i in 1:length(data))
{
 ## El bucle for itera en cada linea, y si encuentra la palabra clave TSS ejecuta
 ## el bloque siguiente. La palabra TSS marca el inicio de una nueva entrada
 ## (Se ha dejado atras la anterior y empieza una nueva de otro gen)
 if (length(grep("TSS", data[i])) == 1)
 {
   print(index)
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
   
   ## Inicializa una lista con las 4 cadenas de caracteres y separa cada una por
   ## el simbolo de igual "=". Quedate con lo que queda despues de = en cada cadena
   ##(el segundo elemento) y finalmente pasalo a numerico
   ## Obtienes un vector numerico que guarda la expresion del mismo TSS en
   ## WT0, WT8, hetR0 y hetR8 por ese orden
   current.tss <- list(current.wt0, current.wt8, current.hetR0, current.hetR8)
   current.tss <- as.numeric(unlist(lapply(current.tss, function(x) strsplit(x, split= "=")))[seq(from = 2, to = 8, by = 2)])
   ids[index] <- current.id
   expression.matrix[index,] <- current.tss
   index <- index + 1
 }
 ## Esto corta el bucle for porque. Cada vez que se encuentra una nueva entrada,
 ## una linea que contiene "TSS", index suma 1. Como sabemos que hay 12797
 ## (grep -c "TSS" fichero .gbk) sabremos que cuando se haya encontrado esa
 ## cantidad no tendra sentido que el bucle continue
 if (index > 12797) break
}
which(ids == "all1873")


exp.matrix <- matrix(nrow = 12797*4, ncol = 3)
colnames(exp.matrix) <- c("TSS", "genotype", "time")
exp.matrix[,1] <- rep(ids, each = 4)
exp.matrix[,2] <- rep(rep(c("WT","hetR"), each = 2), times = 12797)
exp.matrix[,3] <- rep(c(0, 8, 0, 8), times = 12797)


expression.df <- as.data.frame(exp.matrix)

expression.matrix <- t(expression.matrix)
mymatrix <- expression.matrix
expression.df[["score"]] <- as.numeric(mymatrix)
expression.df

all1873 <- expression.df[(expression.df[["TSS"]] == "all1873"),]

all1873[["genotype"]] <- factor(c("WT","WT","hetR","hetR"), levels = c("WT", "hetR"))
all1873[["score"]] <- all1873[["score"]]/1000

p <- ggplot(data = all1873, aes(x = genotype, y = score, fill = time)) + geom_bar(width = 0.7, pos = "dodge", stat = "identity") 
p <- p + scale_fill_discrete(name="Tiempo (h)")
p <- p + labs(x = NULL, y = NULL, title = "Perfil de expresión de all1873")
p <- p + labs(y = "Expresión (kreads)", x = "")
p <- p + scale_fill_discrete(name="Tiempo (h)")
p <- p + theme_set(theme_bw(base_size = 15))
p <- p + theme(plot.title = element_text(hjust = 0)) 
p

subtitle <- "En carencia de nitrógeno, es uno de los transcritos más expresados del WT,\nmientras que en el mutante hetR esa respuesta se pierde"
p <- ggplot_with_subtitle(p, subtitle, bottom_margin = 20, lineheight = 0.9, fontfamily = "")
ggsave("../Documento/TFG/Figures/barplot.pdf", width = 5, height = 7)


# p <- ggplot(data = all1873, aes(x = time, y = score, col = genotype, group = genotype)) + geom_point() + geom_line()
# p <- p + scale_fill_discrete(name="Tiempo (h)")
# #p <- p + ggtitle("Perfil de expresión de all1873")
# p <- p + labs(y = "Expresión (kreads)", x = "")
# p <- p + scale_fill_discrete(name="Tiempo (h)")
# p <- p + theme_set(theme_bw(base_size = 15))
# p



WT.comparison.df <- dplyr::filter(expression.df, genotype == "WT") %>% dplyr::select(-genotype)
nrow(WT.comparison.df)/2

WT0 <- WT.comparison.df[WT.comparison.df[["time"]] == 0, "score"]
WT8 <- WT.comparison.df[WT.comparison.df[["time"]] == 8, "score"]
WT.comparison.df <- data.frame(ids, WT0, WT8) 
WT.comparison.df <- na.omit(WT.comparison.df)

ratiowt8    <- numeric(length=12797)
ratiohetr0  <- numeric(length=12797)
ratiohetr8  <- numeric(length=12797)

index <- 1
for (i in 4*(1:12797)-3)
{
  print(index)
  current.ratiowt8   <- expression.df[i+1,4] / expression.df[i,4]
  current.ratiohetr0 <- expression.df[i+2,4] / expression.df[i,4]
  current.ratiohetr8 <- expression.df[i+3,4] / expression.df[i,4]
  ratiowt8[index] <- current.ratiowt8
  ratiohetr0[index] <- current.ratiohetr0
  ratiohetr8[index] <- current.ratiohetr8
  index <- index +1
}
#unos segundos y listo

ratio.matrix <- matrix(c(ratiowt8, ratiohetr0, ratiohetr8), nrow = 3, byrow = T)
colnames(ratio.matrix) <- ids

ratio.df <- data.frame(rep(ids, each = 3), as.numeric(ratio.matrix))

ggplot(data = WT.comparison.df) + geom_point(aes(x = log10(WT0), y = log10(WT8)))


WT8.fold.change <- WT.comparison.df[["WT8"]] / WT.comparison.df[["WT0"]]
WT.comparison.df[["fold_change"]] <- WT8.fold.change
WT.comparison.df[["significant"]] <- ifelse(WT8.fold.change > 3 | WT8.fold.change < 1/3, "yes", "no" )

rhg_cols <- c("#771C19", "#AA3929", "#E25033", "#F27314", "#F8A31B", 
              "#E2C59F", "#B6C5CC", "#8E9CA3", "#556670", "#000000")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

plot(1:8, pch=16, cex=2, col=cbPalette)

## add the Arial font
font.add("FranklinGothic-Book",
         regular = "/home/antortjim/MEGAsync/Multimedia/Fonts/FranklinGothicFSBook/Fonts/franklingothic-book.otf")


library("extrafont")
View(fonttable()
)setEPS()

p <- ggplot(data = WT.comparison.df)
p <- p + geom_point(aes(x = log10(WT0), y = log10(WT8), colour = significant))
#p <- p + labs(x = NULL, y = NULL, title = "RNA-seq de la transición a ausencia de nitrógeno 8h")
p <- p + labs(y = expression(log[1*0]*"(8h)"), x = expression(log[1*0]*"(0h)"))
p <- p + geom_text(aes(x = log10(2198), y = log10(21648) + 0.15, label = "all1873"))
p <- p + guides(colour = FALSE)
p <- p + scale_color_manual(values = c("#000099", "#E69F00"))
p <- p + theme_set(theme_bw(base_size = 15))
p <- p + theme(plot.title = element_text(hjust = 0)) 
p
ggsave("../Documento/TFG/Figures/scatterplot.pdf")

# subtitle <- "El panorama transcripcional revela que all1873 es uno de los transcritos más expresados"
# p <- ggplot_with_subtitle(p, subtitle, bottom_margin = 20, 
#                      lineheight = 0.9, fontfamily = "")
#grid.gedit("GRID.text",gp=gpar(fontfamily="serif"))

#hetR.comparison.df <- filter(expression.df, genotype == "hetR") %>% select(-genotype)
#hetR0 <- hetR.comparison.df[hetR.comparison.df[["time"]] == 0, "score"]
#hetR8 <- hetR.comparison.df[hetR.comparison.df[["time"]] == 8, "score"]
#hetR.comparison.df <- data.frame(ids, hetR0, hetR8) 
#hetR.comparison.df <- na.omit(hetR.comparison.df)
#
#hetR8.fold.change <- hetR.comparison.df[["hetR8"]] / hetR.comparison.df[["hetR0"]]
#hetR.comparison.df[["fold_change"]] <- hetR8.fold.change
#hetR.comparison.df[["significant"]] <- ifelse(hetR8.fold.change > 3 | hetR8.fold.change < 1/3, "yes", "no" )
#
#
#WT.comparison.df[4215,]
#which(WT.comparison.df[["ids"]] == "all1873")
#
#q <- ggplot(data = hetR.comparison.df) + 
#  geom_point(aes(x = log10(hetR0), y = log10(hetR8), colour = significant)) +
#  geom_text(aes(x = log10(3961), y = log10(7830) + 0.15, label = "all1873")) + 
#  ggtitle("Carencia de N durante 8 h") + guides(colour = FALSE)
#
#q <- q + scale_color_manual(values = c("#000099", "#E69F00"))
#q
#
#