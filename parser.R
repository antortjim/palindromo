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

threshold <- 3
commandArgs <- function() threshold
source("rna-seq-accession.R")


p <- ggplot(data = expression.df, aes(x = log10(WT0), y = log10(WT8), colour = significant != "no")) + geom_point()
#p <- p + labs(x = NULL, y = NULL, title = "RNA-seq de la transición a ausencia de nitrógeno 8h")
p <- p + labs(y = expression(log[1*0]*"(8h)"), x = expression(log[1*0]*"(0h)"))
p <- p + guides(colour = FALSE)
p <- p + scale_color_manual(values = c("#000099", "#E69F00"))
p <- p + theme_set(theme_bw(base_size = 15))
p <- p + theme(plot.title = element_text(hjust = 0))
p <- p + geom_text(aes(x = log10(2198), y = log10(21648) + 0.15, label = "all1873"), color = "black")

p
ggsave("./../Documento/TFG/Figures/scatterplot.pdf")
