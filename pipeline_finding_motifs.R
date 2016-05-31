## Fija el espacio de trabajo
setwd("~/MEGA/CuartoCurso/TFG/Bioinformatica")

##Carga de paquetes necesarios
library("ggplot2")
library("dplyr")
library("grid")
library("gtable")
library("showtext")

threshold <- 3
commandArgs <- function() threshold
source("rna-seq-accession.R")

#remove lines containing NA
expression.df <- expression.df %>% na.omit()
print("Generating scatter plot")
p <- ggplot(data = expression.df, aes(x = log2(WT0), y = log2(WT8), colour = significant != "no")) + geom_point()
#p <- p + labs(x = NULL, y = NULL, title = "RNA-seq de la transición a ausencia de nitrógeno 8h")
p <- p + labs(y = expression(log[2]*"(8h)"), x = expression(log[2]*"(0h)"))
p <- p + guides(colour = FALSE)
p <- p + scale_color_manual(values = c("#000099", "#E69F00"))
p <- p + theme_set(theme_bw(base_size = 15))
p <- p + theme(plot.title = element_text(hjust = 0))

p <- p +  annotate('text', log2(2198), log2(21648),
         label= "italic(all1873)", parse=TRUE,
         size=4, vjust = -1) +
#p <- p + geom_text(aes(x = log2(2198), y = log2(21648) + 0.15, label = "all1873"),
#                   color = "black", nudge_y = 0.5)
  geom_abline(intercept = +log2(threshold), linetype = 2, alpha= 0.5) +
  geom_abline(intercept = -log2(threshold), linetype = 2, alpha= 0.5) +
  geom_text(aes(x = max(log2(WT0)), y = log2(WT8)[which.max(log2(WT0))]),
            nudge_y = -0.75, nudge_x = 0.75, color = "black", label = paste("<1/", threshold, sep = "")) +
  geom_text(aes(y = max(log2(WT8)), x = log2(WT0)[which.max(log2(WT8))]),
            nudge_y = +0.75, nudge_x = -0.75, color = "black", label = paste(">", threshold,   sep = ""))


p
ggsave("./../Documento/TFG/Figures/scatterplot.pdf")


print("Launching database_creator.R")
system("Rscript database_creator.R")