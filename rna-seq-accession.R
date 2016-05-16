## Fija el espacio de trabajo
setwd("~/MEGA/CuartoCurso/TFG/Bioinformatica")

##Carga de paquetes necesarios
#Acceso a tabla excel
library("openxlsx")

#Graficos
library("ggplot2")
library("grid")

threshold <- commandArgs()[1] %>% as.numeric()
## Accede a los datos

print("Reading xlsx file")
expression.df <- read.xlsx("sd01.xlsx", cols = c(1:4, 7:8, 10, 12, 14))
colnames(expression.df) <- expression.df[1,]
expression.df <- expression.df[-1,]

print("Transformic numerical values into numerical data")
expression.df[,4] <- expression.df[,4]  %>% as.numeric()
expression.df[,5] <- expression.df[,5]  %>% as.numeric()
expression.df[,6] <- expression.df[,6]  %>% as.numeric()
expression.df[,7] <- expression.df[,7]  %>% as.numeric()
expression.df[,8] <- expression.df[,8]  %>% as.numeric()


print("Calculating fold change")
expression.df$fold.change.WT <- expression.df[["Reads WT-8"]] / expression.df[["Reads WT-0"]]


print("Classifying TSS based on fold change")
expression.df$significant <- rep("no", times = nrow(expression.df))
expression.df$significant[expression.df$fold.change > threshold]   <- "increase"
expression.df$significant[expression.df$fold.change < 1/threshold] <- "decrease"


colnames(expression.df)[c(1,5:8)] <- c("TSS.class","WT0","WT8","hetR8","hetR0")

print("Generating scatter plot")
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