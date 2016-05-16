## Fija el espacio de trabajo
setwd("~/MEGA/CuartoCurso/TFG/Bioinformatica")

##Carga de paquetes necesarios
#Acceso a tabla excel
library("openxlsx")

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