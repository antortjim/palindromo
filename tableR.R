# if(interactive())
# {
#  tabular.data <- xed.data
#  ntss <- 30
# }

tabular.data <- commandArgs()[[1]]
ntss <- commandArgs()[[2]]
#ntss <- nrow(tabular.data)
colnames(tabular.data)[c(1:7, 17)] <- c("Coordenada","Hebra","Gen","Tipo de TSS","Distancia","FC", "Secuencia"," ")

tabular.data$Secuencia
library("openxlsx")
wb <- createWorkbook()
addWorksheet(wb, "palindromicSequences")
azulOscuro  <- createStyle(fontColour = "#000000", bgFill = "#6464ff", fontSize = 8)
azulClarito <- createStyle(fontColour = "#000000", bgFill = "#9999ff", fontSize = 8)
normal      <- createStyle(fontColour = "#000000", bgFill = "#ffffff", fontSize = 8)


setColWidths(wb, sheet = 1, cols = 2:30, widths = 7)
setColWidths(wb, sheet = 1, cols = 2, widths = 7)
setColWidths(wb, sheet = 1, cols = c(1, 4), widths = 11)
setColWidths(wb, sheet = 1, cols = 7:16, widths = 1)
setColWidths(wb, sheet = 1, cols = 17, widths = 1.6)
setRowHeights(wb, sheet = 1, rows = 2:10000, heights = 14)
setRowHeights(wb, sheet = 1, rows = 1, heights = 25)

modifyBaseFont(wb, fontSize = 12, fontName = "Liberation Serif")


#Escribe los datos
writeData(wb, "palindromicSequences", tabular.data[1:ntss,])

mergeCells(wb, "palindromicSequences", 7:16, 1)
writeData(wb, "palindromicSequences", startCol = 7, startRow = 1, x = "Secuencia")
writeData(wb, "palindromicSequences", startCol = 6, startRow = 1, x = "FC")

# which(tabular.data$n.mismatch == 1) %>% tail()
#15-1
# which(tabular.data$n.mismatch == 2) %>% tail()
#314-1
# which(tabular.data$n.mismatch == 3) %>% tail()
#2888 -1

oscuro <- c(rep(F, 6), T, F, T, T, F, F, T, T, F, T)



for(i in 7:16)
{
  current.equal.rule   <- paste("==$", LETTERS[i], "$2", sep = "")
  current.unequal.rule <- paste("!=$", LETTERS[i], "$2", sep = "")
  if(oscuro[i])
  {
    conditionalFormatting(wb, "palindromicSequences", cols=i, rows=2:10000, rule = current.equal.rule, style = azulOscuro)
  }
  else
  {
    conditionalFormatting(wb, "palindromicSequences", cols=i, rows=2:10000, rule = current.equal.rule, style = azulClarito)
  }
  conditionalFormatting(wb, "palindromicSequences", cols=i, rows=2:10000, rule = current.unequal.rule, style = normal)
}


#Redondea el fold change hasta dos decimales
addStyle(wb, "palindromicSequences", cols=6, rows=1:10000, style = createStyle(numFmt = "0.00"))

#centra el texto en todas las celdas tanto vertical como horizontal
addStyle(wb, "palindromicSequences", gridExpand = T, stack = T, cols = 1:30, rows = 1:10000, style = createStyle(halign = "center", valign = "center"))

#define los estilos de los bordes y la fuente
topBorder    <- createStyle(border="Top", borderColour = "#000000", borderStyle = "thick")
leftBorder   <- createStyle(border="Left", borderColour = "#000000", borderStyle = "thick")
rightBorder  <- createStyle(border="Right", borderColour = "#000000", borderStyle = "thick")
bottomBorder <- createStyle(border="Bottom", borderColour = "#000000", borderStyle = "thick")
dottedBorder  <- createStyle(border="Bottom", borderColour = "#000000", borderStyle = "thick")
pagella      <- createStyle(fontName = "Tex Gyre Pagella Math", fontSize = 1)


#aplica los estilos del borde de la tabla (adorno)
addStyle(wb, sheet = 1, topBorder, stack = T, rows = 1, cols = 1:(ncol(tabular.data)-1), gridExpand = TRUE)
addStyle(wb, sheet = 1, bottomBorder, stack = T, rows = 1, cols = 1:1:(ncol(tabular.data)-1), gridExpand = TRUE)
addStyle(wb, sheet = 1, leftBorder, stack = T, rows = 1:(ntss + 1), cols = 1, gridExpand = TRUE)
addStyle(wb, sheet = 1, rightBorder, stack = T, rows = 1:(ntss + 1), cols = ncol(tabular.data)-1, gridExpand = TRUE)
addStyle(wb, sheet = 1, bottomBorder, stack = T, rows = ntss + 1, cols = 1:(ncol(tabular.data)-1), gridExpand = TRUE)
addStyle(wb, sheet = 1, pagella, stack = T, rows = 2:(ntss + 1), cols = ncol(tabular.data)-1, gridExpand = TRUE)

#Separa la tabla en secciones por el numero de mismatches
addStyle(wb, sheet = 1, dottedBorder, stack = T, rows = c(2, 15, 314, 2888), cols = 1:ncol(tabular.data), gridExpand = TRUE)


## Guarda las configuraciones y exporta un archivo xlsx
saveWorkbook(wb, "palindromicSequences.xlsx", TRUE)
## Conviertelo a pdf
system("unoconv -f pdf -d spreadsheet palindromicSequences.xlsx")
## Elimina el espacio blanco
system("inkscape -f palindromicSequences.pdf -A palindromicSequences.pdf")
system("cp palindromicSequences.pdf ~/MEGA/CuartoCurso/TFG/Documento/TFG/Figures")
system("evince palindromicSequences.pdf")

