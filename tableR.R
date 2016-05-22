# if(interactive())
# {
#  tabular.data <- xed.data
#  ntss <- 30
# }

tabular.data <- commandArgs()[[1]]
ntss <- commandArgs()[[2]]

colnames(tabular.data)[c(1:7, 17)] <- c("Coordenada","Hebra","Gen","Tipo de TSS","Distancia","FC", "Secuencia","N")
library("openxlsx")
wb <- createWorkbook()
addWorksheet(wb, "palindromicSequences")
azulOscuro  <- createStyle(fontColour = "#000000", bgFill = "#6464ff", fontSize = 8)
azulClarito <- createStyle(fontColour = "#000000", bgFill = "#9999ff", fontSize = 8)
normal      <- createStyle(fontColour = "#000000", bgFill = "#ffffff", fontSize = 8)

golden.number <- round(1.6180339, digits = 2)

setColWidths(wb, sheet = 1, cols = 2:30, widths = 7)
setColWidths(wb, sheet = 1, cols = 2, widths = 7)
setColWidths(wb, sheet = 1, cols = c(1, 4), widths = 11)
setColWidths(wb, sheet = 1, cols = 7:16, widths = 0.7)
setColWidths(wb, sheet = 1, cols = 17, widths = 0.8)
setRowHeights(wb, sheet = 1, rows = 2:100, heights = 14)
setRowHeights(wb, sheet = 1, rows = 1, heights = 25)

modifyBaseFont(wb, fontSize = 12, fontName = "Liberation Serif")


#Escribe los datos
writeData(wb, "palindromicSequences", tabular.data[1:ntss,])

#Fusiona los headers de cada letra para que haya solo una celda y ponga secuencias (y no X1, X2...)
mergeCells(wb, "palindromicSequences", 7:16, 1)
mergeCells(wb, "palindromicSequences", cols = 17, rows = 3:15)
mergeCells(wb, "palindromicSequences", cols = 17, rows = 16:(ntss+1))
writeData(wb, "palindromicSequences", startCol = 7, startRow = 1, x = "Secuencia")
writeData(wb, "palindromicSequences", startCol = 6, startRow = 1, x = "FC")


oscuro <- c(rep(F, 6), T, F, T, T, F, F, T, T, F, T)



for(i in 7:16)
{
  current.equal.rule   <- paste("==$", LETTERS[i], "$2", sep = "")
  current.unequal.rule <- paste("!=$", LETTERS[i], "$2", sep = "")
  if(oscuro[i])
  {
    conditionalFormatting(wb, "palindromicSequences", cols=i, rows=2:600, rule = current.equal.rule, style = azulOscuro)
  }
  else
  {
    conditionalFormatting(wb, "palindromicSequences", cols=i, rows=2:600, rule = current.equal.rule, style = azulClarito)
  }
  conditionalFormatting(wb, "palindromicSequences", cols=i, rows=2:600, rule = current.unequal.rule, style = normal)
}


#Redondea el fold change hasta dos decimales
addStyle(wb, "palindromicSequences", cols=6, rows=1:600, style = createStyle(numFmt = "0.00"))

#centra el texto en todas las celdas tanto vertical como horizontal
addStyle(wb, "palindromicSequences", gridExpand = T, stack = T, cols = 1:30, rows = 1:600, style = createStyle(halign = "center", valign = "center"))


topBorder    <- createStyle(border="Top", borderColour = "#000000", borderStyle = "thick")
leftBorder   <- createStyle(border="Left", borderColour = "#000000", borderStyle = "thick")
rightBorder  <- createStyle(border="Right", borderColour = "#000000", borderStyle = "thick")
bottomBorder <- createStyle(border="Bottom", borderColour = "#000000", borderStyle = "thick")
dottedBorder  <- createStyle(border="Bottom", borderColour = "#000000", borderStyle = "dotted")
pagella      <- createStyle(fontName = "Tex Gyre Pagella Math", fontSize = 1)


#borde de la tabla (adorno)
addStyle(wb, sheet = 1, topBorder, stack = T, rows = 1, cols = 1:(ncol(tabular.data)-1), gridExpand = TRUE)
addStyle(wb, sheet = 1, bottomBorder, stack = T, rows = 1, cols = 1:1:(ncol(tabular.data)-1), gridExpand = TRUE)
addStyle(wb, sheet = 1, leftBorder, stack = T, rows = 1:(ntss + 1), cols = 1, gridExpand = TRUE)
addStyle(wb, sheet = 1, rightBorder, stack = T, rows = 1:(ntss + 1), cols = ncol(tabular.data)-1, gridExpand = TRUE)
addStyle(wb, sheet = 1, bottomBorder, stack = T, rows = ntss + 1, cols = 1:(ncol(tabular.data)-1), gridExpand = TRUE)
addStyle(wb, sheet = 1, pagella, stack = T, rows = 2:(ntss + 1), cols = 7:(ncol(tabular.data)-1), gridExpand = TRUE)

#Separa la tabla en ecciones por el numero de mismatches
addStyle(wb, sheet = 1, dottedBorder, stack = T, rows = c(2, 15), cols = 1:ncol(tabular.data), gridExpand = TRUE)


## Guarda las configuraciones y exporta un archivo xlsx
saveWorkbook(wb, "palindromicSequences.xlsx", TRUE)
## Conviertelo a pdf
system("unoconv -f pdf -d spreadsheet palindromicSequences.xlsx")
## Elimina el espacio blanco
system("inkscape -f palindromicSequences.pdf -A palindromicSequences.pdf")
system("cp palindromicSequences.pdf ~/MEGA/CuartoCurso/TFG/Documento/TFG/Figures")
system("evince palindromicSequences.pdf")

