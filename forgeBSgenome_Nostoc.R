library("BSgenome")

## Este bloque de codigo permite generar el objeto BSgenome de mi organismo
## a partir de los fasta de sus moleculas de ADN (cromosomas o plasmidos) y
## un archivo seed
## El seed file tiene que cumplir con los requisitos especificados en la viñeta de BSgenomeForge
## vignette("BSgenomeForge")

## Indica a R donde esta el seed-file
seed_files <- "/home/antortjim/MEGA/CuartoCurso/TFG/Bioinformatica/7120-seed"

## Comprueba su contenido
cat(readLines(seed_files), sep="\n")

## Forge el genoma
forgeBSgenomeDataPkg(seed_files)

## Tendras una carpeta llamada pckge_name (el nombre que le diste en el seed file)
## dentro de la carpeta que tiene el seed file
## Esta carpeta contiene el código necesario para compilar un paquete cargable en R
## Abre una terminal bash en la carpeta y ejecuta en este orden:
system("R CMD build BSgenome.PCC7120.Ensembl.29")            # compila el paquete
system("R CMD INSTALL BSgenome.PCC7120.Ensembl.29*tar.gz")   # instala el paquete

## El paquete ya se puede cargar con un library()

## Carga el BSgenome object generado con el codigo comentado justo encima
## Las instrucciones para generar el objeto estan muy bien descritas en la
## viñeta del paquete