'''
NAME
    01-analisis.R

VERSION
    1.0

AUTHOR
    Daniela Goretti Castillo León <danigore22@gmail.com> <dgoretti@lcg.unam.mx>

DESCRIPTION
    Este programa muestra la expresión diferencial en el proyecto SRP131764,
    relacionado con la transferencia de microbiota perteneciente a personas con
    trastorno del espectro del autismo a ratones sanos.

CATEGORY
    Análisis de expresión diferencial

USAGE
    01-analisis.R [sin opciones]

ARGUMENTS
    No se requieren argumentos.

PACKAGES
    recount3. Es un paquete que ayuda a abalizar datos y democratizar el acceso
              a ellos.
    ggplot2. Crea una gráfica de expresión diferencial a partir de un objeto
              tipo SummarizedExperiment.
    pheatmap. Crea una gráfica de tipo heatmap que ayuda a visualizar la
              expresión de los genes.

INPUT
    Objeto SummarizedExperiment.

OUTPUT
    Gráficas de expresión génica del proyecto.

EXAMPLES
    Se tiene .jrjkfjrjfhejkfhrjkhfjrhfkjherfjhrjkfhjrekhfjrhe

GITHUB LINK
    https://github.com/Danigore25/rnaseq_proyecto/blob/master/R/01-analisis.R
'''

# 1. Cargar librerías.
library(recount3)

# 2. Descargar el estudio deseado de la página de recount3.

# R Code
rse_SRP131764 <- recount3::create_rse_manual(
    project = "SRP131764",
    project_home = "data_sources/sra",
    organism = "mouse",
    annotation = "gencode_v23",
    type = "gene"
)

# 3.
attributes <- expand_sra_attributes(rse_SRP131764)

colData(attributes)[, grepl("^sra_attribute", colnames(colData(attributes)))]

