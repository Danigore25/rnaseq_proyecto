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
library(ggplot2)


# 2. Descargar el estudio deseado de la página de recount3.

# R Code
rse_SRP131764 <- recount3::create_rse_manual(
    project = "SRP131764",
    project_home = "data_sources/sra",
    organism = "mouse",
    annotation = "gencode_v23",
    type = "gene"
)

# 3. Guardar cuentas de genes presentes en muestras.
assay(rse_SRP131764, "counts") <- compute_read_counts(rse_SRP131764)


# 4. Ver los atributos de las muestras.
attributes <- expand_sra_attributes(rse_SRP131764)
colData(attributes)[, grepl("^sra_attribute", colnames(colData(attributes)))]


# 5. Guardar las variables de interés. En este caso serán aquellas que corresponden
# al origen del microbioma humano (proveniente de una persona en el espectro
# autista o un control), la relación del gen con la muestra y el tejido.

attributes$sra_attribute.mouse_status <- as.factor(attributes$sra_attribute.mouse_status)

attributes$sra_attribute.source_name <- as.factor(attributes$sra_attribute.source_name)

attributes$sra_attribute.tissue <- as.factor(attributes$sra_attribute.tissue)


# 6. Analizar las muestras del microbioma, separarlos por un factor.
attributes$sra_attribute.mouse_status <- factor(ifelse(attributes$sra_attribute.mouse_status == "colonized with Human feces (ASD)", "ASD-feces", "control-feces"))

# table (attributes$sra_attribute.mouse_status)

