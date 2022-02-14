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
library(edgeR)
library(ggplot2)
library(limma)
library(pheatmap)
# library(RColorBrewer)


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


# 5. Guardar las variables de interés (origen fecal, relación del gen con la
# muestra y y tejido).

attributes$sra_attribute.mouse_status <- as.factor(attributes$sra_attribute.mouse_status)

attributes$sra_attribute.source_name <- as.factor(attributes$sra_attribute.source_name)

attributes$sra_attribute.tissue <- as.factor(attributes$sra_attribute.tissue)


# 6. Cambiar nombre de tipo de origen de microbioma.
attributes$sra_attribute.mouse_status <- factor(ifelse(attributes$sra_attribute.mouse_status == "colonized with Human feces (ASD)", "ASD-feces", "control-feces"))


# 7. Comparar cuentas de atributos.
attributes$assigned_gene_prop <- attributes$recount_qc.gene_fc_count_all.assigned / attributes$recount_qc.gene_fc_count_all.total
summary(attributes$assigned_gene_prop)


# 8. Plotear atributos.
with(colData(attributes), plot(sra_attribute.mouse_status, assigned_gene_prop))


# 9. Realizar un filtrado de los datos.
attributes_filtered <- attributes
hist(attributes_filtered$assigned_gene_prop)
table(attributes_filtered$assigned_gene_prop < 0.75)
attributes_filtered <- attributes[, attributes$assigned_gene_prop > 0.75]


# 10. Eliminar genes según su media de expresión.
expression_means <- rowMeans(assay(attributes_filtered, "counts"))
attributes_without <- attributes_filtered[expression_means > 0.1, ]


# 11. Normalizar datos.
dge <- DGEList(
  counts = assay(attributes_without, "counts"),
  genes = rowData(attributes_without)
)
dge <- calcNormFactors(dge)


# 12. Expresión diferencial.
ggplot(as.data.frame(colData(attributes_without)), aes(y = assigned_gene_prop, x = sra_attribute.mouse_status)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Microbiome origin")


# 13. Model matrix con atributos de interés.
mod <- model.matrix(~ sra_attribute.mouse_status + sra_attribute.source_name + sra_attribute.tissue + assigned_gene_prop, data = colData(attributes_without))
colnames(mod)


# 14. Hacer gráfica voom con limma.
vGene <- voom(dge, mod, plot = TRUE)


# 15. Ajustar resultados.
eb_results <- eBayes(lmFit(vGene))

de_results <- topTable(eb_results, coef = 2, number = nrow(attributes_without),
                       sort.by = "none"
)
table(de_results$adj.P.Val < 0.05)


# 16. Visualizar resultados estadísticos.
plotMA(eb_results, coef = 2)


# 17. Volcano plot.
volcanoplot(eb_results, coef = 2, highlight = 3, names = de_results$gene_name)


# 18. Ver los 50 genes diferencialmente más expresados.
exprs_heatmap <- vGene$E[rank(de_results$adj.P.Val) <= 50, ]

df <- as.data.frame(colData(attributes_without)[, c("sra_attribute.source_name", "sra_attribute.mouse_status")])
colnames(df) <- c("Source name", "Fecal origin")

# HEATMAP
pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = df
)

# MDS
col.group <- df$Fecal
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
plotMDS(vGene$E, labels = df$Fecal, col = col.group)


