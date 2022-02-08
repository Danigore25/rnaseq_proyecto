# 1. Cargar librería.
library(recount3)

# 2. Descargamos un estudio proveniente de la página de recount3. En este caso
# se utilizará el que tiene el ID SRP131764, con la anotación gencode_v23.

# R Code
rse <- recount3::create_rse_manual(
    project = "SRP131764",
    project_home = "data_sources/sra",
    organism = "mouse",
    annotation = "gencode_v23",
    type = "gene"
)
