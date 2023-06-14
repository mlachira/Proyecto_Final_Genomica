#Melisa Lachira
#Diana Karina Rangel
#Maria del Carmen Ramirez

#Cargamos la libreria ya previamente descargada "phyloseq"
#Si no se tiene descargado "phyloseq" se puede descargar con el siguiente codigo
#if (!require("BiocManager", tranquilamente = VERDADERO))
#install.packages("BiocManager")
#BiocManager::install("phyloseq")
#Phyloseq nos va a permitir importar, almacenar, analizar y mostrar gráficamente datos de secuenciación filogenéticos complejos que ya se han agrupado en unidades taxonómicas operativas (OTU).
library(phyloseq)

#Descarga de la libreria de "MGnifyR"
install.packages("MGnifyR")
devtools::install_github("beadyallen/MGnifyR")
#MGnifyR sirve para buscar y recuperar datos del recurso Metagenomica de EBI.
#Cargamos la libreria "MGnifyR"
library(MGnifyR)


### DATA ###

#En MGnify se encontro una base de datos sobre "Microbioma intestinal de osos en hibernacion"
#Las muestras se recolectaron del colon 17 osos negros (Ursus americanus) en hibernacion en Minessota, EUA.
#Se identificaron las bacterias presentes en la microbioma utilizando 16S aplicon data
#https://www.ebi.ac.uk/metagenomics/studies/MGYS00003951#overview
#Study MGYS00003951
mgclnt <- mgnify_client(usecache = T, cache_dir = '/tmp/MGnify_cache')
accession_list <- mgnify_analyses_from_studies(mgclnt, "MGYS00003951", usecache = T)
meta_dataframe <- mgnify_get_analyses_metadata(mgclnt, accession_list, usecache = T )
oso <- mgnify_get_analyses_phyloseq(mgclnt, meta_dataframe$analysis_accession, usecache = T)


#Al tener el analisis de phyloseq de la base de datos de oso nos dimos cuenta que tiene varias NA en genero y había menos en familia, por lo que se creo un objeto de cada uno 
#Borramos de la base de datos aquellos que en Familia tengan NA
oso_limpio_F_1 <- subset_taxa(oso, !(Family %in% c(NA)))
oso_limpio_F_1
#Con la funcion view me permite ver las muestras y si nivel taxonomico 
View(tax_table(oso_limpio_F_1))
#Con save creamos un objeto llamado oso_limpio_familia para que así otras puedan cargar la base de datos sin la necesidad de hacer todo lo anterior, solo con la función load y nombre del objeto
save(oso_limpio_F_1, file="Data/oso_limpio_1_familia")

#Borramos de la base de datos aquellos que en Genero tengan NA
oso_limpio_G_1 <- subset_taxa(oso, !(Genus %in% c(NA)))
oso_limpio_G_1
#Con la funcion view me permite ver las muestras y si nivel taxonomico 
View(tax_table(oso_limpio_G_1))
#Con save creamos un objeto llamado oso_limpio_familia para que así otras puedan cargar la base de datos sin la necesidad de hacer todo lo anterior, solo con la función load y nombre del objeto
save(oso_limpio_G_1, file="Data/oso_limpio_1_genero")

#En MGnify se encontro una base de datos de "Song Colorado freshwater fish"
#Las muestras son de baba y tripa de diferentes especies de peces de agua dulce del rio Colorado, EUA.
#https://www.ebi.ac.uk/metagenomics/studies/MGYS00003748#overview
#Study MGYS00003748
mgclnt <- mgnify_client(usecache = T, cache_dir = '/tmp/MGnify_cache')
accession_list <- mgnify_analyses_from_studies(mgclnt, "MGYS00003748", usecache = T)
meta_dataframe <- mgnify_get_analyses_metadata(mgclnt, accession_list, usecache = T )
peces <- mgnify_get_analyses_phyloseq(mgclnt, meta_dataframe$analysis_accession, usecache = T)
peces


#Al tener el analisis de phyloseq de la base de datos de peces nos dimos cuenta que tiene varias NA en genero y había menos en familia, por lo que se creo un objeto de cada uno 
#Borramos de la base de datos aquellos que en Familia tengan NA
pez_limpio_F_1 <- subset_taxa(peces, !(Family %in% c(NA)))
pez_limpio_F_1
View(tax_table(pez_limpio_F_1))
#Con save creamos un objeto llamado oso_limpio_familia para que así otras puedan cargar la base de datos sin la necesidad de hacer todo lo anterior, solo con la función load y nombre del objeto
save(pez_limpio_F_1, file="Data/pez_limpio_1_familia")

#Borramos de la base de datos aquellos que en Genero tengan NA
pez_limpio_G_1 <- subset_taxa(peces, !(Genus %in% c(NA)))
pez_limpio_G_1
#Con la función view me permite ver las muestras y si nivel taxonomico 
View(tax_table(pez_limpio_G_1))
#Con save creamos un objeto llamado oso_limpio_familia para que así otras puedan cargar la base de datos sin la necesidad de hacer todo lo anterior, solo con la función load y nombre del objeto
save(pez_limpio_G_1, file="Data/pez_limpio_1_genero")

#Todas las bases de datos que se utilizaron en el proyecto se encuentran en la carpeta "data"

##### Ignora esto ahorita jajaj
otu_table(peces)
sample_data(oso)
tax_table(oso)
View(tax_table(oso))
####Link de phyoseq:
sample_variables(oso)
plot_bar(oso, fill = "Genus")

##### PENDIENTES: #######
##### Crear objetos en DADA2 #####
##### Checar analisis que vienen en el analisis de DADA2 #######
##### Analisis comparativo? #######



# Encontramos algunas bases de datos en zenodo:

# Una es de osos negros y su microbiota en dos organos
  # "Data From: Wild black bears harbor simple gut microbial communities with little difference between the jejunum and colon"

# Y la otra es una base de datos de osos cafes a los que tambien les estudian la microbiota
  # "Correlating gut microbial membership to brown bear health metrics"

###### zen4R #####

# zen4R es una funcion en R que permite descargar los datos desde zenodo a tu computadora
# se necesita del doi que viene en la pagina de zenodo para que identifique que base de datos se descargara

# Se descarga zen4R y se carga:
install.packages("zen4R")
library(zen4R)

getwd()

  # Base de osos negros:
download_zenodo("10.5281/zenodo.4060480", path = "")
    # No me permite descargarla, me pone un error, pero no me lo explica
    # Suponemos que es la base porque la otra si la carga normal y es el mismo codigo (solo cambia el doi)

  # Base de osos cafes:
osos_cafe <- download_zenodo("10.5281/zenodo.5759055", path = "")
    # "path = " ", la tengo que especificar, ahorita no puedo ver mis carpetas ni en donde estoy porque se esta actualizando mi R, pero ahi se especifica en donde quieres los archivos
    # NOTA: cree una nueva carpeta llamada "oso_cafe_zenodo" para que ahi meta los archivos y no se confundan con los demas

# Para nuestra fortuna, dentro de los archivos que se descargan, vienen un script (CH2) que explica como importar los datos, convertirlos en phyloseq y analizarlos
# Como no vamos a usar todo ese codigo, pondremos aqui y editaremos las lineas que nos son utiles



##### CH2 (Visalizar y editar los datos de la base de osos cafes) #####

# Hay que instalar y cargar todo esto:
# No he visto que paquetes son los que vamos a usar y cuales no, asi que puse todos los que dice:

packageVersion(qiime2R)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("microbiome")
library(microbiome)
library("BiocManager")
## data analysis
install.packages("remotes")
remotes::install_github("jbisanz/qiime2R")
library("remotes")
library(qiime2R) # import data
library(phyloseq) # also the basis of data object. Data analysis and visualization
library(vegan) # some utility tools
library(data.table) # alternative to data.frame
library(dplyr) # data handling
library(tidyverse)
library(DT) ## interactive tables
library(ggpubr) ## plotting 
library(ggplot2)
devtools::install_github("leffj/mctoolsr")
library(mctoolsr)
library(picante) ## faith's PD
library(see)
library(Rmisc)## graphing
library(SRS)
library(cowplot)
library(shiny)
library(glue)


# Importacion de los archivos:

# NOTA: cuando mi R termine de actualizarse, pongo bien las ubicaciones de los archivos

# Para leer los archivos tsv descargados:
metadata1 <- read_tsv("brownbearmeta.tsv") 
metadata2 <- read_tsv("metadata.tsv")
head(metadata2) # aqui solo se ven las primeras 10 lineas del objeto

metadata<-full_join(metadata1, metadata2)
head(metadata)
  # full_join: agregan columnas de "y" a "x", haciendo coincidir las observaciones en funcion de las claves.
  # Full_join mantiene todas las observaciones en "x" y "y"

# Para leer los archivos qza:
SVs<-read_qza("clean-brownbear-table-unassigned_Unknown_Arch-rm.qza")
head(SVs)

taxonomy<-read_qza("brownbear-taxonomy_renamed.qza")
taxtable<-taxonomy$data %>% as_tibble() %>% separate(Taxon, sep=";", c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) #convert the table into a tabular split version
tree<-read_qza("filter-rooted-brownbear-tree.qza")

metadata$`Body Fat`<- cut(metadata$`Body Fat (%)`,3)
metadata$NetBodyMass<-cut(metadata$NetBodyMasskgs,3)
metadata$`Fat Mass`<-cut(metadata$`Fat Mass (kg)`,3)
metadata$`Lean Mass`<-cut(metadata$`Lean Mass (kg)`,3)

# Crear un objeto phyloseq:
phy_obj<-phyloseq(
  otu_table(SVs$data, taxa_are_rows = TRUE), 
  phy_tree(tree$data), 
  tax_table(as.data.frame(taxtable) %>% dplyr::select(-Confidence) %>% column_to_rownames("Feature.ID") %>% as.matrix()), #moving the taxonomy to the way phyloseq wants it
  sample_data(metadata %>% as.data.frame() %>% column_to_rownames("SampleID")))


#Visualizar la tabla con los datos
datatable(tax_table(phy_obj))


# NOTA: me parece que lo que nos interesa viene hasta la linea 155 del script CH2 porque despues de esa linea viene lo de las diversidades