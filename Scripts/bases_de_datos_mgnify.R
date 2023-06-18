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

# Una de ellas es una base de datos de osos cafes a los que tambien les estudian la microbiota
  # "Correlating gut microbial membership to brown bear health metrics"
  # Datos disponibles aqui: https://zenodo.org/record/5759055#.YbuUqvHMLVY
  # El doi es: 10.5281/zenodo.5759055


###### zen4R #####

# zen4R es una funcion en R que permite descargar los datos desde zenodo a tu computadora
# se necesita del doi que viene en la pagina de zenodo para que identifique que base de datos se descargara

# Se descarga zen4R y se carga:
install.packages("zen4R")
library(zen4R)


#

# Esta no la vamos a usar (?) #
  # Base de osos negros:
download_zenodo("10.5281/zenodo.4060480", path = "Data/oso_cafe_zenodo/")
    # No me permite descargarla, me pone un error, pero no me lo explica
    # Suponemos que es la base porque la otra si la carga normal y es el mismo codigo (solo cambia el doi)

#


  # Base de osos cafes:
osos_cafe <- download_zenodo("10.5281/zenodo.5759055", path = "Data/oso_cafe_zenodo/")
    # "path = " ", la tengo que especificar, ahorita no puedo ver mis carpetas ni en donde estoy porque se esta actualizando mi R, pero ahi se especifica en donde quieres los archivos
    # NOTA: cree una nueva carpeta llamada "oso_cafe_zenodo" para que ahi meta los archivos y no se confundan con los demas

# Para nuestra fortuna, dentro de los archivos que se descargan, vienen un script (CH2) que explica como importar los datos, convertirlos en phyloseq y analizarlos
# Como no vamos a usar todo ese codigo, pondremos aqui y editaremos las lineas que nos son utiles



##### CH2 (Visalizar y editar los datos de la base de osos cafes) #####

# Hay que instalar y cargar todo esto:
packageVersion(qiime2R)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("microbiome")
#library(microbiome)
library("BiocManager")
## data analysis
library(qiime2R) # import data
library(phyloseq) # also the basis of data object. Data analysis and visualization
library(vegan) # some utility tools
library(data.table) # alternative to data.frame
library(dplyr) # data handling
#install.packages("tidyverse")
library(tidyverse)
#install.packages("DT")
library(DT) ## interactive tables


# Para importar el tsv con datos que nos interesan:
#install.packages("readr")
library(readr)

# Importacion de los archivos:

# NOTA: cuando mi R termine de actualizarse, pongo bien las ubicaciones de los archivos

# Para leer los archivos tsv descargados:
metadata1 <- read_tsv("Data/oso_cafe_zenodo/brownbearmeta.tsv") 
metadata2 <- read_tsv("Data/oso_cafe_zenodo/metadata.tsv")
head(metadata2) # aqui solo se ven las primeras 10 lineas del objeto

metadata<-full_join(metadata1, metadata2)
head(metadata)
  # full_join: agregan columnas de "y" a "x", haciendo coincidir las observaciones en funcion de las claves.
  # Full_join mantiene todas las observaciones en "x" y "y"

# Para leer los archivos qza:
SVs<-read_qza("Data/oso_cafe_zenodo/clean-brownbear-table-unassigned_Unknown_Arch-rm.qza")
head(SVs)

taxonomy<-read_qza("Data/oso_cafe_zenodo/brownbear-taxonomy_renamed.qza")
taxtable<-taxonomy$data %>% as_tibble() %>% separate(Taxon, sep=";", c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) #convert the table into a tabular split version
tree<-read_qza("Data/oso_cafe_zenodo/filter-rooted-brownbear-tree.qza")

metadata$`Body Fat`<- cut(metadata$`Body Fat (%)`,3)
metadata$NetBodyMass<-cut(metadata$NetBodyMasskgs,3)
metadata$`Fat Mass`<-cut(metadata$`Fat Mass (kg)`,3)
metadata$`Lean Mass`<-cut(metadata$`Lean Mass (kg)`,3)

# Crear un objeto phyloseq:
oso_cafe<-phyloseq(
  otu_table(SVs$data, taxa_are_rows = TRUE), 
  phy_tree(tree$data), 
  tax_table(as.data.frame(taxtable) %>% dplyr::select(-Confidence) %>% column_to_rownames("Feature.ID") %>% as.matrix()), #moving the taxonomy to the way phyloseq wants it
  sample_data(metadata %>% as.data.frame() %>% column_to_rownames("SampleID")))


#Visualizar la tabla con los datos
datatable(tax_table(oso_cafe))
  # Esto te da una tabla interactiva.

# NOTA: me parece que lo que nos interesa viene hasta la linea 155 del script CH2 porque despues de esa linea viene lo de las diversidades

#Nuestro objeto phyloseq con la base de datos de zenodo se llama:
oso_cafe

#Con la funcion view me permite ver las muestras y su nivel taxonomico 
View(tax_table(oso_cafe))
#Con save creamos un objeto llamado oso_cafe_limpio para que asi otras puedan cargar la base de datos sin la necesidad de hacer todo lo anterior, solo con la funcionload y nombre del objeto
save(oso_cafe, file="Data/oso_cafe_limpio")


oso_limpio_F_2 <- subset_taxa(oso_cafe, !(Family %in% c(NA)))
oso_limpio_F_2
View(tax_table(oso_limpio_F_2))
#Con save creamos un objeto llamado oso_limpio_familia para que así otras puedan cargar la base de datos sin la necesidad de hacer todo lo anterior, solo con la función load y nombre del objeto
save(oso_limpio_F_2, file="Data/oso_limpio_2_familia")

#Borramos de la base de datos aquellos que en Genero tengan NA
oso_limpio_G_2 <- subset_taxa(oso_cafe, !(Genus %in% c(NA)))
oso_limpio_G_2
#Con la función view me permite ver las muestras y si nivel taxonomico 
View(tax_table(oso_limpio_G_2))
#Con save creamos un objeto llamado oso_limpio_familia para que así otras puedan cargar la base de datos sin la necesidad de hacer todo lo anterior, solo con la función load y nombre del objeto
save(oso_limpio_G_2, file="Data/oso_limpio_2_genero")



# Prueba de cargar las bases de datos con load:
getwd()
setwd("C:/Users/Luis Diego/Documents/Proyecto_Final_Genomica")
load("Data/oso_limpio_1_familia")
oso_limpio_F_1

#Encontre esta función para hacer más chiquito nuestro objeto phyloseq
#El objeot phyloseq tenia 85 samples y cambiamos a que escogiera aleatoriamente 10 samples

# Para osos negros (familia): 
sample_oso <- function(oso_limpio_F_1, FUN = sample,...){
  ids <- sample_names(oso_limpio_F_1)
  sampled_ids <- FUN(ids, ...)
  oso_limpio_F_1 <- prune_samples(sampled_ids, oso_limpio_F_1)
  return(oso_limpio_F_1)
}

oso_corto_F_1<- sample_oso(oso_limpio_F_1, size=10)
oso_corto_F_1
save(oso_corto_F_1, file = "Data/oso_corto_F_1")


# Para osos negros (genero): 

load("Data/oso_limpio_1_genero")
oso_limpio_G_1
#El objeot phyloseq tenia 85 samples y cambiamos a que escogiera aleatoriamente 10 samples
sample_osoo_1 <- function(oso_limpio_G_1, FUN = sample,...){
  ids <- sample_names(oso_limpio_G_1)
  sampled_ids <- FUN(ids, ...)
  oso_limpio_G_1 <- prune_samples(sampled_ids, oso_limpio_G_1)
  return(oso_limpio_G_1)
}

oso_corto_G_1<- sample_osoo_1(oso_limpio_G_1, size=10)
oso_corto_G_1
save(oso_corto_G_1, file = "Data/oso_corto_G_1")


# Para peces (familias):

load("Data/pez_limpio_1_familia")
#Paso de tener 261 samples a 
pez_limpio_F_1
sample_pezz <- function(pez_limpio_F_1, FUN = sample,...){
  ids <- sample_names(pez_limpio_F_1)
  sampled_ids <- FUN(ids, ...)
  pez_limpio_F_1 <- prune_samples(sampled_ids, pez_limpio_F_1)
  return(pez_limpio_F_1)
}

pez_corto_F_1<- sample_pezz(pez_limpio_F_1, size=10)
pez_corto_F_1
save(pez_corto_F_1, file = "Data/pez_corto_F1")


# Para peces (generos):

load("Data/pez_limpio_1_genero")
#Paso de tener 261 samples a 10 samples
pez_limpio_G_1
sample_pezz_1 <- function(pez_limpio_G_1, FUN = sample,...){
  ids <- sample_names(pez_limpio_G_1)
  sampled_ids <- FUN(ids, ...)
  pez_limpio_G_1 <- prune_samples(sampled_ids, pez_limpio_G_1)
  return(pez_limpio_G_1)
}

pez_corto_G_1<- sample_pezz_1(pez_limpio_G_1, size=10)
pez_corto_G_1
save(pez_corto_G_1, file = "Data/pez_corto_1")


# Para osos cafes (familia):

load("Data/oso_cafe_limpio")
load("Data/oso_limpio_2_familia")
#Paso de tener 66 samples a 10
oso_limpio_F_2
sample_oso_2 <- function(oso_limpio_F_2, FUN = sample,...){
  ids <- sample_names(oso_limpio_F_2)
  sampled_ids <- FUN(ids, ...)
  oso_limpio_F_2 <- prune_samples(sampled_ids, oso_limpio_F_2)
  return(oso_limpio_F_2)
}

oso_corto_F_2<- sample_oso_2(oso_limpio_F_2, size=10)
oso_corto_F_2
save(oso_corto_F_2, file = "Data/oso_corto_F_2")


# Para osos cafes (genero): 

load("Data/oso_limpio_2_genero")
#Paso de tener 66 samples a 10 
oso_limpio_G_2
sample_osoo_2 <- function(oso_limpio_G_2, FUN = sample,...){
  ids <- sample_names(oso_limpio_G_2)
  sampled_ids <- FUN(ids, ...)
  oso_limpio_G_2 <- prune_samples(sampled_ids, oso_limpio_G_2)
  return(oso_limpio_G_2)
}

oso_corto_G_2<- sample_osoo_2(oso_limpio_G_2, size=10)
save(oso_corto_G_2,file = "Data/oso_corto_G_2")

# Para cargar las "nuevas" bases de datos:
load("Data/oso_corto_F_1")
load("Data/oso_corto_G_1")
load("Data/oso_corto_F_2")
load("Data/oso_corto_G_2")
load("Data/pez_corto_F1")
load("Data/pez_corto_1")

# Para juntar los phyloseq:
#Se juntaron los phyloseq de familia
data_familia_corto<- merge_phyloseq(oso_corto_F_1, pez_corto_F_1, oso_corto_F_2)
data_familia_corto

#Se juntaron los phyloseq de genero
data_genero_corto<- merge_phyloseq(oso_corto_G_1, pez_corto_G_1, oso_corto_G_2)
data_genero_corto



# Esto hace un plor con las muestras mas abundantes:
plot_bar(data_familia_corto)


# Creo que con esto se pueden ver los mas abundantes (de nuevo):

library(ggplot2)
library(MicrobiotaProcess)

classtaxa <- get_taxadf(obj=data_familia_corto, taxlevel=7, type = "others")

# The 30 most abundant taxonomy will be visualized by default (parameter `topn=30`):
pclass <- ggbartax(obj=classtaxa) +
  xlab(NULL) +
  ylab("relative abundance (%)") +
  scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(31))) +
  guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5))

pclass

# Y con la grafica anterior salen las mismas muestras mas abundantes.


View(otu_table(data_familia_corto))
View(phy_tree((data_familia_corto)))
View(sample_data(data_familia_corto))



# Creo que con esto se pueden ver los mas abundantes (de nuevo):

library(ggplot2)
library(MicrobiotaProcess)

# En taxlevel la verdad no entiendo que cambia, solo se que el maximo es de 7 y las graficas si salen bien diferentes si lo cambias
osofamilia <- get_taxadf(obj=oso_limpio_F_1, taxlevel=3, type = "others")

# The 30 most abundant taxonomy will be visualized by default (parameter `topn=30`):
plot_osofamilia <- ggbartax(obj=osofamilia) +
  xlab(NULL) +
  ylab("relative abundance (%)") +
  scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(31))) +
  guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5))

plot_osofamilia


load("Data/oso_limpio_1_familia")
load("Data/oso_limpio_1_genero")
load("Data/oso_limpio_2_familia")
load("Data/oso_limpio_2_genero")
load("Data/pez_limpio_1_familia")
load("Data/pez_limpio_1_genero")


# Adicion de una columna:
sample_data(oso_limpio_F_1)$nuevacol<-"DATA_OSO1"
sample_data(oso_limpio_F_2)$nuevacol<-"DATA_OSO2"
sample_data(oso_limpio_G_1)$nuevacol<-"DATA_OSO1"
sample_data(oso_limpio_G_2)$nuevacol<-"DATA_OSO2"
sample_data(pez_limpio_F_1)$nuevacol<-"DATA_PEZ1"
sample_data(pez_corto_G_1)$nuevacol<-"DATA_PEZ2"


# Para ver las 6 graficas al mismo tiempo:
par(mfrow=c(3,2))
  # Pero creo que no se puede porque son graficas muy grandes :(

# Para oso_limpio_F_1:

# En taxlevel la verdad no entiendo que cambia, solo se que el maximo es de 7 y las graficas si salen bien diferentes si lo cambias
oso1familia <- get_taxadf(obj=oso_limpio_F_1, taxlevel=7, type = "others")

# The 30 most abundant taxonomy will be visualized by default (parameter `topn=30`):
plot_oso1familia <- ggbartax(obj=oso1familia) +
  xlab(NULL) +
  ylab("relative abundance (%)") +
  scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(31))) +
  guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5))

plot_oso1familia


# Para oso_limpio_G_1:
# En taxlevel la verdad no entiendo que cambia, solo se que el maximo es de 7 y las graficas si salen bien diferentes si lo cambias
oso1genero <- get_taxadf(obj=oso_limpio_G_1, taxlevel=7, type = "others")

# The 30 most abundant taxonomy will be visualized by default (parameter `topn=30`):
plot_oso1genero <- ggbartax(obj=oso1genero) +
  xlab(NULL) +
  ylab("relative abundance (%)") +
  scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(31))) +
  guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5))

plot_oso1genero


# Para oso_limpio_F_2:
# En taxlevel la verdad no entiendo que cambia, solo se que el maximo es de 7 y las graficas si salen bien diferentes si lo cambias
oso2familia <- get_taxadf(obj=oso_limpio_F_2, taxlevel=7, type = "others")

# The 30 most abundant taxonomy will be visualized by default (parameter `topn=30`):
plot_oso2familia <- ggbartax(obj=oso2familia) +
  xlab(NULL) +
  ylab("relative abundance (%)") +
  scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(31))) +
  guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5))

plot_oso2familia



# Para oso_limpio_G_2:
# En taxlevel la verdad no entiendo que cambia, solo se que el maximo es de 7 y las graficas si salen bien diferentes si lo cambias
oso2genero <- get_taxadf(obj=oso_limpio_G_2, taxlevel=7, type = "others")

# The 30 most abundant taxonomy will be visualized by default (parameter `topn=30`):
plot_oso2genero <- ggbartax(obj=oso2genero) +
  xlab(NULL) +
  ylab("relative abundance (%)") +
  scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(31))) +
  guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5))

plot_oso2genero



# Para pez_limpio_f_1
# En taxlevel la verdad no entiendo que cambia, solo se que el maximo es de 7 y las graficas si salen bien diferentes si lo cambias
pez1familia <- get_taxadf(obj=pez_limpio_F_1, taxlevel=7, type = "others")

# The 30 most abundant taxonomy will be visualized by default (parameter `topn=30`):
plot_pez1familia <- ggbartax(obj=pez1familia) +
  xlab(NULL) +
  ylab("relative abundance (%)") +
  scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(31))) +
  guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5))

plot_pez1familia


# Para pez_limpio_g_1
# En taxlevel la verdad no entiendo que cambia, solo se que el maximo es de 7 y las graficas si salen bien diferentes si lo cambias
pez1genero <- get_taxadf(obj=pez_limpio_G_1, taxlevel=7, type = "others")

# The 30 most abundant taxonomy will be visualized by default (parameter `topn=30`):
plot_pez1genero <- ggbartax(obj=pez1genero) +
  xlab(NULL) +
  ylab("relative abundance (%)") +
  scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(31))) +
  guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5))

plot_pez1genero
