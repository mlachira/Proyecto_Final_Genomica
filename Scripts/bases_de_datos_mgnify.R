#Melisa Lachira
#Diana Karina Rangel
#Maria del Carmen Ramírez

#Libreria phyloseq
library(phyloseq)
#Para installar MGnifyR
install.packages("MGnifyR")
devtools::install_github("beadyallen/MGnifyR")
library(MGnifyR)
#Base de datos microbioma de osos
mgclnt <- mgnify_client(usecache = T, cache_dir = '/tmp/MGnify_cache')
accession_list <- mgnify_analyses_from_studies(mgclnt, "MGYS00003951", usecache = T)
meta_dataframe <- mgnify_get_analyses_metadata(mgclnt, accession_list, usecache = T )
oso <- mgnify_get_analyses_phyloseq(mgclnt, meta_dataframe$analysis_accession, usecache = T)
save(oso,file = "")
#Para limpiar la base
oso_limpio <- subset_taxa(oso, !(Genus %in% "NA"))
oso_limpio
otu_table(peces)
sample_data(oso)
tax_table(oso)
View(tax_table(oso))
####Link de phyoseq:
sample_variables(oso)
plot_bar(oso, fill = "Genus")
#Base de datos microbioma peces
mgclnt <- mgnify_client(usecache = T, cache_dir = '/tmp/MGnify_cache')
accession_list <- mgnify_analyses_from_studies(mgclnt, "MGYS00003748", usecache = T)
meta_dataframe <- mgnify_get_analyses_metadata(mgclnt, accession_list, usecache = T )
peces <- mgnify_get_analyses_phyloseq(mgclnt, meta_dataframe$analysis_accession, usecache = T)
peces
#Para limpiar la base
pez_limpio <- subset_taxa(peces, !(Genus %in% "NA"))
pez_limpio
otu_table(peces)
sample_data(peces)
tax_table(peces)
View(tax_table(peces))
sample_variables(oso)
plot_bar(peces, fill = "Genus")



##### PENDIENTES: #######
##### Crear objetos en DADA2 #####
##### Checar analisis que vienen en el analisis de DADA2 #######
##### Analisis comparativo? #######