# ANALISIS A REALIZAR #

#Diversidad alfa
  # Shannon
  # Simpson
  # Chao1 (tal vez)

#Diversidad beta
  # Bray Curtis



library(microbiome)
library("BiocManager")
library(qiime2R) # import data
library(phyloseq) # also the basis of data object. Data analysis and visualization
library(vegan) # some utility tools
library(data.table) # alternative to data.frame
library(dplyr) # data handling
library(tidyverse)
library(DT) ## interactive tables
library(ggplot2)
library(MicrobiotaProcess)
library(microbiome)
#install.packages("eulerr") 
library(eulerr)
library(microbiome)
#devtools::install_github('microsud/microbiomeutilities')
library(microbiomeutilities)



load("Data/oso_limpio_1_familia")
load("Data/oso_limpio_1_genero")
load("Data/oso_limpio_2_familia")
load("Data/oso_limpio_2_genero")
load("Data/pez_limpio_1_familia")
load("Data/pez_limpio_1_genero")


# Adicion de una columna:
sample_data(oso_limpio_F_1)$nuevacol<-"DATA_OSO1"
#sample_data(oso_limpio_F_2)$nuevacol<-"DATA_OSO2"
sample_data(oso_limpio_G_1)$nuevacol<-"DATA_OSO1"
#sample_data(oso_limpio_G_2)$nuevacol<-"DATA_OSO2"
sample_data(pez_limpio_F_1)$nuevacol<-"DATA_PEZ1"
sample_data(pez_limpio_G_1)$nuevacol<-"DATA_PEZ2"

#OSOS
View(sample_data(oso_limpio_F_1))
View(sample_data(oso_limpio_G_1))

#PECES
View(sample_data(pez_limpio_F_1))
View(sample_data(pez_limpio_G_1))



<<<<<<< HEAD

# Para juntar los phyloseq de ambas bases de familias:
intento1<-merge_phyloseq(oso_limpio_F_1, pez_limpio_F_1)
sample_data(intento1)

x1<- tax_glom(intento1, taxrank="Family", NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
x1
save(x1, file="Data/x1")
load("Data/x1")

sample_data(x1)


pseq.rel <- microbiome::transform(x1, "compositional")
disease_states <- unique(as.character(meta(pseq.rel)$nuevacol))
print(disease_states)
list_core <- c() # an empty object to store information

for (n in disease_states){ # for each variable n in DiseaseState
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(pseq.rel, nuevacol == n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in atleast 90% samples 
                         prevalence = 0.00,
                         include.lowest = FALSE)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}
print(list_core)
library(eulerr)
mycols <- c(DATA_OSO1="#d6e2e9", DATA_PEZ1="#fcf5c7") 
plot(venn(list_core),
     fills = mycols)



intento1<- merge_phyloseq(oso_limpio_F_1, pez_limpio_F_1)
library(microbiome)
View(sample_data(intento1))
View(tax_table(intento1))
sample_data(intento1)
=======
  View(sample_data(oso_limpio_F_1))
View(sample_data(oso_limpio_G_1))
View(sample_data(oso_limpio_F_1))


getwd()

# Para juntar los phyloseq de las bases de oso 1 y pez (ambas familias):
intento1<- merge_phyloseq(oso_limpio_F_1, pez_limpio_F_1)

library(microbiome)
View(sample_data(intento1))
View(tax_table(intento1))
View(otu_table(intento1))

>>>>>>> 1e8624a97a891fb40613a921f6df3376d9dfe43c

otu_table(intento1)
otu_table(oso_limpio_F_1)
otu_table(oso_limpio_F_2)
otu_table(pez_limpio_F_1)

data("GlobalPatterns")
otu_table(GlobalPatterns)
View(sample_data(GlobalPatterns))

x1<-tax_glom(intento1, taxrank="Family", NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
x1

View(otu_table(x1))
class(x1) #es un phyloseq ;)


#install.packages("eulerr") 
library(eulerr)
library(microbiome)
#devtools::install_github('microsud/microbiomeutilities')
library(microbiomeutilities)

pseq.rel <- microbiome::transform(x1, "compositional")
pseq.rel

disease_states <- unique(as.character(meta(x1)$nuevacol))
print(disease_states)

View(sample_data(x1))


