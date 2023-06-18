# ANALISIS A REALIZAR #

#Curva de rarefacción
library(vegan)
otu.rare= otu_table(oso_limpio_F_1)
otu.rare= as.data.frame(otu.rare)
sample_names=rownames(otu.rare)
otu_rarecurve<- rarecurve(otu.rare, setp=10000, label=T)
save(otu_rarecurve, file = "Plots/curva_rarefacción_oso")

#Diversidad alfa
  # Shannon
  # Simpson
Tab<- evenness(oso_limpio_F_1, c("pielou","simpson"))
Tab
  # Chao1 (tal vez)

#Diversidad beta
psd5.bray<- ordinate(oso_limpio_F_1, method="MDS", distance="bray")
psd5.bray

par(mfrow=c(2,1))
hist(oso_limpio_F_1$shannon, main="Shannon diversity", xlab="", breaks=10)
hist(pez_limpio_F_1$shannon, main="Shannon diversity", xlab="", breaks=10)


p.shannon<- boxplot_alpha(taxas_juntas,
                          index="Shannon",
                          x_var="nuevacol",
                          fill.color= c(DATA_OSO1="blue", DATA_PEZ1="purple"))

p.shannon<- p.shannon + theme_minimal()+
  labs(x="data", y="Shannon diversity")+
  theme(axis.text=element_text(size=12),
        axis.title = element_text(size=16),
        legend.text = element_text(size=12),
        legend.title = element_text(size = 16))
save(p.shannon, file = "Plots/shannon_boxplot")
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
#load("Data/oso_limpio_2_familia")
#load("Data/oso_limpio_2_genero")
load("Data/pez_limpio_1_familia")
load("Data/pez_limpio_1_genero")


# Adicion de una columna:
sample_data(oso_limpio_F_1)$nuevacol<-"DATA_OSO1"
#sample_data(oso_limpio_F_2)$nuevacol<-"DATA_OSO2"
sample_data(oso_limpio_G_1)$nuevacol<-"DATA_OSO1"
#sample_data(oso_limpio_G_2)$nuevacol<-"DATA_OSO2"
sample_data(pez_limpio_F_1)$nuevacol<-"DATA_PEZ1"
sample_data(pez_limpio_G_1)$nuevacol<-"DATA_PEZ1"

#OSOS
View(sample_data(oso_limpio_F_1))
View(sample_data(oso_limpio_G_1))

#PECES
View(sample_data(pez_limpio_F_1))
View(sample_data(pez_limpio_G_1))




# Para juntar los phyloseq de ambas bases de familias:
base_juntas<-merge_phyloseq(oso_limpio_F_1, pez_limpio_F_1)
View(sample_data(base_juntas))



taxas_juntas<- tax_glom(base_juntas, taxrank="Family", NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
taxas_juntas
save(taxas_juntas, file="Data/taxas_juntas")
load("Data/taxas_juntas")

View(sample_data(taxas_juntas))

#Diagrama de venn de familia
pseq.rel <- microbiome::transform(taxas_juntas, "compositional")
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
venn_familia<-plot(venn(list_core),
     fills = mycols)
save(venn_familia, file="Plots/venn_familia")




####### PARA GENERO ########

# Para juntar los phyloseq de ambas bases de familias:
base_juntas_g<-merge_phyloseq(oso_limpio_G_1, pez_limpio_G_1)
View(sample_data(base_juntas_g))



taxas_juntas_g<- tax_glom(base_juntas_g, taxrank="Genus", NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
taxas_juntas_g
save(taxas_juntas_g, file="Data/taxas_juntas_g")
load("Data/taxas_juntas_g")

View(sample_data(taxas_juntas_g))


pseq.rel_g <- microbiome::transform(taxas_juntas_g, "compositional")
disease_states_g <- unique(as.character(meta(pseq.rel_g)$nuevacol))
print(disease_states_g)
list_core_g <- c() # an empty object to store information

for (n in disease_states_g){ # for each variable n in DiseaseState
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub_g <- subset_samples(pseq.rel_g, nuevacol == n) # Choose sample from DiseaseState by n
  
  core_m_g <- core_members(ps.sub_g, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in atleast 90% samples 
                         prevalence = 0.00,
                         include.lowest = FALSE)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core_g[[n]] <- core_m_g # add to a list core taxa for each group.
  #print(list_core)
}
print(list_core_g)
library(eulerr)
mycols_g <- c(DATA_OSO1="#d6e2e9", DATA_PEZ1="#fcf5c7") 
venn_genero<-plot(venn(list_core_g),
     fills = mycols_g)
save(venn_genero, file="Plots/venn_genero")

getwd()

#####
# Pero creo que no se puede porque son graficas muy grandes :(



# Para oso_limpio_F_1:


# En taxlevel la verdad no entiendo que cambia, solo se que el maximo es de 7 y las graficas si salen bien diferentes si lo cambias
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MicrobiotaProcess")
library(MicrobiotaProcess)
par(mfrow=c(2,1))
oso1familia <- get_taxadf(obj=oso_limpio_F_1, taxlevel=7, type = "others")

# The 30 most abundant taxonomy will be visualized by default (parameter `topn=30`):#
plot_oso1familia <- ggbartax(obj=oso1familia) +
  xlab("nuevacol") +
  ylab("relative abundance (%)") +
  scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(31))) +
  guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5))

plot_oso1familia

pez1familia <- get_taxadf(obj=pez_limpio_F_1, taxlevel=7, type = "others")

# The 30 most abundant taxonomy will be visualized by default (parameter `topn=30`):#
plot_pez1familia <- ggbartax(obj=pez1familia) +
  xlab("nuevacol") +
  ylab("relative abundance (%)") +
  scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(31))) +
  guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5))

plot_pez1familia

#Cambie los valores de detection y prevalence para que saliera bien
p0<- subset_samples(taxas_juntas)
p0<-core(p0, detection= 0.4/100, prevalence=40/100)
plot_taxa_prevalence(p0, "Family", detection = 0.1/100)
