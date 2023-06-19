# ANALISIS A REALIZAR #

#Curva de rarefacción
#Oso limpio Familia
library(vegan)
par(mfrow=c(1,2))
otu.rare= otu_table(oso_limpio_F_1)
otu.rare= as.data.frame(otu.rare)
sample_names=rownames(otu.rare)
otu_rarecurve<- rarecurve(otu.rare, setp=10000, label=T)
save(otu_rarecurve, file = "Plots/curva_rarefacción_oso")
load("Plots/curva_rarefacción_oso")


#Pez limpio Familia
otu.rare= otu_table(pez_limpio_F_1)
otu.rare= as.data.frame(otu.rare)
sample_names=rownames(otu.rare)
otu_rarecurve_p<- rarecurve(otu.rare, setp=10000, label=T)
save(otu_rarecurve_p, file = "Plots/curva_rarefacción_pez")

#install.packages("gtExtras")
library(gtExtras)
library(gt)

#Diversidad alfa
  # Shannon
  # Simpson
library("microbiome")
Tab<- evenness(oso_limpio_F_1, c("pielou","simpson"))
Tab
View(Tab)
class(Tab)

# Esto crea una tabla, pero solo muestra los primeros 6 valores y no en un orden especifico.
head(Tab) %>%
  gt() %>% 
  gt_theme_guardian()


Tab2<- evenness(oso_limpio_G_1, c("pielou","simpson"))
Tab2

Tab3<- evenness(pez_limpio_F_1, c("pielou","simpson"))
Tab3

Tab4<- evenness(pez_limpio_G_1, c("pielou","simpson"))
Tab4

#Diversidad beta
oso.familia.bray<- ordinate(oso_limpio_F_1, method="MDS", distance="bray")
oso.familia.bray
oso_limpio_F_1
View(sample_data(oso_limpio_F_1))
#Por alguna razon ya corre estos
oso.genero.bray<- ordinate(oso_limpio_G_1, method="MDS", distance="bray")
oso.genero.bray
oso_limpio_G_1
View(sample_data(oso_limpio_G_1))
pez.famialia.bray<- ordinate(pez_limpio_F_1, method="MDS", distance="bray")
pez.famialia.bray
pez.genero.bray<- ordinate(pez_limpio_G_1, method="MDS", distance="bray")
pez.genero.bray


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

# Graficas de abundancias:

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


# Para los peces (familia):

pez1familia <- get_taxadf(obj=pez_limpio_F_1, taxlevel=7, type = "others")

# The 30 most abundant taxonomy will be visualized by default (parameter `topn=30`):#
plot_pez1familia <- ggbartax(obj=pez1familia) +
  xlab("nuevacol") +
  ylab("relative abundance (%)") +
  scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(31))) +
  guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5))

plot_pez1familia

oso1genero <- get_taxadf(obj=oso_limpio_G_1, taxlevel=7, type = "others")

# The 30 most abundant taxonomy will be visualized by default (parameter `topn=30`):#
plot_oso1genero <- ggbartax(obj=oso1genero) +
  xlab("nuevacol") +
  ylab("relative abundance (%)") +
  scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(31))) +
  guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5))

plot_oso1genero

pez1genero <- get_taxadf(obj=pez_limpio_G_1, taxlevel=7, type = "others")

# The 30 most abundant taxonomy will be visualized by default (parameter `topn=30`):#
plot_pez1genero <- ggbartax(obj=pez1genero) +
  xlab("nuevacol") +
  ylab("relative abundance (%)") +
  scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(31))) +
  guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5))

plot_pez1genero


#Cambie los valores de detection y prevalence para que saliera bien
p0<- subset_samples(taxas_juntas)
p0<-core(p0, detection= 0.1/100, prevalence=40/100)
prevalencia_familia<-plot_taxa_prevalence(p0, "Family", detection = 0.1/100)
save(prevalencia_familia, file = "Plots/prevalencia_familia")

red<-plot_net(taxas_juntas, type = "taxa", point_label = "Family", point_size = 2, point_alpha = 0.5, maxdist = 0.5, color = "Family", distance = "bray", laymeth = "auto")
save(red, file = "Plots/red")

p1<- subset_samples(taxas_juntas)
p1<-core(p1, detection= 0.1/100, prevalence=40/100)
prevalencia_genero<-plot_taxa_prevalence(p0, "Genus", detection = 0.1/100)
save(prevalencia_familia, file = "Plots/prevalencia_genero")

red_1<-plot_net(taxas_juntas, type = "taxa", point_label = "Genus", point_size = 2, point_alpha = 0.5, maxdist = 0.5, color = "Genus", distance = "bray", laymeth = "auto")
save(red_1, file = "Plots/red_genero")

#####
pseq.rel <- microbiome::transform(taxas_juntas, "compositional")
head(prevalence(pseq.rel, detection = 1/100, sort = TRUE))
head(prevalence(pseq.rel, detection = 1/100, sort = TRUE, count = TRUE))
core.taxa.standard <- core_members(pseq.rel, detection = 0, prevalence = 50/100)
pseq.core <- core(pseq.rel, detection = 0, prevalence = .5)
pseq.core2 <- aggregate_rare(pseq.rel, "Family", detection = 0, prevalence = .5)
core.taxa <- taxa(pseq.core)
det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)
#ggplot(d) + geom_point(aes(x, y)) + scale_x_continuous(trans="log10", limits=c(NA,1))

#Core microbiome
core<- plot_core(pseq.rel, 
          prevalences = prevalences, 
          detections = det, 
          plot.type = "lineplot") + 
  xlab("Relative Abundance (%)")
save(core, file = "Plots/core_microbiome")




###Heatmap
install.packages("RColorBrewer")
library(RColorBrewer)
install.packages("reshape")
library(reshape)

prevalences <- seq(.05, 1, .05)

detections <- round(10^seq(log10(0.01), log10(.2), length = 9), 3)

# Also define gray color palette
gray <- gray(seq(0,1,length=5))

#Added pseq.rel, I thin... must be checked if it was in the the rednred version,; where it is initialized
#pseq.rel<- microbiome::transform(pseq, 'compositional')
#min-prevalence gets the 100th highest prevalence
install.packages("remotes")
remotes::install_github("microbiome/microbiome")
library(microbiome)
p <- plot_core(pseq.rel,
               plot.type = "heatmap", 
               colours = gray,
               prevalences = prevalences, 
               detections = detections, 
               min.prevalence = prevalence(pseq.rel, sort = TRUE)[100]) +
  labs(x = "Detection Threshold\n(Relative Abundance (%))") +
  
  #Adjusts axis text size and legend bar height
  theme(axis.text.y= element_text(size=8, face="italic"),
        axis.text.x.bottom=element_text(size=8),
        axis.title = element_text(size=10),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10))

print(p)
save(p, file = "Plots//heat_map")
detections <- seq(from = 50, to = round(max(abundances(taxas_juntas))/10, -1), by = 100)

p1 <- plot_core(taxas_juntas, plot.type = "heatmap",
               prevalences = prevalences,
               detections = detections,
               colours = rev(brewer.pal(5, "Spectral")),
               min.prevalence = .2, horizontal = TRUE) +
  theme(axis.text.x= element_text(size=8, face="italic", hjust=1),
        axis.text.y= element_text(size=8),
        axis.title = element_text(size=10),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10))

print(p1)
save(p1, file= "Plots/heatmap_mas_cosas")

#Histograma de abundancias
plot_density(pseq, "Dialister") + ggtitle("Absolute abundance")

# Same with log10 compositional abundances
#No se si snos sirve esto, en base a la prevalencia que dio el heat map puse el tax 55775, porque es el que tenía la prevalencia más alta
x <- microbiome::transform(taxas_juntas, "compositional")
tax <- "55775"
plot_density(x, tax, log10 = TRUE) +
  ggtitle("Relative abundance") +
  xlab("Relative abundance (%)")
otu_table(taxas_juntas)

#Boxplot de abundancias
p3 <- boxplot_abundance(taxas_juntas, x = "nuevacol", y = "55775") + scale_y_log10()
print(p3)
save(p3, file = "Plots/abundancia_55775")

#Para ver la variación de la microbiota
#Falta ver como quitar los 0
p4 <- plot_landscape(pseq.rel, method = "NMDS", distance = "bray", col = "nuevacol", size = 3)
print(p4)

