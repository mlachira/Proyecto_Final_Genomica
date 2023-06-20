################ PREPARACION DE LAS BASES ########################

# Paquetes necesarios para los analisis:
library(microbiome)
library("BiocManager")
library(qiime2R) 
library(phyloseq) 
library(vegan) 
library(data.table) 
library(dplyr) 
library(tidyverse)
library(ggplot2)
library(MicrobiotaProcess)
library(microbiome)
#install.packages("eulerr") 
library(eulerr)
library(microbiome)
#devtools::install_github('microsud/microbiomeutilities')
library(microbiomeutilities)

# Para cargar las librerias que limpiamos previamente de las NA (en el script "bases_de_datos_mgnify" esta descrito)
load("Data/oso_limpio_1_familia")
load("Data/oso_limpio_1_genero")
load("Data/pez_limpio_1_familia")
load("Data/pez_limpio_1_genero")

# Adicion de una columna
  # Se adiciono una columna a cada base para facilitar las instrucciones en los codigos de los analisis
sample_data(oso_limpio_F_1)$nuevacol<-"DATA_OSO1"
sample_data(oso_limpio_G_1)$nuevacol<-"DATA_OSO1"
sample_data(pez_limpio_F_1)$nuevacol<-"DATA_PEZ1"
sample_data(pez_limpio_G_1)$nuevacol<-"DATA_PEZ1"

# OSOS
View(sample_data(oso_limpio_F_1)) # observamos la base, que ya contiene una nueva columna
View(sample_data(oso_limpio_G_1)) # observamos la base, que ya contiene una nueva columna

# PECES
View(sample_data(pez_limpio_F_1)) # observamos la base, que ya contiene una nueva columna
View(sample_data(pez_limpio_G_1)) # observamos la base, que ya contiene una nueva columna


####### PARA FAMILIA ########

# Para juntar los phyloseq de ambas bases de familias:

base_juntas<-merge_phyloseq(oso_limpio_F_1, pez_limpio_F_1) # se usa la funcion merge_phyloseq para juntarlos
View(sample_data(base_juntas)) # se visualiza su sample_data

# Para juntar especies con la misma taxonomia, en este caso a nivel de familia.
taxas_juntas<- tax_glom(base_juntas, taxrank="Family", NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
taxas_juntas # se imprime el objeto para asegurarse de que se cargo bien
View(sample_data(taxas_juntas)) # se visualiza su sample_data

save(taxas_juntas, file="Data/taxas_juntas") # se guarda la base porque es muy pesada y para no tener que volver a cargarla
load("Data/taxas_juntas") # con esto los demas integrantes pueden cargar la base sin tener que cagar el objeto desde cero.


####### PARA GENERO ########

# Para juntar los phyloseq de ambas bases de generos:
base_juntas_g<-merge_phyloseq(oso_limpio_G_1, pez_limpio_G_1)  # se usa la funcion merge_phyloseq para juntarlos
View(sample_data(base_juntas_g)) # se visualiza su sample_data

# Para juntar especies con la misma taxonomia, en este caso a nivel de familia.
taxas_juntas_g<- tax_glom(base_juntas_g, taxrank="Genus", NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
taxas_juntas_g # se imprime el objeto para asegurarse de que se cargo bien
View(sample_data(taxas_juntas_g)) # se visualiza su sample_data

save(taxas_juntas_g, file="Data/taxas_juntas_g") # se guarda la base porque es muy pesada y para no tener que volver a cargarla
load("Data/taxas_juntas_g") # con esto los demas integrantes pueden cargar la base sin tener que cagar el objeto desde cero.

#Recortamos las datas porque son muy grandes, lo que nos afceta para ciertos analisis
# Para osos negros (familia): 
sample_oso <- function(oso_limpio_F_1, FUN = sample,...){
  ids <- sample_names(oso_limpio_F_1)
  sampled_ids <- FUN(ids, ...)
  oso_limpio_F_1 <- prune_samples(sampled_ids, oso_limpio_F_1)
  return(oso_limpio_F_1)
}

oso_corto_F_1<- sample_oso(oso_limpio_F_1, size=30)
oso_corto_F_1
#Para osos negro (genero):
sample_osoo_1 <- function(oso_limpio_G_1, FUN = sample,...){
  ids <- sample_names(oso_limpio_G_1)
  sampled_ids <- FUN(ids, ...)
  oso_limpio_G_1 <- prune_samples(sampled_ids, oso_limpio_G_1)
  return(oso_limpio_G_1)
}

oso_corto_G_1<- sample_osoo_1(oso_limpio_G_1, size=30)
oso_corto_G_1
save(oso_corto_G_1, file = "Data/oso_corto_G_1")

#Para peces (familia):
sample_pezz <- function(pez_limpio_F_1, FUN = sample,...){
  ids <- sample_names(pez_limpio_F_1)
  sampled_ids <- FUN(ids, ...)
  pez_limpio_F_1 <- prune_samples(sampled_ids, pez_limpio_F_1)
  return(pez_limpio_F_1)
}
  
pez_corto_F_1<- sample_pezz(pez_limpio_F_1, size=29)
pez_corto_F_1
save(pez_corto_F_1, file = "Data/pez_corto_F1")

#Para peces (genero):
sample_pezz_1 <- function(pez_limpio_G_1, FUN = sample,...){
  ids <- sample_names(pez_limpio_G_1)
  sampled_ids <- FUN(ids, ...)
  pez_limpio_G_1 <- prune_samples(sampled_ids, pez_limpio_G_1)
return(pez_limpio_G_1)
}
  
pez_corto_G_1<- sample_pezz_1(pez_limpio_G_1, size=30)
pez_corto_G_1
save(pez_corto_G_1, file = "Data/pez_corto_1")

################## ANALISIS ################################

# 1. Curvas de rarefaccion
# 2. Diversidad alfa
  # 2.1. Indice de shannon normalizado (pielou) y de simpson
  # 2.2. Boxplot de la diversidad de shannon
# 3. Diagramas de Venn
# 4. Graficas de  abundancia relativa
# 5. Redes de co ocurrencia
# 6. Gráfica de prevalencia de taxones a nivel de familia 
# 7. Heatmaps
# 8. Histograma de abundancias
# 9. Boxplot de abundancias

## Curvas de rarefacción ##

#Oso limpio Familia
library(vegan) # se carga la libreria necesaria
par(mfrow=c(1,2)) # para que aparezcan dos plot en una sola imagen

otu.rare= otu_table(oso_limpio_F_1) # nos referimos al otu_table de la base
otu.rare= as.data.frame(otu.rare) # lo convertimos en data frame
sample_names=rownames(otu.rare) # acomodamos los nombres de las filas 
otu_rarecurve<- rarecurve(otu.rare, setp=10000, label=F, ylab = "Familia") # para obtener la curva de rarefaccion
save(otu_rarecurve, file = "Plots/curva_rarefacción_oso") # para guardar la curva en la carpeta "Plots"
load("Plots/curva_rarefacción_oso") # con esto cargamos lo que salvamos.


#Pez limpio Familia
otu.rare= otu_table(pez_limpio_F_1) # nos referimos al otu_table de la base
otu.rare= as.data.frame(otu.rare) # lo convertimos en data frame
sample_names=rownames(otu.rare) # acomodamos los nombres de las filas 
otu_rarecurve_p<- rarecurve(otu.rare, setp=10000, label=F, ylab = "Famlia") # para obtener la curva de rarefaccion
save(otu_rarecurve_p, file = "Plots/curva_rarefacción_pez") # para guardar la curva en la carpeta "Plots"


# Pez limpio genero
otu.rare= otu_table(pez_limpio_G_1) # nos referimos al otu_table de la base
otu.rare= as.data.frame(otu.rare) # lo convertimos en data frame
sample_names=rownames(otu.rare) # acomodamos los nombres de las filas 
otu_rarecurve_p<- rarecurve(otu.rare, setp=10000, label=F, ylab = "Género") # para obtener la curva de rarefaccion
save(otu_rarecurve_p, file = "Plots/curva_rarefacción_pez") # para guardar la curva en la carpeta "Plots"

#Oso limpio genero
otu.rare= otu_table(oso_limpio_G_1) # nos referimos al otu_table de la base
otu.rare= as.data.frame(otu.rare) # lo convertimos en data frame
sample_names=rownames(otu.rare) # acomodamos los nombres de las filas 
otu_rarecurve_p<- rarecurve(otu.rare, setp=10000, label=F, ylab = "Género") # para obtener la curva de rarefaccion
save(otu_rarecurve_p, file = "Plots/curva_rarefacción_pez") # para guardar la curva en la carpeta "Plots"


# Paquetes necesarios para hacer tablas bonitas en R:
#install.packages("gtExtras")
library(gtExtras)
library(gt)


## Diversidad alfa ##
  # Shannon
  # Simpson

library("microbiome") # se carga la libreria necesaria

# Para osos (familia)
Tab<- evenness(oso_limpio_F_1, c("pielou","simpson")) # con la funcion evenness se puede sacar los indices de shannon y de pielou
Tab # se imprime el objeto para ver que corra y como se ve.
class(Tab) # es un data frame

df1_pielou <- Tab[order(Tab$pielou,decreasing=FALSE),] #Ordena la columna de pielou en el orden que le indico.
df1_pielou # se imprime el objeto para ver que corra y como se ve.

df1_simpson <- Tab[order(Tab$simpson, decreasing = FALSE),] #Ordena la columna de pielou en el orden que le indico.
df1_simpson # se imprime el objeto para ver que corra y como se ve.


# Esto crea una tabla, pero solo muestra los primeros 6 valores y no en un orden especifico
# Por ello, se lo indicamos anteriormente y los imprimimos aqui
head(df1_pielou) %>% # para pielou
  gt() %>% 
  gt_theme_nytimes()

head(df1_simpson) %>% # para simpson
  gt() %>% 
  gt_theme_nytimes()


# Para osos (genero)
Tab2<- evenness(oso_limpio_G_1, c("pielou","simpson")) # con la funcion evenness se puede sacar los indices de shannon y de pielou
Tab2 # se imprime el objeto para ver que corra y como se ve.

df2_pielou <- Tab2[order(Tab2$pielou,decreasing=FALSE),] #Ordena la columna de pielou en el orden que le indico.
df2_pielou # se imprime el objeto para ver que corra y como se ve.

df2_simpson <- Tab2[order(Tab2$simpson, decreasing = FALSE),] #Ordena la columna de pielou en el orden que le indico.
df2_simpson # se imprime el objeto para ver que corra y como se ve.


# Esto crea una tabla, pero solo muestra los primeros 6 valores y no en un orden especifico
# Por ello, se lo indicamos anteriormente y los imprimimos aqui
head(df2_pielou) %>% # para pielou
  gt() %>% 
  gt_theme_nytimes()

head(df2_simpson) %>% # para simpson
  gt() %>% 
  gt_theme_nytimes()


# Para peces (familia)
Tab3<- evenness(pez_limpio_F_1, c("pielou","simpson"))# con la funcion evenness se puede sacar los indices de shannon y de pielou
Tab3 # se imprime el objeto para ver que corra y como se ve.

df3_pielou <- Tab3[order(Tab3$pielou,decreasing=TRUE),] #Ordena la columna de pielou en el orden que le indico.
df3_pielou # se imprime el objeto para ver que corra y como se ve.

df3_simpson <- Tab3[order(Tab3$simpson, decreasing = TRUE),] #Ordena la columna de pielou en el orden que le indico.
df3_simpson # se imprime el objeto para ver que corra y como se ve.


# Esto crea una tabla, pero solo muestra los primeros 6 valores y no en un orden especifico
# Por ello, se lo indicamos anteriormente y los imprimimos aqui
head(df3_pielou) %>% # para pielou
  gt() %>% 
  gt_theme_nytimes()

head(df3_simpson) %>% # para simpson
  gt() %>% 
  gt_theme_nytimes()


# Para peces (genero)
Tab4<- evenness(pez_limpio_G_1, c("pielou","simpson")) # con la funcion evenness se puede sacar los indices de shannon y de pielou
Tab4 # se imprime el objeto para ver que corra y como se ve.

df4_pielou <- Tab4[order(Tab4$pielou,decreasing=FALSE),] #Ordena la columna de pielou en el orden que le indico.
df4_pielou # se imprime el objeto para ver que corra y como se ve.

df4_simpson <- Tab4[order(Tab4$simpson, decreasing = FALSE),] #Ordena la columna de pielou en el orden que le indico.
df4_simpson # se imprime el objeto para ver que corra y como se ve.


# Esto crea una tabla, pero solo muestra los primeros 6 valores y no en un orden especifico
# Por ello, se lo indicamos anteriormente y los imprimimos aqui
head(df4_pielou) %>% # para pielou
  gt() %>% 
  gt_theme_nytimes()

head(df4_simpson) %>% # para simpson
  gt() %>% 
  gt_theme_nytimes()


## Boxplot de la diversidad de Shannon ##

p.shannon<- boxplot_alpha(taxas_juntas,
                          index="Shannon",
                          x_var="nuevacol",
                          fill.color= c(DATA_OSO1="blue", DATA_PEZ1="purple"))

  # Nuestro objeto se llama p.shannon
  # realizamos un boxplot_alpha porque es diversidad alfa
  # usamos el phyloseq que contiene ambas bases juntas
  # indicamos que el indice sea de shannon
  # indicamos que el eje de las x sea con nuestra columna llamada "nuevacol"
  # indicamos los colores para cada boxplot.

# Esto es para estilizar la grafica, indicando los tamanos de letras y textos que se mostraran
p.shannon<- p.shannon + theme_minimal()+
  labs(x="data", y="Shannon diversity")+
  theme(axis.text=element_text(size=12),
        axis.title = element_text(size=16),
        legend.text = element_text(size=12),
        legend.title = element_text(size = 16))
p.shannon # se imprime para observar la grafica.

save(p.shannon, file = "Plots/shannon_boxplot") # salvamos el objeto en la carpeta de Plots


## Diversidad beta ##
 # Bray Curtis #

# *Marca error para genero, por lo que se descarta por completo #

# Para oso (familia)
#oso.familia.bray<- ordinate(oso_limpio_F_1, method="MDS", distance="bray") # con esta funcion se hace la prueba de Bray Curtis
#oso.familia.bray # se imprime el objeto para ver que corra y como se ve.
#oso_limpio_F_1 # se corre la base para ver los cambios
#View(sample_data(oso_limpio_F_1)) # se observa el sample_data de la base 

# Para oso (genero)
#oso.genero.bray<- ordinate(oso_limpio_G_1, method="MDS", distance="bray") # con esta funcion se hace la prueba de Bray Curtis
#oso.genero.bray # se imprime el objeto para ver que corra y como se ve.
#oso_limpio_G_1 # se corre la base para ver los cambios
#View(sample_data(oso_limpio_G_1)) # ses observa el sample_data de la base

# Para pez (familia)
#pez.famialia.bray<- ordinate(pez_limpio_F_1, method="MDS", distance="bray") # con esta funcion se hace la prueba de Bray Curtis
#pez.famialia.bray # se imprime el objeto para ver que corra y como se ve.

# Para pez (genero)
#pez.genero.bray<- ordinate(pez_limpio_G_1, method="MDS", distance="bray") # con esta funcion se hace la prueba de Bray Curtis
#pez.genero.bray # se imprime el objeto para ver que corra y como se ve.



##  Diagrama de venn de familia ##

pseq.rel <- microbiome::transform(taxas_juntas, "compositional")
disease_states <- unique(as.character(meta(pseq.rel)$nuevacol))
print(disease_states) # se imprime el objeto 
list_core <- c() # objeto vacio para llenarlo de informacion
  # Se indica que la base es "taxas_juntas", la cual contiene ambas bases (de pez y oso)
  # Tambien especificamos que tome en cuenta la columna llamada "nuevacol"

for (n in disease_states){ # para cada variable n en "DiseaseState"
  ps.sub <- subset_samples(pseq.rel, nuevacol == n) # Escoge una muestra de nuevacol por n
  core_m <- core_members(ps.sub, # ps.sub es seleccionado por phyloseq con solo muestras de g 
                         detection = 0.001, # 0.001 en al menos 90% de las muestras
                         prevalence = 0.00,
                         include.lowest = FALSE)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # Imprime el core taxa identificado en cada nuevacol.
  list_core[[n]] <- core_m # anade a una lista nucleo "list core" de taxa por cada grupo.
}
print(list_core)

library(eulerr) # se necesita esta libreria para los diagramas

# Creacion del diagrama de Venn 
mycols <- c(DATA_OSO1="#d6e2e9", DATA_PEZ1="#fcf5c7") 
venn_familia<-plot(venn(list_core),
     fills = mycols)
save(venn_familia, file="Plots/venn_familia") # se guarda en la carpeta Plots


## Diagrama de Venn de genero ##

pseq.rel_g <- microbiome::transform(taxas_juntas_g, "compositional")
disease_states_g <- unique(as.character(meta(pseq.rel_g)$nuevacol))
print(disease_states_g) # se imprime el objeto 
list_core_g <- c() # objeto vacio para llenarlo de informacion
  # Se indica que la base es "taxas_juntas", la cual contiene ambas bases (de pez y oso)
  # Tambien especificamos que tome en cuenta la columna llamada "nuevacol"

for (n in disease_states_g){ # para cada variable n en "DiseaseState"
  ps.sub_g <- subset_samples(pseq.rel_g, nuevacol == n) # Escoge una muestra de nuevacol por n
  core_m_g <- core_members(ps.sub_g, # ps.sub es seleccionado por phyloseq con solo muestras de g 
                         detection = 0.001, # 0.001 en al menos 90% de las muestras
                         prevalence = 0.00,
                         include.lowest = FALSE)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # Imprime el core taxa identificado en cada nuevacol.
  list_core_g[[n]] <- core_m_g # anade a una lista nucleo "list core" de taxa por cada grupo.
}
print(list_core_g)

library(eulerr) # se necesita esta libreria para los diagramas

# Creacion del diagrama de venn
mycols_g <- c(DATA_OSO1="#d6e2e9", DATA_PEZ1="#fcf5c7") 
venn_genero<-plot(venn(list_core_g),
     fills = mycols_g)
save(venn_genero, file="Plots/venn_genero") # se guarda en la carpeta Plots



## Graficas de abundancias relativas ##

library(MicrobiotaProcess) # se necesita esta libreria
par(mfrow=c(2,1)) # para obtener dos plots en una imagen

# Pagina de (mucha) ayuda:
# https://yulab-smu.top/MicrobiotaProcessWorkshop/articles/MicrobiotaProcessWorkshop.html


# Para los osos (familia):

oso1familia <- get_taxadf(obj=oso_corto_F_1, taxlevel=6, type = "others") # indico la base que quiero y que tome en cuenta todos los niveles taxonomicos.
  # En taxlevel se pone el nivel taxonomico, hay que cambiarlo para ver familias / generos
  # Taxlevel de 6 es para familias

# Las 30 taxonomias mas abundantes se visualizaron por default
  # pero si queremos cambiar ese valor se agrega topn=# en la primera linea.
plot_oso1familia <- ggbartax(obj=oso1familia) +
  xlab("nuevacol") +
  ylab("relative abundance (%)") +
  scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(31))) +
  guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5))
# Arriba se indican las especificaciones para el plot, como el tamano de letra, los textos, los colores.
plot_oso1familia # se imprime el objeto 


# Para los peces (familia):

pez1familia <- get_taxadf(obj=pez_limpio_F_1, taxlevel=6, type = "others") # indico la base que quiero y que tome en cuenta todos los niveles taxonomicos.

# Las 30 taxonomias mas abundantes se visualizaron por default
plot_pez1familia <- ggbartax(obj=pez1familia) +
  xlab("nuevacol") +
  ylab("relative abundance (%)") +
  scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(31))) +
  guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5))
# Arriba se indican las especificaciones para el plot, como el tamano de letra, los textos, los colores.
plot_pez1familia # se imprime el objeto 


# Para los osos (genero):

oso1genero <- get_taxadf(obj=oso_limpio_G_1, taxlevel=7, type = "others") # indico la base que quiero y que tome en cuenta todos los niveles taxonomicos.

# Las 30 taxonomias mas abundantes se visualizaron por default
plot_oso1genero <- ggbartax(obj=oso1genero) +
  xlab("nuevacol") +
  ylab("relative abundance (%)") +
  scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(31))) +
  guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5))
# Arriba se indican las especificaciones para el plot, como el tamano de letra, los textos, los colores.
plot_oso1genero # se imprime el objeto 


# Para los peces (genero):

pez1genero <- get_taxadf(obj=pez_limpio_G_1, taxlevel=7, type = "others") # indico la base que quiero y que tome en cuenta todos los niveles taxonomicos.

# Las 30 taxonomias mas abundantes se visualizaron por default
plot_pez1genero <- ggbartax(obj=pez1genero) +
  xlab("nuevacol") +
  ylab("relative abundance (%)") +
  scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(31))) +
  guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5))
# Arriba se indican las especificaciones para el plot, como el tamano de letra, los textos, los colores.
plot_pez1genero # se imprime el objeto 



## Redes de co ocurrencia ##

#### Para familia ### 

#Cambie los valores de detection y prevalence para que saliera bien
p0<- subset_samples(taxas_juntas)
p0<-core(p0, detection= 0.1/100, prevalence=40/100)
prevalencia_familia<-plot_taxa_prevalence(p0, "Family", detection = 0.1/100)
save(prevalencia_familia, file = "Plots/prevalencia_familia")

red<-plot_net(taxas_juntas, type = "taxa", point_label = "Family", point_size = 2, point_alpha = 0.5, maxdist = 0.5, color = "Family", distance = "bray", laymeth = "auto")
red
save(red, file = "Plots/red")


### Para genero ###

p1<- subset_samples(taxas_juntas)
p1<-core(p1, detection= 0.1/100, prevalence=40/100)
prevalencia_genero<-plot_taxa_prevalence(p0, "Genus", detection = 0.1/100)
save(prevalencia_familia, file = "Plots/prevalencia_genero")

red_1<-plot_net(taxas_juntas, type = "taxa", point_label = "Genus", point_size = 2, point_alpha = 0.5, maxdist = 0.5, color = "Genus", distance = "bray", laymeth = "auto")
red_1
save(red_1, file = "Plots/red_genero")



## Gráfica de prevalencia de taxones a nivel de familia ##

pseq.rel <- microbiome::transform(taxas_juntas, "compositional") # se indica que queremos trabajar con taxas_juntas, el phyloseq que contiene ambas bases (pex y oso)
head(prevalence(pseq.rel, detection = 1/100, sort = TRUE))
head(prevalence(pseq.rel, detection = 1/100, sort = TRUE, count = TRUE))
core.taxa.standard <- core_members(pseq.rel, detection = 0, prevalence = 50/100)
pseq.core <- core(pseq.rel, detection = 0, prevalence = .5)
pseq.core2 <- aggregate_rare(pseq.rel, "Family", detection = 0, prevalence = .5)
core.taxa <- taxa(pseq.core)
det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)
#ggplot(d) + geom_point(aes(x, y)) + scale_x_continuous(trans="log10", limits=c(NA,1))


# Core microbiome
core<- plot_core(pseq.rel, 
          prevalences = prevalences, 
          detections = det, 
          plot.type = "lineplot") + 
  xlab("Relative Abundance (%)")

save(core, file = "Plots/core_microbiome")



## Heatmaps ##
#Familia osos
total = median(sample_sums(oso_corto_F_1))
carbom_abund <- filter_taxa(oso_corto_F_1, function(x) sum(x > total*0.20) > 0, TRUE)
carbom_abund
otu_table(carbom_abund)[1:7, 1:30]
plot_heatmap(carbom_abund, method = "NMDS", distance = "bray")
plot_heatmap(carbom_abund, method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 
             taxa.label = "Family", taxa.order = "Family", 
             trans=NULL, low="beige", high="red", na.value="beige")
#Genero oso
total = median(sample_sums(oso_corto_G_1))
carbom_abund <- filter_taxa(oso_corto_G_1, function(x) sum(x > total*0.20) > 0, TRUE)
carbom_abund
otu_table(carbom_abund)[1:6, 1:30]
plot_heatmap(carbom_abund, method = "NMDS", distance = "bray")
plot_heatmap(carbom_abund, method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 
             taxa.label = "Genus", taxa.order = "Genus", 
             trans=NULL, low="beige", high="red", na.value="beige")

#Familia peces
total = median(sample_sums(pez_corto_F_1))
total
carbom_abund <- filter_taxa(pez_corto_F_1, function(x) sum(x > total*0.20) > 0, TRUE)
carbom_abund
otu_table(carbom_abund)[1:9, 1:29]
plot_heatmap(carbom_abund, method = "NMDS", distance = "bray")
plot_heatmap(carbom_abund, method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 
             taxa.label = "Family", taxa.order = "Family", 
             trans=NULL, low="beige", high="blue", na.value="beige")

#Genero peces
total = median(sample_sums(pez_corto_G_1))
carbom_abund <- filter_taxa(pez_corto_G_1, function(x) sum(x > total*0.20) > 0, TRUE)
carbom_abund
otu_table(carbom_abund)[1:6, 1:30]
plot_heatmap(carbom_abund, method = "NMDS", distance = "bray")
plot_heatmap(carbom_abund, method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 
             taxa.label = "Genus", taxa.order = "Genus", 
             trans=NULL, low="beige", high="blue", na.value="beige")
