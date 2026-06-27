setwd("/home/batch_correction")
list.files()


library(qiime2R)
library(phyloseq)
library(microeco)
library(ConQuR)
library(doParallel)
library(ggplot2)

# ==============================================================================
# PASO 1: IMPORTAR DATOS DESDE QIIME 2
# ==============================================================================
# Asegúrate de usar la tabla YA FILTRADA (sin mitocondrias/cloroplastos)
physeq <- qza_to_phyloseq(
  features = "concatenated_method/merge_data/final_table_clean.qza", 
  #tree = "manifest_files/manifest_batch/phylogeny/rooted-tree.qza",
  taxonomy = "concatenated_method/merge_data/taxonomy.qza",
  metadata = "metadata/sample-metadata.tsv"
)

# Verificamos que se cargó bien
print(physeq)

# ==============================================================================
# PASO 2: PREPARAR DATOS PARA CONQUR
# ==============================================================================
# ConQuR necesita una matriz numérica y dataframes separados

# A. Extraer Tabla de OTUs y transponerla (Filas=Muestras, Columnas=Taxones)

otu_mat <- as(otu_table(physeq), "matrix")
if(taxa_are_rows(physeq)){
  otu_mat <- t(otu_mat)
}

# B. Extraer Metadatos
meta_df <- as(sample_data(physeq), "data.frame")

batchid <- factor(meta_df$Lot)

covar <- meta_df[, c( "Sx", "agey","BMI", "glu", "SBP", "chol", "Class..SBP_CTRL")]

Plot_PCoA(TAX = otu_mat, factor = batchid, main = "Before Correction, Bray-Curtis")

#Modelo para Conquer
taxa_corrected1 <- ConQuR(
  tax_tab = otu_mat,
  batchid = batchid,
  covariates = covar,
  batch_ref = "3"   # si "3" es tu lote de referencia
)

taxa_corrected1[1:5, 1:5]

#Graficamos nuestros datos
png("concatenated_method/batch_correction/PCoA_comparacion.png", width = 2000, height = 1200, res = 200)

# Divide en 2 filas y 3 columnas (como querías)
par(mfrow = c(2, 1))

# Haz tus gráficas
Plot_PCoA(TAX = otu_mat, factor = batchid, main = "Before Correction, Bray-Curtis")
Plot_PCoA(TAX = taxa_corrected1, factor = batchid, main = "ConQuR (Default), Bray-Curtis")

dev.off()
