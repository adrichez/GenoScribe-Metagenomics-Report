library(dplyr)
library(tibble)
library(magrittr)

# Tipos de muestra únicos
tipos <- unique(metadata$Type)

# Todas las combinaciones posibles
comparaciones <- combn(tipos, 2, simplify = FALSE)

for (comparacion in comparaciones) {
  grupo1 <- comparacion[1]
  grupo2 <- comparacion[2]
  
  # Filtrar los resultados para esa comparativa
  df_comp <- MetaCyc_daa_annotated_results_df %>%
    filter((group1 == grupo1 & group2 == grupo2) |
             (group1 == grupo2 & group2 == grupo1))
  
  if (nrow(df_comp) == 0) {
    message(paste("❌ No hay datos para:", grupo1, "vs", grupo2))
    next
  }
  
  # Convertir p-values a numéricos
  df_comp$p_values <- as.numeric(df_comp$p_values)
  df_comp$p_adjust <- as.numeric(df_comp$p_adjust)
  
  # Ordenar y extraer el umbral de significancia (30º valor)
  df_comp <- df_comp[order(df_comp$p_adjust), ]
  p_threshold <- if (nrow(df_comp) >= 30) df_comp$p_adjust[30] else max(df_comp$p_adjust, na.rm = TRUE)
  
  # Filtrar metadata solo con muestras de esos grupos
  met <- metadata %>% filter(Type %in% c(grupo1, grupo2))
  muestras_incluidas <- met$SampleID
  
  # Filtrar abundance_data para conservar solo columnas de esas muestras
  abundance_bis <- abundance_data[, c("pathway", muestras_incluidas), drop = FALSE]
  
  # Filtrar solo pathways con abundancia no nula en al menos una muestra
  abundance_bis_filtrada <- abundance_bis[rowSums(abundance_bis[, -1] != 0) > 0, ]
  
  # Crear nombre seguro para los archivos
  nombre_seguro <- paste0(gsub("\\+", "plus", grupo1), "_vs_", gsub("\\+", "plus", grupo2))
  
  # === ERRORBAR ===
  
  p<-pathway_errorbar(
    abundance = abundance_bis %>% column_to_rownames("pathway"),
    daa_results_df = df_comp,
    order = "group",
    Group = met$Type,
    p_values_threshold = p_threshold,
    select = NULL,
    p_value_bar = TRUE,
    colors = NULL,
    x_lab = "description"
  )
  
  png(paste0("pathway_MetagenomeSeq_", nombre_seguro, ".png"), units="in", width=15, height=15, res=300)
  print(p)
  dev.off()
  
  # === PCA ===
  #png(paste0("pathway_pca_MetagenomeSeq_", nombre_seguro, ".png"), units="in", width=15, height=15, res=300)
  #pathway_pca(
  #  abundance = abundance_bis_filtrada %>% column_to_rownames("pathway"),
  #  metadata = met,
  #  group = "Type"
  #)
  #dev.off()
  
  message(paste("✅ Completado:", grupo1, "vs", grupo2))
}

