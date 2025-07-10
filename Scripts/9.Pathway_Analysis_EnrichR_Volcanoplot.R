library(dplyr)
library(enrichR)
library(openxlsx)
library(ggplot2)

# -----------------------------
# Enrichr databases
# -----------------------------
databases <- c(
  "TRRUST_Transcription_Factors_2019", "ChEA_2022", "TRANSFAC_and_JASPAR_PWMs",
  "KEGG_2021_Human", "WikiPathways_2024_Human", "GO_Biological_Process_2023",
  "MSigDB_Hallmark_2020", "Panther_2016", "Reactome_2022", "BioPlanet_2019"
)

tf_databases <- c("TRRUST_Transcription_Factors_2019", "ChEA_2022", "TRANSFAC_and_JASPAR_PWMs")
pathway_databases <- setdiff(databases, tf_databases)

# -----------------------------
# Input directories for DGE CSVs
# -----------------------------
input_dirs <- list(
  HEIvsHEU = "/home/akshay-iyer/Documents/CD8_Longitudinal/Differential_Expression/TARA_ALL/HEIvsHEU",
  HEUvsHUU = "/home/akshay-iyer/Documents/CD8_Longitudinal/Differential_Expression/TARA_ALL/HEUvsHUU",
  PostARTvsPreART = "/home/akshay-iyer/Documents/CD8_Longitudinal/Differential_Expression/TARA_ALL/PostARTvsPreART",
  HighvsLowVL = "/home/akshay-iyer/Documents/CD8_Longitudinal/Differential_Expression/TARA_Entry/HighvsLowVL"
)

# -----------------------------
# Output base directory
# -----------------------------
base_output <- "/home/akshay-iyer/Documents/CD8_Longitudinal/Pathway_Analysis_EnrichR"

# -----------------------------
# Run Enrichr on significant genes (p_adj < 0.05)
# -----------------------------
run_enrichment <- function(gene_df, label, base_output) {
  gene_list <- rownames(gene_df)
  if (length(gene_list) == 0) {
    message("No genes found for ", label)
    return(NULL)
  }
  
  enrichment <- enrichr(gene_list, databases)
  
  # Create output subfolders
  dir.create(file.path(base_output, "CSVs"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(base_output, "Plots"), recursive = TRUE, showWarnings = FALSE)
  
  # Save Excel results
  excel_out <- file.path(base_output, "CSVs", paste0(label, "_Enrichment.xlsx"))
  wb <- createWorkbook()
  for (db in names(enrichment)) {
    addWorksheet(wb, substr(db, 1, 31))
    writeData(wb, substr(db, 1, 31), enrichment[[db]])
  }
  saveWorkbook(wb, excel_out, overwrite = TRUE)
  
  # Process top TF results
  top_tf_list <- list()
  top_pathway_list <- list()
  
  for (db_name in names(enrichment)) {
    db_results <- enrichment[[db_name]]
    
    # Rename Combined Score if needed
    if ("Combined Score" %in% colnames(db_results)) {
      db_results <- db_results %>% rename(Combined.Score = `Combined Score`)
    }
    
    if (!"Combined.Score" %in% colnames(db_results)) {
      message("Skipping ", db_name, " for ", label, " — no Combined.Score.")
      next
    }
    
    sig_results <- db_results %>% filter(Adjusted.P.value < 0.05)
    
    if (nrow(sig_results) > 0) {
      top_terms <- sig_results %>%
        arrange(desc(Combined.Score)) %>%
        slice_head(n = 10) %>%
        mutate(Database = db_name)
      
      if (db_name %in% tf_databases) {
        top_tf_list[[db_name]] <- top_terms
      } else if (db_name %in% pathway_databases) {
        top_pathway_list[[db_name]] <- top_terms
      }
    }
  }
  
  # Plot TFs
  tf_df <- bind_rows(top_tf_list)
  if ("Combined.Score" %in% colnames(tf_df) && nrow(tf_df) > 0) {
    tf_df <- tf_df %>% arrange(desc(Combined.Score)) %>% slice_head(n = 20)
    
    p_tf <- ggplot(tf_df, aes(x = reorder(Term, Combined.Score), y = Combined.Score, fill = Database)) +
      geom_bar(stat = "identity") +
      scale_y_log10() +  # <-- log scale here
      coord_flip() +
      labs(
        title = paste("Top Transcription Factors -", label),
        x = "TF Term", y = "log10(Combined Score)"
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5)
      )
    
    ggsave(
      filename = file.path(base_output, "Plots", paste0(label, "_Transcription_Factors.png")),
      plot = p_tf, width = 12, height = 10, dpi = 300, bg = "white"
    )
  }
  
  # Plot Pathways
  pathway_df <- bind_rows(top_pathway_list)
  if ("Combined.Score" %in% colnames(pathway_df) && nrow(pathway_df) > 0) {
    pathway_df <- pathway_df %>% arrange(desc(Combined.Score)) %>% slice_head(n = 20)
    
    p_path <- ggplot(pathway_df, aes(x = reorder(Term, Combined.Score), y = Combined.Score, fill = Database)) +
      geom_bar(stat = "identity") +
      scale_y_log10() +  # <-- log scale here
      coord_flip() +
      labs(
        title = paste("Top Pathways -", label),
        x = "Pathway Term", y = "log10(Combined Score)"
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5)
      )
    
    ggsave(
      filename = file.path(base_output, "Plots", paste0(label, "_Pathways.png")),
      plot = p_path, width = 12, height = 10, dpi = 300, bg = "white"
    )
  }
}

# -----------------------------
# Loop over all comparisons and clusters
# -----------------------------
for (comparison in names(input_dirs)) {
  input_dir <- input_dirs[[comparison]]
  comparison_output <- file.path(base_output, comparison)
  
  csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
  
  for (csv in csv_files) {
    dge <- read.csv(csv, row.names = 1)
    sig_genes <- dge %>% filter(p_val_adj < 0.05)
    
    if (nrow(sig_genes) == 0) {
      message("No significant genes in ", basename(csv))
      next
    }
    
    label <- tools::file_path_sans_ext(basename(csv))  # e.g., 0__CD4_T_cell_HEIvsHEU
    message("Running Enrichr for: ", label)
    run_enrichment(sig_genes, label, comparison_output)
  }
}
