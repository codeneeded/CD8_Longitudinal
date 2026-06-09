################################################################################
## CMV-predicted clones vs other expanding clones - DGE
##
## Memory-conscious: TARA is processed end-to-end, freed with gc(), THEN EARTH.
## Each cohort loads its .RData into a temp env and keeps ONLY the *_TRB_<ed>
## object (base object + other ED objects are dropped immediately).
##
## Groups (per cell):
##   CMV             = expanding clone (clonalFrequency > 1) whose TRB epitope
##                     species (Trex annotateDB) matches a CMV pattern
##   Other_Expanding = expanding clone (clonalFrequency > 1) NOT matching CMV
################################################################################


#### 0. Setup ##################################################################

library(Seurat)
library(tidyverse)
library(ggrepel)
library(qs2)

## Parameters
ed_choice    <- "1"                                    # Trex ED object: "0","1","2"
cmv_pattern  <- "CMV|cytomegalo|herpesvirus 5|HHV-?5"  # case-insensitive
species_col  <- "TRB_Epitope.species"
freq_col     <- "clonalFrequency"
assay        <- "RNA"                                  # log-norm 'data' slot used
other_must_be_annotated <- FALSE                       # TRUE = vs other KNOWN specificities only
min_pct         <- 0.10
logfc_threshold <- 0.25

## Paths
rdata_dir <- "/home/akshay-iyer/Documents/CD8_Longitudinal/saved_R_data"                             # inputs + qs2 outputs
fig_dir   <- "/home/akshay-iyer/Documents/CD8_Longitudinal/Analysis/VDJ/TCR/DGE_CMV_vs_Expanding"    # CSV + PNG
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## RUNTIME ASSUMPTIONS:
##  - RNA 'data' slot is log-normalized upstream (NormalizeData already run).
##  - CMV specificity is read from TRB_Epitope.species.


################################################################################
#### 1. TARA ###################################################################
################################################################################

## --- 1a. Load (keep only the ED object, drop the rest) ---
e <- new.env()
load(file.path(rdata_dir, "TARA_TCR_Combined.RData"), envir = e)   # TARA_ALL, TARA_ALL_TRB_0/1/2
tara <- get(paste0("TARA_ALL_TRB_", ed_choice), envir = e)
rm(e); gc()

## --- 1b. Classify cells: CMV vs Other_Expanding ---
sp           <- tara@meta.data[[species_col]]
is_expanding <- !is.na(tara@meta.data[[freq_col]]) & tara@meta.data[[freq_col]] > 1
is_cmv       <- !is.na(sp) & grepl(cmv_pattern, sp, ignore.case = TRUE)

tara$CMV_DGE_Group <- NA_character_
tara$CMV_DGE_Group[is_expanding & is_cmv] <- "CMV"
if (other_must_be_annotated) {
  tara$CMV_DGE_Group[is_expanding & !is_cmv & !is.na(sp) & sp != ""] <- "Other_Expanding"
} else {
  tara$CMV_DGE_Group[is_expanding & !is_cmv] <- "Other_Expanding"
}
table(tara$CMV_DGE_Group, useNA = "ifany")

## --- 1c. DGE (FindMarkers; +avg_log2FC = enriched in CMV clones) ---
DefaultAssay(tara) <- assay
if (inherits(tara[[assay]], "Assay5")) tara[[assay]] <- JoinLayers(tara[[assay]])
Idents(tara) <- "CMV_DGE_Group"
tara_markers <- FindMarkers(tara,
                            ident.1 = "CMV", ident.2 = "Other_Expanding",
                            min.pct = min_pct, logfc.threshold = logfc_threshold)
tara_markers$gene <- rownames(tara_markers)
tara_markers <- tara_markers %>% arrange(p_val_adj, desc(abs(avg_log2FC)))

## --- 1d. Save (qs2 -> saved_R_data ; CSV -> analysis folder) ---
qs_save(tara,         file.path(rdata_dir, paste0("TARA_TRB_ED", ed_choice, "_CMVgroup.qs2")))
qs_save(tara_markers, file.path(rdata_dir, paste0("TARA_CMV_vs_OtherExpanding_ED", ed_choice, "_markers.qs2")))
write_csv(tara_markers, file.path(fig_dir, paste0("TARA_CMV_vs_OtherExpanding_ED", ed_choice, "_DGE.csv")))

## --- 1e. Volcano ---
n       <- table(tara$CMV_DGE_Group)
floor_p <- min(tara_markers$p_val_adj[tara_markers$p_val_adj > 0], na.rm = TRUE)
volc <- tara_markers %>%
  mutate(neg_log10_padj = -log10(ifelse(p_val_adj == 0, floor_p, p_val_adj)),
         direction = case_when(
           p_val_adj < 0.05 & avg_log2FC >  logfc_threshold ~ "Up in CMV",
           p_val_adj < 0.05 & avg_log2FC < -logfc_threshold ~ "Down in CMV",
           TRUE ~ "n.s."))
ggplot(volc, aes(avg_log2FC, neg_log10_padj, color = direction)) +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_text_repel(data = slice_head(filter(volc, direction != "n.s."), n = 20),
                  aes(label = gene), size = 3.5, max.overlaps = 30, show.legend = FALSE) +
  scale_color_manual(values = c("Up in CMV" = "#E41A1C", "Down in CMV" = "#377EB8", "n.s." = "grey70")) +
  labs(title = paste0("TARA: CMV-predicted vs other expanding clones (ED", ed_choice, ")"),
       subtitle = paste0("CMV: ", n["CMV"], " cells  |  Other_Expanding: ", n["Other_Expanding"], " cells"),
       x = "avg log2FC  (>0 = enriched in CMV clones)", y = "-log10(adjusted p-value)", color = NULL) +
  theme_minimal(base_size = 14) + theme(plot.title = element_text(face = "bold"))
ggsave(file.path(fig_dir, paste0("TARA_CMV_vs_OtherExpanding_ED", ed_choice, "_volcano.png")),
       width = 11, height = 9, dpi = 300, bg = "white")

## --- 1f. Free memory before EARTH ---
rm(tara, tara_markers, volc, n, floor_p, sp, is_expanding, is_cmv); gc()


################################################################################
#### 2. EARTH ##################################################################
################################################################################

## --- 2a. Load (keep only the ED object, drop the rest) ---
e <- new.env()
load(file.path(rdata_dir, "EARTH_TCR_Combined.RData"), envir = e)   # EARTH, EARTH_TRB_0/1/2
earth <- get(paste0("EARTH_TRB_", ed_choice), envir = e)
rm(e); gc()

## --- 2b. Classify cells: CMV vs Other_Expanding ---
sp           <- earth@meta.data[[species_col]]
is_expanding <- !is.na(earth@meta.data[[freq_col]]) & earth@meta.data[[freq_col]] > 1
is_cmv       <- !is.na(sp) & grepl(cmv_pattern, sp, ignore.case = TRUE)

earth$CMV_DGE_Group <- NA_character_
earth$CMV_DGE_Group[is_expanding & is_cmv] <- "CMV"
if (other_must_be_annotated) {
  earth$CMV_DGE_Group[is_expanding & !is_cmv & !is.na(sp) & sp != ""] <- "Other_Expanding"
} else {
  earth$CMV_DGE_Group[is_expanding & !is_cmv] <- "Other_Expanding"
}
table(earth$CMV_DGE_Group, useNA = "ifany")

## --- 2c. DGE (FindMarkers; +avg_log2FC = enriched in CMV clones) ---
DefaultAssay(earth) <- assay
if (inherits(earth[[assay]], "Assay5")) earth[[assay]] <- JoinLayers(earth[[assay]])
Idents(earth) <- "CMV_DGE_Group"
earth_markers <- FindMarkers(earth,
                             ident.1 = "CMV", ident.2 = "Other_Expanding",
                             min.pct = min_pct, logfc.threshold = logfc_threshold)
earth_markers$gene <- rownames(earth_markers)
earth_markers <- earth_markers %>% arrange(p_val_adj, desc(abs(avg_log2FC)))

## --- 2d. Save (qs2 -> saved_R_data ; CSV -> analysis folder) ---
qs_save(earth,         file.path(rdata_dir, paste0("EARTH_TRB_ED", ed_choice, "_CMVgroup.qs2")))
qs_save(earth_markers, file.path(rdata_dir, paste0("EARTH_CMV_vs_OtherExpanding_ED", ed_choice, "_markers.qs2")))
write_csv(earth_markers, file.path(fig_dir, paste0("EARTH_CMV_vs_OtherExpanding_ED", ed_choice, "_DGE.csv")))

## --- 2e. Volcano ---
n       <- table(earth$CMV_DGE_Group)
floor_p <- min(earth_markers$p_val_adj[earth_markers$p_val_adj > 0], na.rm = TRUE)
volc <- earth_markers %>%
  mutate(neg_log10_padj = -log10(ifelse(p_val_adj == 0, floor_p, p_val_adj)),
         direction = case_when(
           p_val_adj < 0.05 & avg_log2FC >  logfc_threshold ~ "Up in CMV",
           p_val_adj < 0.05 & avg_log2FC < -logfc_threshold ~ "Down in CMV",
           TRUE ~ "n.s."))
ggplot(volc, aes(avg_log2FC, neg_log10_padj, color = direction)) +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_text_repel(data = slice_head(filter(volc, direction != "n.s."), n = 20),
                  aes(label = gene), size = 3.5, max.overlaps = 30, show.legend = FALSE) +
  scale_color_manual(values = c("Up in CMV" = "#E41A1C", "Down in CMV" = "#377EB8", "n.s." = "grey70")) +
  labs(title = paste0("EARTH: CMV-predicted vs other expanding clones (ED", ed_choice, ")"),
       subtitle = paste0("CMV: ", n["CMV"], " cells  |  Other_Expanding: ", n["Other_Expanding"], " cells"),
       x = "avg log2FC  (>0 = enriched in CMV clones)", y = "-log10(adjusted p-value)", color = NULL) +
  theme_minimal(base_size = 14) + theme(plot.title = element_text(face = "bold"))
ggsave(file.path(fig_dir, paste0("EARTH_CMV_vs_OtherExpanding_ED", ed_choice, "_volcano.png")),
       width = 11, height = 9, dpi = 300, bg = "white")

## --- 2f. Free memory ---
rm(earth, earth_markers, volc, n, floor_p, sp, is_expanding, is_cmv); gc()

