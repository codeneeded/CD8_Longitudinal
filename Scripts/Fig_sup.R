################################################################################
# SUPPLEMENTARY FIGURES — CD8 sub-cluster supporting data
#
# Structure:
#   S1/ — Annotation heatmaps (ADT + RNA)             [supports Fig 4A]
#   S2/ — Per-cluster volcanos + ADT protein dot plot  [supports Fig 4E/F]
#   S3/ — Cluster composition heatmaps                 [supports Fig 4 text]
#   S4/ — Module + gene expression along pseudotime    [supports Fig 5C/D]
#   S5/ — Naïve pairwise volcanos (all 3 comparisons)  [supports Fig 5A text]
#
# Prerequisite: Run Fig4_5_Prep_updated.R
################################################################################

# ── Libraries ─────────────────────────────────────────────────────────────────
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(pheatmap)
library(grid)
library(qs2)
library(ggrepel)
library(ggpubr)
library(SeuratExtend)

# ── Paths ─────────────────────────────────────────────────────────────────────
base_dir     <- "~/Documents/CD8_Longitudinal"
saved_dir    <- file.path(base_dir, "saved_R_data")
manuscript   <- "/home/akshay-iyer/Documents/CD8_Longitudinal/Manuscript/Fig 4-5"
supp_dir     <- file.path(manuscript, "Supplementary")
analysis_dir <- file.path(manuscript, "analysis")

s1_dir <- file.path(supp_dir, "S1_Annotation_Heatmaps")
s2_dir <- file.path(supp_dir, "S2_Effector_Volcanos_and_ADT")
s3_dir <- file.path(supp_dir, "S3_Cluster_Composition")
s4_dir <- file.path(supp_dir, "S4_Module_Pseudotime")
s5_dir <- file.path(supp_dir, "S5_Naive_Volcanos")

for (d in c(s1_dir, s2_dir, s3_dir, s4_dir, s5_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ── Load object ──────────────────────────────────────────────────────────────
TARA_cd8 <- qs_read(file.path(saved_dir, "TARA_cd8_HEI_annotated_final.qs2"))
cat("Loaded:", ncol(TARA_cd8), "cells\n")

# ── Shared aesthetics ────────────────────────────────────────────────────────
art_colors <- c(
  "PreART_Entry"           = "#4A90D9",
  "PostART_Suppressed"     = "#52B788",
  "PostART_Unsuppressed"   = "#E76F51"
)

col_order_cd8 <- c(
  "Naïve CD8", "Naïve CD8 2", "Naïve CD8 3",
  "Transitional Tem CD8", "TEM CD8", "TEMRA CD8",
  "KIR+ innate-like CD8", "MAIT-like Trm",
  "γδ1 T cell", "Naïve γδ1 T cell", "γδ2 T cell"
)

effector_clusters <- c("TEM CD8", "TEMRA CD8", "Transitional Tem CD8")

################################################################################
# ═══════════════ S1: ANNOTATION HEATMAPS ═══════════════════════════════════
################################################################################
cat("\n=== S1: Annotation heatmaps ===\n")

adt_marker_groups <- c(
  "CD45RA"="Naïve", "SELL"="Naïve", "CD7"="Naïve",
  "FAS"="Memory/Tscm", "IL2RB"="Memory/Tscm", "CD45RO"="Memory/Tscm",
  "IL7R"="Memory/Tscm", "CD44"="Memory/Tscm",
  "CD27"="Co-stimulation", "CD28"="Co-stimulation", "ICOS"="Co-stimulation",
  "B3GAT1"="Effector/TEMRA", "KLRG1"="Effector/TEMRA", "KIR3DL1"="Effector/TEMRA",
  "CX3CR1"="Effector/TEMRA",
  "TIGIT"="Exhaustion", "PDCD1"="Exhaustion", "LAG3"="Exhaustion",
  "CD38"="Activation", "CD69"="Activation", "ENTPD1"="Activation", "NT5E"="Activation",
  "CXCR3"="Homing", "ITGB7"="Homing", "ITGA1"="Homing",
  "NCAM1"="NK-like/Innate", "FCGR3A"="NK-like/Innate", "SIGLEC7"="NK-like/Innate",
  "KLRD1"="NK-like/Innate", "KLRB1"="NK-like/Innate",
  "TCR-AB"="TCR Identity", "TCR-vA7.2"="TCR Identity", "TCR-vD2"="TCR Identity"
)

rna_marker_groups <- c(
  "CCR7"="Naïve/Stemness", "SELL"="Naïve/Stemness", "TCF7"="Naïve/Stemness",
  "LEF1"="Naïve/Stemness", "IL7R"="Naïve/Stemness", "KLF2"="Naïve/Stemness",
  "S1PR1"="Naïve/Stemness", "FOXP1"="Naïve/Stemness",
  "BACH2"="Memory/Survival", "BCL2"="Memory/Survival", "CD27"="Memory/Survival",
  "CD28"="Memory/Survival", "ID3"="Memory/Survival",
  "FAS"="Transitional", "IL2RB"="Transitional", "CD44"="Transitional",
  "CXCR3"="Transitional",
  "GZMK"="Effector Memory", "EOMES"="Effector Memory",
  "GZMB"="TEMRA/Cytotoxic", "GNLY"="TEMRA/Cytotoxic", "PRF1"="TEMRA/Cytotoxic",
  "NKG7"="TEMRA/Cytotoxic", "TBX21"="TEMRA/Cytotoxic", "CX3CR1"="TEMRA/Cytotoxic",
  "FGFBP2"="TEMRA/Cytotoxic", "GZMA"="TEMRA/Cytotoxic", "GZMH"="TEMRA/Cytotoxic",
  "GZMM"="TEMRA/Cytotoxic",
  "RUNX3"="Effector TFs", "ZEB2"="Effector TFs", "PRDM1"="Effector TFs", "ID2"="Effector TFs",
  "TOX"="Exhaustion", "PDCD1"="Exhaustion", "TIGIT"="Exhaustion",
  "HAVCR2"="Exhaustion", "LAG3"="Exhaustion",
  "MKI67"="Proliferation", "TOP2A"="Proliferation",
  "ITGAE"="Tissue Residency", "CXCR6"="Tissue Residency", "CD69"="Tissue Residency",
  "IFNG"="Cytokines", "TNF"="Cytokines",
  "TRDV1"="γδ TCR", "TRDV2"="γδ TCR", "TRGV9"="γδ TCR", "TRDC"="γδ TCR",
  "TYROBP"="NK-like/Innate", "KLRD1"="NK-like/Innate", "KLRB1"="NK-like/Innate",
  "FCGR3A"="NK-like/Innate",
  "SLC4A10"="MAIT", "ZBTB16"="MAIT", "NCR3"="MAIT",
  "HLA-DRA"="Activation", "CD38"="Activation", "S1PR5"="Activation"
)

group_colors <- c(
  "Naïve"="#74C2E1", "Naïve/Stemness"="#74C2E1", "Memory/Tscm"="#66BB6A",
  "Memory/Survival"="#66BB6A", "Co-stimulation"="#AED581", "Transitional"="#FFCC80",
  "Effector Memory"="#FFA726", "Effector/TEMRA"="#EF5350", "TEMRA/Cytotoxic"="#EF5350",
  "Effector TFs"="#FF8A65", "Exhaustion"="#A1887F", "Proliferation"="#CE93D8",
  "Tissue Residency"="#80DEEA", "Cytokines"="#FFE082", "Activation"="#FFD54F",
  "Homing"="#B2DFDB", "γδ TCR"="#9C6FD6", "NK-like/Innate"="#F06292",
  "MAIT"="#DCE775", "TCR Identity"="#B0BEC5"
)

adt_present <- names(adt_marker_groups)[names(adt_marker_groups) %in% rownames(TARA_cd8[["ADT"]])]
rna_present <- names(rna_marker_groups)[names(rna_marker_groups) %in% rownames(TARA_cd8[["RNA"]])]

DefaultAssay(TARA_cd8) <- "ADT"
avg_adt <- AverageExpression(TARA_cd8, assays = "ADT", features = adt_present,
                             group.by = "CD8_Annotation", slot = "data")$ADT
DefaultAssay(TARA_cd8) <- "RNA"
avg_rna <- AverageExpression(TARA_cd8, assays = "RNA", features = rna_present,
                             group.by = "CD8_Annotation", slot = "data")$RNA

scale_01 <- function(mat) {
  t(apply(mat, 1, function(x) { r <- max(x)-min(x); if(r==0) rep(0.5,length(x)) else (x-min(x))/r }))
}

avg_adt_s <- scale_01(log10(avg_adt+1)); avg_rna_s <- scale_01(avg_rna)
colnames(avg_adt_s) <- gsub("^g ","",colnames(avg_adt_s))
colnames(avg_rna_s) <- gsub("^g ","",colnames(avg_rna_s))
co <- col_order_cd8[col_order_cd8 %in% colnames(avg_adt_s)]
avg_adt_s <- avg_adt_s[,co]; avg_rna_s <- avg_rna_s[,col_order_cd8[col_order_cd8 %in% colnames(avg_rna_s)]]

make_annot_gaps <- function(mat, marker_groups, grp_colors) {
  annot <- data.frame(Group=factor(marker_groups[rownames(mat)],levels=names(grp_colors)),row.names=rownames(mat))
  groups <- levels(droplevels(annot$Group)); row_order <- character(0)
  for(grp in groups){rows<-rownames(annot)[annot$Group==grp]
  if(length(rows)==1){row_order<-c(row_order,rows)}else{
    hc<-hclust(dist(mat[rows,,drop=FALSE]),method="ward.D2");row_order<-c(row_order,rows[hc$order])}}
  ao<-annot[row_order,,drop=FALSE]; grps<-as.character(ao$Group); gaps<-which(grps[-length(grps)]!=grps[-1])
  list(mat=mat[row_order,],annot=ao,gaps=gaps,colors=list(Group=grp_colors[levels(droplevels(ao$Group))]))
}

ag_adt <- make_annot_gaps(avg_adt_s, adt_marker_groups, group_colors)
ag_rna <- make_annot_gaps(avg_rna_s, rna_marker_groups, group_colors)
hm_cols <- colorRampPalette(c("#F7FCF5","#C7E9C0","#74C476","#31A354","#006D2C"))(100)

p_adt <- pheatmap(ag_adt$mat,cluster_rows=F,cluster_cols=F,scale="none",color=hm_cols,border_color="white",
                  annotation_row=ag_adt$annot,annotation_colors=ag_adt$colors,annotation_names_row=FALSE,
                  gaps_row=ag_adt$gaps,cellwidth=60,cellheight=22,fontsize=16,fontsize_row=14,fontsize_col=14,
                  angle_col=45,main="Protein (ADT)",silent=TRUE)
p_rna <- pheatmap(ag_rna$mat,cluster_rows=F,cluster_cols=F,scale="none",color=hm_cols,border_color="white",
                  annotation_row=ag_rna$annot,annotation_colors=ag_rna$colors,annotation_names_row=FALSE,
                  gaps_row=ag_rna$gaps,cellwidth=60,cellheight=22,fontsize=16,fontsize_row=14,fontsize_col=14,
                  angle_col=45,main="mRNA (RNA)",silent=TRUE)

png(file.path(s1_dir,"S1_combined_ADT_RNA_heatmap.png"),width=34,height=22,units="in",res=300,bg="white")
grid.newpage(); pushViewport(viewport(layout=grid.layout(1,2)))
pushViewport(viewport(layout.pos.col=1)); grid.draw(p_adt$gtable); popViewport()
pushViewport(viewport(layout.pos.col=2)); grid.draw(p_rna$gtable); popViewport()
dev.off()
png(file.path(s1_dir,"S1a_ADT_heatmap.png"),width=17,height=16,units="in",res=300,bg="white")
grid.draw(p_adt$gtable); dev.off()
png(file.path(s1_dir,"S1b_RNA_heatmap.png"),width=17,height=24,units="in",res=300,bg="white")
grid.draw(p_rna$gtable); dev.off()
cat("  S1 done.\n")

################################################################################
# ═══════════════ S2: EFFECTOR VOLCANOS + ADT DOT PLOT ══════════════════════
################################################################################
cat("\n=== S2: Per-cluster volcanos + ADT protein ===\n")

# ── Shared volcano gene categories ──────────────────────────────────────────
exhaustion_genes   <- c("TOX","PDCD1","HAVCR2","LAG3","TIGIT","CTLA4","ENTPD1")
stemness_genes     <- c("TCF7","SELL","CCR7","LEF1","BACH2","BCL2","IL7R","S1PR1","KLF2")
cytotoxic_genes    <- c("GZMB","GNLY","PRF1","NKG7","GZMA","GZMH","FGFBP2")
stress_genes       <- c("HSPA1A","HSPA1B","DNAJB1","IFI27","IFI44L")
ifn_memory_genes   <- c("IFIT1","IFIT3","ISG15","MX1")
chemokine_genes    <- c("CCL3","CCL3L3","CCL4","CCL4L2","CCL5")
terminal_genes     <- c("ZEB2","PRDM1","TBX21","CX3CR1","S1PR5","ID2","EOMES")
activation_genes   <- c("HLA-DRA","CD38","FAS","CD69")

highlight_cols <- c(
  "Exhaustion"="#A32D2D","Stemness"="#185FA5","Cytotoxicity"="#3B6D11",
  "Stress response"="#BA7517","Type I IFN memory"="#0F6E56","Chemokines"="#993556",
  "Terminal diff."="#534AB7","Activation"="#D85A30","Other"="grey80"
)

genes_to_label_all <- c(exhaustion_genes,stemness_genes,cytotoxic_genes,stress_genes,
                        ifn_memory_genes,chemokine_genes,terminal_genes,activation_genes)

# ── S2a–c: Per-cluster volcanos ──────────────────────────────────────────────
cluster_volcanos <- list(
  list(label="S2a",title="Transitional Tem CD8",
       candidates=c("DGE_MAST_Expanding_TransitionalTem_Sup_vs_Unsup.csv")),
  list(label="S2b",title="TEM CD8",
       candidates=c("DGE_MAST_Expanding_TEMCD8_Sup_vs_Unsup.csv",
                    "DGE_MAST_Expanding_EffectorCD8_Sup_vs_Unsup.csv")),
  list(label="S2c",title="TEMRA CD8",
       candidates=c("DGE_MAST_Expanding_TEMRACD8_Sup_vs_Unsup.csv",
                    "DGE_MAST_Expanding_TerminalTEMRA_Sup_vs_Unsup.csv"))
)

for (cv in cluster_volcanos) {
  dge_path <- NULL
  for (cand in cv$candidates) {
    fp <- file.path(analysis_dir,"03_DGE",cand)
    if (file.exists(fp)) { dge_path <- fp; break }
  }
  if (is.null(dge_path)) { cat("    SKIP",cv$label,"\n"); next }
  cat("    ",cv$label,":",cv$title,"...\n")
  
  dge <- read.csv(dge_path,row.names=1)
  dge$gene <- rownames(dge); dge$neg_log10_padj <- -log10(dge$p_val_adj+1e-300)
  dge$highlight <- "Other"
  dge$highlight[dge$gene %in% exhaustion_genes]  <- "Exhaustion"
  dge$highlight[dge$gene %in% stemness_genes]    <- "Stemness"
  dge$highlight[dge$gene %in% cytotoxic_genes]   <- "Cytotoxicity"
  dge$highlight[dge$gene %in% stress_genes]      <- "Stress response"
  dge$highlight[dge$gene %in% ifn_memory_genes]  <- "Type I IFN memory"
  dge$highlight[dge$gene %in% chemokine_genes]   <- "Chemokines"
  dge$highlight[dge$gene %in% terminal_genes]    <- "Terminal diff."
  dge$highlight[dge$gene %in% activation_genes]  <- "Activation"
  dge$highlight <- factor(dge$highlight, levels=names(highlight_cols))
  
  gp <- genes_to_label_all[genes_to_label_all %in% dge$gene]
  dge$label <- ""; dge$label[dge$gene %in% gp & dge$p_val_adj<0.05] <- dge$gene[dge$gene %in% gp & dge$p_val_adj<0.05]
  
  x_lim <- quantile(abs(dge$avg_log2FC[is.finite(dge$avg_log2FC)]),0.995)*1.15
  y_cap <- quantile(dge$neg_log10_padj[is.finite(dge$neg_log10_padj)],0.995)*1.10
  dge$x_plot <- pmax(pmin(dge$avg_log2FC,x_lim*0.98),-x_lim*0.98)
  dge$y_plot <- pmin(dge$neg_log10_padj,y_cap)
  
  p <- ggplot(dge,aes(x=x_plot,y=y_plot,color=highlight)) +
    geom_point(data=dge%>%filter(highlight=="Other"),size=2,alpha=0.3) +
    geom_point(data=dge%>%filter(highlight!="Other"),size=7,alpha=0.85) +
    geom_label_repel(data=dge%>%filter(label!=""),aes(label=label,fill=highlight),color="white",
                     size=8,fontface="bold.italic",max.overlaps=50,segment.size=0.4,segment.color="grey40",
                     min.segment.length=0.1,box.padding=0.5,label.size=0.2,show.legend=FALSE) +
    geom_hline(yintercept=-log10(0.05),linetype="dashed",color="grey40",linewidth=0.4) +
    geom_vline(xintercept=c(-0.25,0.25),linetype="dashed",color="grey40",linewidth=0.4) +
    scale_color_manual(values=highlight_cols,name="Gene category",drop=FALSE) +
    scale_fill_manual(values=highlight_cols,guide="none") +
    scale_y_continuous(expand=expansion(mult=c(0.02,0.02))) +
    coord_cartesian(xlim=c(-x_lim,x_lim),ylim=c(0,y_cap),clip="off") +
    annotate("segment",x=0.15,xend=x_lim*0.85,y=-y_cap*0.12,yend=-y_cap*0.12,
             arrow=arrow(length=unit(0.3,"cm"),type="closed"),color="#52B788",linewidth=1.3) +
    annotate("text",x=x_lim*0.52,y=-y_cap*0.17,label="Higher in suppressed",
             color="#52B788",size=6.5,fontface="bold") +
    annotate("segment",x=-0.15,xend=-x_lim*0.85,y=-y_cap*0.12,yend=-y_cap*0.12,
             arrow=arrow(length=unit(0.3,"cm"),type="closed"),color="#E76F51",linewidth=1.3) +
    annotate("text",x=-x_lim*0.52,y=-y_cap*0.17,label="Higher in unsuppressed",
             color="#E76F51",size=6.5,fontface="bold") +
    labs(x=expression(log[2]~fold~change),y=expression(-log[10]~adjusted~italic(p)),
         title=paste0(cv$title," — expanding clones: suppressed vs unsuppressed")) +
    theme_cowplot(font_size=20) +
    theme(legend.position="right",legend.text=element_text(size=22),
          legend.title=element_text(size=23,face="bold"),legend.key.size=unit(1.5,"cm"),
          axis.text=element_text(size=18),axis.title=element_text(size=20),
          plot.title=element_text(size=18,face="bold"),
          plot.background=element_rect(fill="white",color=NA),
          panel.background=element_rect(fill="white",color=NA),
          plot.margin=margin(10,10,75,10)) +
    guides(color=guide_legend(override.aes=list(size=9,alpha=1)))
  
  ggsave(file.path(s2_dir,paste0(cv$label,"_Volcano_",gsub(" ","_",cv$title),".png")),
         plot=p,width=14,height=11,dpi=300,bg="white")
}

# ── S2d: ADT protein dot plot — key surface markers by ART status ────────────
cat("    S2d: ADT protein dot plot...\n")

# Select proteins that showed significant DPE and are biologically meaningful
adt_proteins_to_plot <- c(
  # Exhaustion / inhibitory
  "TIGIT", "PDCD1", "LAG3", "HAVCR2",
  # Activation / costimulation
  "FAS", "CD38", "ICOS", "CD27", "CD28",
  # Effector / differentiation
  "B3GAT1", "KLRG1", "CX3CR1", "CD45RA", "CD45RO",
  # Homing / adhesion
  "SELL", "ITGB1", "CD44",
  # NK-like
  "KIR3DL1", "SIGLEC7", "KLRD1"
)

# Filter to proteins present in the object
adt_proteins_present <- adt_proteins_to_plot[adt_proteins_to_plot %in% rownames(TARA_cd8[["ADT"]])]

# Subset to expanding effector clones only (matches what Cohen's d was computed on)
expand_eff <- subset(TARA_cd8,
                     CD8_Annotation %in% effector_clusters &
                       !is.na(cloneSize) & cloneSize != "Single (0 < X <= 1)")

DefaultAssay(expand_eff) <- "ADT"

# Compute mean expression and % expressing per protein × condition
adt_plot_data <- list()
for (prot in adt_proteins_present) {
  expr <- FetchData(expand_eff, vars = prot)[, 1]
  df <- data.frame(
    protein = prot,
    Timepoint_Group = expand_eff$Timepoint_Group,
    expression = expr
  )
  summary_df <- df %>%
    group_by(protein, Timepoint_Group) %>%
    summarise(
      mean_expr = mean(expression, na.rm = TRUE),
      pct_pos   = mean(expression > 0, na.rm = TRUE) * 100,
      .groups = "drop"
    )
  adt_plot_data[[prot]] <- summary_df
}
adt_df <- bind_rows(adt_plot_data)

adt_df$Timepoint_Group <- factor(adt_df$Timepoint_Group,
                                 levels = c("PreART_Entry", "PostART_Suppressed", "PostART_Unsuppressed"))
adt_df$protein <- factor(adt_df$protein, levels = rev(adt_proteins_present))

# Clean labels
adt_df$Condition <- dplyr::recode(adt_df$Timepoint_Group,
                                  "PreART_Entry" = "Pre-ART", "PostART_Suppressed" = "Suppressed",
                                  "PostART_Unsuppressed" = "Unsuppressed")
adt_df$Condition <- factor(adt_df$Condition, levels = c("Pre-ART", "Suppressed", "Unsuppressed"))

p_adt_dot <- ggplot(adt_df, aes(x = Condition, y = protein)) +
  geom_point(aes(size = pct_pos, color = mean_expr)) +
  scale_size_continuous(name = "% expressing", range = c(2, 12), limits = c(0, 100)) +
  scale_color_viridis_c(option = "magma", name = "Mean expression\n(DSB-normalized)",
                        direction = -1) +
  labs(x = NULL, y = NULL,
       title = "Surface protein expression — expanding effector clones") +
  theme_cowplot(font_size = 18) +
  theme(
    axis.text.x      = element_text(size = 16, angle = 30, hjust = 1),
    axis.text.y      = element_text(size = 14),
    plot.title       = element_text(size = 18, face = "bold"),
    legend.position  = "right",
    legend.text      = element_text(size = 12),
    legend.title     = element_text(size = 13),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3)
  )

ggsave(file.path(s2_dir, "S2d_ADT_protein_dotplot.png"),
       plot = p_adt_dot, width = 10, height = 12, dpi = 300, bg = "white")

cat("  S2 done.\n")

################################################################################
# ═══════════════ S3: CLUSTER COMPOSITION HEATMAPS ══════════════════════════
################################################################################
cat("\n=== S3: Cluster composition ===\n")

# Row-normalized: what % of each cluster comes from each condition
comp_mat <- as.data.frame.matrix(
  prop.table(table(TARA_cd8$CD8_Annotation, TARA_cd8$Timepoint_Group), margin = 1)
)
comp_mat <- comp_mat[col_order_cd8[col_order_cd8 %in% rownames(comp_mat)], , drop = FALSE]
colnames(comp_mat) <- c("Pre-ART", "Suppressed", "Unsuppressed")[
  match(colnames(comp_mat), c("PreART_Entry", "PostART_Suppressed", "PostART_Unsuppressed"))]

count_mat <- as.data.frame.matrix(table(TARA_cd8$CD8_Annotation, TARA_cd8$Timepoint_Group))
count_mat <- count_mat[rownames(comp_mat), , drop = FALSE]

display_mat <- matrix(paste0(round(as.matrix(comp_mat)*100,1),"%"),
                      nrow=nrow(comp_mat), dimnames=dimnames(comp_mat))

prop_colors <- colorRampPalette(c("#FFFFFF","#D6E8F7","#7FB3D8","#2171B5","#08306B"))(100)

col_annot <- data.frame(Condition=colnames(comp_mat), row.names=colnames(comp_mat))
col_annot_colors <- list(Condition=c("Pre-ART"="#4A90D9","Suppressed"="#52B788","Unsuppressed"="#E76F51"))

# FIXED: wider canvas + larger margins for labels
png(file.path(s3_dir, "S3a_Cluster_composition_by_ART.png"),
    width = 12, height = 14, units = "in", res = 300, bg = "white")
pheatmap(
  as.matrix(comp_mat),
  cluster_rows = FALSE, cluster_cols = FALSE,
  color = prop_colors, border_color = "white",
  cellwidth = 100, cellheight = 40,
  fontsize = 16, fontsize_row = 14, fontsize_col = 16,
  angle_col = 0,
  display_numbers = display_mat,
  fontsize_number = 13,
  number_color = ifelse(as.matrix(comp_mat) > 0.45, "white", "black"),
  annotation_col = col_annot,
  annotation_colors = col_annot_colors,
  annotation_names_col = FALSE,
  main = "Cluster composition by ART status\n(% of each cluster from each condition)"
)
dev.off()

# Column-normalized: where do cells from each condition distribute?
comp_col <- as.data.frame.matrix(
  prop.table(table(TARA_cd8$CD8_Annotation, TARA_cd8$Timepoint_Group), margin = 2)
)
comp_col <- comp_col[col_order_cd8[col_order_cd8 %in% rownames(comp_col)], , drop = FALSE]
colnames(comp_col) <- c("Pre-ART", "Suppressed", "Unsuppressed")[
  match(colnames(comp_col), c("PreART_Entry", "PostART_Suppressed", "PostART_Unsuppressed"))]

display_col <- matrix(paste0(round(as.matrix(comp_col)*100,1),"%"),
                      nrow=nrow(comp_col), dimnames=dimnames(comp_col))

png(file.path(s3_dir, "S3b_Condition_distribution_across_clusters.png"),
    width = 12, height = 14, units = "in", res = 300, bg = "white")
pheatmap(
  as.matrix(comp_col),
  cluster_rows = FALSE, cluster_cols = FALSE,
  color = prop_colors, border_color = "white",
  cellwidth = 100, cellheight = 40,
  fontsize = 16, fontsize_row = 14, fontsize_col = 16,
  angle_col = 0,
  display_numbers = display_col,
  fontsize_number = 13,
  number_color = ifelse(as.matrix(comp_col) > 0.15, "white", "black"),
  annotation_col = col_annot,
  annotation_colors = col_annot_colors,
  annotation_names_col = FALSE,
  main = "Distribution of each condition across clusters\n(% of each condition in each cluster)"
)
dev.off()

write.csv(data.frame(Cluster=rownames(comp_mat),comp_mat,Total=rowSums(count_mat)),
          file.path(s3_dir,"S3_composition_proportions.csv"),row.names=FALSE)
write.csv(data.frame(Cluster=rownames(count_mat),count_mat,Total=rowSums(count_mat)),
          file.path(s3_dir,"S3_composition_counts.csv"),row.names=FALSE)
cat("  S3 done.\n")

################################################################################
# ═══════════════ S4: MODULE + GENE PSEUDOTIME ══════════════════════════════
################################################################################
cat("\n=== S4: Module & gene pseudotime panels ===\n")

m3_dir <- file.path(analysis_dir, "07_trajectory", "monocle3")

# Copy individual module plots (not the combined AllModules/ExhaustionStemness in main figs)
m3_files <- list.files(m3_dir, pattern="^Monocle3_[A-Z].*_along_pseudotime\\.png$", full.names=TRUE)
m3_files <- m3_files[!grepl("AllModules|ExhaustionStemness", m3_files)]
for (f in m3_files) file.copy(f, file.path(s4_dir, paste0("S4_",basename(f))), overwrite=TRUE)

# Copy gene expression panels + trajectory UMAPs
m3_gene <- list.files(m3_dir, pattern="Genes_in_pseudotime|NarrativeGenes|TopVariable|TopGenes_UMAP|UMAP_clusters|UMAP_ART",
                      full.names=TRUE)
for (f in m3_gene) file.copy(f, file.path(s4_dir, paste0("S4_",basename(f))), overwrite=TRUE)

cat("  S4: Copied", length(m3_files)+length(m3_gene), "panels.\n")

################################################################################
# ═══════════════ S5: NAÏVE PAIRWISE VOLCANOS ═══════════════════════════════
################################################################################
cat("\n=== S5: Naïve pairwise volcanos ===\n")

quiescence_genes  <- c("KLF2","S1PR1","SELL","BACH2","BCL2","IL7R","FOXP1","CCR7","LEF1","TCF7")
hypoxia_genes     <- c("IGFBP2","HILPDA","PFKFB4","EGLN3","GBP5","AK4","HIF1A","MIR210HG","DARS","FGF11")
ifn_genes_naive   <- c("MX1","IFIT1","ISG15","IFI44L","STAT1","IRF7","IFIT3")
survival_genes    <- c("BACH2","BCL2","FOXP1","SATB1","ID3")
activation_naive  <- c("HLA-DRA","CD38","ICOS","CD69","FAS")

naive_cols <- c(
  "Quiescence / homing"="#185FA5","Hypoxia / metabolism"="#A32D2D",
  "Type I IFN"="#0F6E56",
  "Survival / quiescence"="#7B68EE","Activation"="#D85A30","Other"="grey80"
)

naive_label_genes <- unique(c(quiescence_genes,hypoxia_genes,
                              ifn_genes_naive,survival_genes,activation_naive))

naive_comps <- list(
  list(label="S5a",title="Naïve CD8 vs Naïve CD8 2",
       id1="Naïve CD8",id2="Naïve CD8 2",
       files=c("DGE_MAST_RNA_NaveCD8_vs_NaveCD82.csv","DGE_MAST_RNA_NaveCD8_vs_NaveCD8innatelike.csv")),
  list(label="S5b",title="Naïve CD8 vs Naïve CD8 3",
       id1="Naïve CD8",id2="Naïve CD8 3",
       files=c("DGE_MAST_RNA_NaveCD8_vs_NaveCD83.csv","DGE_MAST_RNA_NaveCD8_vs_NaveCD8postART.csv")),
  list(label="S5c",title="Naïve CD8 2 vs Naïve CD8 3",
       id1="Naïve CD8 2",id2="Naïve CD8 3",
       files=c("DGE_MAST_RNA_NaveCD82_vs_NaveCD83.csv","DGE_MAST_RNA_NaveCD8innatelike_vs_NaveCD8postART.csv"))
)

for (nc in naive_comps) {
  dge_path <- NULL
  for (f in nc$files) { fp <- file.path(analysis_dir,"03_DGE",f); if(file.exists(fp)){dge_path<-fp;break} }
  if (is.null(dge_path)) { cat("    SKIP",nc$label,"\n"); next }
  cat("    ",nc$label,":",nc$title,"...\n")
  
  dge <- read.csv(dge_path,row.names=1)
  dge$gene <- rownames(dge); dge$neg_log10_padj <- -log10(dge$p_val_adj+1e-300)
  
  dge$highlight <- "Other"
  dge$highlight[dge$gene %in% quiescence_genes] <- "Quiescence / homing"
  dge$highlight[dge$gene %in% hypoxia_genes]    <- "Hypoxia / metabolism"
  dge$highlight[dge$gene %in% ifn_genes_naive]  <- "Type I IFN"
  dge$highlight[dge$gene %in% survival_genes]   <- "Survival / quiescence"
  dge$highlight[dge$gene %in% activation_naive] <- "Activation"
  dge$highlight <- factor(dge$highlight, levels=names(naive_cols))
  
  gp <- naive_label_genes[naive_label_genes %in% dge$gene]
  dge$label <- ""; dge$label[dge$gene%in%gp & dge$p_val_adj<0.05] <- dge$gene[dge$gene%in%gp & dge$p_val_adj<0.05]
  dge$label_ns <- ""; dge$label_ns[dge$gene%in%gp & dge$p_val_adj>=0.05] <- dge$gene[dge$gene%in%gp & dge$p_val_adj>=0.05]
  
  x_lim <- quantile(abs(dge$avg_log2FC[is.finite(dge$avg_log2FC)]),0.995)*1.15
  y_cap <- quantile(dge$neg_log10_padj[is.finite(dge$neg_log10_padj)],0.995)*1.10
  dge$x_plot <- pmax(pmin(dge$avg_log2FC,x_lim*0.98),-x_lim*0.98)
  dge$y_plot <- pmin(dge$neg_log10_padj,y_cap)
  
  p <- ggplot(dge,aes(x=x_plot,y=y_plot,color=highlight)) +
    geom_point(data=dge%>%filter(highlight=="Other"),size=2,alpha=0.3) +
    geom_point(data=dge%>%filter(highlight!="Other"),size=7,alpha=0.85) +
    geom_label_repel(data=dge%>%filter(label!=""),aes(label=label,fill=highlight),color="white",
                     size=8,fontface="bold.italic",max.overlaps=50,segment.size=0.4,segment.color="grey40",
                     min.segment.length=0.1,box.padding=0.5,label.size=0.2,show.legend=FALSE) +
    ggrepel::geom_text_repel(data=dge%>%filter(label_ns!=""),aes(label=label_ns),size=5,
                             fontface="italic",color="grey40",max.overlaps=20,segment.size=0.2,segment.color="grey60",
                             segment.linetype="dashed",min.segment.length=0.1,box.padding=0.3,show.legend=FALSE) +
    geom_hline(yintercept=-log10(0.05),linetype="dashed",color="grey40",linewidth=0.4) +
    geom_vline(xintercept=c(-0.25,0.25),linetype="dashed",color="grey40",linewidth=0.4) +
    scale_color_manual(values=naive_cols,name="Gene category",drop=FALSE) +
    scale_fill_manual(values=naive_cols,guide="none") +
    scale_y_continuous(expand=expansion(mult=c(0.02,0.02))) +
    coord_cartesian(xlim=c(-x_lim,x_lim),ylim=c(0,y_cap),clip="off") +
    annotate("segment",x=0.15,xend=x_lim*0.85,y=-y_cap*0.12,yend=-y_cap*0.12,
             arrow=arrow(length=unit(0.3,"cm"),type="closed"),color="#185FA5",linewidth=1.3) +
    annotate("text",x=x_lim*0.52,y=-y_cap*0.17,label=paste0("Higher in ",nc$id1),
             color="#185FA5",size=5.5,fontface="bold") +
    annotate("segment",x=-0.15,xend=-x_lim*0.85,y=-y_cap*0.12,yend=-y_cap*0.12,
             arrow=arrow(length=unit(0.3,"cm"),type="closed"),color="#A32D2D",linewidth=1.3) +
    annotate("text",x=-x_lim*0.52,y=-y_cap*0.17,label=paste0("Higher in ",nc$id2),
             color="#A32D2D",size=5.5,fontface="bold") +
    labs(x=expression(log[2]~fold~change),y=expression(-log[10]~adjusted~italic(p)),title=nc$title) +
    theme_cowplot(font_size=20) +
    theme(legend.position="right",legend.text=element_text(size=22),
          legend.title=element_text(size=23,face="bold"),legend.key.size=unit(1.5,"cm"),
          axis.text=element_text(size=18),axis.title=element_text(size=20),
          plot.title=element_text(size=20,face="bold"),
          plot.background=element_rect(fill="white",color=NA),
          panel.background=element_rect(fill="white",color=NA),
          plot.margin=margin(10,10,75,10)) +
    guides(color=guide_legend(override.aes=list(size=9,alpha=1)))
  
  ggsave(file.path(s5_dir,paste0(nc$label,"_Volcano_",gsub(" ","_",nc$title),".png")),
         plot=p,width=14,height=11,dpi=300,bg="white")
}
cat("  S5 done.\n")

################################################################################
# SUMMARY + INTERPRETATION EXPORTS
################################################################################
cat("\n", paste(rep("=",60),collapse=""), "\n")
cat("SUPPLEMENTARY COMPLETE\n")
for (d in list(c("S1",s1_dir),c("S2",s2_dir),c("S3",s3_dir),c("S4",s4_dir),c("S5",s5_dir))) {
  cat(" ",d[1],":",length(list.files(d[2],"\\.png$|\\.csv$")),"files\n")
}

# ── Interpretation summaries ─────────────────────────────────────────────────
summary_dir <- file.path(supp_dir, "interpretation_summaries")
dir.create(summary_dir, recursive=TRUE, showWarnings=FALSE)

dge_dir <- file.path(analysis_dir,"03_DGE"); dpe_dir <- file.path(analysis_dir,"04_DPE")

# Top DGE per cluster
for (cl in list(
  list(n="TEM CD8",f=c("DGE_MAST_Expanding_TEMCD8_Sup_vs_Unsup.csv","DGE_MAST_Expanding_EffectorCD8_Sup_vs_Unsup.csv")),
  list(n="TEMRA CD8",f=c("DGE_MAST_Expanding_TEMRACD8_Sup_vs_Unsup.csv","DGE_MAST_Expanding_TerminalTEMRA_Sup_vs_Unsup.csv")),
  list(n="Transitional Tem",f=c("DGE_MAST_Expanding_TransitionalTem_Sup_vs_Unsup.csv")),
  list(n="All",f=c("DGE_MAST_Expanding_AllClusters_Sup_vs_Unsup.csv"))
)) { for(fi in cl$f){fp<-file.path(dge_dir,fi);if(file.exists(fp)){
  file.copy(fp,file.path(summary_dir,paste0("FullDGE_",gsub(" ","_",cl$n),".csv")),overwrite=TRUE);break}}}

# ADT DPE
for (cl in list(
  list(n="TEM CD8",f=c("DPE_MAST_ADT_Expanding_TEMCD8_Sup_vs_Unsup.csv","DPE_MAST_ADT_Expanding_EffectorCD8_Sup_vs_Unsup.csv")),
  list(n="TEMRA CD8",f=c("DPE_MAST_ADT_Expanding_TEMRACD8_Sup_vs_Unsup.csv","DPE_MAST_ADT_Expanding_TerminalTEMRA_Sup_vs_Unsup.csv")),
  list(n="All",f=c("DPE_MAST_ADT_Expanding_AllClusters_Sup_vs_Unsup.csv"))
)) { for(fi in cl$f){fp<-file.path(dpe_dir,fi);if(file.exists(fp)){
  file.copy(fp,file.path(summary_dir,paste0("FullDPE_ADT_",gsub(" ","_",cl$n),".csv")),overwrite=TRUE);break}}}

# Module scores + Cohen's d + VL
for (src in c("10_module_scores/ModuleScores_per_cell.csv",
              "09_viral_load/ViralLoad_Spearman_correlations.csv",
              "09_viral_load/PerSample_ModuleScores_and_ViralLoad.csv")) {
  fp <- file.path(analysis_dir, src)
  if (file.exists(fp)) file.copy(fp, file.path(summary_dir, basename(fp)), overwrite=TRUE)
}
cd_fp <- file.path(manuscript,"Figure4_panels","Fig4F_CohenD_stats.csv")
if (file.exists(cd_fp)) file.copy(cd_fp, file.path(summary_dir,"CohenD_all_comparisons.csv"), overwrite=TRUE)

cat("\n  Interpretation summaries saved to:", summary_dir, "\n")

