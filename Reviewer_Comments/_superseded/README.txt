Superseded outputs. Kept for traceability, NOT for use.

RevFig_4B_cluster_crosswalk.png (+ source data, + CD8_cluster_crosswalk.csv)
  WRONG AXIS LABEL. It labelled the x-axis "Cluster in the whole-PBMC
  clustering (Fig 1)" and its caption claimed the off-diagonal structure showed
  substructure the whole-PBMC clustering had merged.

  In fact $source_cluster is the DEDICATED CD8 WNN clustering at resolution 0.4
  (4.CD8_Subseting.R line 602, set from $seurat_clusters AFTER
  FindMultiModalNeighbors). Verified: purity of CD8_Annotation against
  wsnn_res.0.4 is 0.859, identical to source_cluster.

  The figure's real content is the ADT sub-gating step redistributing the four
  naive WNN clusters into Naive Intermediate (CD45RO+) and Tscm (CD45RO-/FAS+).
  Replaced by RevFig_4B_ADT_subgating_crosswalk.png.

Do not cite anything in this directory in the response letter.
