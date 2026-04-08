############################################################
# Lab 7 Analysis
# Evaluating Which Human Cell Model Best Reflects
# COVID-19 Lung Biology Using RNA-seq Data
#
# Dataset: GSE147507
# Analysis by: Ashley Batugo
############################################################


############################################################
# Load packages
############################################################

library(tidyverse)
library(edgeR)
library(matrixStats)
library(cowplot)
library(DescTools)
library(tximport)
library(reshape2)
library(plotly)
library(patchwork)
library(stringr)
library(tibble)


############################################################
# Load study data
############################################################

studyDesign <- read.table(
  "~/MBMI/CAMB7140/labs/lab7/studyDesign.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
) %>%
  as.tibble()

sampleLabels <- studyDesign$sample

human_raw_reads <- read_tsv(
  "~/MBMI/CAMB7140/labs/lab7/GSE147507_RawReadCounts_Human.tsv"
)

names(human_raw_reads)[1] <- "geneID"

human_raw_reads <- human_raw_reads %>%
  column_to_rownames("geneID") %>%
  as.matrix()

human_DGEList <- DGEList(human_raw_reads)


############################################################
# Convert raw counts to CPM
############################################################

cpm <- edgeR::cpm(human_DGEList)
colSums(cpm)

log2.cpm <- edgeR::cpm(human_DGEList, log = TRUE)


############################################################
# Convert matrix to dataframe
############################################################

log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")

log2.cpm.df.pivot <- pivot_longer(
  log2.cpm.df,
  cols = Series1_NHBE_Mock_1:`Series16_A549-ACE2_SARS-CoV-2_Rux_3`,
  names_to = "samples",
  values_to = "expression"
)


############################################################
# Filtering low expression genes
############################################################

cpm <- cpm(human_DGEList)

keepers <- rowSums(cpm > 1) >= 5

human_DGEList.filtered <- human_DGEList[keepers, ]

log2.cpm.filtered <- cpm(human_DGEList.filtered, log = TRUE)

log2.cpm.filtered.df <- as_tibble(
  log2.cpm.filtered,
  rownames = "geneID"
)

colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)


############################################################
# TMM normalization
############################################################

human_DGEList.filtered.norm <- calcNormFactors(
  human_DGEList.filtered,
  method = "TMM"
)

log2.cpm.filtered.norm <- cpm(
  human_DGEList.filtered.norm,
  log = TRUE
)

log2.cpm.filtered.norm.df <- as_tibble(
  log2.cpm.filtered.norm,
  rownames = "geneID"
)

colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)

log2.cpm.filtered.norm.df.pivot <- pivot_longer(
  log2.cpm.filtered.norm.df,
  cols = -1,
  names_to = "samples",
  values_to = "expression"
)


############################################################
# Differential expression function
############################################################

run_dge_analysis <- function(cell_line){
  
  if (cell_line %in% c("primary NHBE","A549","A549-ACE2","Calu-3")) {
    
    studyDesign_local <- studyDesign %>%
      filter(host_common_name == "human") %>%
      filter(treatment %in% c(
        "Mock treatment",
        "SARS-CoV-2 infected (MOI 2)"
      )) %>%
      mutate(group = case_when(
        treatment == "Mock treatment" ~ "control",
        treatment == "SARS-CoV-2 infected (MOI 2)" ~ "infected"
      ))
    
  } else if (cell_line == "Lung Biopsy") {
    
    studyDesign_local <- studyDesign %>%
      filter(host_common_name == "human") %>%
      filter(tissue_cell_type == "Lung Biopsy") %>%
      mutate(group = case_when(
        str_detect(tolower(subject_status),"healthy") ~ "control",
        str_detect(tolower(subject_status),"covid") ~ "infected"
      ))
  }
  
  if (cell_line == "primary NHBE") {
    
    targets <- studyDesign_local %>%
      filter(cell_line == "primary NHBE")
    
  } else if (cell_line == "A549") {
    
    targets <- studyDesign_local %>%
      filter(cell_line == "A549")
    
  } else if (cell_line == "A549-ACE2") {
    
    targets <- studyDesign_local %>%
      filter(cell_line == "A549-ACE2")
    
  } else if (cell_line == "Calu-3") {
    
    targets <- studyDesign_local %>%
      filter(cell_line == "Calu-3")
    
  } else if (cell_line == "Lung Biopsy") {
    
    targets <- studyDesign_local %>%
      filter(tissue_cell_type == "Lung Biopsy")
  }
  
  samples <- targets$sample
  group <- factor(targets$group)
  
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  
  v.DEGList.filtered.norm <- voom(
    human_DGEList.filtered.norm[, samples],
    design,
    plot = FALSE
  )
  
  fit <- lmFit(v.DEGList.filtered.norm, design)
  
  contrast.matrix <- makeContrasts(
    infection = infected - control,
    levels = design
  )
  
  fits <- contrasts.fit(fit, contrast.matrix)
  ebFit <- eBayes(fits)
  
  myTopHits <- topTable(
    ebFit,
    adjust ="BH",
    coef = 1,
    number = 40000,
    sort.by = "logFC"
  )
  
  myTopHits.df <- myTopHits %>%
    as_tibble(rownames = "geneID")
  
  results <- decideTests(
    ebFit,
    method="global",
    adjust.method="BH",
    p.value=0.01,
    lfc=2
  )
  
  colnames(v.DEGList.filtered.norm$E) <- samples
  
  diffGenes <- v.DEGList.filtered.norm$E[
    results[,1] !=0,
  ]
  
  diffGenes.df <- as_tibble(diffGenes, rownames="geneID")
  
  return(list(
    DESeq2_results = myTopHits.df,
    DEGs = diffGenes.df
  ))
}


############################################################
# Compare cell line DEGs to lung biopsy DEGs
############################################################

run_upregulated_share <- function(
    cell_lines = c("primary NHBE","A549","A549-ACE2","Calu-3")
){
  
  lung_upregulated_genes <- run_dge_analysis("Lung Biopsy")$DESeq2_results %>%
    filter(logFC >= 0) %>%
    filter(adj.P.Val <= 0.05)
  
  cell_line_upregulated_list <- list()
  lung_cell_line_shared_upregulated_list <- list()
  
  for (cell_line in cell_lines) {
    
    cell_line_upregulated_genes <- run_dge_analysis(cell_line)$DESeq2_results %>%
      filter(logFC >= 0) %>%
      filter(adj.P.Val <= 0.05)
    
    cell_line_upregulated_list[[cell_line]] <- cell_line_upregulated_genes
    
    shared_upregulated_genes <- lung_upregulated_genes %>%
      select(geneID) %>%
      inner_join(
        cell_line_upregulated_genes %>% select(geneID),
        by="geneID"
      ) %>%
      mutate(shared_upregulated = n()) %>%
      distinct(shared_upregulated) %>%
      mutate(
        cell_line = cell_line,
        pct_of_all_lung_upregulated =
          (shared_upregulated /
             nrow(lung_upregulated_genes)) * 100
      ) %>%
      select(cell_line,shared_upregulated,pct_of_all_lung_upregulated)
    
    lung_cell_line_shared_upregulated_list[[cell_line]] <-
      shared_upregulated_genes
  }
  
  combined_results <- bind_rows(
    lung_cell_line_shared_upregulated_list
  )
  
  return(list(
    cell_line_upregulated = cell_line_upregulated_list,
    shared_counts = lung_cell_line_shared_upregulated_list,
    combined_results = combined_results
  ))
}


############################################################
# Bar plot of shared DEG percentages
############################################################

dge_results <- run_upregulated_share()$combined_results

dge_plot <- ggplot(
  dge_results %>% mutate(cell_line = factor(cell_line)),
  aes(y = cell_line,
      x = pct_of_all_lung_upregulated,
      fill = cell_line)
) +
  geom_col() +
  xlim(0,100) +
  labs(
    title = "Percentage of Lung Upregulated Genes Shared with Cell Lines",
    x = "Percentage of Lung Upregulated Genes",
    y = "Cell Line"
  ) +
  theme_light() +
  theme(legend.position="none")


############################################################
# PCA analysis
############################################################

filtered_studyDesign <- studyDesign %>%
  filter(host_common_name=="human") %>%
  filter(cell_line %in% c(
    "primary NHBE","A549","A549-ACE2","Calu-3"
  ) | tissue_cell_type=="Lung Biopsy") %>%
  filter(
    treatment %in% c(
      "Mock treatment",
      "SARS-CoV-2 infected (MOI 2)"
    ) |
      str_detect(tolower(subject_status),"healthy") |
      str_detect(tolower(subject_status),"covid")
  ) %>%
  mutate(
    group = case_when(
      treatment=="Mock treatment" ~ "control",
      treatment=="SARS-CoV-2 infected (MOI 2)" ~ "infected",
      str_detect(tolower(subject_status),"healthy") ~ "control",
      str_detect(tolower(subject_status),"covid") ~ "infected"
    )
  ) %>%
  mutate(
    cell_line = case_when(
      tissue_cell_type=="Lung Biopsy" ~ "Lung Biopsy",
      TRUE ~ cell_line
    )
  )

selected_samples <- filtered_studyDesign$sample

expr_matrix <- log2.cpm.filtered.norm[, selected_samples]


############################################################
# Remove housekeeping genes via variance filtering
############################################################

gene_var <- apply(expr_matrix,1,var)

top_genes <- names(sort(gene_var,decreasing=TRUE))[1:1000]

pca_matrix <- expr_matrix[top_genes, ]


############################################################
# Run PCA
############################################################

pca.res <- prcomp(
  t(pca_matrix),
  scale.=TRUE,
  retx=TRUE
)

pc.var <- pca.res$sdev^2
pc.per <- round(pc.var/sum(pc.var)*100,1)

pca.res.df <- as_tibble(
  pca.res$x,
  rownames="sample"
) %>%
  left_join(
    filtered_studyDesign %>%
      mutate(model_condition=paste(cell_line,group,sep="_")) %>%
      select(sample,cell_line,group,model_condition),
    by="sample"
  )


############################################################
# PCA plot
############################################################

pca.plot <- ggplot(pca.res.df) +
  aes(
    x=PC1,
    y=PC2,
    color=cell_line,
    shape=group
  ) +
  geom_point(size=4) +
  stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%)")) +
  ylab(paste0("PC2 (",pc.per[2],"%)")) +
  labs(
    title="PCA plot (Top 1000 most variable genes)"
  ) +
  coord_fixed() +
  theme_bw() +
  guides(
    color = guide_legend(nrow=1),
    shape = guide_legend(nrow=1)
  ) +
  theme(
    legend.position="bottom",
    legend.box="vertical"
  )


############################################################
# Combine plots
############################################################

pca_sub <- pca.plot + labs(title=NULL)
dge_sub <- dge_plot + labs(title=NULL)

combined_plot <-
  (pca_sub | dge_sub) +
  plot_layout(widths=c(2,1)) +
  plot_annotation(
    title="Which Human Cell Model Best Reflects Infected COVID-19 Lung Tissue?",
    subtitle="Comparison of SARS-CoV-2 infected cell lines and human lung tissue",
    caption=paste(
      "PCA computed using top 1000 most variable genes.",
      "",
      "Conclusion: Calu-3 infected with SARS-CoV-2 most closely resembles COVID lung tissue.",
      "",
      "Dataset: GSE147507",
      "Analysis by: Ashley Batugo",
      sep="\n"
    )
  )

combined_plot


############################################################
# Save final figure
############################################################

ggsave(
  "lab7_results.png",
  combined_plot,
  width=14,
  height=7,
  dpi=300
)