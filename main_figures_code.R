
# ============================================================
# Main figures for data/code upload
# Figures included: Fig2, Fig3, Fig4
# Notes:
#   1) This script reorganizes the uploaded working code into
#      figure-wise blocks for direct code deposition.
#   2) Fig2 and Fig4 are data-driven.
#   3) Fig3 mixes data-driven panels with the manually curated
#      mechanism panels used in the final figure version.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(vegan)
  library(randomForest)
  library(pROC)
  library(pheatmap)
  library(patchwork)
  library(cowplot)
  library(ggplotify)
  library(limma)
  library(readr)
  library(janitor)
  library(RColorBrewer)
  library(scales)
  library(stringr)
})

# ---------------------------
# Global settings
# ---------------------------
dir.create("figures/main", recursive = TRUE, showWarnings = FALSE)

pal_group <- c("HC" = "#7FA0D0", "NC" = "#7FA0D0", "SLE" = "#CF3D3E")
pal_direction <- c("SLE Enriched" = "#CF3D3E", "SLE Depleted" = "#7FA0D0", "NS" = "grey80")
pal_updown <- c("Up_in_SLE" = "#C43E3D", "Down_in_SLE" = "#7E9FD3")

theme_clean <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "black", face = "bold"),
      axis.title = element_text(color = "black", face = "bold"),
      plot.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold")
    )
}

# ---------------------------
# Shared data loaders
# ---------------------------
species <- read.csv("species.csv", row.names = 1, check.names = FALSE)
genus   <- read.csv("genus.csv",   row.names = 1, check.names = FALSE)
ko      <- read.csv("ko.csv",      row.names = 1, check.names = FALSE)

species[species < 0] <- 0

if (file.exists("group.csv")) {
  grp_raw <- read.csv("group.csv", stringsAsFactors = FALSE, check.names = FALSE)
  g_names <- names(grp_raw)
  sample_col <- g_names[grepl("sample|id|name", g_names, ignore.case = TRUE) &
                          !grepl("group|grp|condition", g_names, ignore.case = TRUE)][1]
  group_col  <- g_names[grepl("group|grp|condition", g_names, ignore.case = TRUE)][1]

  metadata <- grp_raw %>%
    transmute(
      SampleID = as.character(.data[[sample_col]]),
      Group    = as.character(.data[[group_col]])
    )
  metadata$Group[metadata$Group == "NC"] <- "HC"
  metadata$Group <- factor(metadata$Group, levels = c("HC", "SLE"))
  rownames(metadata) <- metadata$SampleID
} else {
  metadata <- data.frame(
    SampleID = colnames(species),
    Group = ifelse(grepl("SLE", colnames(species), ignore.case = TRUE), "SLE", "HC")
  )
  metadata$Group <- factor(metadata$Group, levels = c("HC", "SLE"))
  rownames(metadata) <- metadata$SampleID
}

common_samples <- intersect(colnames(ko), rownames(metadata))
metadata2 <- metadata[common_samples, , drop = FALSE]

# differential genera for RF feature construction
diff_genus <- apply(genus, 1, function(x) {
  g1 <- x[metadata$Group == "HC"]
  g2 <- x[metadata$Group == "SLE"]
  pval <- tryCatch(wilcox.test(g1, g2)$p.value, error = function(e) 1)
  log2fc <- log2((mean(g2, na.rm = TRUE) + 1e-9) / (mean(g1, na.rm = TRUE) + 1e-9))
  c(p.value = pval, log2fc = log2fc)
}) %>% t() %>% as.data.frame()

diff_genus$Genus <- rownames(diff_genus)
top_6_genera <- diff_genus %>%
  arrange(p.value) %>%
  slice_head(n = 6) %>%
  pull(Genus)

# ============================================================
# Fig2. Functional separation and classifier performance
# ============================================================

# ---- Fig2a: KO PCA ----
ko_clean <- ko[, common_samples, drop = FALSE]
ko_clean <- ko_clean[apply(ko_clean, 1, var, na.rm = TRUE) > 0, , drop = FALSE]

pca_res <- prcomp(t(ko_clean), scale. = TRUE)
var_exp <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)

pca_df <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  Group = metadata2$Group
)

permanova <- adonis2(vegdist(t(ko_clean), method = "bray") ~ Group, data = metadata2)

p_fig2a <- ggplot(pca_df, aes(PC1, PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  scale_color_manual(values = pal_group[c("HC", "SLE")]) +
  annotate(
    "text",
    x = max(pca_df$PC1) * 0.95,
    y = max(pca_df$PC2) * 0.95,
    hjust = 1, vjust = 1,
    label = sprintf("PERMANOVA P=%.3f, R²=%.3f",
                    permanova$`Pr(>F)`[1], permanova$R2[1]),
    fontface = "bold.italic",
    size = 5
  ) +
  labs(x = paste0("PC1 (", var_exp[1], "%)"),
       y = paste0("PC2 (", var_exp[2], "%)"),
       title = "a") +
  theme_clean(12) +
  theme(legend.position = c(0.86, 0.15))

# ---- Fig2b: volcano ----
res <- apply(ko_clean, 1, function(x) {
  g1 <- x[metadata2$Group == "HC"]
  g2 <- x[metadata2$Group == "SLE"]
  pval <- tryCatch(wilcox.test(g1, g2)$p.value, error = function(e) 1)
  log2fc <- log2((mean(g2, na.rm = TRUE) + 1e-9) / (mean(g1, na.rm = TRUE) + 1e-9))
  c(p.value = pval, log2FC = log2fc)
}) %>% t() %>% as.data.frame()

res$FDR <- p.adjust(res$p.value, method = "fdr")
fc_cutoff <- 0.58
p_cutoff  <- 0.05

res$label <- "NS"
res$label[res$log2FC >  fc_cutoff & res$FDR < p_cutoff] <- "SLE Enriched"
res$label[res$log2FC < -fc_cutoff & res$FDR < p_cutoff] <- "SLE Depleted"
sig_genes_list <- rownames(res[res$label != "NS", , drop = FALSE])

label_data <- res %>%
  tibble::rownames_to_column("KO") %>%
  filter(label != "NS") %>%
  arrange(FDR)

p_fig2b <- ggplot(res %>% tibble::rownames_to_column("KO"),
                  aes(log2FC, -log10(FDR))) +
  geom_point(aes(color = label), alpha = 0.7, size = 2.5) +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed") +
  geom_text_repel(
    data = label_data,
    aes(label = KO, color = label),
    max.overlaps = Inf,
    size = 4.5,
    box.padding = 0.35,
    point.padding = 0.2,
    segment.color = "grey50",
    show.legend = FALSE
  ) +
  scale_color_manual(values = pal_direction) +
  labs(x = "log2FC", y = "-log10(FDR)", title = "b") +
  theme_clean(12) +
  theme(legend.position = c(0.82, 0.86))

# ---- Fig2c / Fig2d: repeated RF CV + importance ----
common_rf <- intersect(common_samples, colnames(genus))
feat_genus <- genus[top_6_genera, common_rf, drop = FALSE]
feat_ko    <- ko[sig_genes_list, common_rf, drop = FALSE]

rf_input <- cbind(t(feat_genus), t(feat_ko)) %>% as.data.frame()
rf_input$Group <- factor(metadata[common_rf, "Group"], levels = c("HC", "SLE"))
colnames(rf_input) <- make.names(colnames(rf_input))

one_cv_run <- function(data, seed, k = 5) {
  set.seed(seed)
  y <- data$Group
  idx_by_class <- split(seq_len(nrow(data)), y)
  fold_id <- integer(nrow(data))
  for (cls in names(idx_by_class)) {
    ids <- sample(idx_by_class[[cls]])
    fold_id[ids] <- rep(seq_len(k), length.out = length(ids))
  }

  oof_prob <- numeric(nrow(data))
  for (f in seq_len(k)) {
    train_idx <- which(fold_id != f)
    test_idx  <- which(fold_id == f)
    rf_fold <- randomForest(
      Group ~ .,
      data = data[train_idx, , drop = FALSE],
      ntree = 1000,
      mtry = floor(sqrt(ncol(data) - 1))
    )
    oof_prob[test_idx] <- predict(rf_fold, data[test_idx, , drop = FALSE], type = "prob")[, "SLE"]
  }

  roc_obj <- suppressMessages(pROC::roc(y, oof_prob, levels = c("HC", "SLE"), direction = "<"))
  list(
    auc = as.numeric(pROC::auc(roc_obj)),
    prob = oof_prob,
    roc = roc_obj
  )
}

set.seed(2026)
n_repeats <- 100
cv_runs <- lapply(seq_len(n_repeats), function(i) one_cv_run(rf_input, seed = i, k = 5))
auc_vec <- vapply(cv_runs, `[[`, numeric(1), "auc")
auc_mean <- mean(auc_vec)
auc_ci   <- quantile(auc_vec, c(0.025, 0.975))

median_id <- order(abs(auc_vec - median(auc_vec)))[1]
roc_rep <- cv_runs[[median_id]]$roc

p_fig2c <- ggroc(roc_rep, linewidth = 1.2, colour = "#CF3D3E") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60") +
  annotate(
    "text",
    x = 0.47, y = 0.22,
    hjust = 0,
    label = sprintf("Representative AUC = %.3f\nMean AUC = %.3f\n95%% CI %.3f-%.3f",
                    as.numeric(pROC::auc(roc_rep)), auc_mean, auc_ci[1], auc_ci[2]),
    size = 5,
    fontface = "bold"
  ) +
  labs(x = "1 - Specificity", y = "Sensitivity", title = "c") +
  theme_clean(12)

set.seed(123)
rf_full <- randomForest(Group ~ ., data = rf_input, importance = TRUE, ntree = 1000)
imp <- importance(rf_full)
imp_df <- data.frame(
  Feature = rownames(imp),
  Accuracy = imp[, "MeanDecreaseAccuracy"],
  row.names = NULL
) %>%
  arrange(desc(Accuracy)) %>%
  slice_head(n = 15)

p_fig2d <- ggplot(imp_df, aes(reorder(Feature, Accuracy), Accuracy)) +
  geom_col(fill = "#F4A300") +
  coord_flip() +
  labs(x = NULL, y = "Mean decrease in accuracy", title = "d") +
  theme_clean(12)

fig2 <- (p_fig2a | p_fig2b) / (p_fig2c | p_fig2d)
ggsave("figures/main/Fig2.png", fig2, width = 14, height = 10, dpi = 300)
ggsave("figures/main/Fig2.pdf", fig2, width = 14, height = 10, device = cairo_pdf)

# ============================================================
# Fig3. Metabolic mechanism and multi-omics linkage
# ============================================================

# ---- proteomics shared inputs ----
read_clean <- function(path) {
  read_csv(path, show_col_types = FALSE, na = c("NA", "NaN", "", "null", "NULL")) %>%
    rename_with(~ gsub("\ufeff", "", .x)) %>%
    janitor::clean_names()
}

prot_raw <- read_clean("dia-proteinSummary-all.csv")
annot    <- read_clean("annotation_allprotein.csv")
grp_raw  <- read_clean("group.csv")

g_names <- names(grp_raw)
sample_col <- g_names[grepl("sample|id|name", g_names) & !grepl("group|grp|condition", g_names)][1]
group_col  <- g_names[grepl("group|grp|condition", g_names)][1]

grp <- grp_raw %>%
  transmute(
    Sample = janitor::make_clean_names(.data[[sample_col]]),
    Group  = factor(ifelse(.data[[group_col]] == "NC", "HC", as.character(.data[[group_col]])),
                    levels = c("HC", "SLE"))
  )

prot_id_col <- names(prot_raw)[grepl("^protein|^accession", names(prot_raw))][1]
prot_m <- prot_raw %>%
  rename(ProteinID = all_of(prot_id_col)) %>%
  mutate(ProteinID = gsub("\\[|\\]", "", ProteinID)) %>%
  select(ProteinID, any_of(grp$Sample))

annot_id_col <- names(annot)[grepl("^protein|^accession", names(annot))][1]
annot <- annot %>%
  rename(ProteinID = all_of(annot_id_col)) %>%
  mutate(ProteinID = gsub("\\[|\\]", "", ProteinID))

ko_col <- names(annot)[grepl("kegg|ko_id|ko_code", names(annot), ignore.case = TRUE)][1]
genus_col_anno <- names(annot)[grepl("genus|taxonomy", names(annot), ignore.case = TRUE)][1]

# ---- relaxed limma for proteomics ----
expr <- as.matrix(prot_m %>% column_to_rownames("ProteinID"))
mode(expr) <- "numeric"
expr <- expr[rowSums(!is.na(expr)) > 0, , drop = FALSE]
expr[is.na(expr)] <- min(expr, na.rm = TRUE)

design <- model.matrix(~ 0 + grp$Group)
colnames(design) <- levels(grp$Group)

fit <- lmFit(expr, design)
contr <- makeContrasts(SLEvsHC = SLE - HC, levels = design)
fit2 <- eBayes(contrasts.fit(fit, contr))

deg <- topTable(fit2, coef = "SLEvsHC", number = Inf) %>%
  rownames_to_column("ProteinID") %>%
  as_tibble() %>%
  mutate(Direction = case_when(
    P.Value < 0.05 & logFC >  0.58 ~ "Up_in_SLE",
    P.Value < 0.05 & logFC < -0.58 ~ "Down_in_SLE",
    TRUE ~ "NS"
  ))

# ---- Fig3a: KEGG enrichment ----
p_fig3a <- ggplot() + theme_void() + ggtitle("a")
if (requireNamespace("clusterProfiler", quietly = TRUE) && !is.na(ko_col)) {
  suppressPackageStartupMessages(library(clusterProfiler))
  deg_ko_df <- deg %>%
    filter(Direction != "NS") %>%
    left_join(annot, by = "ProteinID") %>%
    mutate(KO = str_extract(.data[[ko_col]], "K\\d+")) %>%
    filter(!is.na(KO)) %>%
    distinct(KO, .keep_all = TRUE)

  try({
    ekegg <- enrichKEGG(
      gene = deg_ko_df$KO,
      organism = "ko",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
    if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
      ekegg_df <- as.data.frame(ekegg) %>% slice_head(n = 12)
      ekegg_df$Description <- factor(ekegg_df$Description, levels = rev(ekegg_df$Description))
      p_fig3a <- ggplot(ekegg_df, aes(GeneRatio, Description, size = Count, colour = p.adjust)) +
        geom_point() +
        scale_colour_gradient(low = "#E64B35", high = "#3C7DC4", trans = "reverse") +
        labs(x = "GeneRatio", y = NULL, title = "a", colour = "p.adjust") +
        theme_clean(11)
    }
  }, silent = TRUE)
}

# ---- Fig3b / Fig3d / Fig3f: curated mechanism panels ----
mechanism_tbl <- tibble::tribble(
  ~Gene,   ~logFC,    ~P.Value,   ~Category,
  "IMPDH",  664.291,  6.26e-8,    "Purine",
  "GuaA",   281.655,  0.00182,    "Purine",
  "ArcA",   150.000,  0.00500,    "ADI",
  "OTC",    200.000,  0.00800,    "ADI",
  "ArcC",   180.000,  0.00300,    "ADI",
  "ArgG",  -253.775,  0.00596,    "Arg_Synthesis",
  "DnaE",  -190.444,  0.00520,    "Replication",
  "NrdA",  -248.447,  0.00289,    "Replication",
  "MalP",  -350.000,  0.00100,    "Energy",
  "TktA",  -200.000,  0.01000,    "Energy"
) %>%
  mutate(
    Direction = case_when(logFC > 0 ~ "Up", TRUE ~ "Down"),
    Significance = case_when(
      P.Value < 0.001 ~ "***",
      P.Value < 0.01  ~ "**",
      P.Value < 0.05  ~ "*",
      TRUE ~ ""
    ),
    Gene = factor(Gene, levels = rev(Gene))
  )

xmax_mech <- max(abs(mechanism_tbl$logFC)) * 1.25
band_df <- tibble::tribble(
  ~ymin, ~ymax, ~fill,      ~lab_x, ~lab_y, ~label,                         ~col,
   0.5,   2.5,  "#E8F0FE",  -xmax_mech*0.90, 1.5, "Energy\nMetabolism",     "#3C5488",
   2.5,   4.5,  "#E8F0FE",  -xmax_mech*0.90, 3.5, "DNA\nReplication",       "#3C5488",
   4.5,   5.5,  "#E8F0FE",  -xmax_mech*0.90, 5.0, "Arg Synthesis\n(Blocked)","#3C5488",
   5.5,   8.5,  "#FDECEA",   xmax_mech*0.90, 7.0, "ADI Pathway\n(Arg Catabolism)", "#E64B35",
   8.5,  10.5,  "#FDECEA",   xmax_mech*0.90, 9.5, "Purine/GTP\nSynthesis",  "#E64B35"
)

p_fig3b <- ggplot(mechanism_tbl, aes(logFC, Gene)) +
  geom_rect(data = band_df, inherit.aes = FALSE,
            aes(xmin = -xmax_mech, xmax = xmax_mech, ymin = ymin, ymax = ymax),
            fill = band_df$fill, alpha = 0.5) +
  geom_vline(xintercept = 0, colour = "grey40") +
  geom_segment(aes(x = 0, xend = logFC, y = Gene, yend = Gene, colour = Direction), linewidth = 1.5) +
  geom_point(aes(colour = Direction), size = 4) +
  geom_text(aes(label = Significance,
                x = ifelse(logFC > 0, logFC + xmax_mech * 0.05, logFC - xmax_mech * 0.05)),
            fontface = "bold", size = 4.5) +
  geom_text(data = band_df, inherit.aes = FALSE,
            aes(x = lab_x, y = lab_y, label = label), colour = band_df$col,
            fontface = "bold", size = 3.8) +
  scale_colour_manual(values = c("Up" = "#E64B35", "Down" = "#3C5488")) +
  labs(x = expression(bold(Delta ~ "Abundance (SLE - HC)")), y = NULL, title = "b") +
  theme_clean(11) +
  theme(legend.position = "none",
        axis.text.y = element_text(face = "bold.italic"))

arg_balance <- data.frame(
  Process = c("Synthesis\n(ArgG)", "Catabolism\n(ADI pathway)"),
  Value   = c(-253.775, 150 + 200 + 180),
  Type    = c("Down", "Up")
)
p_fig3d <- ggplot(arg_balance, aes(Process, Value, fill = Type)) +
  geom_hline(yintercept = 0) +
  geom_col(width = 0.6, colour = "black") +
  scale_fill_manual(values = c("Down" = "#3C5488", "Up" = "#FF4D33")) +
  labs(x = NULL, y = expression(Delta ~ "Abundance (SLE - HC)"), title = "d") +
  theme_clean(11) +
  theme(legend.position = "none")

nodes <- tibble::tribble(
  ~name, ~x, ~y, ~size, ~col,
  "Arginine",          2, 6, 20, "#3C5488",
  "Citrulline",        5, 6, 16, "grey60",
  "Ornithine",         8, 6, 16, "grey60",
  "ATP",              11, 6, 18, "#FF6248",
  "Aspartate",         2, 4, 16, "grey60",
  "Argininosuccinate", 5, 4, 16, "grey60",
  "IMP",               2, 2, 16, "grey60",
  "XMP",               5, 2, 16, "grey60",
  "GMP",               8, 2, 16, "grey60",
  "GTP",              11, 2, 20, "#FF6248",
  "DNA",              14, 2, 18, "#3C5488",
  "c-di-GMP",         11, 0.4, 18, "#FF6248"
)
edges <- tibble::tribble(
  ~x, ~xend, ~y, ~yend, ~enzyme,   ~col,
   3,   4,    6,   6,    "ArcA",   "#E64B35",
   6,   7,    6,   6,    "OTC",    "#E64B35",
   9,  10,    6,   6,    "ArcC",   "#E64B35",
   2,   4,    4,   4,    "ArgG",   "#3C5488",
   3,   4,    2,   2,    "IMPDH",  "#E64B35",
   6,   7,    2,   2,    "GuaA",   "#E64B35",
   9,  10,    2,   2,    "Kinases","grey40",
  11,  13,    2,   2,    "NrdA",   "#3C5488",
  11,  11,  1.5, 0.7,    "DGC/PDE","#E64B35"
)

p_fig3f <- ggplot() +
  annotate("rect", xmin = 0.5, xmax = 12.5, ymin = 5, ymax = 7, fill = "#FDECEA", alpha = 0.25, color = NA) +
  annotate("rect", xmin = 0.5, xmax = 7.0, ymin = 3, ymax = 5, fill = "#E8F0FE", alpha = 0.28, color = NA) +
  annotate("rect", xmin = 0.5, xmax = 15.5, ymin = -0.5, ymax = 3, fill = "#FDECEA", alpha = 0.12, color = NA) +
  annotate("text", x = 6.5, y = 6.85, label = "ADI Pathway  (Arginine Catabolism → ATP)",
           fontface = "bold", size = 4, colour = "#E64B35") +
  annotate("text", x = 4.0, y = 4.85, label = "Arginine Biosynthesis  (Blocked)",
           fontface = "bold", size = 4, colour = "#3C5488") +
  annotate("text", x = 8.2, y = 2.85, label = "Purine / GTP Biosynthesis",
           fontface = "bold", size = 4, colour = "#E64B35") +
  geom_segment(data = edges, aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               linewidth = 1.1, colour = edges$col) +
  annotate("text", x = 3.0, y = 4.0, label = "\u2716", size = 8, colour = "#3C5488") +
  annotate("text", x = 12.0, y = 2.0, label = "\u2716", size = 8, colour = "#3C5488") +
  geom_point(data = nodes, aes(x = x, y = y), colour = nodes$col, size = nodes$size) +
  geom_text(data = nodes, aes(x = x, y = y, label = name),
            colour = "white", fontface = "bold", size = 3.3) +
  geom_text(data = edges %>% filter(y == yend), aes(x = (x + xend)/2, y = y + 0.38, label = enzyme),
            colour = edges$col[edges$y == edges$yend], fontface = "bold.italic", size = 3.6) +
  geom_text(data = edges %>% filter(y != yend), aes(x = x + 0.75, y = (y + yend)/2, label = enzyme),
            colour = edges$col[edges$y != edges$yend], fontface = "bold.italic", size = 3.6) +
  annotate("segment", x = 5, xend = 8.2, y = 6.1, yend = 4.4,
           linetype = "dashed", colour = "grey55") +
  annotate("segment", x = 5, xend = 2.5, y = 6.0, yend = 4.3,
           linetype = "dashed", colour = "grey55") +
  annotate("segment", x = 11.0, xend = 11.0, y = 5.8, yend = 2.3,
           linetype = "dotted", colour = "#7E6148") +
  annotate("text", x = 12.2, y = 4.0, label = "NH3 / Energy",
           fontface = "bold.italic", size = 3.3, colour = "#7E6148") +
  coord_cartesian(xlim = c(-0.5, 16), ylim = c(-0.8, 7.4)) +
  theme_void() +
  labs(title = "f")

# ---- Fig3c: KO-protein genus source ----
p_fig3c <- ggplot() + theme_void() + ggtitle("c")
if (!is.na(ko_col) && !is.na(genus_col_anno)) {
  contribution_raw <- deg %>%
    left_join(annot, by = "ProteinID") %>%
    mutate(KO = str_extract(.data[[ko_col]], "K\\d+")) %>%
    filter(KO %in% sig_genes_list) %>%
    filter(!is.na(.data[[genus_col_anno]]) & .data[[genus_col_anno]] != "" & .data[[genus_col_anno]] != "-") %>%
    filter(.data[[genus_col_anno]] != "Bacteria") %>%
    rename(Genus = !!sym(genus_col_anno))

  if (nrow(contribution_raw) > 0) {
    top_genera <- contribution_raw %>% count(Genus, sort = TRUE) %>% slice_head(n = 12) %>% pull(Genus)
    plot_data <- contribution_raw %>%
      mutate(Genus_Label = ifelse(Genus %in% top_genera, Genus, "Others (Minor)")) %>%
      count(KO, Genus_Label)

    n_cols <- length(unique(plot_data$Genus_Label))
    pal <- colorRampPalette(brewer.pal(min(12, n_cols), "Paired"))(n_cols)
    if ("Others (Minor)" %in% unique(plot_data$Genus_Label)) {
      plot_data$Genus_Label <- forcats::fct_relevel(plot_data$Genus_Label, "Others (Minor)", after = Inf)
      pal[length(pal)] <- "grey70"
    }

    p_fig3c <- ggplot(plot_data, aes(KO, n, fill = Genus_Label)) +
      geom_col(color = "white", linewidth = 0.15) +
      scale_fill_manual(values = pal, name = "Genus Source") +
      labs(x = NULL, y = "Protein Count", title = "c") +
      theme_clean(10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "right",
            legend.text = element_text(face = "italic", size = 8))
  }
}

# ---- Fig3e: enzyme-contributing genera ----
enzyme_genus_map <- tibble::tribble(
  ~Enzyme,  ~KO,       ~Genus,
  "OTC",    "K00611",  "Anaerostipes",
  "OTC",    "K00611",  "Roseburia",
  "ArgG",   "K01940",  "Ruminococcus",
  "ArcC",   "K00926",  "Anaerobutyricum",
  "ArcC",   "K00926",  "Roseburia",
  "IMPDH",  "K00088",  "Phascolarctobacterium",
  "IMPDH",  "K00088",  "Parasutterella",
  "IMPDH",  "K00088",  "Anaerotignum",
  "IMPDH",  "K00088",  "Bifidobacterium"
)

group_for_genus <- read.csv("group.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
  select(1, 2)
names(group_for_genus) <- c("SampleID", "Group")
group_for_genus$Group[group_for_genus$Group == "NC"] <- "HC"

genus_long <- genus %>%
  rownames_to_column("Genus") %>%
  filter(Genus %in% unique(enzyme_genus_map$Genus)) %>%
  pivot_longer(-Genus, names_to = "SampleID", values_to = "Abundance") %>%
  left_join(group_for_genus, by = "SampleID") %>%
  left_join(enzyme_genus_map, by = "Genus") %>%
  mutate(Group = factor(Group, levels = c("HC", "SLE")))

stat_results <- genus_long %>%
  group_by(Genus, Enzyme, KO) %>%
  summarise(
    HC_mean  = mean(Abundance[Group == "HC"],  na.rm = TRUE),
    SLE_mean = mean(Abundance[Group == "SLE"], na.rm = TRUE),
    log2FC = log2((SLE_mean + 1e-6) / (HC_mean + 1e-6)),
    p_value = tryCatch(wilcox.test(Abundance ~ Group)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ ""
    ),
    bar_color = case_when(
      !is.na(p_value) & p_value < 0.05 & log2FC < 0 ~ "Decreased*",
      !is.na(p_value) & p_value >= 0.05 & log2FC < 0 ~ "Decreased (ns)",
      !is.na(p_value) & p_value < 0.05 & log2FC > 0 ~ "Increased*",
      TRUE ~ "Increased (ns)"
    ),
    Enzyme = factor(Enzyme, levels = c("OTC", "ArgG", "ArcC", "IMPDH"))
  )

p_fig3e <- ggplot(stat_results, aes(Genus, log2FC, fill = bar_color)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = significance), vjust = -0.4, size = 4) +
  facet_wrap(~Enzyme, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = c(
    "Decreased*" = "#3498DB",
    "Decreased (ns)" = "#AED6F1",
    "Increased*" = "#E74C3C",
    "Increased (ns)" = "#F5B7B1"
  )) +
  labs(x = NULL, y = "log2 Fold Change (SLE / HC)", title = "e") +
  theme_clean(10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
        legend.position = "right",
        strip.background = element_rect(fill = "grey95"))

top_fig3 <- p_fig3a + p_fig3b + p_fig3d + plot_layout(widths = c(1.0, 1.5, 1.0))
mid_fig3 <- p_fig3c + p_fig3e + plot_layout(widths = c(1.2, 1.8))
fig3 <- (top_fig3 / mid_fig3 / p_fig3f) + plot_layout(heights = c(1.0, 1.1, 0.95))
ggsave("figures/main/Fig3.png", fig3, width = 16, height = 14, dpi = 300)
ggsave("figures/main/Fig3.pdf", fig3, width = 16, height = 14, device = cairo_pdf)

# ============================================================
# Fig4. Metaproteomic COG architecture
# ============================================================

deg_strict <- topTable(fit2, coef = "SLEvsHC", number = Inf) %>%
  rownames_to_column("ProteinID") %>%
  as_tibble() %>%
  mutate(Direction = case_when(
    adj.P.Val < 0.05 & logFC >  1 ~ "Up_in_SLE",
    adj.P.Val < 0.05 & logFC < -1 ~ "Down_in_SLE",
    TRUE ~ "NS"
  ))

cog_col <- names(annot)[grepl("cog_cat|cog_class|nog_cat|functional_cat|cog_code", names(annot), ignore.case = TRUE)][1]
cog_desc_col <- names(annot)[grepl("cog_desc|nog_desc|function|description|def", names(annot), ignore.case = TRUE)][1]

deg_cog <- deg_strict %>%
  filter(Direction != "NS") %>%
  inner_join(
    annot %>% select(ProteinID, COG_Code = all_of(cog_col), COG_Desc = all_of(cog_desc_col)),
    by = "ProteinID"
  ) %>%
  mutate(COG_Code = str_extract(COG_Code, "[A-Z]")) %>%
  filter(!is.na(COG_Code))

# ---- Fig4a: mirrored counts by COG ----
cog_order <- c("C","E","F","G","H","I","J","K","L","M","O","Q","R","S","T","V")
panel_a_df <- deg_cog %>%
  count(COG_Code, Direction) %>%
  mutate(
    Count = ifelse(Direction == "Down_in_SLE", -n, n),
    COG_Code = factor(COG_Code, levels = cog_order)
  ) %>%
  filter(!is.na(COG_Code))

p_fig4a <- ggplot(panel_a_df, aes(COG_Code, Count, fill = Direction)) +
  geom_col(width = 0.72) +
  geom_hline(yintercept = 0, linewidth = 0.6) +
  scale_fill_manual(values = c("Down_in_SLE" = "#7FA0D0", "Up_in_SLE" = "#EA3B37"),
                    labels = c("Down in SLE", "Up in SLE")) +
  scale_y_continuous(labels = abs) +
  labs(x = "COG Category", y = "Number of Differential Proteins", title = "a", fill = NULL) +
  theme_clean(12) +
  theme(legend.position = c(0.78, 0.12))

# ---- Fig4b / Fig4c: COG G and J top features ----
target_df <- deg_cog %>%
  filter(COG_Code %in% c("G", "J")) %>%
  left_join(annot, by = "ProteinID")

desc_col <- names(target_df)[grepl("swissprot|description|product", names(target_df), ignore.case = TRUE)][1]
target_df$Display_Name <- target_df$ProteinID
if (!is.na(desc_col)) {
  clean_desc <- target_df[[desc_col]] %>%
    str_remove(" OS=.*") %>%
    str_remove(" OX=.*") %>%
    str_remove(" GN=.*") %>%
    str_remove("^sp\\|.*\\|") %>%
    str_remove("RecName: Full=") %>%
    str_remove(";.*") %>%
    str_trim()
  idx <- !is.na(clean_desc) & clean_desc != ""
  target_df$Display_Name[idx] <- clean_desc[idx]
}
target_df$Display_Name <- str_trunc(target_df$Display_Name, 18, ellipsis = "")

plot_cog_panel <- function(code, title_text, tag_letter, n_show = 12) {
  df <- target_df %>%
    filter(COG_Code == code) %>%
    group_by(Display_Name, Direction) %>%
    summarise(logFC = logFC[which.max(abs(logFC))], .groups = "drop") %>%
    arrange(logFC) %>%
    slice_tail(n = n_show)

  ggplot(df, aes(logFC, reorder(Display_Name, logFC), fill = Direction)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, linewidth = 0.6) +
    scale_fill_manual(values = c("Down_in_SLE" = "#6D9ACB", "Up_in_SLE" = "#EA3B37"),
                      guide = "none") +
    labs(x = "Log2 Fold Change (SLE / HC)", y = NULL, title = tag_letter) +
    annotate("label", x = 0, y = Inf, label = title_text, vjust = 1.6, fontface = "bold.italic",
             fill = "grey92") +
    theme_clean(11)
}

p_fig4b <- plot_cog_panel("G", "G (Carbohydrate transport and metabolism)", "b", n_show = 12)
p_fig4c <- plot_cog_panel("J", "J (Translation, ribosomal structure)", "c", n_show = 15)

fig4 <- p_fig4a | (p_fig4b / p_fig4c)
ggsave("figures/main/Fig4.png", fig4, width = 16, height = 10, dpi = 300)
ggsave("figures/main/Fig4.pdf", fig4, width = 16, height = 10, device = cairo_pdf)

message("Main figures finished: Fig2, Fig3, Fig4")
