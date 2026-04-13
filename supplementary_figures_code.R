
# ============================================================
# Supplementary figures for data/code upload
# Figures included: FigS2, FigS3, FigS4
# FigS1 note:
#   The final FigS1 image was provided in this conversation, but
#   its source R code was not included among the uploaded scripts.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(vegan)
  library(randomForest)
  library(pROC)
  library(patchwork)
  library(data.table)
  library(scales)
  library(grid)
})

dir.create("figures/supplementary", recursive = TRUE, showWarnings = FALSE)

pal_group <- c("HC" = "#7FA0D0", "NC" = "#7FA0D0", "QC" = "#28B89A", "SLE" = "#CF3D3E")
pal_module <- c(
  "cdiGMP" = "#00A087",
  "biofilm_EPS" = "#3C5488",
  "adhesion" = "#F39B7F",
  "motility" = "#8491B4"
)
pal_direction <- c("Up_in_SLE" = "#E64B35", "Down_in_SLE" = "#3C5488")

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

# shared microbiome tables
species <- read.csv("species.csv", row.names = 1, check.names = FALSE)
genus   <- read.csv("genus.csv", row.names = 1, check.names = FALSE)
ko      <- read.csv("ko.csv", row.names = 1, check.names = FALSE)

if (file.exists("group.csv")) {
  group_raw <- read.csv("group.csv", stringsAsFactors = FALSE, check.names = FALSE)
  names(group_raw)[1:2] <- c("SampleID", "Group")
  group_raw$Group[group_raw$Group == "NC"] <- "HC"
  metadata <- group_raw %>%
    distinct(SampleID, Group) %>%
    mutate(Group = factor(Group, levels = c("HC", "SLE")))
  rownames(metadata) <- metadata$SampleID
} else {
  metadata <- data.frame(
    SampleID = colnames(species),
    Group = ifelse(grepl("SLE", colnames(species), ignore.case = TRUE), "SLE", "HC")
  )
  metadata$Group <- factor(metadata$Group, levels = c("HC", "SLE"))
  rownames(metadata) <- metadata$SampleID
}

# ============================================================
# FigS2. Community structure and top differential genera
# ============================================================

# ---- panel a: PCoA ----
bray <- vegdist(t(species), method = "bray")
pcoa_res <- cmdscale(bray, eig = TRUE, k = 2)
var_exp <- round(pcoa_res$eig / sum(pcoa_res$eig) * 100, 1)

pcoa_df <- data.frame(
  PC1 = pcoa_res$points[, 1],
  PC2 = pcoa_res$points[, 2],
  Group = metadata[colnames(species), "Group"]
)
permanova <- adonis2(bray ~ Group, data = metadata[colnames(species), , drop = FALSE])

p_figS2a <- ggplot(pcoa_df, aes(PC1, PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.85) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  scale_color_manual(values = pal_group[c("HC", "SLE")]) +
  annotate("text",
           x = min(pcoa_df$PC1) + 0.05 * diff(range(pcoa_df$PC1)),
           y = min(pcoa_df$PC2) + 0.10 * diff(range(pcoa_df$PC2)),
           hjust = 0,
           label = sprintf("PERMANOVA: P = %.3f, R² = %.3f",
                           permanova$`Pr(>F)`[1], permanova$R2[1]),
           fontface = "bold.italic",
           size = 4) +
  labs(x = paste0("PC1 (", var_exp[1], "%)"),
       y = paste0("PC2 (", var_exp[2], "%)"),
       title = "a") +
  theme_clean(11)

# ---- panel b/c: differential genera ----
diff_res <- apply(genus, 1, function(x) {
  g1 <- x[metadata$Group == "HC"]
  g2 <- x[metadata$Group == "SLE"]
  pval <- tryCatch(wilcox.test(g1, g2)$p.value, error = function(e) 1)
  log2fc <- log2((mean(g2, na.rm = TRUE) + 1e-9) / (mean(g1, na.rm = TRUE) + 1e-9))
  c(p.value = pval, log2fc = log2fc)
}) %>% t() %>% as.data.frame()

diff_res$Genus <- rownames(diff_res)
top_6_genera <- diff_res %>% arrange(p.value) %>% slice_head(n = 6) %>% pull(Genus)

bar_df <- diff_res %>%
  filter(Genus %in% top_6_genera) %>%
  mutate(Direction = ifelse(log2fc > 0, "Enriched in SLE", "Enriched in HC"),
         Genus = factor(Genus, levels = Genus[order(log2fc)]))

p_figS2b <- ggplot(bar_df, aes(Genus, log2fc, fill = Direction)) +
  geom_col(width = 0.72) +
  coord_flip() +
  scale_fill_manual(values = c("Enriched in HC" = "#7FA0D0", "Enriched in SLE" = "#C94443")) +
  labs(x = NULL, y = "log2fc", title = "b") +
  theme_clean(11)

box_df <- genus[top_6_genera, , drop = FALSE] %>%
  rownames_to_column("Genus") %>%
  pivot_longer(-Genus, names_to = "SampleID", values_to = "Abundance") %>%
  left_join(metadata %>% rownames_to_column("SampleID"), by = "SampleID") %>%
  mutate(Group = factor(Group, levels = c("HC", "SLE")))

p_figS2c <- ggplot(box_df, aes(Group, Abundance)) +
  geom_jitter(aes(color = Group), width = 0.18, alpha = 0.7, size = 1.8) +
  geom_boxplot(aes(fill = Group, color = Group), alpha = 0.35, width = 0.56, outlier.shape = NA) +
  stat_compare_means(method = "wilcox.test", label.y.npc = 0.95, fontface = "bold.italic") +
  scale_color_manual(values = pal_group[c("HC", "SLE")]) +
  scale_fill_manual(values = pal_group[c("HC", "SLE")]) +
  facet_wrap(~Genus, scales = "free_y", nrow = 2) +
  labs(x = NULL, y = "Relative abundance (14%)", title = "c") +
  theme_clean(10) +
  theme(legend.position = "none")

figS2 <- (p_figS2a | p_figS2b) / p_figS2c
ggsave("figures/supplementary/FigS2.png", figS2, width = 16, height = 10, dpi = 300)
ggsave("figures/supplementary/FigS2.pdf", figS2, width = 16, height = 10, device = cairo_pdf)

# ============================================================
# FigS3. Biofilm-related module analysis
# ============================================================

prot <- fread("dia-proteinSummary-all.csv")
anno <- fread("annotation_allprotein.csv")
group <- fread("group.csv")

margin <- ggplot2::margin
theme_nature <- theme_classic(base_size = 12) +
  theme(
    text = element_text(family = "sans", face = "bold", color = "black"),
    plot.title = element_text(face = "bold", size = 16, margin = margin(b = 10)),
    plot.subtitle = element_text(face = "bold", size = 12, color = "grey30"),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11, face = "bold", color = "black"),
    axis.line = element_line(linewidth = 0.8, color = "black"),
    axis.ticks = element_line(linewidth = 0.8, color = "black"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    legend.position = "top",
    strip.background = element_rect(fill = NA, color = NA),
    strip.text = element_text(size = 13, face = "bold", color = "black"),
    plot.margin = margin(15, 15, 15, 15)
  )

anno2 <- anno %>%
  rowwise() %>%
  mutate(
    anno_text = str_to_lower(
      paste(
        c_across(any_of(c(
          "Swissprot_description",
          "NR_subject",
          "NOG_orthologous_id_description",
          "COG_function_description"
        ))),
        collapse = " "
      )
    )
  ) %>%
  ungroup()

kw <- list(
  cdiGMP = c("ggdef","eal","hd-gyp","diguanylate cyclase","phosphodiesterase","pilz"),
  biofilm_EPS = c("biofilm","exopolysaccharide","eps","capsule","wza","wzb","wzc","wzx","wzy","wzz","pnag","pgaa","pgab","pgac","pgad","cellulose synthase","bcsa","bcsb","bcsz"),
  adhesion = c("adhesin","fimbr","pili","pilin","usher","autotransporter","surface protein","lpxtg"),
  motility = c("flagell","motility","chemotaxis","fli","flg","mota","motb","chea","chey")
)

biofilm_all <- purrr::imap_dfr(
  kw,
  function(words, module) {
    pat <- paste(words, collapse = "|")
    anno2 %>%
      filter(str_detect(anno_text, pat)) %>%
      mutate(module = module)
  }
) %>%
  distinct(Protein_id, module, .keep_all = TRUE)

prot <- prot %>% rename(Protein_id = colnames(prot)[1])

group_col <- colnames(group)[sapply(group, function(x) any(x %in% c("NC", "SLE")))]
sample_col <- setdiff(colnames(group), group_col)

group <- group %>%
  rename(sample_id = !!sample_col, group = !!group_col) %>%
  mutate(group = ifelse(group == "NC", "HC", group))

prot_long <- prot %>%
  pivot_longer(-Protein_id, names_to = "sample_id", values_to = "abundance") %>%
  left_join(group, by = "sample_id") %>%
  filter(!is.na(group))

diff_all <- prot_long %>%
  semi_join(biofilm_all %>% select(Protein_id), by = "Protein_id") %>%
  group_by(Protein_id) %>%
  summarise(
    mean_SLE = mean(abundance[group == "SLE"], na.rm = TRUE),
    mean_HC  = mean(abundance[group == "HC"],  na.rm = TRUE),
    log2FC = ifelse(is.na(mean_SLE) | is.na(mean_HC), NA, log2((mean_SLE + 1e-9)/(mean_HC + 1e-9))),
    p_value = tryCatch(wilcox.test(abundance[group == "SLE"], abundance[group == "HC"])$p.value,
                       error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    FDR = p.adjust(p_value, "BH"),
    direction = case_when(
      log2FC > 0 ~ "Up_in_SLE",
      log2FC < 0 ~ "Down_in_SLE",
      TRUE ~ NA_character_
    )
  ) %>%
  left_join(biofilm_all %>% select(Protein_id, module), by = "Protein_id", relationship = "many-to-many")

abundance_data <- prot_long %>%
  semi_join(biofilm_all, by = "Protein_id") %>%
  left_join(biofilm_all %>% select(Protein_id, module), by = "Protein_id", relationship = "many-to-many") %>%
  mutate(
    module = factor(module, levels = c("cdiGMP","biofilm_EPS","adhesion","motility")),
    group = factor(group, levels = c("HC","SLE"))
  )

stat_df <- abundance_data %>%
  group_by(module) %>%
  summarise(
    max_y = max(abundance, na.rm = TRUE),
    p_val = tryCatch(wilcox.test(abundance ~ group)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    p_label = case_when(
      p_val < 0.001 ~ "***",
      p_val < 0.01  ~ "**",
      p_val < 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    y_position = max_y * 1.05,
    y_text = max_y * 1.12
  )

panel_a <- biofilm_all %>%
  count(module) %>%
  mutate(module = factor(module, levels = names(pal_module))) %>%
  ggplot(aes(module, n, fill = module)) +
  geom_col(width = 0.7, color = "black", linewidth = 1) +
  geom_text(aes(label = n), vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = pal_module) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "a", x = NULL, y = "Count") +
  theme_nature +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

panel_b <- diff_all %>%
  ggplot(aes(log2FC, -log10(p_value))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, color = "grey60") +
  geom_point(data = diff_all %>% filter(p_value >= 0.05), color = "grey80", size = 2, alpha = 0.5) +
  geom_point(data = diff_all %>% filter(p_value < 0.05), aes(fill = module), shape = 21, size = 3, color = "black") +
  scale_fill_manual(values = pal_module) +
  labs(title = "b", x = expression(log[2]*"FC"), y = expression(-log[10]*"P")) +
  theme_nature +
  theme(legend.position = "right")

panel_c <- diff_all %>%
  filter(!is.na(direction)) %>%
  count(module, direction) %>%
  group_by(module) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup() %>%
  ggplot(aes(module, pct, fill = direction)) +
  geom_col(color = "black") +
  geom_text(aes(label = paste0(round(pct), "%")), position = position_stack(vjust = 0.5),
            color = "white", size = 4) +
  scale_fill_manual(values = pal_direction) +
  labs(title = "c", x = NULL, y = "Percentage") +
  theme_nature

panel_d <- ggplot(abundance_data, aes(group, abundance)) +
  geom_jitter(aes(color = group), width = 0.2, size = 1.5, alpha = 0.6) +
  geom_boxplot(width = 0.5, fill = NA, color = "black", outlier.shape = NA) +
  geom_segment(data = stat_df, aes(x = 1, xend = 2, y = y_position, yend = y_position), inherit.aes = FALSE) +
  geom_text(data = stat_df, aes(x = 1.5, y = y_text, label = p_label), inherit.aes = FALSE, size = 6) +
  facet_wrap(~module, scales = "free_y", nrow = 1) +
  scale_color_manual(values = pal_group[c("HC","SLE")]) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  labs(title = "d", x = NULL, y = "Intensity") +
  theme_nature +
  theme(legend.position = "none")

figS3 <- (panel_a | panel_b) / (panel_c | panel_d)
ggsave("figures/supplementary/FigS3.png", figS3, width = 12, height = 8, dpi = 300)
ggsave("figures/supplementary/FigS3.pdf", figS3, width = 12, height = 8, device = cairo_pdf)

# ============================================================
# FigS4. Repeated cross-validation diagnostics
# ============================================================

common_samples <- intersect(colnames(ko), rownames(metadata))
ko_clean <- ko[, common_samples, drop = FALSE]
ko_clean <- ko_clean[apply(ko_clean, 1, var, na.rm = TRUE) > 0, , drop = FALSE]

diff_genus <- apply(genus, 1, function(x) {
  g1 <- x[metadata$Group == "HC"]
  g2 <- x[metadata$Group == "SLE"]
  pval <- tryCatch(wilcox.test(g1, g2)$p.value, error = function(e) 1)
  log2fc <- log2((mean(g2, na.rm = TRUE) + 1e-9) / (mean(g1, na.rm = TRUE) + 1e-9))
  c(p.value = pval, log2fc = log2fc)
}) %>% t() %>% as.data.frame()
diff_genus$Genus <- rownames(diff_genus)
top_6_genera <- diff_genus %>% arrange(p.value) %>% slice_head(n = 6) %>% pull(Genus)

res <- apply(ko_clean, 1, function(x) {
  g1 <- x[metadata[common_samples, "Group"] == "HC"]
  g2 <- x[metadata[common_samples, "Group"] == "SLE"]
  pval <- tryCatch(wilcox.test(g1, g2)$p.value, error = function(e) 1)
  log2fc <- log2((mean(g2, na.rm = TRUE) + 1e-9) / (mean(g1, na.rm = TRUE) + 1e-9))
  c(p.value = pval, log2FC = log2fc)
}) %>% t() %>% as.data.frame()
res$FDR <- p.adjust(res$p.value, "fdr")
sig_genes_list <- rownames(res[abs(res$log2FC) > 0.58 & res$FDR < 0.05, , drop = FALSE])

rf_input <- cbind(
  t(genus[top_6_genera, common_samples, drop = FALSE]),
  t(ko[sig_genes_list, common_samples, drop = FALSE])
) %>% as.data.frame()
rf_input$Group <- factor(metadata[common_samples, "Group"], levels = c("HC", "SLE"))
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
    rf_fold <- randomForest(Group ~ ., data = data[train_idx, ], ntree = 1000,
                            mtry = floor(sqrt(ncol(data) - 1)))
    oof_prob[test_idx] <- predict(rf_fold, data[test_idx, ], type = "prob")[, "SLE"]
  }
  roc_obj <- suppressMessages(pROC::roc(y, oof_prob, levels = c("HC", "SLE"), direction = "<"))
  list(auc = as.numeric(pROC::auc(roc_obj)), prob = oof_prob)
}

cv_runs <- lapply(seq_len(100), function(i) one_cv_run(rf_input, seed = i, k = 5))
auc_vec <- vapply(cv_runs, `[[`, numeric(1), "auc")
auc_mean <- mean(auc_vec)
auc_ci <- quantile(auc_vec, c(0.025, 0.975))
median_id <- order(abs(auc_vec - median(auc_vec)))[1]

pred_df <- data.frame(
  Sample = seq_len(nrow(rf_input)),
  Probability = cv_runs[[median_id]]$prob,
  Group = rf_input$Group
)

p_figS4a <- ggplot(pred_df, aes(Sample, Probability, color = Group)) +
  geom_point(size = 2.4, alpha = 0.9) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  scale_color_manual(values = pal_group[c("HC", "SLE")]) +
  labs(x = "Sample", y = "Predicted probability of SLE", title = "a") +
  theme_clean(11)

p_figS4b <- ggplot(data.frame(AUC = auc_vec), aes(AUC)) +
  geom_histogram(bins = 20, fill = "#CF5C5C", color = "white", alpha = 0.9) +
  geom_vline(xintercept = auc_mean, linewidth = 0.9) +
  geom_vline(xintercept = auc_ci, linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = auc_mean, y = Inf,
           label = sprintf("Mean AUC = %.3f\n95%% CI: %.3f - %.3f", auc_mean, auc_ci[1], auc_ci[2]),
           vjust = 1.6, size = 4.2, fontface = "bold.italic") +
  labs(x = "Out-of-fold AUC (100 × 5-fold CV)", y = "Count (out of 100 repeats)", title = "b") +
  theme_clean(11)

figS4 <- p_figS4a | p_figS4b
ggsave("figures/supplementary/FigS4.png", figS4, width = 12, height = 5.5, dpi = 300)
ggsave("figures/supplementary/FigS4.pdf", figS4, width = 12, height = 5.5, device = cairo_pdf)

message("Supplementary figures finished: FigS2, FigS3, FigS4. FigS1 code not available in the uploaded scripts.")
