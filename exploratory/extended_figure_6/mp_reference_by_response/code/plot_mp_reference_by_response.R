suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(grid)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- args_all[grepl("^--file=", args_all)]
script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else getwd()
script_dir <- dirname(normalizePath(script_path))
project_root <- normalizePath(file.path(script_dir, "..", "..", "..", ".."))

source_path <- file.path(project_root, "final", "extended_figure_6", "mp7_nonccr_association_candidates", "outputs", "extended_figure6_mp_abundance_by_sample.csv")
out_dir <- file.path(project_root, "exploratory", "extended_figure_6", "mp_reference_by_response", "outputs")
legend_dir <- file.path(project_root, "exploratory", "extended_figure_6", "mp_reference_by_response", "legends")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(legend_dir, recursive = TRUE, showWarnings = FALSE)

mp_levels <- paste0("MP_", 1:7)
mp_palette <- c(
  MP_1 = "#1F77B4",
  MP_2 = "#6FA4C8",
  MP_3 = "#E1812C",
  MP_4 = "#F0C08A",
  MP_5 = "#2CA02C",
  MP_6 = "#A1D99B",
  MP_7 = "#C52B2F"
)
outcome_cols <- c("cCR" = "#7BB6A4", "non-cCR" = "#F2A38A")

theme_nc <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      axis.line = element_line(linewidth = 0.5, colour = "black"),
      axis.ticks = element_line(linewidth = 0.5, colour = "black"),
      axis.text = element_text(colour = "black", size = 6),
      axis.title = element_text(colour = "black", size = 7, face = "bold"),
      plot.title = element_text(size = 7, face = "bold", hjust = 0),
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 7, face = "bold"),
      plot.margin = margin(4, 4, 4, 4, "pt")
    )
}

mm_pdf <- function(filename, width_mm, height_mm) {
  grDevices::pdf(
    file.path(out_dir, filename),
    width = width_mm / 25.4,
    height = height_mm / 25.4,
    useDingbats = FALSE,
    family = "Helvetica"
  )
}

dt <- fread(source_path)
dt[, MP := factor(MP, levels = mp_levels)]
dt[, cCR := factor(cCR, levels = c("cCR", "non-cCR"))]

stats <- rbindlist(lapply(mp_levels, function(mp) {
  sub <- dt[MP == mp]
  p <- tryCatch(wilcox.test(mp_percent ~ cCR, data = sub, exact = FALSE)$p.value, error = function(e) NA_real_)
  data.table(
    MP = mp,
    n_cCR = sum(sub$cCR == "cCR", na.rm = TRUE),
    n_non_cCR = sum(sub$cCR == "non-cCR", na.rm = TRUE),
    median_cCR = median(sub[cCR == "cCR", mp_percent], na.rm = TRUE),
    median_non_cCR = median(sub[cCR == "non-cCR", mp_percent], na.rm = TRUE),
    p_value = p
  )
}))
stats[, p_adjusted_BH := p.adjust(p_value, method = "BH")]
fwrite(stats, file.path(out_dir, "hard_mp_fraction_response_reference_stats.csv"))

all_mp <- ggplot(dt, aes(cCR, mp_percent, colour = cCR)) +
  geom_boxplot(width = 0.52, outlier.shape = NA, linewidth = 0.5, colour = "#7FAEC8", fill = "white") +
  geom_point(size = 0.65, alpha = 1, position = position_jitter(width = 0.08, height = 0), show.legend = FALSE) +
  facet_wrap(~ MP, nrow = 2, scales = "free_y") +
  scale_colour_manual(values = outcome_cols) +
  labs(
    title = "Hard MP fractions by clinical outcome",
    x = NULL,
    y = "Cells among malignant MPs (%)"
  ) +
  theme_nc()
mm_pdf("reference_all_hard_mp_fraction_by_response.pdf", 120, 78)
print(all_mp)
dev.off()

composition <- copy(dt)
sample_order <- dt[MP == "MP_7"][order(cCR, mp_percent), sampleID]
composition[, sampleID := factor(sampleID, levels = sample_order)]
composition[, MP := factor(MP, levels = rev(mp_levels))]
stacked <- ggplot(composition, aes(sampleID, mp_percent, fill = MP)) +
  geom_col(width = 0.74, colour = NA) +
  scale_fill_manual(values = mp_palette, name = "Hard MP") +
  labs(
    title = "Hard MP composition by sample",
    x = NULL,
    y = "Cells among malignant MPs (%)"
  ) +
  theme_nc() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
mm_pdf("reference_hard_mp_composition_by_sample.pdf", 120, 62)
print(stacked)
dev.off()

legend_text <- c(
  "Reference hard MP assignment plots",
  "",
  "These plots use hard MP assignment rather than continuous MP scores. They are provided as reference only because continuous scores for MP1-MP6 are not present in the local exported score table."
)
writeLines(legend_text, file.path(legend_dir, "legend_draft.md"))

message("Saved hard MP reference plots to: ", out_dir)
