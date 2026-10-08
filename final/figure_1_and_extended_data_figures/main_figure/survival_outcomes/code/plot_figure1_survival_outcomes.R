suppressPackageStartupMessages({
  library(readxl)
  library(survival)
  library(ggplot2)
  library(grid)
  library(gridExtra)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg))
part_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
project_root <- normalizePath(file.path(part_dir, "..", "..", "..", ".."), mustWork = TRUE)
data_path <- file.path(project_root, "data", "_updated_survival", "20261008TENERGY初発症例データ一覧.xlsx")
endpoint_arg <- grep("^--endpoint=", commandArgs(trailingOnly = TRUE), value = TRUE)
endpoint_selection <- if (length(endpoint_arg)) sub("^--endpoint=", "", endpoint_arg[[1]]) else "all"
stopifnot(endpoint_selection %in% c("all", "pfs", "os"))
output_dir <- file.path(part_dir, "outputs")
legend_dir <- file.path(part_dir, "legends")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(legend_dir, recursive = TRUE, showWarnings = FALSE)

# Build grobs on Quartz so Arial metrics match the final PDF device.
metrics_file <- tempfile(fileext = ".pdf")
grDevices::quartz(type = "pdf", file = metrics_file, family = "Arial")

# Nature Cancer final-artwork minimum: 0.5 pt lines. ggplot2 line widths are mm.
line_width <- 0.5 / ggplot2::.pt
clinical_colors <- c(cCR = "#7BB6A4", `non-cCR` = "#F2A38A")

dat <- read_excel(data_path, sheet = "TENERGY長期追跡")
stopifnot(nrow(dat) == 40, !anyDuplicated(dat$caseno),
          all(dat$CR %in% c(0, 1)), all(dat$event_pfs_obs %in% c(0, 1)),
          all(is.finite(dat$leng_pfs_obsm)), all(dat$leng_pfs_obsm >= 0))
dat$CR_group <- factor(ifelse(dat$CR == 1, "cCR", "non-cCR"),
                       levels = c("cCR", "non-cCR"))

format_p <- function(p) {
  if (p < 0.001) "P < 0.001" else paste0("P = ", formatC(p, format = "f", digits = 3))
}

build_km_panel <- function(data, time_col, event_col, endpoint, by_response = FALSE) {
  if (by_response) {
    fit_formula <- as.formula(paste0("Surv(`", time_col, "`, `", event_col, "`) ~ CR_group"))
    fit <- survfit(fit_formula, data = data, conf.type = "log-log")
    lr <- survdiff(fit_formula, data = data)
    p_value <- pchisq(lr$chisq, df = length(lr$n) - 1, lower.tail = FALSE)
    cox_fit <- coxph(fit_formula, data = data)
    cox_summary <- summary(cox_fit)
    hr <- unname(cox_summary$coefficients["CR_groupnon-cCR", "exp(coef)"])
    ci_low <- unname(cox_summary$conf.int["CR_groupnon-cCR", "lower .95"])
    ci_high <- unname(cox_summary$conf.int["CR_groupnon-cCR", "upper .95"])
    group_levels <- c("cCR", "non-cCR")
    group_colors <- clinical_colors
    group_title <- "cCR status"
  } else {
    fit_formula <- as.formula(paste0("Surv(`", time_col, "`, `", event_col, "`) ~ 1"))
    fit <- survfit(fit_formula, data = data, conf.type = "log-log")
    p_value <- hr <- ci_low <- ci_high <- NA_real_
    group_levels <- "All patients"
    group_colors <- c(`All patients` = "#000000")
    group_title <- "All patients"
  }

  # censored = TRUE is essential: it retains censor-only times after the final event.
  curve_summary <- summary(fit, censored = TRUE)
  if (by_response) {
    curve_df <- data.frame(
      time = curve_summary$time,
      survival = curve_summary$surv,
      group = sub("CR_group=", "", curve_summary$strata),
      n_censor = curve_summary$n.censor
    )
  } else {
    curve_df <- data.frame(
      time = curve_summary$time,
      survival = curve_summary$surv,
      group = "All patients",
      n_censor = curve_summary$n.censor
    )
  }
  curve_df$group <- factor(curve_df$group, levels = group_levels)
  baseline_df <- data.frame(
    time = 0, survival = 1,
    group = factor(group_levels, levels = group_levels),
    n_censor = 0
  )
  curve_df <- rbind(baseline_df, curve_df)

  risk_times <- c(0, 12, 24, 36, 48, 60)
  risk_summary <- summary(fit, times = risk_times, extend = TRUE)
  if (by_response) {
    risk_df <- data.frame(
      time = risk_summary$time,
      n_risk = risk_summary$n.risk,
      group = sub("CR_group=", "", as.character(risk_summary$strata))
    )
  } else {
    risk_df <- data.frame(
      time = risk_summary$time,
      n_risk = risk_summary$n.risk,
      group = "All patients"
    )
  }
  risk_df$group <- factor(risk_df$group, levels = group_levels)

  statistic_annotations <- NULL
  if (by_response) {
    statistic_annotations <- data.frame(
      x = 2,
      y = c(0.16, 0.09),
      label = c(
        paste0("log-rank ", format_p(p_value)),
        paste0("HR (fixed group) = ", formatC(hr, format = "f", digits = 2),
               " (95% CI ", formatC(ci_low, format = "f", digits = 2),
               "-", formatC(ci_high, format = "f", digits = 2), ")")
      )
    )
  }

  main_plot <- ggplot(curve_df, aes(time, survival, colour = group, group = group)) +
    geom_step(linewidth = line_width, direction = "hv") +
    geom_point(
      data = subset(curve_df, n_censor > 0),
      shape = 3, size = 1.5, stroke = line_width
    ) +
    scale_colour_manual(values = group_colors, name = NULL) +
    scale_x_continuous(
      limits = c(0, 60), breaks = risk_times,
      expand = expansion(mult = c(0, 0.025))
    ) +
    scale_y_continuous(
      labels = function(x) paste0(round(100 * x), "%"),
      breaks = seq(0, 1, 0.2), limits = c(0, 1.08), expand = c(0, 0)
    ) +
    labs(title = group_title, x = "Months after registration", y = endpoint) +
    theme_classic(base_family = "Arial", base_size = 7) +
    theme(
      axis.line = element_line(linewidth = line_width, colour = "black"),
      axis.ticks = element_line(linewidth = line_width, colour = "black"),
      axis.ticks.length = unit(1.5, "mm"),
      axis.text = element_text(size = 7, colour = "black"),
      axis.title = element_text(size = 7, colour = "black"),
      plot.title = element_text(size = 7, face = "plain", hjust = 0.5),
      legend.position = if (by_response) c(0.78, 0.84) else "none",
      legend.text = element_text(size = 7),
      legend.key.width = unit(6, "mm"),
      legend.key.height = unit(3, "mm"),
      plot.margin = margin(3, 4, 1, 4, "mm")
    )
  if (by_response) {
    main_plot <- main_plot + geom_text(
      data = statistic_annotations,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE, family = "Arial", size = 2.45, hjust = 0
    )
  }

  risk_plot <- ggplot(risk_df, aes(time, group, label = n_risk, colour = group)) +
    geom_text(size = 2.25, family = "Arial", show.legend = FALSE) +
    scale_colour_manual(values = group_colors, guide = "none") +
    scale_x_continuous(
      limits = c(0, 60), breaks = risk_times, labels = risk_times,
      expand = expansion(mult = c(0, 0.025))
    ) +
    scale_y_discrete(limits = rev(group_levels), expand = expansion(add = 0.38)) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_classic(base_family = "Arial", base_size = 7) +
    theme(
      axis.line = element_blank(), axis.ticks = element_blank(),
      axis.text.x = element_text(size = 6, colour = "black"),
      axis.text.y = element_text(size = 6, colour = "black", margin = margin(r = 8)),
      plot.margin = margin(0, 4, 1, 4, "mm")
    )

  main_grob <- ggplotGrob(main_plot)
  risk_grob <- ggplotGrob(risk_plot)
  aligned_widths <- unit.pmax(main_grob$widths, risk_grob$widths)
  main_grob$widths <- risk_grob$widths <- aligned_widths
  list(
    grob = arrangeGrob(main_grob, risk_grob, ncol = 1, heights = c(4.6, 1)),
    p_value = p_value, hr = hr, ci_low = ci_low, ci_high = ci_high,
    fit = fit, curves = curve_df, risk = risk_df
  )
}

# The macOS Quartz device registers the installed Arial TrueType font reliably.
save_arial_pdf <- function(grob, filename, width = 6.7, height = 2.6) {
  grDevices::quartz(type = "pdf", file = filename, width = width, height = height,
                    family = "Arial")
  grid.newpage()
  grid.draw(grob)
  grDevices::dev.off()
}

if (endpoint_selection %in% c("all", "os")) {
  os_all <- build_km_panel(dat, "leng_os_obsm", "event_os_obs", "Overall survival", by_response = FALSE)
  os_cr <- build_km_panel(dat, "leng_os_obsm", "event_os_obs", "Overall survival", by_response = TRUE)
  save_arial_pdf(arrangeGrob(os_all$grob, os_cr$grob, ncol = 2),
                 file.path(output_dir, "figure1b_overall_survival.pdf"))
  message("OS log-rank ", format_p(os_cr$p_value), "; HR = ", formatC(os_cr$hr, format = "f", digits = 4))
}
if (endpoint_selection %in% c("all", "pfs")) {
  pfs_all <- build_km_panel(dat, "leng_pfs_obsm", "event_pfs_obs", "Progression-free survival", by_response = FALSE)
  pfs_cr <- build_km_panel(dat, "leng_pfs_obsm", "event_pfs_obs", "Progression-free survival", by_response = TRUE)
  save_arial_pdf(arrangeGrob(pfs_all$grob, pfs_cr$grob, ncol = 2),
                 file.path(output_dir, "figure1c_progression_free_survival.pdf"))
  save_arial_pdf(pfs_all$grob, file.path(output_dir, "figure1c_pfs_all_patients.pdf"), 3.35, 2.6)
  save_arial_pdf(pfs_cr$grob, file.path(output_dir, "figure1c_pfs_by_response.pdf"), 3.35, 2.6)
  write.csv(pfs_all$curves, file.path(output_dir, "pfs_all_curve.csv"), row.names = FALSE)
  write.csv(pfs_cr$curves, file.path(output_dir, "pfs_by_response_curve.csv"), row.names = FALSE)
  write.csv(pfs_cr$risk, file.path(output_dir, "pfs_by_response_at_risk.csv"), row.names = FALSE)
  median_tables <- rbind(
    data.frame(group = "All patients", t(summary(pfs_all$fit)$table), check.names = FALSE),
    data.frame(group = c("cCR", "non-cCR"), summary(pfs_cr$fit)$table, check.names = FALSE)
  )
  write.csv(median_tables, file.path(output_dir, "pfs_km_summary.csv"), row.names = FALSE)
  write.csv(data.frame(source = basename(data_path), model = "Fixed-response-group Cox, non-cCR versus cCR; not time-dependent",
                       logrank_p = pfs_cr$p_value, HR = pfs_cr$hr,
                       CI_low = pfs_cr$ci_low, CI_high = pfs_cr$ci_high),
            file.path(output_dir, "pfs_statistics.csv"), row.names = FALSE)
  message("PFS log-rank p = ", format(pfs_cr$p_value, digits = 16),
          "; fixed-group HR = ", formatC(pfs_cr$hr, format = "f", digits = 4))
}

writeLines(c(
  "# Figure 1b-c survival outcomes", "",
  "**Figure 1b | Overall survival in the updated TENERGY initial cohort.** Kaplan-Meier estimates of overall survival from clinical-trial registration for all patients (left) and by clinical response (right; cCR, n = 16; non-cCR, n = 24).", "",
  "**Figure 1c | Progression-free survival in the updated TENERGY initial cohort.** Kaplan-Meier estimates of progression-free survival from clinical-trial registration for all patients (left) and by clinical response (right; cCR, n = 16; non-cCR, n = 24).", "",
  "Crosses denote censored observations. Numbers below each curve are patients at risk. Between-response comparisons use two-sided log-rank tests; hazard ratios (HRs; non-cCR versus cCR) and 95% confidence intervals are from univariable Cox proportional-hazards models with eventual response treated as a fixed group. These are descriptive post-treatment group comparisons, not time-dependent Cox estimates, and are susceptible to guarantee-time bias. OS uses `leng_os_obsm` and `event_os_obs`; PFS uses `leng_pfs_obsm` and `event_pfs_obs`. Median survival confidence intervals use the log-log transformation.", "",
  "PFS source updated to `20261008TENERGY初発症例データ一覧.xlsx`. OS artwork was not regenerated in the October 8 PFS-only run. The statistical department's time-dependent Cox estimate for PFS is a different analysis: cCR versus non-cCR HR 0.39 (95% CI 0.15-1.05), P = 0.063, as reported in the October 8 presentation; this script does not reconstruct that model without response-onset dates."
), file.path(legend_dir, "figure1b_c_survival_legends.md"))

grDevices::dev.off()
unlink(metrics_file)
