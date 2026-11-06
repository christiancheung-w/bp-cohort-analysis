# 5_2_clinical_bins.R
# SAP §5.2: clinical-threshold proportions.
#
# Per-user 3-bin classification of ΔBP (MO6 − BASELINE), separately for SBP
# and DBP, tabulated by baseline AHA stratum.
#
#   improved:  Δ ≤ −5 mmHg
#   stable:    −5 < Δ < +5 mmHg
#   worsened:  Δ ≥ +5 mmHg
#
# No formal test (per SAP). Wilson 95% CIs added for the proportions because
# they're free and useful for the manuscript text.
#
# Cohort: §2.2 hypertensive (Elevated + Stage 1 + Stage 2). Users without MO6
# for a given outcome are excluded from that outcome's tabulation.
#
# Input:  Build/derived/qc/user_analytic_frame_1X.rds
#
# Outputs (Analyses/5_2/derived/):
#   results/5_2_bin_counts.csv          — long: outcome × stratum × bin → n
#   results/5_2_bin_percentages.csv     — wide: stratum × outcome → n_improved, n_stable, n_worsened + Wilson CIs
#   results/5_2_paper.md                — paper-ready text with percentages
#   results/5_2_clinical_bins.log
#   figures/5_2_figure_1b.png           — faux-3D (perspective) stacked bar chart, 2-panel (SBP, DBP)

library(tidyverse)

project_dir <- "./PRIMARY"
input_path  <- file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")
out_dir     <- file.path(project_dir, "Analyses", "5_2", "derived")
results_dir <- file.path(out_dir, "results")
figures_dir <- file.path(out_dir, "figures")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

log_path <- file.path(results_dir, "5_2_clinical_bins.log")
sink(log_path, split = TRUE)
cat("§5.2 clinical-threshold proportions —", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Read frame and filter to §2.2 cohort ---
user_frame <- read_rds(input_path)
cat("User frame:", nrow(user_frame), "rows ×", ncol(user_frame), "cols\n")

cohort <- user_frame %>%
  filter(BASELINE_AHA %in% c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"))
cat("§2.2 cohort:", nrow(cohort), "users\n\n")

# --- Derive deltas + bins (re-derived for safety; not relying on stored columns) ---
binned <- cohort %>%
  mutate(
    DELTA_SBP_MO6 = MO6_SBP - BASELINE_SBP,
    DELTA_DBP_MO6 = MO6_DBP - BASELINE_DBP,
    bin_SBP = case_when(
      is.na(DELTA_SBP_MO6) ~ NA_character_,
      DELTA_SBP_MO6 <= -5  ~ "improved",
      DELTA_SBP_MO6 >=  5  ~ "worsened",
      TRUE                 ~ "stable"
    ),
    bin_DBP = case_when(
      is.na(DELTA_DBP_MO6) ~ NA_character_,
      DELTA_DBP_MO6 <= -5  ~ "improved",
      DELTA_DBP_MO6 >=  5  ~ "worsened",
      TRUE                 ~ "stable"
    ),
    stratum = factor(BASELINE_AHA,
                     levels = c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
                     labels = c("Elevated", "Stage 1", "Stage 2"))
  )

# --- Long-form bin counts (one row per outcome × stratum × bin) ---
bin_long <- bind_rows(
  binned %>% filter(!is.na(bin_SBP)) %>%
    count(stratum, bin = bin_SBP) %>% mutate(outcome = "SBP"),
  binned %>% filter(!is.na(bin_DBP)) %>%
    count(stratum, bin = bin_DBP) %>% mutate(outcome = "DBP")
) %>%
  mutate(bin = factor(bin, levels = c("improved", "stable", "worsened"))) %>%
  group_by(outcome, stratum) %>%
  mutate(n_total = sum(n),
         prop    = n / n_total) %>%
  ungroup() %>%
  arrange(outcome, stratum, bin) %>%
  select(outcome, stratum, bin, n, n_total, prop)

# Wilson 95% CI for binomial proportion
wilson_ci <- function(x, n, z = 1.959964) {
  p <- x / n
  denom  <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denom
  half   <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  list(lcl = center - half, ucl = center + half)
}

bin_long <- bin_long %>%
  rowwise() %>%
  mutate(
    ci = list(wilson_ci(n, n_total)),
    lcl = ci$lcl,
    ucl = ci$ucl
  ) %>%
  ungroup() %>%
  select(-ci)

write_csv(bin_long, file.path(results_dir, "5_2_bin_counts.csv"))

# --- Wide-format percentages for the manuscript table ---
bin_wide <- bin_long %>%
  mutate(
    cell = sprintf("%d / %d (%.1f%% [%.1f, %.1f])",
                   n, n_total, 100*prop, 100*lcl, 100*ucl)
  ) %>%
  select(outcome, stratum, bin, cell) %>%
  pivot_wider(names_from = bin, values_from = cell) %>%
  arrange(outcome, stratum)

write_csv(bin_wide, file.path(results_dir, "5_2_bin_percentages.csv"))

cat("\n--- Bin tabulation ---\n")
print(bin_long, n = Inf, width = Inf)
cat("\n--- Wide percentages (for paper text) ---\n")
print(bin_wide, n = Inf, width = Inf)

# --- Figure 1B: faux-3D (perspective) stacked bar chart ---
# NOTE: the 3D perspective distorts the very proportions this chart conveys; it
# is the deliberately chosen presentation figure. Depth is faked with offset
# polygon faces: front (base colour) + darker right side + lighter top.
bin_levels <- c("worsened", "stable", "improved")   # bottom-to-top stacking
base_cols  <- c("improved" = "#478A43",  # muted forest green — favorable
                "stable"   = "#bdbdbd",  # gray — neutral
                "worsened" = "#C0392B")  # muted red — adverse

# Shade a hex toward white (>1) or black (<1) for the lit top / shadowed side.
shade <- function(hex, f) {
  v <- pmin(pmax(col2rgb(hex)[, 1] / 255 * f, 0), 1)
  grDevices::rgb(v[1], v[2], v[3])
}
top_col  <- vapply(base_cols, shade, character(1), f = 1.15)  # lighter (lit from above)
side_col <- vapply(base_cols, shade, character(1), f = 0.70)  # darker  (in shadow)

# Geometry of the fake-3D box.
hw <- 0.31     # front-face half-width
dx <- 0.16     # depth offset in x
dy <- 0.035    # depth offset in y (proportion units)

bars <- bin_long %>%
  mutate(
    outcome = factor(outcome, levels = c("SBP", "DBP")),
    stratum = factor(stratum, levels = c("Elevated", "Stage 1", "Stage 2")),
    bin     = factor(bin, levels = bin_levels)
  ) %>%
  arrange(outcome, stratum, bin) %>%
  group_by(outcome, stratum) %>%
  mutate(x = as.integer(stratum),
         ymax = cumsum(prop),
         ymin = ymax - prop) %>%
  ungroup()

flat <- bars %>%
  transmute(outcome = as.character(outcome), x, ymin, ymax, prop,
            bin = as.character(bin))

# Right-side faces (one per segment): darker shade.
side <- pmap_dfr(flat, function(outcome, x, ymin, ymax, prop, bin) {
  tibble(outcome = outcome, id = paste(outcome, x, bin, "s"),
         fill = unname(side_col[bin]),
         px = c(x + hw, x + hw, x + hw + dx, x + hw + dx),
         py = c(ymin,   ymax,   ymax + dy,   ymin + dy))
})

# Top faces (only the topmost segment of each bar): lighter shade.
topdf <- bars %>%
  group_by(outcome, stratum) %>% slice_max(ymax, n = 1, with_ties = FALSE) %>% ungroup() %>%
  transmute(outcome = as.character(outcome), x, ymax, bin = as.character(bin))
topf <- pmap_dfr(topdf, function(outcome, x, ymax, bin) {
  tibble(outcome = outcome, id = paste(outcome, x, "t"),
         fill = unname(top_col[bin]),
         px = c(x - hw, x + hw, x + hw + dx, x - hw + dx),
         py = c(ymax,   ymax,   ymax + dy,   ymax + dy))
})

# Front faces: base colour.
front <- flat %>%
  mutate(fill = unname(base_cols[bin]), xmin = x - hw, xmax = x + hw)

# Restore facet order on the polygon frames.
fac <- function(d) mutate(d, outcome = factor(outcome, levels = c("SBP", "DBP")))
side <- fac(side); topf <- fac(topf); front <- fac(front)

p <- ggplot() +
  geom_polygon(data = side, aes(px, py, group = id, fill = fill)) +
  geom_polygon(data = topf, aes(px, py, group = id, fill = fill)) +
  geom_rect(data = front, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            color = NA) +
  geom_text(data = filter(front, prop > 0.02),
            aes(x = x, y = (ymin + ymax) / 2, label = sprintf("%.1f%%", 100 * prop)),
            color = "black", fontface = "bold", size = 3.3) +
  facet_wrap(~ outcome, ncol = 2,
             labeller = as_labeller(c("SBP" = "Systolic Blood Pressure",
                                      "DBP" = "Diastolic Blood Pressure"))) +
  scale_x_continuous(breaks = 1:3, labels = c("Elevated", "Stage 1", "Stage 2")) +
  scale_y_continuous(labels = scales::label_number(scale = 100),
                     expand = expansion(mult = c(0, 0.06))) +
  scale_fill_identity(name = "Δ Blood Pressure", guide = "legend",
                      breaks = unname(base_cols[c("improved", "stable", "worsened")]),
                      labels = c("Improved (Δ ≤ −5 mmHg)",
                                 "Stable (−5 < Δ < +5)",
                                 "Worsened (Δ ≥ +5)")) +
  labs(x = "Baseline AHA Stratum", y = "Proportion of Stratum (%)",
       title = "6-Month Blood Pressure Change by Baseline AHA Stratum") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid       = element_blank(),
    panel.spacing.x  = grid::unit(2.2, "lines"),
    axis.line        = element_line(color = "black", linewidth = 0.4),
    axis.ticks       = element_line(color = "black", linewidth = 0.4),
    axis.text.x      = element_text(face = "bold"),
    axis.title.x     = element_text(face = "bold", size = 12, margin = margin(t = 16)),
    axis.title.y     = element_text(face = "bold", size = 12),
    strip.text       = element_text(face = "bold", size = 12),
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold", size = 12),
    legend.text      = element_text(size = 11),
    plot.title       = element_text(face = "bold", size = 12.5)
  )

out_png <- file.path(figures_dir, "5_2_figure_1b.png")
ggsave(out_png, p, width = 9, height = 6, dpi = 300)
cat("\nWrote figure (faux-3D):\n  ", out_png, "\n")

# --- Paper-ready markdown text ---
narrative_row <- function(outcome_label, stratum_label) {
  r <- bin_long %>%
    filter(outcome == outcome_label, stratum == stratum_label)
  if (nrow(r) == 0) return("")
  improved <- r %>% filter(bin == "improved")
  stable   <- r %>% filter(bin == "stable")
  worsened <- r %>% filter(bin == "worsened")
  sprintf(
    "%s × %s: improved %.1f%% (n = %d), stable %.1f%% (n = %d), worsened %.1f%% (n = %d). Total %s = %d.",
    outcome_label, stratum_label,
    100*improved$prop, improved$n,
    100*stable$prop,   stable$n,
    100*worsened$prop, worsened$n,
    "n", improved$n_total
  )
}

md_lines <- c(
  "# §5.2 — Clinical-threshold bin proportions",
  "",
  paste0("_Generated ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "_"),
  "",
  "## Per-stratum bin percentages (paper-quote-ready)",
  ""
)
for (otc in c("SBP", "DBP")) {
  md_lines <- c(md_lines, paste0("### ", otc))
  for (st in c("Elevated", "Stage 1", "Stage 2")) {
    md_lines <- c(md_lines, paste0("- ", narrative_row(otc, st)))
  }
  md_lines <- c(md_lines, "")
}
md_lines <- c(md_lines,
  "## Notes",
  "",
  "- Bin thresholds at ±5 mmHg are the conventional clinical-meaningfulness cutoff for BP change (per SAP §5.2).",
  "- Proportions are within-stratum (rows sum to 100% within each outcome × stratum cell).",
  "- Per SAP §5.2, no formal hypothesis test is performed. Wilson 95% CIs are reported for transparency.",
  "- Users with missing MO6 BP for the relevant outcome are excluded from that outcome's tabulation; per-stratum n shown above reflects the analyzable subset."
)
writeLines(md_lines, file.path(results_dir, "5_2_paper.md"))
cat("\nWrote markdown:\n  ", file.path(results_dir, "5_2_paper.md"), "\n")

cat("\n=== §5.2 clinical bins done ===\n")
sink()
