# 8_figure_5_combined.R
# SAP §8 Figure 5 (combined): engagement-metric distributions by BP-responder
# status, SBP and DBP on a SINGLE 7-panel figure using split (half) violins.
#
# Per-panel layout: two split-violins side by side — an SBP violin and a DBP
# violin. Within each, the LEFT half is the responder (SBP-Responder blue /
# DBP-Responder green) and the RIGHT half is the non-responder for that BP
# definition (gray). A slim white box (median + IQR) sits on each half.
#
# Half-violins are built manually from kernel densities. The 7 panels are built as
# separate plots and composed with patchwork on a centered 3-over-4 design, so all
# panels render at equal width (a single facet_manual grid can't center 3 over 4
# without overlapping columns and squishing the rightmost panels).
#
# Standalone figure script: reads the frozen user analytic frame (read-only) and
# writes only the PNG; the §8 generator and all §8 result CSVs stay frozen.
#
# Input:  Build/derived/qc/user_analytic_frame_1X.rds
# Output: Analyses/8_engagement/derived/figures/8_figure_5_violin_combined.png

library(tidyverse)
library(patchwork)

project_dir <- "./PRIMARY"
input_path  <- file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")
figures_dir <- file.path(project_dir, "Analyses", "8_engagement", "derived", "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

cohort <- read_rds(input_path) %>%
  filter(BASELINE_AHA %in% c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
         SEX %in% c("male", "female")) %>%
  mutate(DBP_RESPONDER_6MO = as.integer((MO6_DBP - BASELINE_DBP) <= -5))

engagement_preds_violin <- c(
  "PCT_DAYS_JOURNALED", "LOG_AVG_DAILY_TIME_IN_APP_MIN", "PCT_DAYS_ENGAGED_NOTIFICATIONS",
  "PCT_DAYS_USED_COACH", "AVG_SLEEP_NEED_PCT", "PCT_DAYS_MET_OPTIMAL_STRAIN",
  "LOG_N_BP_READINGS_STUDY")
# Order = facet panel order (A–G) for the design below: top row Sleep / Physical
# Activity / BP Readings; bottom row Daily App Time / Journal Completion /
# Notification / AI Coach (Daily App Time and Journal Completion swapped per request).
pred_labels <- c(
  AVG_SLEEP_NEED_PCT             = "Average Sleep Need Achieved (%)",
  PCT_DAYS_MET_OPTIMAL_STRAIN    = "Physical Activity Recommendation (% of Days)",
  LOG_N_BP_READINGS_STUDY        = "Blood Pressure Readings (log Count)",
  LOG_AVG_DAILY_TIME_IN_APP_MIN  = "Daily App Time (log min)",
  PCT_DAYS_JOURNALED             = "Journal Completion (% of Days)",
  PCT_DAYS_ENGAGED_NOTIFICATIONS = "Notification Engagement (% of Days)",
  PCT_DAYS_USED_COACH            = "AI Coach Engagement (% of Days)")

# Metrics definitionally bounded at 0–100% — y-axis fixed at 0–100.
pct_bounded <- c("PCT_DAYS_JOURNALED", "PCT_DAYS_ENGAGED_NOTIFICATIONS",
                 "PCT_DAYS_USED_COACH", "PCT_DAYS_MET_OPTIMAL_STRAIN")

# --- Long frame: one row per (user, metric) for each BP definition × status ---
long_one <- function(tag) {
  resp_col <- paste0(tag, "_RESPONDER_6MO")
  cohort %>%
    filter(!is.na(.data[[resp_col]])) %>%
    transmute(bp_type = tag,
              status  = ifelse(.data[[resp_col]] == 1L, "Responder", "Non-Responder"),
              across(all_of(engagement_preds_violin))) %>%
    pivot_longer(all_of(engagement_preds_violin), names_to = "predictor", values_to = "value") %>%
    filter(!is.na(value))
}
vd <- bind_rows(long_one("SBP"), long_one("DBP"))


# --- Geometry constants ---
x_of      <- c(SBP = 1, DBP = 2.4)   # the two violin centers within each panel
halfwidth <- 0.40                    # max half-violin extent from center
box_off   <- 0.11                    # box offset from center toward each half
box_w     <- 0.09                    # box width
sign_of   <- c(Responder = -1, `Non-Responder` = 1)  # responder left, non-responder right
fill_of   <- function(bp, st) dplyr::case_when(
  st == "Non-Responder" ~ "Non-Responder",
  bp == "SBP"           ~ "SBP-Responder",
  TRUE                  ~ "DBP-Responder")

# --- Build half-violin polygons + box stats per (metric, bp_type, status) ---
grid <- vd %>% distinct(predictor, bp_type, status)
poly_rows <- list(); box_rows <- list()
for (i in seq_len(nrow(grid))) {
  g  <- grid[i, ]
  v  <- vd$value[vd$predictor == g$predictor & vd$bp_type == g$bp_type & vd$status == g$status]
  v  <- v[is.finite(v)]
  if (length(v) < 5 || sd(v) == 0) next
  x0 <- x_of[[g$bp_type]]; sgn <- sign_of[[g$status]]
  d    <- density(v, n = 256)
  keep <- d$x >= min(v) & d$x <= max(v)
  dx   <- d$x[keep]; dy <- d$y[keep]
  if (length(dx) < 3) { dx <- d$x; dy <- d$y }
  w   <- dy / max(dy) * halfwidth
  n  <- length(dx)
  fillgrp <- fill_of(g$bp_type, g$status)
  poly_rows[[i]] <- tibble(
    predictor = g$predictor, bp_type = g$bp_type, status = g$status,
    poly_id = paste(g$predictor, g$bp_type, g$status, sep = "|"),
    fillgrp = fillgrp,
    x = c(x0 + sgn * w, x0, x0),
    y = c(dx, dx[n], dx[1]))
  qs <- quantile(v, c(0.25, 0.5, 0.75), names = FALSE)
  bx <- x0 + sgn * box_off
  box_rows[[i]] <- tibble(
    predictor = g$predictor, bp_type = g$bp_type, status = g$status,
    xmin = bx - box_w/2, xmax = bx + box_w/2, q1 = qs[1], med = qs[2], q3 = qs[3])
}
poly_df <- bind_rows(poly_rows)
box_df  <- bind_rows(box_rows)

lab_fac <- function(p) factor(pred_labels[p], levels = pred_labels)
poly_df$predictor_label <- lab_fac(poly_df$predictor)
box_df$predictor_label  <- lab_fac(box_df$predictor)
poly_df$fillgrp <- factor(poly_df$fillgrp,
                          levels = c("SBP-Responder", "DBP-Responder", "Non-Responder"))

fill_cols <- c("SBP-Responder" = "#1f78b4", "DBP-Responder" = "#33a02c",
               "Non-Responder" = "#9e9e9e")
# Significance (responder vs non-responder association) is reported in the §8
# forest plot, not encoded on this descriptive distribution figure.

panel_theme <- theme_minimal(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    axis.line       = element_line(color = "black", linewidth = 0.4),
    axis.ticks      = element_line(color = "black", linewidth = 0.4),
    axis.text.x     = element_text(face = "bold"),
    axis.title.y    = element_text(face = "bold", size = 9),   # metric name = rotated left title
    legend.position = "bottom",
    plot.margin     = margin(4, 6, 4, 4)
  )

# One self-contained plot per metric (its two split-violins + median/IQR boxes).
# Composed with patchwork below so every panel renders at equal width AND the top
# three center over the bottom four — a single facet_manual grid can't do both:
# centering 3 over 4 overlaps columns and squishes the rightmost bottom panels.
make_panel <- function(metric) {
  pd <- poly_df %>% filter(predictor == metric)
  bd <- box_df  %>% filter(predictor == metric)
  yscale <- if (metric %in% pct_bounded)
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                       expand = expansion(mult = c(0.03, 0.05)))
  else if (metric == "AVG_SLEEP_NEED_PCT")
    scale_y_continuous(limits = c(50, 100), expand = expansion(mult = c(0.03, 0.05)))
  else if (metric == "LOG_N_BP_READINGS_STUDY")
    scale_y_continuous(limits = c(1, NA), breaks = 1:5,
                       expand = expansion(mult = c(0, 0.05)))
  else
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.05)))
  ggplot() +
    geom_polygon(data = pd, aes(x = x, y = y, group = poly_id, fill = fillgrp),
                 color = "gray30", linewidth = 0.25) +
    geom_rect(data = bd, aes(xmin = xmin, xmax = xmax, ymin = q1, ymax = q3),
              fill = "white", color = "gray30", linewidth = 0.25) +
    geom_segment(data = bd, aes(x = xmin, xend = xmax, y = med, yend = med),
                 color = "gray15", linewidth = 0.5) +
    scale_x_continuous(breaks = unname(x_of), labels = names(x_of), limits = c(0.45, 2.95)) +
    yscale +
    scale_fill_manual(values = fill_cols, name = NULL,
                      breaks = c("SBP-Responder", "DBP-Responder", "Non-Responder")) +
    labs(x = NULL, y = pred_labels[[metric]]) +
    panel_theme
}

# Panel order = names(pred_labels): top row Sleep / Physical Activity / BP Readings;
# bottom row Daily App Time / Journal / Notification / AI Coach.
# The rows are built INDEPENDENTLY and stacked. A single spanning design (3 over 4)
# makes the top and bottom panels share grid columns, which squeezes the rightmost
# bottom panel; nesting two equal-outer-width rows avoids that. The top row is
# centered with half-panel spacers each side (0.5,1,1,1,0.5 → each panel = 1/4 of
# the width, matching the bottom row's four 1/4-width panels).
plots <- lapply(names(pred_labels), make_panel)
top_row    <- wrap_plots(plot_spacer(), plots[[1]], plots[[2]], plots[[3]], plot_spacer(),
                         nrow = 1, widths = c(0.5, 1, 1, 1, 0.5))
bottom_row <- wrap_plots(plots[[4]], plots[[5]], plots[[6]], plots[[7]], nrow = 1)
combined <- (top_row / bottom_row) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Six-Month WHOOP Engagement Metrics by Responder Status",
    subtitle = "Per-BP split violins: responder vs non-responder distributions; inner box = median and IQR",
    theme    = theme(plot.title    = element_text(face = "bold", size = 13),
                     plot.subtitle = element_text(size = 11))
  ) &
  theme(legend.position = "bottom")

out_png <- file.path(figures_dir, "8_figure_5_violin_combined.png")
ggsave(out_png, combined, width = 11, height = 8.5, dpi = 300)
cat("Wrote figure:\n  ", out_png, "\n")
