# 5_4_aha_category_shift.R
# SAP §5.4: AHA category proportions at baseline vs 6 months.
#
# Cohort-overall (§2.2 hypertensive) AHA distribution at baseline and at MO6,
# compared side-by-side. Hypothesis (descriptive): proportion of lower BP
# categories (Normal, Elevated) increases by 6 months.
#
# Beyond SAP spec (added because nearly free and useful for narrative):
#   - Transition matrix: baseline AHA × MO6 AHA per-user crosstab.
#
# Per SAP §5.4: no formal test.
#
# Input:  Build/derived/qc/user_analytic_frame_1X.rds
#
# Outputs (Analyses/5_4/derived/):
#   results/5_4_distribution.csv         — baseline vs MO6 marginal: category × {n, %}
#   results/5_4_transition_matrix.csv    — baseline AHA × MO6 AHA crosstab (counts + row %)
#   results/5_4_paper.md                 — paper-ready narrative
#   results/5_4_aha_category_shift.log
#   figures/5_4_figure_1d.png            — gradient-ribbon alluvial (baseline → MO6 reclassification)

library(tidyverse)

project_dir <- "./PRIMARY"
input_path  <- file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")
out_dir     <- file.path(project_dir, "Analyses", "5_4", "derived")
results_dir <- file.path(out_dir, "results")
figures_dir <- file.path(out_dir, "figures")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

log_path <- file.path(results_dir, "5_4_aha_category_shift.log")
sink(log_path, split = TRUE)
cat("§5.4 AHA category shift —", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Read frame ---
user_frame <- read_rds(input_path)
cat("User frame:", nrow(user_frame), "rows ×", ncol(user_frame), "cols\n")

stopifnot("BASELINE_AHA" %in% names(user_frame))
stopifnot("MO6_AHA"      %in% names(user_frame))

# --- §2.2 cohort ---
cohort <- user_frame %>%
  filter(BASELINE_AHA %in% c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"))
cat("§2.2 cohort:", nrow(cohort), "users\n")

# Restrict to users with non-missing MO6_AHA (paired analysis)
paired_users <- cohort %>% filter(!is.na(MO6_AHA))
cat("Paired (have both BASELINE_AHA and MO6_AHA):", nrow(paired_users), "users\n\n")

# --- Inspect MO6_AHA category levels to set factor order defensively ---
mo6_levels <- sort(unique(c(paired_users$BASELINE_AHA, paired_users$MO6_AHA)))
cat("All AHA levels observed (baseline ∪ MO6):\n  ", paste(mo6_levels, collapse = ", "), "\n\n")

# Conventional AHA ordering: Normal → Elevated → Stage 1 → Stage 2 (→ Crisis if present)
preferred_order <- c("Normal", "Elevated",
                     "Stage 1 Hypertension", "Stage 2 Hypertension",
                     "Hypertensive Crisis")
aha_order <- preferred_order[preferred_order %in% mo6_levels]
# Append anything unexpected at the end (defensive)
extras <- setdiff(mo6_levels, aha_order)
aha_order <- c(aha_order, extras)
cat("AHA factor order used:\n  ", paste(aha_order, collapse = " → "), "\n\n")

paired_users <- paired_users %>%
  mutate(
    BASELINE_AHA_f = factor(BASELINE_AHA, levels = aha_order),
    MO6_AHA_f      = factor(MO6_AHA,      levels = aha_order)
  )

# --- Marginal distribution: baseline vs MO6 ---
n_total <- nrow(paired_users)
dist <- bind_rows(
  paired_users %>% count(category = BASELINE_AHA_f, .drop = FALSE) %>%
    mutate(timepoint = "Baseline"),
  paired_users %>% count(category = MO6_AHA_f, .drop = FALSE) %>%
    mutate(timepoint = "6 months")
) %>%
  group_by(timepoint) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(timepoint = factor(timepoint, levels = c("Baseline", "6 months"))) %>%
  arrange(timepoint, category)

write_csv(dist, file.path(results_dir, "5_4_distribution.csv"))

cat("--- Marginal distribution (baseline vs MO6) ---\n")
print(dist %>%
        select(timepoint, category, n, prop) %>%
        mutate(prop = scales::percent(prop, accuracy = 0.1)) %>%
        pivot_wider(names_from = timepoint, values_from = c(n, prop)))

# --- Transition matrix: baseline × MO6 ---
trans <- paired_users %>%
  count(BASELINE_AHA_f, MO6_AHA_f, .drop = FALSE, name = "n") %>%
  group_by(BASELINE_AHA_f) %>%
  mutate(row_total = sum(n),
         row_pct   = ifelse(row_total > 0, n / row_total, 0)) %>%
  ungroup() %>%
  arrange(BASELINE_AHA_f, MO6_AHA_f)

write_csv(trans, file.path(results_dir, "5_4_transition_matrix.csv"))

cat("\n--- Transition matrix (baseline → MO6), row percentages within each baseline stratum ---\n")
print(trans %>%
        mutate(cell = sprintf("%d (%.1f%%)", n, 100*row_pct)) %>%
        select(BASELINE_AHA_f, MO6_AHA_f, cell) %>%
        pivot_wider(names_from = MO6_AHA_f, values_from = cell, values_fill = "0 (0.0%)"))

# --- Figure 1D: gradient-ribbon alluvial (baseline → MO6 reclassification) ---
# ggalluvial fills each ribbon with one solid colour, so the origin→destination
# gradient is built by hand: cosine-eased sankey curves sliced into K thin bands
# whose fill interpolates (Lab space) from the baseline colour to the 6-month
# colour. Solid-colour slices => renders on any device (no ggalluvial/ragg).
# Stratum colours matched to Figs 5_1/5_2: AHA chart hues + 5_2 forest green.
aha_palette <- c(
  "Normal"               = "#478A43",  # forest green  (Fig 5_2 "improved")
  "Elevated"             = "#EDCE44",  # AHA yellow, muted + warmed for readability
  "Stage 1 Hypertension" = "#FFB600",  # AHA amber      (Fig 5_1)
  "Stage 2 Hypertension" = "#BB3A03",  # AHA red-brown  (Fig 5_1)
  "Hypertensive Crisis"  = "#990712"   # AHA crisis red (per house style)
)
aha_palette <- aha_palette[names(aha_palette) %in% aha_order]

# Flow counts: every paired user on each baseline → MO6 path.
flows <- paired_users %>%
  transmute(origin = BASELINE_AHA, dest = MO6_AHA) %>%
  count(origin, dest, name = "n")

# Axis stratum y-extents (stacked top-down in aha_order, full height = n_total).
left_strata <- flows %>%
  group_by(origin) %>% summarise(tot = sum(n), .groups = "drop") %>%
  mutate(o = factor(origin, levels = aha_order)) %>% arrange(o) %>%
  mutate(top = n_total - cumsum(lag(tot, default = 0)), bot = top - tot)
right_strata <- flows %>%
  group_by(dest) %>% summarise(tot = sum(n), .groups = "drop") %>%
  mutate(d = factor(dest, levels = aha_order)) %>% arrange(d) %>%
  mutate(top = n_total - cumsum(lag(tot, default = 0)), bot = top - tot)

# Flow band extents at each axis (left sub-ordered by dest, right by origin).
flows_left <- flows %>%
  left_join(left_strata %>% select(origin, otop = top), by = "origin") %>%
  mutate(d = factor(dest, levels = aha_order)) %>%
  group_by(origin) %>% arrange(d, .by_group = TRUE) %>%
  mutate(yL_top = otop - cumsum(lag(n, default = 0)), yL_bot = yL_top - n) %>% ungroup()
flows_right <- flows %>%
  left_join(right_strata %>% select(dest, dtop = top), by = "dest") %>%
  mutate(o = factor(origin, levels = aha_order)) %>%
  group_by(dest) %>% arrange(o, .by_group = TRUE) %>%
  mutate(yR_top = dtop - cumsum(lag(n, default = 0)), yR_bot = yR_top - n) %>% ungroup()
flows2 <- flows_left %>%
  select(origin, dest, n, yL_top, yL_bot) %>%
  left_join(flows_right %>% select(origin, dest, yR_top, yR_bot), by = c("origin", "dest")) %>%
  arrange(desc(n)) %>%                 # big ribbons drawn first (behind)
  mutate(flow_id = row_number())

# Geometry of the manual sankey ribbons.
hw <- 0.08                              # stratum half-width
xL <- 1 + hw                            # left ribbon edge
xR <- 2 - hw                            # right ribbon edge
K  <- 100                               # slices per ribbon (gradient resolution)
ease <- function(t) 0.5 - 0.5 * cos(pi * t)   # cosine S-curve

make_slices <- function(flow_id, origin, dest, yL_top, yL_bot, yR_top, yR_bot, ...) {
  x <- seq(xL, xR, length.out = K + 1)
  ytop_at <- function(xx) yL_top + (yR_top - yL_top) * ease((xx - xL) / (xR - xL))
  ybot_at <- function(xx) yL_bot + (yR_bot - yL_bot) * ease((xx - xL) / (xR - xL))
  ramp <- grDevices::colorRampPalette(c(aha_palette[[origin]], aha_palette[[dest]]),
                                      space = "Lab")(K)
  xeps <- (xR - xL) / K * 0.6           # curve-following overlap to hide seams
  purrr::map_dfr(seq_len(K), function(k) {
    xa <- x[k]; xb <- min(x[k + 1] + xeps, xR)
    tibble(poly = paste(flow_id, k, sep = "_"), fill = ramp[k],
           px = c(xa, xb, xb, xa),
           py = c(ybot_at(xa), ybot_at(xb), ytop_at(xb), ytop_at(xa)))
  })
}
slices <- purrr::pmap_dfr(flows2, make_slices)

# Strata rectangles + labels ("Stage 1/2 Hypertension" wraps to fit inside bar).
strata_rects <- bind_rows(
  left_strata  %>% transmute(cat = origin, ymin = bot, ymax = top, xmin = 1 - hw, xmax = 1 + hw),
  right_strata %>% transmute(cat = dest,   ymin = bot, ymax = top, xmin = 2 - hw, xmax = 2 + hw)
) %>%
  mutate(fill = unname(aha_palette[cat]),
         xmid = (xmin + xmax) / 2, ymid = (ymin + ymax) / 2,
         cat_lab = str_replace(cat, " Hypertension", "\nHypertension"),
         lab     = sprintf("%s\n%.1f%%", cat_lab, 100 * (ymax - ymin) / n_total))

p <- ggplot() +
  geom_polygon(data = slices, aes(px, py, group = poly, fill = fill), alpha = 1) +
  geom_rect(data = strata_rects,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            color = "white", linewidth = 0.4) +
  geom_text(data = strata_rects, aes(xmid, ymid, label = lab),
            size = 2.7, fontface = "bold", color = "black", lineheight = 0.9) +
  scale_fill_identity() +
  scale_x_continuous(breaks = c(1, 2), labels = c("Baseline", "6 months"),
                     expand = expansion(mult = c(0.10, 0.10))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = NULL,
       title = "Six-Month Reclassification Across AHA Blood Pressure Categories") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid      = element_blank(),
    axis.text.y     = element_blank(),
    axis.ticks      = element_blank(),
    axis.text.x     = element_text(face = "bold", size = 12),
    legend.position = "none",
    plot.title      = element_text(face = "bold")
  )

out_png <- file.path(figures_dir, "5_4_figure_1d.png")
ggsave(out_png, p, width = 7.5, height = 6.5, dpi = 300)
cat("\nWrote figure (gradient alluvial):\n  ", out_png, "\n")

# --- Paper-ready markdown ---
pct <- function(x) scales::percent(x, accuracy = 0.1)

bl <- dist %>% filter(timepoint == "Baseline")
mo <- dist %>% filter(timepoint == "6 months")

# Numbers for narrative
bl_line <- paste(sprintf("%s %s (n = %s)", bl$category, pct(bl$prop),
                         format(bl$n, big.mark = ",")), collapse = "; ")
mo_line <- paste(sprintf("%s %s (n = %s)", mo$category, pct(mo$prop),
                         format(mo$n, big.mark = ",")), collapse = "; ")

# Per-stratum "down-shift" rates (moved to a lower-severity category)
trans_down <- trans %>%
  mutate(
    bl_idx  = as.integer(BASELINE_AHA_f),
    mo_idx  = as.integer(MO6_AHA_f),
    is_down = mo_idx < bl_idx,
    is_up   = mo_idx > bl_idx,
    is_same = mo_idx == bl_idx
  ) %>%
  group_by(BASELINE_AHA_f) %>%
  summarise(
    n_total   = sum(n),
    n_down    = sum(n[is_down]),
    n_same    = sum(n[is_same]),
    n_up      = sum(n[is_up]),
    pct_down  = n_down / n_total,
    pct_same  = n_same / n_total,
    pct_up    = n_up   / n_total,
    .groups = "drop"
  )

cat("\n--- Per-stratum down/same/up movement ---\n")
print(trans_down %>%
        mutate(across(starts_with("pct_"), ~ scales::percent(., accuracy = 0.1))))

md_lines <- c(
  "# §5.4 — AHA category proportions, baseline vs 6 months",
  "",
  paste0("_Generated ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "_"),
  "",
  "## Cohort-overall AHA distribution",
  "",
  paste0("- **Baseline:** ", bl_line, "."),
  paste0("- **6 months:** ", mo_line, "."),
  paste0("- Paired cohort n = ", format(n_total, big.mark = ","),
         " users with both baseline AHA and MO6 AHA in the §2.2 cohort."),
  "",
  "## Per-stratum movement (baseline → MO6)",
  ""
)
for (i in seq_len(nrow(trans_down))) {
  r <- trans_down[i, ]
  md_lines <- c(md_lines, sprintf(
    "- **%s (baseline n = %s):** moved to a lower-severity category in %s (n = %s); remained in same category in %s (n = %s); moved to a higher-severity category in %s (n = %s).",
    as.character(r$BASELINE_AHA_f), format(r$n_total, big.mark = ","),
    pct(r$pct_down), format(r$n_down, big.mark = ","),
    pct(r$pct_same), format(r$n_same, big.mark = ","),
    pct(r$pct_up),   format(r$n_up,   big.mark = ",")
  ))
}
md_lines <- c(md_lines,
  "",
  "## Notes",
  "",
  "- Cohort = §2.2 hypertensive (Elevated + Stage 1 + Stage 2 at baseline) restricted to users with non-missing MO6_AHA.",
  "- Baseline distribution is by construction limited to Elevated / Stage 1 / Stage 2 (Normal is 0% at baseline by cohort definition).",
  "- MO6 distribution can populate any of the five AHA categories (Normal, Elevated, Stage 1, Stage 2, Hypertensive Crisis if observed).",
  "- Per SAP §5.4, no formal hypothesis test is performed. Movement is reported as observed counts and percentages.",
  "- Per-stratum down/same/up movement is added beyond SAP spec because it directly answers 'how many users crossed a clinical boundary' — useful narrative material that the marginal distribution alone obscures.",
  "- RTM caveat: AHA category at a single timepoint is a noisy classification, particularly near category boundaries. A user with true SBP near 130 mmHg may classify into Elevated or Stage 1 depending on measurement noise. The category-shift counts therefore inflate apparent transitions; this is structural and not corrected here.",
  "- Pairs with: §5.1 trajectory (continuous Δ), §5.2 clinical-threshold bins (per-user ±5 mmHg), Figure 1D (this artifact)."
)
writeLines(md_lines, file.path(results_dir, "5_4_paper.md"))
cat("\nWrote markdown:\n  ", file.path(results_dir, "5_4_paper.md"), "\n")

cat("\n=== §5.4 AHA category shift done ===\n")
sink()
