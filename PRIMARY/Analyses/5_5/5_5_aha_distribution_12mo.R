# 5_5_aha_distribution_12mo.R
# SAP §5.5 companion (parallels §5.4 at MO12).
#
# Cohort-overall AHA category distribution at baseline / 6mo / 12mo, reported
# side-by-side. Per-stratum down/same/up movement tabulated baseline → MO12.
#
# Cohort: §2.2 hypertensive, restricted to MO12 subset (users with non-missing
# MO12_AHA). Baseline and MO6 distributions reported on the SAME MO12 subset
# (for apples-to-apples comparison across timepoints — not the full-cohort
# MO6 distribution).
#
# Input:  Build/derived/qc/user_analytic_frame_1X.rds
#
# Outputs (Analyses/5_5/derived/):
#   results/5_5_aha_distribution_12mo.csv                                  — marginal: timepoint × category × {n, %}
#   results/5_5_aha_transition_12mo.csv                                    — baseline → MO12 crosstab
#   results/5_5_aha_transition_baseline_to_mo6_within_mo12_subset.csv      — baseline → MO6 crosstab on MO12 subset (added 2026-05-28 per coherence audit)
#   results/5_5_aha_movement_baseline_to_mo6_within_mo12_subset.csv        — per-stratum down/same/up baseline → MO6 on MO12 subset
#   results/5_5_aha_movement_mo6_vs_mo12_comparison.csv                    — pct_down comparison MO12 vs MO6 within MO12 subset
#   results/5_5_aha_distribution_12mo_paper.md
#   results/5_5_aha_distribution_12mo.log
#   figures/5_5_figure_aha_12mo.png            — 3-axis gradient alluvial: Baseline → MO6 → MO12

library(tidyverse)

project_dir <- "./PRIMARY"
input_path  <- file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")
out_dir     <- file.path(project_dir, "Analyses", "5_5", "derived")
results_dir <- file.path(out_dir, "results")
figures_dir <- file.path(out_dir, "figures")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

log_path <- file.path(results_dir, "5_5_aha_distribution_12mo.log")
sink(log_path, split = TRUE)
cat("§5.5 AHA distribution (MO12) —", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

user_frame <- read_rds(input_path)
stopifnot(all(c("BASELINE_AHA", "MO6_AHA", "MO12_AHA") %in% names(user_frame)))

cohort <- user_frame %>%
  filter(BASELINE_AHA %in% c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"))
cat("§2.2 cohort:", nrow(cohort), "users\n")

paired_users <- cohort %>% filter(!is.na(MO12_AHA))
cat("MO12 paired subset:", nrow(paired_users), "users\n\n")

# --- AHA factor order ---
all_levels <- sort(unique(c(paired_users$BASELINE_AHA,
                            paired_users$MO6_AHA,
                            paired_users$MO12_AHA)))
preferred_order <- c("Normal", "Elevated",
                     "Stage 1 Hypertension", "Stage 2 Hypertension",
                     "Hypertensive Crisis")
aha_order <- preferred_order[preferred_order %in% all_levels]
aha_order <- c(aha_order, setdiff(all_levels, aha_order))
cat("AHA factor order:\n  ", paste(aha_order, collapse = " → "), "\n\n")

paired_users <- paired_users %>%
  mutate(
    BASELINE_AHA_f = factor(BASELINE_AHA, levels = aha_order),
    MO6_AHA_f      = factor(MO6_AHA,      levels = aha_order),
    MO12_AHA_f     = factor(MO12_AHA,     levels = aha_order)
  )

n_total <- nrow(paired_users)

# --- Marginal distributions at all three timepoints (MO12 subset only) ---
dist <- bind_rows(
  paired_users %>% count(category = BASELINE_AHA_f, .drop = FALSE) %>%
    mutate(timepoint = "Baseline"),
  paired_users %>% filter(!is.na(MO6_AHA_f)) %>%
    count(category = MO6_AHA_f, .drop = FALSE) %>%
    mutate(timepoint = "6 months"),
  paired_users %>% count(category = MO12_AHA_f, .drop = FALSE) %>%
    mutate(timepoint = "12 months")
) %>%
  group_by(timepoint) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(timepoint = factor(timepoint,
                            levels = c("Baseline", "6 months", "12 months"))) %>%
  arrange(timepoint, category)

write_csv(dist, file.path(results_dir, "5_5_aha_distribution_12mo.csv"))

cat("--- Marginal distribution (Baseline / MO6 / MO12, MO12 subset only) ---\n")
print(dist %>%
        select(timepoint, category, n, prop) %>%
        mutate(prop = scales::percent(prop, accuracy = 0.1)) %>%
        pivot_wider(names_from = timepoint, values_from = c(n, prop)))

# --- Baseline → MO12 transition matrix ---
trans <- paired_users %>%
  count(BASELINE_AHA_f, MO12_AHA_f, .drop = FALSE, name = "n") %>%
  group_by(BASELINE_AHA_f) %>%
  mutate(row_total = sum(n),
         row_pct   = ifelse(row_total > 0, n / row_total, 0)) %>%
  ungroup() %>%
  arrange(BASELINE_AHA_f, MO12_AHA_f)

write_csv(trans, file.path(results_dir, "5_5_aha_transition_12mo.csv"))

cat("\n--- Transition matrix (baseline → MO12), row % within each baseline stratum ---\n")
print(trans %>%
        mutate(cell = sprintf("%d (%.1f%%)", n, 100*row_pct)) %>%
        select(BASELINE_AHA_f, MO12_AHA_f, cell) %>%
        pivot_wider(names_from = MO12_AHA_f, values_from = cell, values_fill = "0 (0.0%)"))

# --- Per-stratum down/same/up movement to MO12 ---
trans_down <- trans %>%
  mutate(
    bl_idx  = as.integer(BASELINE_AHA_f),
    mo_idx  = as.integer(MO12_AHA_f),
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
  ) %>%
  filter(n_total > 0)

cat("\n--- Per-stratum down/same/up movement (baseline → MO12) ---\n")
print(trans_down %>%
        mutate(across(starts_with("pct_"), ~ scales::percent(., accuracy = 0.1))))

# --- Within-MO12-subset baseline → MO6 transitions ---
# Added 2026-05-28 per coherence audit (chapter_5_5/triage.md) CANNOT-VERIFY #1.
# Excerpt claims per-stratum MO12 migration is "each 7 to 10 percentage points
# higher than the corresponding migration at 6 months in the same subset" —
# verification requires baseline→MO6 movement rates computed ON THE MO12
# SUBSET (not on the full §2.2 cohort) for apples-to-apples comparison.

paired_users_mo6 <- paired_users %>% filter(!is.na(MO6_AHA))

trans_mo6 <- paired_users_mo6 %>%
  count(BASELINE_AHA_f, MO6_AHA_f, .drop = FALSE, name = "n") %>%
  group_by(BASELINE_AHA_f) %>%
  mutate(row_total = sum(n),
         row_pct   = ifelse(row_total > 0, n / row_total, 0)) %>%
  ungroup() %>%
  arrange(BASELINE_AHA_f, MO6_AHA_f)

write_csv(trans_mo6, file.path(results_dir, "5_5_aha_transition_baseline_to_mo6_within_mo12_subset.csv"))

trans_down_mo6 <- trans_mo6 %>%
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
  ) %>%
  filter(n_total > 0)

write_csv(trans_down_mo6, file.path(results_dir, "5_5_aha_movement_baseline_to_mo6_within_mo12_subset.csv"))

cat("\n--- Per-stratum down/same/up movement (baseline → MO6, MO12 subset only) ---\n")
print(trans_down_mo6 %>%
        mutate(across(starts_with("pct_"), ~ scales::percent(., accuracy = 0.1))))

# --- Cross-comparison: pct_down at MO12 vs at MO6 (same subset) ---
movement_comparison <- trans_down %>%
  select(BASELINE_AHA_f, pct_down_mo12 = pct_down) %>%
  left_join(
    trans_down_mo6 %>% select(BASELINE_AHA_f, pct_down_mo6 = pct_down),
    by = "BASELINE_AHA_f"
  ) %>%
  mutate(
    diff_pp_mo12_minus_mo6 = (pct_down_mo12 - pct_down_mo6) * 100
  )

write_csv(movement_comparison,
          file.path(results_dir, "5_5_aha_movement_mo6_vs_mo12_comparison.csv"))

cat("\n--- MO6 → MO12 movement comparison (within MO12 subset) ---\n")
cat("    Verifies the excerpt claim that MO12 down-migration is ~7-10 pp higher than MO6 within the same users.\n")
print(movement_comparison %>%
        mutate(
          pct_down_mo6  = scales::percent(pct_down_mo6,  accuracy = 0.1),
          pct_down_mo12 = scales::percent(pct_down_mo12, accuracy = 0.1),
          diff_pp_mo12_minus_mo6 = sprintf("%+0.1f pp", diff_pp_mo12_minus_mo6)
        ))

# --- Figure 2_b: 3-axis gradient-ribbon alluvial (Baseline → MO6 → MO12) ---
# Parallels §5.4 Figure 1D, extended to three timepoints. ggalluvial fills each
# ribbon with one solid colour, so origin→destination gradients are built by
# hand: cosine-eased sankey curves sliced into K thin Lab-interpolated bands.
# Computed on the MO12 paired subset (all have BL+MO6+MO12 non-missing), so the
# three bars are equal height by construction — apples-to-apples within the same
# users. Stratum colours matched to Figs 5_1/5_2.
aha_palette <- c(
  "Normal"               = "#478A43",  # forest green  (Fig 5_2 "improved")
  "Elevated"             = "#EDCE44",  # AHA yellow, muted + warmed for readability
  "Stage 1 Hypertension" = "#FFB600",  # AHA amber      (Fig 5_1)
  "Stage 2 Hypertension" = "#BB3A03",  # AHA red-brown  (Fig 5_1)
  "Hypertensive Crisis"  = "#990712"   # AHA crisis red (per house style)
)
aha_palette <- aha_palette[names(aha_palette) %in% aha_order]

# Per-leg flow counts.
flows1 <- paired_users %>% transmute(origin = BASELINE_AHA, dest = MO6_AHA)  %>% count(origin, dest, name = "n")
flows2 <- paired_users %>% transmute(origin = MO6_AHA,      dest = MO12_AHA) %>% count(origin, dest, name = "n")

# Axis stratum y-extents (stacked top-down in aha_order, full height = n_total).
axis_extents <- function(tab) {
  tab %>% mutate(o = factor(cat, levels = aha_order)) %>% arrange(o) %>%
    mutate(top = n_total - cumsum(lag(tot, default = 0)), bot = top - tot) %>%
    select(cat, tot, top, bot)
}
base_strata <- axis_extents(flows1 %>% group_by(cat = origin) %>% summarise(tot = sum(n), .groups = "drop"))
mo6_strata  <- axis_extents(flows1 %>% group_by(cat = dest)   %>% summarise(tot = sum(n), .groups = "drop"))
mo12_strata <- axis_extents(flows2 %>% group_by(cat = dest)   %>% summarise(tot = sum(n), .groups = "drop"))

# Ribbon geometry. Wider bars (hw) so the stratum labels fit within bar width.
hw <- 0.14; K <- 100
ease <- function(t) 0.5 - 0.5 * cos(pi * t)   # cosine S-curve

# One leg of gradient ribbons between origin axis (right edge xL) and destination
# axis (left edge xR). Left bands stack by dest order within each origin stratum;
# right bands stack by origin order within each dest stratum (so the middle MO6
# bar conserves: incoming bands and outgoing bands each tile its full height).
build_leg <- function(flows_ab, left_ext, right_ext, xL, xR, tag) {
  fl <- flows_ab %>%
    left_join(left_ext %>% select(origin = cat, otop = top), by = "origin") %>%
    mutate(d = factor(dest, levels = aha_order)) %>%
    group_by(origin) %>% arrange(d, .by_group = TRUE) %>%
    mutate(yL_top = otop - cumsum(lag(n, default = 0)), yL_bot = yL_top - n) %>% ungroup()
  fr <- flows_ab %>%
    left_join(right_ext %>% select(dest = cat, dtop = top), by = "dest") %>%
    mutate(o = factor(origin, levels = aha_order)) %>%
    group_by(dest) %>% arrange(o, .by_group = TRUE) %>%
    mutate(yR_top = dtop - cumsum(lag(n, default = 0)), yR_bot = yR_top - n) %>% ungroup()
  f2 <- fl %>% select(origin, dest, n, yL_top, yL_bot) %>%
    left_join(fr %>% select(origin, dest, yR_top, yR_bot), by = c("origin", "dest")) %>%
    arrange(desc(n)) %>% mutate(flow_id = paste(tag, row_number(), sep = "_"))
  xeps <- (xR - xL) / K * 0.6
  pmap_dfr(f2, function(flow_id, origin, dest, n, yL_top, yL_bot, yR_top, yR_bot, ...) {
    x <- seq(xL, xR, length.out = K + 1)
    ytop_at <- function(xx) yL_top + (yR_top - yL_top) * ease((xx - xL) / (xR - xL))
    ybot_at <- function(xx) yL_bot + (yR_bot - yL_bot) * ease((xx - xL) / (xR - xL))
    ramp <- grDevices::colorRampPalette(c(aha_palette[[origin]], aha_palette[[dest]]),
                                        space = "Lab")(K)
    map_dfr(seq_len(K), function(k) {
      xa <- x[k]; xb <- min(x[k + 1] + xeps, xR)
      tibble(poly = paste(flow_id, k, sep = "_"), fill = ramp[k],
             px = c(xa, xb, xb, xa),
             py = c(ybot_at(xa), ybot_at(xb), ytop_at(xb), ytop_at(xa)))
    })
  })
}

slices <- bind_rows(
  build_leg(flows1, base_strata, mo6_strata,  1 + hw, 2 - hw, "L1"),
  build_leg(flows2, mo6_strata,  mo12_strata, 2 + hw, 3 - hw, "L2")
)

# Strata rectangles + labels ("Stage 1/2 Hypertension" wraps to fit inside bar).
strata_rects <- bind_rows(
  base_strata %>% mutate(x = 1),
  mo6_strata  %>% mutate(x = 2),
  mo12_strata %>% mutate(x = 3)
) %>%
  mutate(xmin = x - hw, xmax = x + hw, ymin = bot, ymax = top,
         fill = unname(aha_palette[cat]), xmid = x, ymid = (bot + top) / 2,
         cat_lab = str_replace(cat, " Hypertension", "\nHypertension"),
         lab     = sprintf("%s\n%.1f%%", cat_lab, 100 * tot / n_total))

p <- ggplot() +
  geom_polygon(data = slices, aes(px, py, group = poly, fill = fill), alpha = 1) +
  geom_rect(data = strata_rects,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            color = "white", linewidth = 0.4) +
  geom_text(data = strata_rects, aes(xmid, ymid, label = lab),
            size = 2.6, fontface = "bold", color = "black", lineheight = 0.9) +
  scale_fill_identity() +
  scale_x_continuous(breaks = c(1, 2, 3),
                     labels = c("Baseline", "6 months", "12 months"),
                     expand = expansion(mult = c(0.10, 0.10))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = NULL,
       title = "Reclassification Across AHA Blood Pressure Categories Through 12 Months") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid      = element_blank(),
    axis.text.y     = element_blank(),
    axis.ticks      = element_blank(),
    axis.text.x     = element_text(face = "bold", size = 12),
    legend.position = "none",
    plot.title      = element_text(face = "bold")
  )

out_png <- file.path(figures_dir, "5_5_figure_aha_12mo.png")
ggsave(out_png, p, width = 10, height = 6.5, dpi = 300)
cat("\nWrote figure (3-axis gradient alluvial):\n  ", out_png, "\n")

# --- Paper-ready markdown ---
pct <- function(x) scales::percent(x, accuracy = 0.1)

build_line <- function(tp) {
  rows <- dist %>% filter(timepoint == tp)
  paste(sprintf("%s %s (n = %s)", rows$category, pct(rows$prop),
                format(rows$n, big.mark = ",")), collapse = "; ")
}

md_lines <- c(
  "# §5.5 — AHA category distribution at 12 months",
  "",
  paste0("_Generated ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "_"),
  "",
  "## Cohort-overall AHA distribution (MO12 paired subset)",
  "",
  paste0("- **Baseline:** ", build_line("Baseline"), "."),
  paste0("- **6 months:** ", build_line("6 months"), "."),
  paste0("- **12 months:** ", build_line("12 months"), "."),
  paste0("- MO12 paired subset n = ", format(n_total, big.mark = ","), "."),
  "",
  "## Per-stratum movement baseline → MO12",
  ""
)
for (i in seq_len(nrow(trans_down))) {
  r <- trans_down[i, ]
  md_lines <- c(md_lines, sprintf(
    "- **%s (baseline n = %s):** moved down %s (n = %s); same %s (n = %s); moved up %s (n = %s).",
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
  "- All three timepoint distributions computed on the MO12 paired subset only — apples-to-apples within the same users across baseline / 6mo / 12mo.",
  "- The Baseline + MO6 distributions reported here are NOT the same as §5.4's full-cohort distributions; they reflect only the ~7% MO12 paired subset. Comparing within this subset isolates trajectory from selection composition.",
  "- MO12 subset is engagement-selected. Higher 'Normal' proportions at MO12 than at MO6 may reflect continued improvement, sustained engagement, or both — not separable from the data.",
  "- Per SAP, no formal test. Descriptive only.",
  "- Pairs with: §5.4 (full-cohort MO6 AHA distribution), §5.5 LMM (trajectory analysis), §5.5 clinical bins (MO12 per-user proportions)."
)
writeLines(md_lines, file.path(results_dir, "5_5_aha_distribution_12mo_paper.md"))

cat("\nWrote markdown:\n  ", file.path(results_dir, "5_5_aha_distribution_12mo_paper.md"), "\n")
cat("\n=== §5.5 12mo AHA distribution done ===\n")
sink()
