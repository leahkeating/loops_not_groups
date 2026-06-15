import os

DST = "/Users/leah/Documents/loops_not_groups/regular_networks"
models = [
    "loops_no_groups",
    "groups_no_loops",
    "loops_and_groups",
    "simple_contagion",
]

TEMPLATE = r'''library(tidyverse)
library(ggplot2)
library(cowplot)
library(patchwork)
library(RColorBrewer)
theme_set(theme_cowplot())

setwd("/Users/leah/Documents/loops_not_groups/regular_networks")

# Regular networks: curves are parametrised by num_tri (triangles per node),
# deg = 6. Theory only -- there are no regular-network simulations to overlay.
read_theory <- function(path, num_tri_val) {{
  read_delim(path,
             col_names = c("T", "p"),
             delim = " ",
             na = c("", "NA", "nan", "NaN"),
             col_types = "dd") |>
    mutate(num_tri = factor(num_tri_val))
}}

theory_size_lower_3 <- read_theory("data/size_{m}_num_tri_3_lower_branch.txt", 3)
theory_size_lower_2 <- read_theory("data/size_{m}_num_tri_2_lower_branch.txt", 2)
theory_size_lower_1 <- read_theory("data/size_{m}_num_tri_1_lower_branch.txt", 1)

theory_size_upper_3 <- read_theory("data/size_{m}_num_tri_3_upper_branch.txt", 3)
theory_size_upper_2 <- read_theory("data/size_{m}_num_tri_2_upper_branch.txt", 2)
theory_size_upper_1 <- read_theory("data/size_{m}_num_tri_1_upper_branch.txt", 1)

theory_rho_boundary_3 <- read_theory("data/{m}_rho_boundary_num_tri_3.txt", 3)
theory_rho_boundary_2 <- read_theory("data/{m}_rho_boundary_num_tri_2.txt", 2)
theory_rho_boundary_1 <- read_theory("data/{m}_rho_boundary_num_tri_1.txt", 1)

theory_prob_3 <- read_theory("data/prob_{m}_num_tri_3.txt", 3)
theory_prob_2 <- read_theory("data/prob_{m}_num_tri_2.txt", 2)
theory_prob_1 <- read_theory("data/prob_{m}_num_tri_1.txt", 1)

# A tibble with all three values of num_tri for the gcc size (lower and upper branches)
theory_gcc_size_subcrit <- theory_size_lower_1 |> add_row(theory_size_lower_2) |> add_row(theory_size_lower_3)
theory_gcc_size_supcrit <- theory_size_upper_1 |> add_row(theory_size_upper_2) |> add_row(theory_size_upper_3)
connecting_pts <- theory_gcc_size_supcrit |>
  group_by(num_tri) |>
  slice_min(T, n = 1) |>
  ungroup()

theory_gcc_size_rho_boundary <- theory_rho_boundary_1 |> add_row(theory_rho_boundary_2) |> add_row(theory_rho_boundary_3) |>
  filter(!is.na(p))                    # keep only the bistable region (drop NaN tails)
# pin the unstable branch to the lower stable branch (Size 0) at the lower fold,
# so the dashed line meets the flat lower line instead of stopping just above it
bottom_pts <- theory_gcc_size_rho_boundary |>
  group_by(num_tri) |> slice_max(T, n = 1) |> ungroup() |> mutate(p = 0)
theory_gcc_size_rho_boundary <- theory_gcc_size_rho_boundary |>
  bind_rows(connecting_pts) |>         # upper end: meet the upper stable tip
  bind_rows(bottom_pts) |>             # lower end: meet the lower stable branch
  arrange(num_tri, T, desc(p))

# Continuous transitions (upper branch emanates from ~0) have no bistable region,
# so the lower and upper stable branches are really one curve. Merge them into a
# single line per num_tri so they join cleanly. Bistable branches (upper tip well
# above 0) keep their separate lower/upper pieces, joined by the dashed branch.
continuous_nt <- theory_gcc_size_supcrit |>
  group_by(num_tri) |> slice_min(T, n = 1) |> ungroup() |>
  filter(p < 0.1) |> pull(num_tri)
theory_gcc_size_continuous <- bind_rows(
  theory_gcc_size_subcrit |> filter(num_tri %in% continuous_nt),
  theory_gcc_size_supcrit |> filter(num_tri %in% continuous_nt)
) |> arrange(num_tri, T)
theory_gcc_size_subcrit <- theory_gcc_size_subcrit |> filter(!(num_tri %in% continuous_nt))
theory_gcc_size_supcrit <- theory_gcc_size_supcrit |> filter(!(num_tri %in% continuous_nt))

# A tibble with all three values of num_tri for the gcc prob
theory_gcc_prob <- theory_prob_1 |> add_row(theory_prob_2) |> add_row(theory_prob_3)

size_plot <- ggplot() +
  geom_line(data = theory_gcc_size_continuous, aes(x = T, y = p, color = num_tri), linewidth = 1) +
  geom_line(data = theory_gcc_size_subcrit, aes(x = T, y = p, color = num_tri), linewidth = 1) +
  geom_line(data = theory_gcc_size_supcrit, aes(x = T, y = p, color = num_tri), linewidth = 1) +
  geom_line(data = theory_gcc_size_rho_boundary, aes(x = T, y = p, color = num_tri), linetype = "dashed", linewidth = 1) +
  coord_cartesian(xlim = c(0.01, 0.4),
                  ylim = c(0,1)) +
  labs(
    x = expression(T),
    y = expression("Size"),
    color = expression(k[Delta])
  ) +
  theme_minimal() +
  scale_color_brewer(palette = "Set2") +  theme(
    text = element_text(family = "Times"),
    legend.position = c(0.85, 0.35),  # adjust position as needed
    legend.title = element_text(hjust = 0.5),
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA),  # optional: semi-transparent bg
    legend.key = element_blank(),  # optional: remove legend key border
    legend.text = element_text(size = 13),
    panel.grid.major = element_blank(),  # remove major gridlines
    panel.grid.minor = element_blank(),  # remove minor gridlines
    panel.background = element_blank(),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 13),
    axis.line = element_line(color = "black", size = 0.5)
  )

prob_plot <- theory_gcc_prob |>
  ggplot() +
  geom_line(aes(x = T, y = p, color = num_tri), linewidth = 1) +
  coord_cartesian(xlim = c(0.01, 0.4),
                  ylim = c(0,1)) +
  labs(
    x = expression(T),
    y = expression("Probability"),
    color = expression(k[Delta])
  ) +
  theme_minimal() +
  scale_color_brewer(palette = "Set2") +  theme(
    text = element_text(family = "Times"),
    legend.position = c(0.85, 0.35),  # adjust position as needed
    legend.title = element_text(hjust = 0.5),
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA),  # optional: semi-transparent bg
    legend.key = element_blank(),  # optional: remove legend key border
    legend.text = element_text(size = 13),
    panel.grid.major = element_blank(),  # remove major gridlines
    panel.grid.minor = element_blank(),  # remove minor gridlines
    panel.background = element_blank(),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 13),
    axis.line = element_line(color = "black", size = 0.5)
  )

size_plot | prob_plot
'''

for m in models:
    out = os.path.join(DST, f"{m}_plotting.R")
    with open(out, "w") as f:
        f.write(TEMPLATE.format(m=m))
    print("wrote", out)
