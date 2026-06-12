library(tidyverse)
library(ggplot2)
library(cowplot)
library(patchwork)
library(RColorBrewer)
theme_set(theme_cowplot())

setwd("/Users/leah/Documents/loops_not_groups")

theory_size_lower_3 <- read_delim("data/size_simple_contagion_nu_3_lower_branch.txt",
                            col_names = c("T", "p"),
                            delim = " ") |> mutate(nu = factor(3))

theory_size_lower_2 <- read_delim("data/size_simple_contagion_nu_2_lower_branch.txt",
                            col_names = c("T", "p"),
                            delim = " ") |> mutate(nu = factor(2))

theory_size_lower_1 <- read_delim("data/size_simple_contagion_nu_1_lower_branch.txt",
                            col_names = c("T", "p"),
                            delim = " ") |> mutate(nu = factor(1))

theory_size_upper_3 <- read_delim("data/size_simple_contagion_nu_3_upper_branch.txt",
                            col_names = c("T", "p"),
                            delim = " ") |> mutate(nu = factor(3))

theory_size_upper_2 <- read_delim("data/size_simple_contagion_nu_2_upper_branch.txt",
                            col_names = c("T", "p"),
                            delim = " ") |> mutate(nu = factor(2))

theory_size_upper_1 <- read_delim("data/size_simple_contagion_nu_1_upper_branch.txt",
                            col_names = c("T", "p"),
                            delim = " ") |> mutate(nu = factor(1))

theory_rho_boundary_3 <- read_delim("data/simple_contagion_rho_boundary_nu_3.txt",
                            col_names = c("T", "p"),
                            delim = " ") |> mutate(nu = factor(3))

theory_rho_boundary_2 <- read_delim("data/simple_contagion_rho_boundary_nu_2.txt",
                            col_names = c("T", "p"),
                            delim = " ") |> mutate(nu = factor(2))

theory_rho_boundary_1 <- read_delim("data/simple_contagion_rho_boundary_nu_1.txt",
                            col_names = c("T", "p"),
                            delim = " ") |> mutate(nu = factor(1))

theory_prob_3 <- read_delim("data/prob_simple_contagion_nu_3.txt",
                            col_names = c("T", "p"),
                            delim = " ") |> mutate(nu = factor(3))

theory_prob_2 <- read_delim("data/prob_simple_contagion_nu_2.txt",
                            col_names = c("T", "p"),
                            delim = " ") |> mutate(nu = factor(2))

theory_prob_1 <- read_delim("data/prob_simple_contagion_nu_1.txt",
                            col_names = c("T", "p"),
                            delim = " ") |> mutate(nu = factor(1))

sims_3 <- read_delim("data/sc_nu_3_sims.txt",
                     col_names = c("T", "cascade_size"),
                     delim = " ") |> mutate(nu = factor(3)) |>
  mutate(rho = cascade_size/100000)

sims_3_near_boundary <- read_delim("data/logs_2000000/sc_nu_3_sims.txt",
                     col_names = c("T", "cascade_size"),
                     delim = " ") |> mutate(nu = factor(3)) |>
  mutate(rho = cascade_size/2000000)

# thin the near-boundary data to match the T spacing of the 100000-node grid
T_spacing <- 0.4/30
near_boundary_T <- sims_3_near_boundary |> distinct(T) |>
  mutate(bin = round(T / T_spacing)) |>
  group_by(bin) |>
  slice_min(abs(T - bin * T_spacing), n = 1) |>
  pull(T)

sims_3 <- sims_3 |> filter((T < 0.1) | (T > 0.15)) |>
  add_row(sims_3_near_boundary |> filter(T >= 0.1, T %in% near_boundary_T))

sims_2 <- read_delim("data/sc_nu_2_sims.txt",
                     col_names = c("T", "cascade_size"),
                     delim = " ") |> mutate(nu = factor(2)) |>
  mutate(rho = cascade_size/100000)

sims_1 <- read_delim("data/sc_nu_1_sims.txt",
                     col_names = c("T", "cascade_size"),
                     delim = " ") |> mutate(nu = factor(1)) |>
  mutate(rho = cascade_size/100000)

# A tibble with all three values for nu for the gcc size (lower and upper branches)
theory_gcc_size_subcrit <- theory_size_lower_1 |> add_row(theory_size_lower_2) |> add_row(theory_size_lower_3)
theory_gcc_size_supcrit <- theory_size_upper_1 |> add_row(theory_size_upper_2) |> add_row(theory_size_upper_3)
connecting_pts <- theory_gcc_size_supcrit |>
  group_by(nu) |>
  slice_min(T, n = 1) |>
  ungroup()

theory_gcc_size_rho_boundary <- theory_rho_boundary_1 |> add_row(theory_rho_boundary_2) |> add_row(theory_rho_boundary_3) |>
  bind_rows(connecting_pts) |>
  arrange(nu, T)

# A tibble with all three values for nu for the gcc prob
theory_gcc_prob <- theory_prob_1 |> add_row(theory_prob_2) |> add_row(theory_prob_3)

sims <- sims_1 |> add_row(sims_2) |> add_row(sims_3)

size_sims <- sims |>
  mutate(significant_frac = rho > 0.01) |>
  group_by(nu, T, significant_frac) |>
  summarise(giant_component_size = median(rho)) |>
  group_by(nu, T) |>
  summarise(giant_component_size = max(giant_component_size)) |> ungroup()

prob_sims <- sims |>
  mutate(significant_frac = rho > 0.01) |>
  group_by(nu, T) |>
  summarise(p = sum(significant_frac) / n()) |> ungroup()

size_plot <- ggplot() +
  geom_line(data = theory_gcc_size_subcrit, aes(x = T, y = p, color = nu)) +
  geom_line(data = theory_gcc_size_supcrit, aes(x = T, y = p, color = nu)) +
  geom_line(data = theory_gcc_size_rho_boundary, aes(x = T, y = p, color = nu), linetype = "dashed") +
  geom_point(data = size_sims, aes(x = T, y = giant_component_size, color = nu, shape = nu)) +
  coord_cartesian(xlim = c(0.01, 0.4),
                  ylim = c(0,1)) +
  labs(
    x = expression(T),
    y = expression("Size"),
    color = expression(nu),
    shape = expression(nu)
  ) +
  theme_minimal() +
  scale_color_manual(
    values = (brewer.pal(9, "Blues")[5:10])
  ) +  theme(
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
  geom_line(aes(x = T, y = p, color = nu)) +
  geom_point(data = prob_sims, aes(x = T, y = p, color = nu, shape = nu)) +
  coord_cartesian(xlim = c(0.01, 0.4),
                  ylim = c(0,1)) +
  labs(
    x = expression(T),
    y = expression("Probability"),
    color = expression(nu),
    shape = expression(nu)
  ) +
  theme_minimal() +
  scale_color_manual(
    values = (brewer.pal(9, "Blues")[5:10])
  ) +  theme(
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
