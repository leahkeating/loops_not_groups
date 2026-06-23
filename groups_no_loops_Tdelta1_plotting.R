library(tidyverse)
library(ggplot2)
library(cowplot)
library(patchwork)
library(RColorBrewer)
theme_set(theme_cowplot())

setwd("/Users/leah/Documents/loops_not_groups")

theory_size_lower_3 <- read_delim("data/size_groups_no_loops_Tdelta1_nu_3_lower_branch.txt",
                            col_names = c("T", "p"),
                            delim = " ") |> mutate(nu = factor(3))

theory_size_lower_2 <- read_delim("data/size_groups_no_loops_Tdelta1_nu_2_lower_branch.txt",
                            col_names = c("T", "p"),
                            delim = " ") |> mutate(nu = factor(2))

theory_size_lower_1 <- read_delim("data/size_groups_no_loops_Tdelta1_nu_1_lower_branch.txt",
                            col_names = c("T", "p"),
                            delim = " ") |> mutate(nu = factor(1))

theory_size_upper_3 <- read_delim("data/size_groups_no_loops_Tdelta1_nu_3_upper_branch.txt",
                            col_names = c("T", "p"),
                            delim = " ") |> mutate(nu = factor(3))

theory_size_upper_2 <- read_delim("data/size_groups_no_loops_Tdelta1_nu_2_upper_branch.txt",
                            col_names = c("T", "p"),
                            delim = " ") |> mutate(nu = factor(2))

theory_size_upper_1 <- read_delim("data/size_groups_no_loops_Tdelta1_nu_1_upper_branch.txt",
                            col_names = c("T", "p"),
                            delim = " ") |> mutate(nu = factor(1))

theory_rho_boundary_3 <- read_delim("data/groups_no_loops_Tdelta1_rho_boundary_nu_3.txt",
                            col_names = c("T", "p"),
                            delim = " ") |> mutate(nu = factor(3))

theory_rho_boundary_2 <- read_delim("data/groups_no_loops_Tdelta1_rho_boundary_nu_2.txt",
                            col_names = c("T", "p"),
                            delim = " ") |> mutate(nu = factor(2))

theory_rho_boundary_1 <- read_delim("data/groups_no_loops_Tdelta1_rho_boundary_nu_1.txt",
                            col_names = c("T", "p"),
                            delim = " ") |> mutate(nu = factor(1))

theory_prob_3 <- read_delim("data/prob_groups_no_loops_Tdelta1_nu_3.txt",
                            col_names = c("T", "p"),
                            delim = " ") |> mutate(nu = factor(3))

theory_prob_2 <- read_delim("data/prob_groups_no_loops_Tdelta1_nu_2.txt",
                            col_names = c("T", "p"),
                            delim = " ") |> mutate(nu = factor(2))

theory_prob_1 <- read_delim("data/prob_groups_no_loops_Tdelta1_nu_1.txt",
                            col_names = c("T", "p"),
                            delim = " ") |> mutate(nu = factor(1))

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

size_plot <- ggplot() +
  geom_line(data = theory_gcc_size_subcrit, aes(x = T, y = p, color = nu)) +
  geom_line(data = theory_gcc_size_supcrit, aes(x = T, y = p, color = nu)) +
  geom_line(data = theory_gcc_size_rho_boundary, aes(x = T, y = p, color = nu), linetype = "dashed") +
  coord_cartesian(xlim = c(0.01, 0.4),
                  ylim = c(0,1)) +
  labs(
    x = expression(T),
    y = expression("Size"),
    color = expression(nu),
    shape = expression(nu)
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
  geom_line(aes(x = T, y = p, color = nu)) +
  coord_cartesian(xlim = c(0.01, 0.4),
                  ylim = c(0,1)) +
  labs(
    x = expression(T),
    y = expression("Probability"),
    color = expression(nu),
    shape = expression(nu)
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
