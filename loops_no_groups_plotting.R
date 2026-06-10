library(tidyverse)
library(ggplot2)
library(cowplot)
library(patchwork)
library(RColorBrewer)
theme_set(theme_cowplot())

setwd("/Users/leah/Documents/loops_not_groups")

theory_size_3 <- read_delim("data/size_loops_not_groups_nu_3.txt",
                            col_names = c("T", "p"),
                            delim = " ") |> mutate(nu = factor(3))

theory_size_2 <- read_delim("data/size_loops_not_groups_nu_2.txt",
                            col_names = c("T", "p"),
                            delim = " ") |> mutate(nu = factor(2))

theory_size_1 <- read_delim("data/size_loops_not_groups_nu_1.txt",
                            col_names = c("T", "p"),
                            delim = " ") |> mutate(nu = factor(1))

theory_prob_3 <- read_delim("data/prob_loops_not_groups_nu_3.txt",
                            col_names = c("T", "p"),
                            delim = " ") |> mutate(nu = factor(3))

theory_prob_2 <- read_delim("data/prob_loops_not_groups_nu_2.txt",
                            col_names = c("T", "p"),
                            delim = " ") |> mutate(nu = factor(2))

theory_prob_1 <- read_delim("data/prob_loops_not_groups_nu_1.txt",
                            col_names = c("T", "p"),
                            delim = " ") |> mutate(nu = factor(1))

# A tibble with all three values for nu for the gcc size
theory_gcc_size <- theory_size_1 |> add_row(theory_size_2) |> add_row(theory_size_3)

# A tibble with all three values for nu for the gcc prob
theory_gcc_prob <- theory_prob_1 |> add_row(theory_prob_2) |> add_row(theory_prob_3)

size_plot <- theory_gcc_size |>
  ggplot() +
  geom_line(aes(x = T, y = p, color = nu)) +
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

size_plot | prob_plot
