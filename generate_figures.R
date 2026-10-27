library(ggplot2)

alpha_star = exp(-1)

# SIS 1d ---------------------------------------------------------------

load("SIS_high1_final_df.RData")

ggplot(df, aes(x = log10(delta), y = M, z = ESS_min / M)) +
  geom_contour_filled(bins = 20) +
  geom_contour(color = "black") +
  geom_contour(
    aes(z = alpha_ref1),
    breaks = alpha_star,
    color = "red",
    linetype = "dashed",
    linewidth = 1.2
  ) +
  scale_fill_viridis_d(
    guide = guide_coloursteps(
      show.limits = FALSE,
      ticks = TRUE
    ),
    labels = function(x) {
      keep <- round(seq(1, length(x), length.out = 5))
      ifelse(seq_along(x) %in% keep, x, "")
    }
  ) +
  labs(
    x = expression(log[10]~delta),
    y = "M",
    fill = "sample efficiency"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  )

# SIS 3d ---------------------------------------------------------------

load("SIS_high3_final_df.RData")

ggplot(df, aes(x = log10(delta), y = M, z = ESS_min / M)) +
  geom_contour_filled(bins = 20) +
  geom_contour(color = "black") +
  geom_contour(
    aes(z = alpha_ref1),
    breaks = alpha_star,
    color = "red",
    linetype = "dashed",
    linewidth = 1.2
  ) +
  scale_fill_viridis_d(
    guide = guide_coloursteps(
      show.limits = FALSE,
      ticks = TRUE
    ),
    labels = function(x) {
      keep <- round(seq(1, length(x), length.out = 5))
      ifelse(seq_along(x) %in% keep, x, "")
    }
  ) +
  labs(
    x = expression(log[10]~delta),
    y = "M",
    fill = "sample efficiency"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  )

# SIS 5d ---------------------------------------------------------------

load("SIS_high5_final_df.RData")

ggplot(df, aes(x = log10(delta), y = M, z = ESS_min / M)) +
  geom_contour_filled(bins = 20) +
  geom_contour(color = "black") +
  geom_contour(
    aes(z = alpha_ref1),
    breaks = alpha_star,
    color = "red",
    linetype = "dashed",
    linewidth = 1.2
  ) +
  scale_fill_viridis_d(
    guide = guide_coloursteps(
      show.limits = FALSE,
      ticks = TRUE
    ),
    labels = function(x) {
      keep <- round(seq(1, length(x), length.out = 5))
      ifelse(seq_along(x) %in% keep, x, "")
    }
  ) +
  labs(
    x = expression(log[10]~delta),
    y = "M",
    fill = "sample efficiency"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  )

# SIS 7d ---------------------------------------------------------------

load("SIS_high7_final_df.RData")

ggplot(df, aes(x = log10(delta), y = M, z = ESS_min / M)) +
  geom_contour_filled(bins = 20) +
  geom_contour(color = "black") +
  geom_contour(
    aes(z = alpha_ref1),
    breaks = alpha_star,
    color = "red",
    linetype = "dashed",
    linewidth = 1.2
  ) +
  scale_fill_viridis_d(
    guide = guide_coloursteps(
      show.limits = FALSE,
      ticks = TRUE
    ),
    labels = function(x) {
      keep <- round(seq(1, length(x), length.out = 5))
      ifelse(seq_along(x) %in% keep, x, "")
    }
  ) +
  labs(
    x = expression(log[10]~delta),
    y = "M",
    fill = "sample efficiency"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  )

# four coupled SIR ---------------------------------------------------------------

load("SIR_high4_final_df.RData")

ggplot(df, aes(x = log10(delta), y = M, z = ESS_min / M)) +
  geom_contour_filled(bins = 20) +
  geom_contour(color = "black") +
  scale_fill_viridis_d(
    guide = guide_coloursteps(
      show.limits = FALSE,
      ticks = TRUE
    ),
    labels = function(x) {
      keep <- round(seq(1, length(x), length.out = 5))
      ifelse(seq_along(x) %in% keep, x, "")
    }
  ) +
  labs(
    x = expression(log[10]~delta),
    y = "M",
    fill = "sample efficiency"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 20),
    plot.title = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  )

ggplot(df, aes(x = log10(delta), y = M, z = ESS_I_min / M)) +
  geom_contour_filled(bins = 20) +
  geom_contour(color = "black") +
  geom_contour(
    aes(z = alpha_ref1),
    breaks = alpha_star,
    color = "red",
    linetype = "dashed",
    linewidth = 1.2
  ) +
  scale_fill_viridis_d(
    guide = guide_coloursteps(
      show.limits = FALSE,
      ticks = TRUE
    ),
    labels = function(x) {
      keep <- round(seq(1, length(x), length.out = 5))
      ifelse(seq_along(x) %in% keep, x, "")
    }
  ) +
  labs(
    x = expression(log[10]~delta),
    y = "M",
    fill = "sample efficiency"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 20),
    plot.title = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  )

ggplot(df, aes(x = log10(delta), y = M, z = ESS_t1_min / M)) +
  geom_contour_filled(bins = 20) +
  geom_contour(color = "black") +
  geom_contour(
    aes(z = alpha_ref1),
    breaks = alpha_star,
    color = "red",
    linetype = "dashed",
    linewidth = 1.2
  ) +
  scale_fill_viridis_d(
    guide = guide_coloursteps(
      show.limits = FALSE,
      ticks = TRUE
    ),
    labels = function(x) {
      keep <- round(seq(1, length(x), length.out = 5))
      ifelse(seq_along(x) %in% keep, x, "")
    }
  ) +
  labs(
    x = expression(log[10]~delta),
    y = "M",
    fill = "sample efficiency"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 20),
    plot.title = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  )


# four coupled SIR with BS ------------------------------------------------

load("SIR_BS_high4_final_df.RData")

ggplot(df, aes(x = log10(delta), y = M, z = ESS_min / M)) +
  geom_contour_filled(bins = 20) +
  geom_contour(color = "black") +
  scale_fill_viridis_d(
    guide = guide_coloursteps(
      show.limits = FALSE,
      ticks = TRUE
    ),
    labels = function(x) {
      keep <- round(seq(1, length(x), length.out = 5))
      ifelse(seq_along(x) %in% keep, x, "")
    }
  ) +
  labs(
    x = expression(log[10]~delta),
    y = "M",
    fill = "sample efficiency"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 20),
    plot.title = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  )

# four coupled SEIR ---------------------------------------------------------------

load("SEIR_high4_final_df.RData")

# since the xPG chains when M=10 and delta=0.01 exhibit poor mixing and 
# do not provide sufficiently reliable approximations to the posterior 
# distribution within the simulation length considered,
# we change the ESS_min at this point to a smaller value

df <- df %>%
  mutate(
    ESS_min = if_else(
      M == 10 & abs(delta - 0.01) < 1e-10,
      100,
      ESS_min
    )
  )

ggplot(df, aes(x = log10(delta), y = M, z = ESS_min / M)) +
  geom_contour_filled(bins = 20) +
  geom_contour(color = "black") +
  scale_fill_viridis_d(
    guide = guide_coloursteps(
      show.limits = FALSE,
      ticks = TRUE
    ),
    labels = function(x) {
      keep <- round(seq(1, length(x), length.out = 5))
      ifelse(seq_along(x) %in% keep, x, "")
    }
  ) +
  labs(
    x = expression(log[10]~delta),
    y = "M",
    fill = "sample efficiency"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 20),
    plot.title = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  )

ggplot(df, aes(x = log10(delta), y = M, z = ESS_I_min / M)) +
  geom_contour_filled(bins = 20) +
  geom_contour(color = "black") +
  geom_contour(
    aes(z = alpha_ref1),
    breaks = alpha_star,
    color = "red",
    linetype = "dashed",
    linewidth = 1.2
  ) +
  scale_fill_viridis_d(
    guide = guide_coloursteps(
      show.limits = FALSE,
      ticks = TRUE
    ),
    labels = function(x) {
      keep <- round(seq(1, length(x), length.out = 5))
      ifelse(seq_along(x) %in% keep, x, "")
    }
  ) +
  labs(
    x = expression(log[10]~delta),
    y = "M",
    fill = "sample efficiency"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 20),
    plot.title = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  )

ggplot(df, aes(x = log10(delta), y = M, z = ESS_t1_min / M)) +
  geom_contour_filled(bins = 20) +
  geom_contour(color = "black") +
  geom_contour(
    aes(z = alpha_ref1),
    breaks = alpha_star,
    color = "red",
    linetype = "dashed",
    linewidth = 1.2
  ) +
  scale_fill_viridis_d(
    guide = guide_coloursteps(
      show.limits = FALSE,
      ticks = TRUE
    ),
    labels = function(x) {
      keep <- round(seq(1, length(x), length.out = 5))
      ifelse(seq_along(x) %in% keep, x, "")
    }
  ) +
  labs(
    x = expression(log[10]~delta),
    y = "M",
    fill = "sample efficiency"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 20),
    plot.title = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  )


# four coupled SEIR with BS -----------------------------------------------

load("SEIR_BS_high4_final_df.RData")

# since the Markov chains when (M, delta)=(1, 0.01) and (M, delta)=(5, 0.01) are not mixing well, 
# we change the ESS_min at this point to a smaller value 

df <- df %>%
  mutate(
    ESS_min = if_else(
      M == 1 & abs(delta - 0.01) < 1e-10,
      30,
      ESS_min
    )
  )

df <- df %>%
  mutate(
    ESS_min = if_else(
      M == 5 & abs(delta - 0.01) < 1e-10,
      180,
      ESS_min
    )
  )

ggplot(df, aes(x = log10(delta), y = M, z = ESS_min / M)) +
  geom_contour_filled(bins = 20) +
  geom_contour(color = "black") +
  scale_fill_viridis_d(
    guide = guide_coloursteps(
      show.limits = FALSE,
      ticks = TRUE
    ),
    labels = function(x) {
      keep <- round(seq(1, length(x), length.out = 5))
      ifelse(seq_along(x) %in% keep, x, "")
    }
  ) +
  labs(
    x = expression(log[10]~delta),
    y = "M",
    fill = "sample efficiency"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 20),
    plot.title = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  )


# Autoregulatory model ----------------------------------------------------

load("Autoreg_df.RData")

ggplot(df, aes(x = log10(delta), y = ESS_X4_t0/M, color = factor(M))) +
  geom_point(size = 3) +  # Scatter plot points
  geom_line() +  # Connect points with lines
  # geom_text(aes(label = M), vjust = -0.5, size = 4.5) +   # bigger labels
  labs(
    # title = "sample efficiency vs. log10(delta)",
    x = expression(log[10]~delta),
    y = "sample efficiency",
    color = "M"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 22),   # bigger axis labels
    axis.text  = element_text(size = 20),                  # bigger tick labels
    legend.title = element_text(size = 20), # bigger legend title
    legend.text  = element_text(size = 20),                # bigger legend text
    plot.title   = element_text(size = 22),  # bigger title
    theme(legend.position = "right")
  )


# posteriors of model parameters for four coupled SEIR --------------------

load("SEIR_high4_final_seed3_Thetas_delta_0.02_M30_double.RData")

load("SEIR_high4_final_seed1_Thetas_delta_0.02_M30_double.RData")

load("SEIR_high4_BS_final_seed3_Thetas_delta_0.02_M10_double.RData")

load("SEIR_high4_BS_final_seed1_Thetas_delta_0.02_M10_double.RData")

# Iterations to use
iter <- (10^2):(2*10^5)

par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

theta_labels <- c(
  expression(beta),
  expression(eta),
  expression(gamma)
)

for (i in 1:3) {
  # Seed 3
  xpg_delta_0.01_M10_seed3 <-
    Thetas_delta_0.01_M10_seed3[i, iter]
  
  xpgbs_delta_0.02_M10_seed3 <-
    Thetas_delta_0.02_M10_seed3_AS[i, iter]
  
  # Seed 1
  xpg_delta_0.01_M10_seed1 <-
    Thetas_delta_0.01_M10_seed1[i, iter]
  
  xpgbs_delta_0.02_M10_seed1 <-
    Thetas_delta_0.02_M10_seed1_AS[i, iter]
  
  # ESS
  
  ess_xpg_seed3 <- as.numeric(
    effectiveSize(xpg_delta_0.01_M10_seed3)
  )
  
  ess_xpgbs_seed3 <- as.numeric(
    effectiveSize(xpgbs_delta_0.02_M10_seed3)
  )
  
  ess_xpg_seed1 <- as.numeric(
    effectiveSize(xpg_delta_0.01_M10_seed1)
  )
  
  ess_xpgbs_seed1 <- as.numeric(
    effectiveSize(xpgbs_delta_0.02_M10_seed1)
  )
  
  # Density estimates
  
  d_xpg_seed3 <- density(
    xpg_delta_0.01_M10_seed3,
    na.rm = TRUE
  )
  
  d_xpgbs_seed3 <- density(
    xpgbs_delta_0.02_M10_seed3,
    na.rm = TRUE
  )
  
  d_xpg_seed1 <- density(
    xpg_delta_0.01_M10_seed1,
    na.rm = TRUE
  )
  
  d_xpgbs_seed1 <- density(
    xpgbs_delta_0.02_M10_seed1,
    na.rm = TRUE
  )
  
  # Common x- and y-axis ranges
  
  xlim_range <- range(
    d_xpg_seed3$x,
    d_xpgbs_seed3$x,
    d_xpg_seed1$x,
    d_xpgbs_seed1$x
  )
  
  ylim_range <- c(
    0,
    1.1 * max(
      d_xpg_seed3$y,
      d_xpgbs_seed3$y,
      d_xpg_seed1$y,
      d_xpgbs_seed1$y
    )
  )

  # Plot seed 3 xPG
  
  plot(
    d_xpg_seed3,
    col = "blue",
    lwd = 2,
    lty = 1,
    main = theta_labels[[i]],
    xlab = theta_labels[[i]],
    ylab = "Density",
    xlim = xlim_range,
    ylim = ylim_range
  )
  
  
  # PG, seed 3
  lines(
    d_xpgbs_seed3,
    col = "red",
    lwd = 2,
    lty = 1
  )
  
  
  # xPG, seed 1
  lines(
    d_xpg_seed1,
    col = "blue",
    lwd = 2,
    lty = 2
  )
  
  
  # PG, seed 1
  lines(
    d_xpgbs_seed1,
    col = "red",
    lwd = 2,
    lty = 2
  )
  
  # Legend
  
  legend(
    "topright",
    legend = c(
      paste0(
        "xPG: delta = 0.01, seed 3, ESS = ",
        round(ess_xpg_seed3, 1)
      ),
      paste0(
        "xPGBS: delta = 0.02, seed 3, ESS = ",
        round(ess_xpgbs_seed3, 1)
      ),
      paste0(
        "xPG: delta = 0.01, seed 1, ESS = ",
        round(ess_xpg_seed1, 1)
      ),
      paste0(
        "xPGBS: delta = 0.02, seed 1, ESS = ",
        round(ess_xpgbs_seed1, 1)
      )
    ),
    col = c(
      "blue",
      "red",
      "blue",
      "red"
    ),
    lwd = 2,
    lty = c(
      1,
      1,
      2,
      2
    ),
    bty = "n",
    cex = 0.75
  )
}

par(mfrow = c(1, 1))

# Inferred infectious state trace -----------------------------------------

load("SEIR_high4_eta_0.6_tau_0.1_seed3_Paths.RData")

load("SEIR_high4_BS_eta_0.6_tau_0.1_seed3_Paths.RData")


time_grid <- (0:200) * 0.1

make_plot <- function(d) {
  
  Path <- Paths_delta_0.05_M30_seed3[d, , ] # comment it and uncomment the next line to generate the plot for xPGBS
  # Path <- Paths_delta_0.02_M10_seed3_BS[d, , ]
  
  q_path <- apply(
    Path, 1, quantile,
    probs = c(0.05, 0.5, 0.95),
    na.rm = TRUE
  )
  
  df_q <- data.frame(
    time = time_grid,
    q05 = q_path[1, ],
    q50 = q_path[2, ],
    q95 = q_path[3, ]
  )
  
  df_ref <- data.frame(
    time = time_grid,
    x_ref = X.ref[d, ]
  )
  
  df_obs <- data.frame(
    time = 1:20,
    y = ys[d - 8, ]
  )
  
  ggplot(df_q, aes(x = time)) +
    geom_ribbon(
      aes(ymin = q05, ymax = q95),
      alpha = 0.3
    ) +
    geom_line(
      aes(
        y = q05,
        linetype = "5% and 95% quantiles",
        color = "5% and 95% quantiles"
      ),
      linewidth = 0.8
    ) +
    geom_line(
      aes(
        y = q95,
        linetype = "5% and 95% quantiles",
        color = "5% and 95% quantiles"
      ),
      linewidth = 0.8
    ) +
    geom_line(
      aes(
        y = q50,
        linetype = "50% quantile",
        color = "50% quantile"
      ),
      linewidth = 0.8
    ) +
    geom_line(
      data = df_ref,
      aes(
        x = time,
        y = x_ref,
        linetype = "Truth",
        color = "Truth"
      ),
      inherit.aes = FALSE,
      linewidth = 1
    ) +
    geom_point(
      data = df_obs,
      aes(x = time, y = y),
      inherit.aes = FALSE,
      shape = 4,
      size = 3,
      stroke = 1.2
    ) +
    scale_color_manual(
      name = NULL,
      values = c(
        "Truth" = "blue",
        "5% and 95% quantiles" = "black",
        "50% quantile" = "red"
      )
    ) +
    scale_linetype_manual(
      name = NULL,
      values = c(
        "Truth" = "solid",
        "5% and 95% quantiles" = "dashed",
        "50% quantile" = "dashed"
      )
    ) +
    labs(
      title = paste0("region ", d-8),
      x = "Time",
      y = "Infectious state"
    ) +
    theme_minimal()
}

p9  <- make_plot(9)
p10 <- make_plot(10)
p11 <- make_plot(11)
p12 <- make_plot(12)


p_all <- (p9 + p10) / (p11 + p12) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 13),
    legend.key.size = unit(1.2, "cm")
  ) &
  guides(
    color = guide_legend(
      override.aes = list(linewidth = 1.5, size = 4)
    ),
    linetype = guide_legend(
      override.aes = list(linewidth = 1.5)
    )
  )

p_all


# SIR high4 acc -----------------------------------------------------------

load("SIR_high4_final_df.RData")

df_sub1 <- df[df$delta<=0.1, ]

delta_levels <- sort(unique(df_sub1$delta))

delta_cols <- setNames(
  scales::hue_pal()(length(delta_levels)),
  delta_levels
)

df_sub1$delta <- factor(df_sub1$delta, levels = delta_levels)

p <- ggplot(df_sub1, aes(x = alpha_ref1, y = ESS_min/M, color = delta)) +
  geom_point(size = 3) +  # Scatter plot points
  geom_line() +  # Connect points with lines
  # geom_text(aes(label = delta), vjust = -0.5, size = 3) + 
  geom_vline(xintercept = exp(-1), linetype = "dashed", color = "red")+
  scale_color_manual(values = delta_cols)+
  labs(
    # title = "sample efficiency vs. acc",
    x = expression(alpha),
    y = "sample efficiency",
    color = expression(delta)
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  )
theme(legend.position = "right")
ggsave("/home/yangl35/h-drive/xPG_MJP/SIS_SIR_high/high4_SIR_acc.png", plot = p, width = 13, height = 6, dpi = 300)

# SEIR high4 acc ----------------------------------------------------------

load("SEIR_high4_final_df.RData")

df_sub1 <- df[df$delta<=0.1, ]

delta_levels <- sort(unique(df_sub1$delta))

delta_cols <- setNames(
  scales::hue_pal()(length(delta_levels)),
  delta_levels
)

df_sub1$delta <- factor(df_sub1$delta, levels = delta_levels)

p <- ggplot(df_sub1, aes(x = alpha_ref1, y = ESS_min/M, color = delta)) +
  geom_point(size = 3) +  # Scatter plot points
  geom_line() +  # Connect points with lines
  # geom_text(aes(label = delta), vjust = -0.5, size = 3) + 
  geom_vline(xintercept = exp(-1), linetype = "dashed", color = "red")+
  scale_color_manual(values = delta_cols)+
  labs(
    # title = "sample efficiency vs. acc",
    x = expression(alpha),
    y = "sample efficiency",
    color = expression(delta)
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  )+
  theme(legend.position = "right")

