library(ggplot2)


# ----------------------------------------------------------------------
# 0. Colour scale matching the book's red -> orange -> yellow -> green
#    -> blue -> purple "Speed [km/h]" colour bar (0 .. 120 km/h)
# ----------------------------------------------------------------------
SPEED_LIMITS <- c(0, 120)
SPEED_COLORS <- c(
  "#B2182B",  # 0   km/h  red
  "#EF6C25",  # 20  km/h  orange
  "#F9D423",  # 40  km/h  yellow
  "#3E9F3E",  # 60  km/h  green
  "#2E6FCC",  # 80  km/h  blue
  "#7B2FA3",  # 100 km/h  purple
  "#7B2FA3"   # 120 km/h  purple (saturate)
)

speed_scale_color <- function(name = "Speed [km/h]") {
  scale_color_gradientn(colors = SPEED_COLORS, limits = SPEED_LIMITS, name = name)
}

speed_scale_fill <- function(name = "Speed [km/h]") {
  scale_fill_gradientn(colors = SPEED_COLORS, limits = SPEED_LIMITS, name = name)
}

speed_scale_fill_binned <- function(breaks, name = "Speed [km/h]") {
  n <- length(breaks) - 1
  ramp <- colorRampPalette(SPEED_COLORS)(n)
  scale_fill_manual(values = ramp, name = name)
}

# ----------------------------------------------------------------------
# 1. Synthetic "ground truth" stop-and-go field + detector sampling
# ----------------------------------------------------------------------
DEFAULT_DETECTORS <- seq(472, 480, by = 1)          # 472..480 km, 1 km spacing
DEFAULT_T_RANGE <- c(7 * 60 + 40, 8 * 60 + 55)       # minutes after midnight

true_speed_field <- function(x, t,
                              x_top = 477.0, x_bottom = 472.0,
                              wave_speed = -18.5, period = 9.5,
                              t_first_onset = 8 * 60 - 3, n_waves = 8,
                              v_free = 100.0, v_jam = 8.0,
                              pulse_half_width = 0.55) {
  V <- rep(v_free, length(x))
  V <- V + 3.0 * sin(t / 11.0 + x * 0.7)

  for (k in 0:(n_waves - 1)) {
    t0 <- t_first_onset + k * period
    x_center <- x_top + (wave_speed / 60.0) * (t - t0)
    depth_from_top <- pmin(pmax((x_top - x_center) / 1.5, 0.0), 1.0)
    in_range <- (x_center <= x_top + 0.3) & (x_center >= x_bottom - 0.3)
    pulse <- exp(-0.5 * ((x - x_center) / pulse_half_width)^2)
    drop <- (v_free - v_jam) * depth_from_top * pulse * in_range
    V <- V - drop
  }

  pmin(pmax(V, v_jam * 0.8), v_free + 6.0)
}

generate_detector_data <- function(detectors = DEFAULT_DETECTORS,
                                    t_range = DEFAULT_T_RANGE,
                                    dt = 1.0, noise_sd = 2.5, seed = 7,
                                    ...) {
  set.seed(seed)
  t_grid <- seq(t_range[1], t_range[2], by = dt)
  grid <- expand.grid(x = detectors, t = t_grid)
  grid$speed <- true_speed_field(grid$x, grid$t, ...) +
    rnorm(nrow(grid), mean = 0, sd = noise_sd)
  grid$speed <- pmin(pmax(grid$speed, 0), 130)
  grid[order(grid$x, grid$t), ]
}

wave_guide_lines <- function(x_top = 477.0, x_bottom = 472.0,
                              wave_speed = -18.5, period = 9.5,
                              t_first_onset = 8 * 60 - 3, n_waves = 8) {
  rows <- lapply(0:(n_waves - 1), function(k) {
    t0 <- t_first_onset + k * period
    t1 <- t0 + (x_bottom - x_top) / (wave_speed / 60.0)
    data.frame(wave = k, x = x_top, t = t0, xend = x_bottom, tend = t1)
  })
  do.call(rbind, rows)
}

# ----------------------------------------------------------------------
# 2. Isotropic smoothing -- Eqs. (6.1)-(6.3)
# ----------------------------------------------------------------------
#' Naive isotropic spatiotemporal interpolation.
#'
#'   phi_0(x - x_i, t - t_i) = exp[-(|x-x_i|/sigma + |t-t_i|/tau)]   (6.2)
#'   V(x,t) = sum_i phi_0 * v_i  /  sum_i phi_0                     (6.1),(6.3)
#'
#' Returns a tidy data.frame with columns x, t, speed on the full grid
#' (every combination of x_grid and t_grid).
isotropic_smoothing <- function(df, x_grid, t_grid, sigma, tau) {
  grid <- expand.grid(x = x_grid, t = t_grid)
  num <- numeric(nrow(grid))
  den <- numeric(nrow(grid))

  for (k in seq_len(nrow(df))) {
    w <- exp(-(abs(grid$x - df$x[k]) / sigma + abs(grid$t - df$t[k]) / tau))
    num <- num + w * df$speed[k]
    den <- den + w
  }

  grid$speed <- num / pmax(den, 1e-12)
  grid
}

# ----------------------------------------------------------------------
# 3. Adaptive smoothing method (ASM) -- Eqs. (6.4)-(6.7)
# ----------------------------------------------------------------------
#' One directional branch of the ASM: Eq. (6.4) if c=c_free (free
#' traffic, perturbations advected downstream), or Eq. (6.5) if c=w
#' (congested traffic, jams propagate upstream, so c is negative).
#'
#'   V_c(x,t) = sum_i phi_0(x-x_i, t-t_i-(x-x_i)/c) v_i / N_c(x,t)
#'
#' NOTE: `c` (c_free or w) is given in km/h, but `t`/`tau` are in
#' minutes -- c is converted to km/min before use.
directional_branch <- function(df, grid, sigma, tau, c) {
  c_km_per_min <- c / 60.0
  num <- numeric(nrow(grid))
  den <- numeric(nrow(grid))

  for (k in seq_len(nrow(df))) {
    dx <- grid$x - df$x[k]
    shifted_dt <- grid$t - df$t[k] - dx / c_km_per_min
    w <- exp(-(abs(dx) / sigma + abs(shifted_dt) / tau))
    num <- num + w * df$speed[k]
    den <- den + w
  }

  num / pmax(den, 1e-12)
}

#' Adaptive Smoothing Method (ASM), Eqs. (6.4)-(6.7).
#'
#' sigma, tau   : smoothing widths [km], [min]
#' c_free       : propagation speed of perturbations in free traffic
#'                [km/h] (positive, downstream)
#' w_speed      : propagation speed of jam fronts in congested traffic
#'                [km/h] (negative, upstream)
#' V_c, delta_V : center and width of the free/congested transition,
#'                Eq. (6.7): w_cong = 1/2 [1 + tanh((V_c - V*)/deltaV)]
#'
#' Returns a tidy data.frame with columns x, t, speed, V_free, V_cong,
#' w_cong on the full grid.
adaptive_smoothing <- function(df, x_grid, t_grid,
                                sigma = 0.5, tau = 0.75,
                                c_free = 70.0, w_speed = -16.0,
                                V_c = 60.0, delta_V = 10.0) {
  grid <- expand.grid(x = x_grid, t = t_grid)

  V_free <- directional_branch(df, grid, sigma, tau, c_free)
  V_cong <- directional_branch(df, grid, sigma, tau, w_speed)

  # Eq. (6.7): w_cong(x,t) = 1/2 [1 + tanh((V_c - V*)/deltaV)], V* = min(V_free, V_cong)
  V_star <- pmin(V_free, V_cong)
  w_cong <- 0.5 * (1.0 + tanh((V_c - V_star) / delta_V))

  grid$speed <- w_cong * V_cong + (1.0 - w_cong) * V_free
  grid$V_free <- V_free
  grid$V_cong <- V_cong
  grid$w_cong <- w_cong
  grid
}

# ----------------------------------------------------------------------
# 4. Plotting -- Fig. 6.1 top (raw stacked detector traces)
# ----------------------------------------------------------------------
minute_to_hhmm <- function(m) {
  m <- round(m) %% (24 * 60)
  sprintf("%02d:%02d", m %/% 60, m %% 60)
}

time_axis <- function(t_range, step = 20) {
  breaks <- seq(ceiling(t_range[1] / step) * step, t_range[2], by = step)
  scale_x_continuous(breaks = breaks, labels = minute_to_hhmm(breaks))
}

# Recreates the top panel of Fig. 6.1: one near-horizontal strip per
# detector, colour-coded by speed, with a small vertical wiggle inside
# each detector's own strip that also encodes speed (matching the look
# of the original figure), plus optional diagonal wave guide lines.
plot_raw_detectors <- function(df, guides_df = NULL, strip_height = 0.42,
                                v_free = 100.0, v_jam = 8.0,
                                t_range = DEFAULT_T_RANGE,
                                title = "A5 South (synthetic), reconstructed like Fig. 6.1 (top)") {
  d <- df[order(df$x, df$t), ]
  norm <- pmin(pmax((d$speed - v_jam) / (v_free - v_jam), 0), 1)
  d$y <- d$x + (norm - 1.0) * strip_height

  seg_list <- lapply(split(d, d$x), function(g) {
    g <- g[order(g$t), ]
    n <- nrow(g)
    data.frame(
      x0 = g$t[-n], x1 = g$t[-1],
      y0 = g$y[-n], y1 = g$y[-1],
      speed = (g$speed[-n] + g$speed[-1]) / 2.0,
      detector = g$x[1]
    )
  })
  seg_df <- do.call(rbind, seg_list)

  p <- ggplot() +
    geom_hline(data = data.frame(x = unique(d$x)),
               mapping = aes(yintercept = x), color = "#888888", linewidth = 0.3) +
    geom_segment(data = seg_df,
                 mapping = aes(x = x0, xend = x1, y = y0, yend = y1, color = speed),
                 linewidth = 1.1)

  if (!is.null(guides_df) && nrow(guides_df) > 0) {
    p <- p + geom_segment(data = guides_df,
                           mapping = aes(x = t, xend = tend, y = x, yend = xend),
                           color = "black", linewidth = 0.5)
  }

  p +
    speed_scale_color() +
    time_axis(t_range) +
    scale_y_continuous(breaks = sort(unique(as.integer(d$x)))) +
    coord_cartesian(xlim = t_range) +
    labs(x = "Time", y = "Detector position [km]", title = title) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())
}

# ----------------------------------------------------------------------
# 5. Plotting -- Fig. 6.1 bottom (interpolated speed field)
# ----------------------------------------------------------------------
#' Recreates the bottom panel of Fig. 6.1 using geom_raster(interpolate =
#' TRUE), ggplot2's native smooth/interpolated raster -- no workaround
#' needed here, unlike plotnine.
plot_speed_field <- function(field_df, t_range = DEFAULT_T_RANGE,
                              title = "Reconstructed speed field",
                              interpolate = TRUE) {
  ggplot(field_df, aes(x = t, y = x, fill = speed)) +
    geom_raster(interpolate = interpolate) +
    speed_scale_fill() +
    time_axis(t_range) +
    scale_y_continuous(breaks = seq(floor(min(field_df$x)), ceiling(max(field_df$x)))) +
    coord_cartesian(xlim = t_range) +
    labs(x = "Time", y = "Location [km]", title = title) +
    theme_minimal() +
    theme(panel.grid = element_blank())
}

# Filled-contour version, using ggplot2's *native* geom_contour_filled()
# -- unlike plotnine (which has no contour geom at all and needs a
# matplotlib-polygon workaround), ggplot2 supports this directly.
plot_speed_field_contour <- function(field_df, t_range = DEFAULT_T_RANGE,
                                      title = "Reconstructed speed field (filled contour)",
                                      breaks = seq(0, 130, by = 10)) {
  ggplot(field_df, aes(x = t, y = x, z = speed)) +
    geom_contour_filled(breaks = breaks) +
    speed_scale_fill_binned(breaks) +
    time_axis(t_range) +
    scale_y_continuous(breaks = seq(floor(min(field_df$x)), ceiling(max(field_df$x)))) +
    coord_cartesian(xlim = t_range) +
    labs(x = "Time", y = "Location [km]", title = title) +
    theme_minimal() +
    theme(panel.grid = element_blank())
}

# ----------------------------------------------------------------------
# 6. Convenience: default grid for the reconstruction
# ----------------------------------------------------------------------
default_grid <- function(detectors = DEFAULT_DETECTORS, t_range = DEFAULT_T_RANGE,
                          dx = 0.1, dt = 1.0) {
  list(
    x_grid = seq(min(detectors), max(detectors), by = dx),
    t_grid = seq(t_range[1], t_range[2], by = dt)
  )
}
