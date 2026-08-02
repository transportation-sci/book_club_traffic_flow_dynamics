source("scripts/chapter06/common_functions.R")

dir.create("scripts/chapter06/output", showWarnings = FALSE)

# ------------------------------------------------------------------
# 1. Synthetic data (stand-in for the real A5-South detector data)
# ------------------------------------------------------------------
detector_df <- generate_detector_data(
  detectors = DEFAULT_DETECTORS,
  t_range = DEFAULT_T_RANGE,
  dt = 1.0,
  noise_sd = 2.5
)
guides_df <- wave_guide_lines()


# ------------------------------------------------------------------
# 2. Fig. 6.1 (top): raw stacked traces
# ------------------------------------------------------------------
p_top <- plot_raw_detectors(
  detector_df, guides_df = guides_df, t_range = DEFAULT_T_RANGE,
  title = "Synthetic A5-South-like data -- cf. Fig. 6.1 (top)"
)
ggsave("scripts/chapter06/output/fig_6_1_top.png", p_top, width = 9, height = 6, dpi = 150)


# ------------------------------------------------------------------
# 3. Reconstructions: adaptive smoothing (Eqs. 6.4-6.7) and isotropic
#    smoothing (Eqs. 6.1-6.3), both on the same grid
# ------------------------------------------------------------------
grid <- default_grid(dx = 0.08, dt = 1.0)

adaptive_df <- adaptive_smoothing(
  detector_df, grid$x_grid, grid$t_grid,
  sigma = 0.5, tau = 0.75, c_free = 70.0, w_speed = -16.0,
  V_c = 60.0, delta_V = 10.0
)
# p_bottom_adaptive <- plot_speed_field(
#   adaptive_df, t_range = DEFAULT_T_RANGE,
#   title = "Adaptive smoothing (ASM) -- cf. Fig. 6.1 (bottom)"
# )
# ggsave("scripts/chapter06/output/fig_6_1_bottom_adaptive.png", p_bottom_adaptive, width = 9, height = 6, dpi = 150)


iso_df <- isotropic_smoothing(
  detector_df, grid$x_grid, grid$t_grid, sigma = 0.5, tau = 0.75
)
p_bottom_isotropic <- plot_speed_field(
  iso_df, t_range = DEFAULT_T_RANGE,
  title = "Naive isotropic smoothing -- cf. Fig. 6.3 (right / egg-carton)"
)
ggsave("scripts/chapter06/output/fig_6_1_bottom_isotropic.png", p_bottom_isotropic, width = 9, height = 6, dpi = 150)


# ------------------------------------------------------------------
# 4. Filled-contour rendering (native geom_contour_filled in ggplot2)
# ------------------------------------------------------------------
p_contour <- plot_speed_field_contour(
  adaptive_df, t_range = DEFAULT_T_RANGE,
  title = "Adaptive smoothing -- filled contour (geom_contour_filled)"
)
ggsave("scripts/chapter06/output/fig_6_1_bottom_adaptive_contour.png", p_contour, width = 9, height = 6, dpi = 150)
