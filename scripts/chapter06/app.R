source("common_functions.R")
library(shiny)

# ----------------------------------------------------------------------
# UI
# ----------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Spatiotemporal Traffic State: Isotropic vs Adaptive Smoothing"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h5("Synthetic data"),
      sliderInput("det_spacing", "Detector spacing [km]",
                  min = 0.5, max = 3.0, value = 1.0, step = 0.25),
      sliderInput("noise_sd", "Measurement noise sd [km/h]",
                  min = 0, max = 10, value = 2.5, step = 0.5),
      actionButton("reseed", "Resample noise"),
      hr(),

      h5("Smoothing method"),
      radioButtons(
        "method", NULL,
        choices = c(
          "Naive isotropic  (Eq. 6.1-6.3)" = "isotropic",
          "Adaptive / ASM  (Eq. 6.4-6.7)"  = "adaptive"
        ),
        selected = "adaptive"
      ),

      conditionalPanel(
        condition = "input.method == 'isotropic'",
        sliderInput("iso_sigma", "\u03c3 \u2013 spatial width [km]",
                    min = 0.1, max = 3.0, value = 0.5, step = 0.05),
        sliderInput("iso_tau", "\u03c4 \u2013 temporal width [min]",
                    min = 0.1, max = 5.0, value = 0.75, step = 0.05)
      ),
      conditionalPanel(
        condition = "input.method == 'adaptive'",
        sliderInput("ad_sigma", "\u03c3 \u2013 spatial width [km]",
                    min = 0.1, max = 3.0, value = 0.5, step = 0.05),
        sliderInput("ad_tau", "\u03c4 \u2013 temporal width [min]",
                    min = 0.1, max = 5.0, value = 0.75, step = 0.05),
        sliderInput("c_free", "c_free \u2013 free-flow propagation speed [km/h]",
                    min = 20, max = 120, value = 70, step = 5),
        sliderInput("w_speed", "w \u2013 congested (jam) propagation speed [km/h]",
                    min = -30, max = -5, value = -16, step = 1),
        sliderInput("Vc", "V_c \u2013 free/congested transition center [km/h]",
                    min = 20, max = 100, value = 60, step = 5),
        sliderInput("dV", "\u0394V \u2013 transition width [km/h]",
                    min = 1, max = 30, value = 10, step = 1)
      ),

      hr(),
      radioButtons(
        "render_style", "Bottom-panel rendering",
        choices = c(
          "Raster  (geom_raster, interpolate=TRUE)" = "raster",
          "Filled contour  (geom_contour_filled)"   = "contour"
        ),
        selected = "raster"
      )
    ),

    mainPanel(
      width = 9,
      h4("Raw detector data -- cf. Fig. 6.1 (top)"),
      plotOutput("plot_top", height = "420px"),
      h4("Reconstructed speed field"),
      plotOutput("plot_bottom", height = "420px")
    )
  )
)

# ----------------------------------------------------------------------
# Server
# ----------------------------------------------------------------------
server <- function(input, output, session) {

  detectors <- reactive({
    seq(472, 480, by = input$det_spacing)
  })

  seed_val <- reactiveVal(7)
  observeEvent(input$reseed, {
    seed_val(sample.int(1e6, 1))
  })

  detector_data <- reactive({
    generate_detector_data(
      detectors = detectors(),
      noise_sd = input$noise_sd,
      seed = seed_val()
    )
  })

  guides_df_r <- reactive({ wave_guide_lines() })

  grid <- reactive({
    default_grid(detectors = detectors(), dx = 0.08, dt = 1.0)
  })

  output$plot_top <- renderPlot({
    plot_raw_detectors(
      detector_data(), guides_df = guides_df_r(), t_range = DEFAULT_T_RANGE,
      title = sprintf("%d detectors, %.2f km spacing, noise sd=%.1f km/h",
                       length(detectors()), input$det_spacing, input$noise_sd)
    )
  })

  output$plot_bottom <- renderPlot({
    g <- grid()

    if (input$method == "isotropic") {
      field_df <- isotropic_smoothing(
        detector_data(), g$x_grid, g$t_grid,
        sigma = input$iso_sigma, tau = input$iso_tau
      )
      title <- sprintf("Naive isotropic  |  \u03c3=%.2f km, \u03c4=%.2f min",
                        input$iso_sigma, input$iso_tau)
    } else {
      field_df <- adaptive_smoothing(
        detector_data(), g$x_grid, g$t_grid,
        sigma = input$ad_sigma, tau = input$ad_tau,
        c_free = input$c_free, w_speed = input$w_speed,
        V_c = input$Vc, delta_V = input$dV
      )
      title <- sprintf(
        "Adaptive (ASM)  |  \u03c3=%.2f km, \u03c4=%.2f min, c_free=%d km/h, w=%d km/h",
        input$ad_sigma, input$ad_tau, input$c_free, input$w_speed
      )
    }

    if (input$render_style == "raster") {
      plot_speed_field(field_df, t_range = DEFAULT_T_RANGE, title = title)
    } else {
      plot_speed_field_contour(field_df, t_range = DEFAULT_T_RANGE, title = title)
    }
  })
}

shinyApp(ui, server)