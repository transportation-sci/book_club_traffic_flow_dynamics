library(tidyverse)

# Q1-----------------------------------------------------
df <- data.frame(
  Time = c(2, 7, 7, 10, 12, 18, 20, 21, 25, 29),
  Speed = c(26, 24, 32, 32, 29, 28, 34, 22, 26, 38),
  Lane = c(1, 1, 2, 2, 1, 1, 2, 1, 1, 2),
  VehicleLength = c(5, 12, 4, 5, 4, 4, 5, 15, 3, 5)
)

View(df)

res <- df |> 
  group_by(Lane) |> 
  summarise(
    V = mean(Speed),
    Q = n() * (60*60)/30
  ) |> 
  mutate(
    rho = Q/(V*3.6)
  )


rho_tot <- sum(res$rho)
Q_tot <- sum(res$Q)
V_s = Q_tot/rho_tot


# truck percentage
2/6 *100
2/10 * 100


# Q3---------------------------------------------------------
L = 5
V_0 = 33.33
s_0 = 2
T = 1.6
rho_max = 1/(L + s_0)
rho_critical = 1/((V_0 * T) + (L + s_0))
capacity <- (1/T) * (1/(1 + ((L + s_0)/(V_0 * T))))

rho_vals_vpk <- seq(0, 120, 20)
rho_vals_vpm = rho_vals_vpk/1000

## free traffic
Q_free <- function(rho, desired_speed = V_0) {
 rho * desired_speed
}

Q_free_vals <- Q_free(rho_vals_vpm[1:2]) * 60 * 60

## congested traffic
V_cong <- function(rho, min_gap = s_0, time_gap = T, len = L){
  (1/time_gap) * ((1/rho) - (len + min_gap))
}

Q_cong <- function(rho, min_gap = s_0, time_gap = T, len = L){
  (1/time_gap) * (1 - (rho * (len + min_gap)))
}

Q_cong_vals <- Q_cong(rho_vals_vpm[3:7]) * 60 * 60

Q_final <- c(Q_free_vals, Q_cong_vals)

V_cong_vals <- V_cong(rho_vals_vpm[3:7]) * 3.6

library(ggplot2)
theme_set(theme_classic())

flow_density <- ggplot(mapping = aes(rho_vals_vpk, Q_final)) +
  geom_point() +
  geom_smooth(
    method = "lm", 
    aes(rho_vals_vpk[2:7], Q_final[2:7]),
    se = FALSE
  ) +
  geom_smooth(
    method = "lm", 
    aes(rho_vals_vpk[1:2], Q_final[1:2]),
    se = FALSE
  ) +
  labs(x = "Density (veh/km)", y = "Flow (veh/hr)",
title = "Flow-Density Plot\n& Fundamental Diagram")


speed_density <- ggplot(mapping = aes(rho_vals_vpk, c(V_0, V_0, V_cong_vals)*3.6)) +
  geom_point() +
  geom_smooth(
    se = FALSE
  ) +
  labs(x = "Density (veh/km)", y = "Speed (km/h)",
title = "Speed-Density Plot")


library(patchwork)

flow_density | speed_density

# plot(rho_vals_vpk, Q_final)
# lines(rho_vals_vpk, Q_final)

# plot(rho_vals_vpk, c(V_0, V_0, V_cong_vals))
# lines(rho_vals_vpk, c(V_0, V_0, V_cong_vals))
