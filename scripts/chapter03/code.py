
#Data: https://data.transportation.gov/Automobiles/Next-Generation-Simulation-NGSIM-Vehicle-Trajector/8ect-6jqj/about_data
import pandas as pd
from plotnine import *

# 1. Load and Preprocess
traj = pd.read_csv("data/I80/0400pm-0415pm/trajectories-0400-0415.csv")

# Filter for Lane 5 and convert units
# Local_Y (ft) -> Meters: 0.3048
# v_Vel (ft/s) -> km/h: 1.09728
traj5 = traj[traj['Lane_ID'] == 5].copy()
traj5['Location'] = traj5['Local_Y'] * 0.3048
traj5['Speed_kmh'] = traj5['v_Vel'] * 1.09728

# Normalize Time to "Seconds after 4:00 PM"
t_min = traj5['Global_Time'].min()
traj5['Time_s'] = (traj5['Global_Time'] - t_min) / 1000.0

# 2. Build the Plot
st_plot = (
    ggplot(traj5, aes(x='Time_s', y='Location', group='Vehicle_ID', color='Speed_kmh'))
    + geom_path(size=0.4) 
    + scale_color_cmap(cmap_name='nipy_spectral', limits=(0, 40))
    + scale_x_continuous(limits=(350, 950), breaks=range(480, 960, 120))
    + scale_y_continuous(limits=(50, 480))
    + labs(
        x='Time [s] after 4:00 pm',
        y='Location [m]',
        color='speed [km/h]'
    )
    + annotate("rect", xmin=370, xmax=550, ymin=380, ymax=460, 
               fill="white", color="black", size=0.5)
    + annotate("text", x=460, y=420, label="NGSIM I-80\nLane 5", 
               ha='center', va='center', size=10)
    + theme_bw()
    + theme(
        figure_size=(10, 6),
        panel_grid_major=element_blank(),
        panel_grid_minor=element_blank(),
        legend_position=(0.92, 0.5),
        legend_direction='vertical'
    )
)

st_plot
