"""main."""
import time
import numpy as np
from numpy import linalg as LA

from constants import SUN_MASS, MERCURY_MASS, MERCURY_PEREHELION_DISTANCE, MERCURY_PEREHELION_SPEED, YEAR, MERCURY_YEAR
from functions import main
import functions_plotting as fp

scale_factor     = 50
central_mass     = SUN_MASS
orbiting_mass    = MERCURY_MASS
initial_position = np.array([MERCURY_PEREHELION_DISTANCE, 0.])  # start at perihelion
initial_velocity = np.array([0., MERCURY_PEREHELION_SPEED])  # speed at perihelion
earth_years      = 10**5
number_of_orbits = int(np.ceil(earth_years * (YEAR / MERCURY_YEAR)))  # i.e. mercury_years

start = time.time()

periapsides, periapsis_angle_list, grav_wave_freq_list = main(
    central_mass, orbiting_mass,
    initial_position, initial_velocity,
    number_of_orbits, scale_factor
)

end = time.time()
print("Simulation runtime: " + str(np.around((end - start) / 3600, 2)) + " hours")  # seconds -> hours

fp.grav_wave_freq_plot(grav_wave_freq_list)
fp.trajectory_plot(periapsides)
fp.periapsis_angle_plot(periapsis_angle_list)
fp.distance_plot(LA.norm(periapsides, axis=1))
