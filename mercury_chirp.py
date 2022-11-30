"""main."""
import numpy as np
from numpy import linalg as LA

from particle_class import Particle
from constants import SUN_MASS, MERCURY_MASS, MERCURY_PEREHELION_DISTANCE, MERCURY_PEREHELION_SPEED, YEAR, MERCURY_YEAR
from functions import main
from functions_plotting import grav_wave_freq_plot, trajectory_plot, distance_plot

mercury = Particle(MERCURY_MASS, np.array([MERCURY_PEREHELION_DISTANCE, 0.]), np.array([0., MERCURY_PEREHELION_SPEED]))

modifier         = "GR"  # Can be GR or Newtonian
scale_factor     = 50
central_mass     = SUN_MASS
orbiting_mass    = MERCURY_MASS
initial_position = np.array([MERCURY_PEREHELION_DISTANCE, 0.])  # start at perihelion
initial_velocity = np.array([0., MERCURY_PEREHELION_SPEED])  # speed at perihelion
earth_years      = 10.**3
number_of_orbits = int(np.ceil(earth_years * (YEAR / MERCURY_YEAR)))  # i.e. mercury_years
orbits           = np.arange(0, earth_years)

position, periapsides, grav_wave_freq_list = main(
    central_mass, orbiting_mass,
    initial_position, initial_velocity,
    number_of_orbits, modifier, scale_factor
)

grav_wave_freq_plot(grav_wave_freq_list[1:])
trajectory_plot(position, str(number_of_orbits))
distance_plot(LA.norm(periapsides, axis=1))
