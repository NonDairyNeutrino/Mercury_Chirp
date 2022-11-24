"""FUNCTIONS."""

import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
# import matplotlib
# from mpl_toolkits.mplot3d import Axes3D
# import multiprocessing as mp

from constants import *

# Schwarzchild radius
# Input: scalar mass of central body e.g. the Sun
# Output: Schwarzchild radius for given mass


def schwarz_radius(mass):
    r_s = 2 * NEWTON_G * mass / (LIGHT_SPEED**2)
    return r_s

# r_L from paper
# Input: 2D vectors of position and velocity
# Output: Scalar (? this might not be right) angular momentum radius


def ang_mom_radius_squared(position, velocity):
    r_L2 = np.sum(np.cross(position, velocity)**2) / (LIGHT_SPEED**2)
    return r_L2

# Modified gravity from paper
# Input: scalar mass of central body, vector position of planet, vector velocity of planet
# Output: vector acceleration due to modified gravity of the central body at that position


def gravity(mass, position, velocity, modifier):
    distance = LA.norm(position)  # This implicitly defines the sun as the origin
    if modifier == "GR":
        acceleration = -0.5 * LIGHT_SPEED**2 * \
            schwarz_radius(mass) * (1 + 3 * ang_mom_radius_squared(position, velocity) / distance**2) * position / (distance**3)
        # print(acceleration*YEAR**2/AU)
    elif modifier == "Newtonian":
        acceleration = -1*NEWTON_G*mass*position / distance**3
        # print(acceleration*YEAR**2/AU)
    else:
        print('no')
        quit()
    return acceleration

# Verlet method
# Input: current point, previous point, current second derivative, time step
# Output: the next point


def verlet(y_i, y_im1, current_second_derivative, time_step):
    y_ip1 = 2.*y_i - y_im1 + current_second_derivative*time_step**2
    return y_ip1

# Function that steps the simulation forward by one iteration
# Input: mass of central body, position array with at least 2 elements, velocity array
# Output: next position, next velocity


def step(mass, position, velocity, time_step, modifier):
    a_i = gravity(mass, position[-1], velocity[-1], modifier)  # acceleration due to modified gravity
    r_ip1 = verlet(position[-1], position[-2], a_i, time_step)  # verlet iteration of position
    v_ip1 = 0.5*(r_ip1 - position[-2])/time_step  # velocity iteration from centered difference
    return r_ip1, v_ip1

# Function to determine the first point after the initial value via Euler-Cromer method
# Input: mass of central body, scalar initial position,
# Output: first position, first velocity


def first_step(central_mass, initial_position, initial_velocity, modifier, time_scale):
    acceleration = gravity(central_mass, initial_position, initial_velocity, modifier)
    time_step = 2*LA.norm(initial_velocity)/LA.norm(acceleration)/time_scale  # to make sure the time step is small enough for the evolution below
    print("dt = " + str(np.around(time_step/60, 1)) + " mins")
    next_velocity = initial_velocity + acceleration * time_step
    next_position = initial_position + next_velocity * time_step
    return next_position, next_velocity, time_step

# Find the angle of a vector in Cartesian coordinates
# Input: a 2D vector
# Output: the principal angle the vector makes with the x-axis


def vector_angle(vector):
    angle = np.arctan2(vector[1], vector[0])
    # This gives the counterclockwise angle from the x-axis
    # if vector[1] >= 0:
    # 	angle = np.arctan2(vector[1], vector[0])
    # else:
    # 	angle = 2*np.pi + np.arctan2(vector[1], vector[0])
    return angle

# Loop step through one orbit/till it gets back to periapsis
# Input: mass of central body, array of previous two positions, array of previous two velocities, eccentricity, time step
# Output: array of positions for this orbit, array of velocities for this orbit


def orbit(central_mass, orbit_mass, position_array, velocity_array, time_step, modifier):
    angular_momentum = orbit_mass * np.cross(position_array[-1], velocity_array[-1])  # ASSUMING CENTRAL BODY IS FIXED
    reduced_mass = central_mass * orbit_mass / (central_mass + orbit_mass)
    # Principal Orbital Parameter (c in Taylor), NOT SPEED OF LIGHT
    pop = LA.norm(angular_momentum)**2 / (NEWTON_G * central_mass * orbit_mass * reduced_mass)
    eccentricity = (pop / LA.norm(position_array[0])) - 1
    semi_major_axis = pop/(1 - eccentricity**2)
    period = np.sqrt(4*np.pi**2 * semi_major_axis**3 * reduced_mass / (NEWTON_G * central_mass * orbit_mass))  # Orbital period
    grav_wave_freq = 2/period
    # print(period/60/60/24) #Earth days
    position, velocity = (position_array, velocity_array)  # Initialize position and velocity

    for time in np.arange(0., period, time_step):
        next_position, next_velocity = step(central_mass, position, velocity, time_step, modifier)
        position = np.append(position, [next_position], axis=0)
        velocity = np.append(velocity, [next_velocity], axis=0)

    # Now position is the array of positions for this orbit plus the last two positions from the previous orbit

    # Find the periapsis of the orbit
    distances = LA.norm(position[2:], axis=1)  # array of distances of this orbit
    index_of_periapsis = np.argmin(distances)  # index_of_periapsis of this orbit
    # print(index_of_periapsis)
    periapsis = (position[2:])[index_of_periapsis]  # periapsis of this orbit
    periapsis_angle = vector_angle(periapsis)*(180/np.pi)*3600  # convert from radians to arcseconds #angle of periapsis of this orbit

    return position[2:], velocity[2:], periapsis, periapsis_angle, eccentricity, grav_wave_freq

# Loop orbit
# Input: mass of central body, array of previous two positions, array of previous two velocities, eccentricity, time step, number of orbits
# Output: array of positions for whole time, array of velocities for whole time, array of positions of periapsides for whole time, array of angles of periapsides for whole time, array of distances for whole time


def multiple_orbits(central_mass, orbit_mass, position_array, velocity_array, time_step, number_of_orbits, modifier):
    position, velocity = (position_array, velocity_array)  # Initialize multiple orbit position and velocity arrays
    periapsides = np.zeros(len(position_array[0])*number_of_orbits).reshape(number_of_orbits,
                                                                            len(position_array[0]))  # pre-allocate periapsides array
    periapsides_angle = np.zeros(number_of_orbits)  # pre-allocate periapsides array
    grav_wave_freq_list = np.zeros(number_of_orbits)

    # Force the first periapsis to be the first element of the position array becAUse WE'RE STARTING AT PERIAPSIS
    periapsides[0] = position[0]
    periapsides_angle[0] = vector_angle(position[0])

    # CHANGE FOR LOOP TO WHILE ORBITAL DECAY TIMESCALE != ORBITAL PERIOD (B&S, eq 1.73 pg 40)

    for i in range(1, number_of_orbits):
        orbit_position, orbit_velocity, orbit_periapsis, orbit_periapsis_angle, orbit_eccentricity, orbit_grav_wave_freq = orbit(
            central_mass, orbit_mass, position[-2:], velocity[-2:], time_step, modifier)
        if i%100 == 0: 
            print("orbit #" + str(i))
            print("f_GW = " + str(np.around(orbit_grav_wave_freq * 10**6, 4)) + " micro-Hz")
        
        grav_wave_freq_list[i] = orbit_grav_wave_freq
        if i == 1:
            position = np.append(position, orbit_position, axis=0)
            velocity = np.append(velocity, orbit_velocity, axis=0)
            periapsides[i] = orbit_periapsis
            periapsides_angle[i] = orbit_periapsis_angle
            eccentricity = orbit_eccentricity
            #distances_array = np.append(LA.norm(position_array, axis = 1), orbit_distances)
        else:
            position = np.append(position, orbit_position, axis=0)
            velocity = np.append(velocity, orbit_velocity, axis=0)
            periapsides[i] = orbit_periapsis
            periapsides_angle[i] = orbit_periapsis_angle
            #distances_array = np.append(distances_array, orbit_distances)

    return position, velocity, periapsides, periapsides_angle, orbit_eccentricity, grav_wave_freq_list

# Best fit line for periapsis angles
# Input: array of times, array of periapsides angles
# Output: slope/rate of angular precession


def precession_rate(times, angles):
    rate = np.polyfit(times, angles, 1)[0]
    return rate

# MAIN ROUTINE
# Input:
# Output:


def main(central_mass, orbit_mass, initial_position, initial_velocity, number_of_orbits, modifier, time_scale):
    first_position, first_velocity, time_step = first_step(central_mass, initial_position, initial_velocity, modifier, time_scale)
    position, velocity = (np.array([initial_position, first_position]), np.array([initial_velocity, first_velocity]))
    position, velocity, periapsides, periapsides_angle, eccentricity, grav_wave_freq_list = multiple_orbits(
        central_mass, orbit_mass, position, velocity, time_step, number_of_orbits, modifier)

    return position, velocity, periapsides, periapsides_angle, eccentricity, grav_wave_freq_list

# Trajectory plot
# Input:
# Output:


def trajectory_plot(positions, mercury_orbit_number):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(positions[:, 0]/AU, positions[:, 1]/AU, linewidth=2)
    ax.plot([0], [0], 'yo', markersize=15)
    ax.set_xlabel('x [AU]')
    ax.set_ylabel('y [AU]')
    ax.grid(True)
    ax.set_title(mercury_orbit_number + " Mercury orbits")
    ax.set_aspect('equal')
    # fig.savefig(vis_dir+'Trajectory_'+str(int(number_of_orbits*MERCURY_YEAR/YEAR))+'.pdf')
    plt.show()

#
# Input:
# Output:


def periapsis_angle_plot(times, periapsides_angles):
    rate = precession_rate(times, periapsides_angles)
    plt.plot(times, periapsides_angles, times, rate*times + periapsides_angles[0])
    plt.xlabel('Time (Earth YEARs)')
    plt.ylabel('Periapsis angle (arcseconds)')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    # 'Precession of Mercury over '+str(int(number_of_orbits*MERCURY_YEAR/YEAR))+" Earth YEARs")
    plt.title('Precession rate: '+str(round(rate*100, 0))+' arcseconds/century')
    plt.tight_layout()
    plt.savefig(vis_dir+'Precession_Angle_'+str(int(number_of_orbits*MERCURY_YEAR/YEAR))+'.pdf')
    plt.show()

# Distance plot
# Input:
# Output:


def distance_plot(times, distances):
    plt.plot(times, distances/AU)  # *MERCURY_YEAR/YEAR
    plt.xlabel('Time (Earth YEARs)')
    plt.ylabel('Distance from Sun (AU)')
    plt.tight_layout()
    #plt.title(str(int(number_of_orbits))+" Earth YEARs")
    # plt.savefig(vis_dir+'Distance_'+str(int(number_of_orbits))+'.pdf')
    plt.show()

# Comparing Newtonian to GR
# Input: distance lower bound, distance upper bound, mass of central object, velocity vector of object
# Output: plot of difference between Newtonian gravity and GR


def gravity_comparison(lower_distance, upper_distance, mass, velocity):
    distances = np.linspace(lower_distance, upper_distance, 100)*AU
    newton = np.array([LA.norm(gravity(mass, np.array([r, 0.]), velocity, "Newtonian")) for r in distances])  # np.array([0, 2*np.pi*AU/YEAR])
    GR = np.array([LA.norm(gravity(mass, np.array([r, 0.]), velocity, "GR")) for r in distances])

    plt.plot(distances/AU, np.abs(newton - GR))
    plt.plot([0.31, 0.31], [np.abs(newton - GR)[-1], np.abs(newton - GR)[0]], 'k--')  # Mercury's perihelion
    plt.plot([0.47, 0.47], [np.abs(newton - GR)[-1], np.abs(newton - GR)[0]], 'k--')  # Mercury's aphelion
    plt.xlabel("Distance (AU)")
    plt.ylabel("Acceleration (m/s^2)")
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    #plt.title("Newtonian-GR Comparison")
    plt.tight_layout()
    # plt.savefig(vis_dir+'NEWTON_GR_Comp_Mercury.pdf')
    plt.show()


def grav_wave_freq_plot(orbit_list, grav_wave_freq_list):
    plt.plot(orbit_list, grav_wave_freq_list)
    plt.xlabel('Orbits')
    plt.ylabel('f_GW [micro-Hz]')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    # 'Precession of Mercury over '+str(int(number_of_orbits*MERCURY_YEAR/YEAR))+" Earth YEARs")
    # plt.title('Precession rate: '+str(round(rate*100, 0))+' arcseconds/century')
    plt.tight_layout()
    # plt.savefig(vis_dir+'Precession_Angle_'+str(int(number_of_orbits*MERCURY_YEAR/YEAR))+'.pdf')
    plt.show()
