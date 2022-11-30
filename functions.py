"""FUNCTIONS."""

import numpy as np
from numpy import linalg as LA

from particle_class import Particle


def time_step(particle, time_scale):
    """
    Generate time step based on initial velocity and acceleration.

    Parameters
    ----------
    central_mass : float
        Mass around which the body is orbiting.
    initial_position : [float, float]
        Starting position of the body relative to the central mass.
    initial_velocity : [float, float]
        Starting velocity of the body relative to the central mass
    modifier : string
        Either "GR" or "Newtonian".
    time_scale : int or float
        By which to inversely scale the operation.

    Returns
    -------
    dt : float
        time step.

    """
    speed = LA.norm(particle.velocity)
    acceleration = LA.norm(particle.acceleration)
    dt = 2 * speed / acceleration / time_scale  # to make sure the time step is small enough for the evolution below
    print("dt = " + str(np.around(dt / 60)) + " mins")
    return dt


def step(particle, dt):
    """
    Perform leapfrog integration using the 4-th order Yoshida method.

    Parameters
    ----------
    particle : Particle
        Object of the Particle class.
    dt : float
        Time step.

    Returns
    -------
    None.

    """
    # integration method parameters
    w_1 = 1 / (2 - 2**(1 / 3))
    w_0 = - 2**(1 / 3) * w_1
    # c_1 = c_4 = w_1 / 2
    # c_2 = c_3 = 0.5 * (w_0 + w_1)
    c = np.dot(0.5 * np.array([[0, 1], [1, 1], [1, 1], [0, 1]]), [w_0, w_1])
    # d_1 = d_3 = w_1
    # d_2 = w_0
    # d_4 = 0
    d = np.dot([[0, 1], [1, 0], [0, 1], [0, 0]], [w_0, w_1])

    # In-place updating of the given particle object
    for k in range(4):
        x_k = particle.get_position()
        v_k = particle.get_velocity()

        x_k = x_k + c[k] * v_k * dt
        particle.set_position(x_k)

        particle.set_acceleration()
        a_k = particle.get_acceleration()

        v_k = v_k + d[k] * a_k * dt
        particle.set_velocity(v_k)
        # Including v_k in the acceleration is a naive modification of the real method
    # Final state: x_k = x_{i + 1}, v_k = {i + 1}


def vector_angle(vector):
    """
    Find the angle of a vector with respect to the vector [1,0].

    Parameters
    ----------
    vector : [float, float]
        Vector.

    Returns
    -------
    angle : float
        Angle.

    """
    # This gives the counterclockwise angle from the x-axis
    if vector[1] >= 0:
        angle = np.arctan2(vector[1], vector[0])
    else:
        angle = 2 * np.pi + np.arctan2(vector[1], vector[0])
    return angle


def orbit(particle, dt):
    """
    Step until the particle completes a single orbit.

    Parameters
    ----------
    particle : Particle
        Object of the Particle class.
    dt : float
        Time step.

    Returns
    -------
    periapsis : [float, float]
        Position of the periapsis of the orbit.
    periapsis_angle : float
        Angle the periapsis makes with [1, 0].
    grav_wave_freq : float
        Leading order frequency of the gravitational wave produced by the orbit of the particle.

    """
    grav_wave_freq = 2 / particle.get_period()
    periapsis = particle.get_position()
    for time in np.arange(0, particle.get_period(), dt):
        step(particle, dt)

        if LA.norm(particle.get_position()) < LA.norm(periapsis):
            periapsis = particle.get_position()

    periapsis_angle = np.rad2deg(vector_angle(periapsis)) * 3600  # convert from radians to arcseconds #angle of periapsis of this orbit

    return periapsis, periapsis_angle, grav_wave_freq


def multiple_orbits(particle, dt, number_of_orbits):
    """
    Make the particle orbit several times.

    Parameters
    ----------
    particle : Particle
        Object of the Particle class.
    dt : float
        Time step.
    number_of_orbits : int
        Number of orbits desired.

    Returns
    -------
    periapsides : [[float, float], ...]
        Positions of the periapsides.
    periapsis_angle_list : [float, ...]
        List of periapsis angles.
    grav_wave_freq_list : [float, ...]
        List of the gravitational wave frequencies.

    """
    periapsides = np.zeros((number_of_orbits, 2))
    periapsis_angle_list = np.zeros(number_of_orbits)
    grav_wave_freq_list = np.zeros(number_of_orbits)

    # CHANGE FOR LOOP TO WHILE ORBITAL DECAY TIMESCALE != ORBITAL PERIOD (B&S, eq 1.73 pg 40)
    for i in range(number_of_orbits):
        periapsis, periapsis_angle, grav_wave_freq = orbit(particle, dt)
        periapsides[i] = periapsis
        periapsis_angle_list[i] = periapsis_angle
        grav_wave_freq_list[i] = grav_wave_freq

        if i % 100 == 0:
            print("orbit #" + str(i))
            print("f_GW = " + str(np.around(grav_wave_freq * 10**6, 3)) + " micro-Hz")

    return periapsides, periapsis_angle_list, grav_wave_freq_list


def precession_rate(times, angles):
    """
    Get the rate of precession.

    Parameters
    ----------
    times : [float, ...]
        List of times.
    angles : [float, ...]
        List of angles of periapsides.

    Returns
    -------
    rate : float
        Rate of periapsis precession.

    """
    rate = np.polyfit(times, angles, 1)[0]
    return rate


def main(orbit_mass, mass, initial_position, initial_velocity, number_of_orbits, time_scale):
    """
    Orbit a particle with given properties.

    Parameters
    ----------
    orbit_mass : float
        Mass of the body being orbited.
    mass : float
        Mass of the particle.
    initial_position : [float, float]
        Initial position of the particle.
    initial_velocity : [float, float]
        Initial velocity of the particle.
    number_of_orbits : int
        Number of orbits to do.
    time_scale : int or float
        Parameter with which to inversely scale the time step.

    Returns
    -------
    periapsides : [[float, float], ...]
        Positions of the periapsides.
    periapsis_angle_list : [float, ...]
        List of periapsis angles.
    grav_wave_freq_list : [float, ...]
        List of the gravitational wave frequencies.

    """
    particle = Particle(orbit_mass, mass, initial_position, initial_velocity)
    dt = time_step(particle, time_scale)
    periapsides, periapsis_angle_list, grav_wave_freq_list = multiple_orbits(particle, dt, number_of_orbits)

    return periapsides, periapsis_angle_list, grav_wave_freq_list
