"""PLOTTING FUNCTIONS."""
import numpy as np
import matplotlib.pyplot as plt
from constants import AU, MERCURY_PEREHELION_DISTANCE
from functions import precession_rate


def trajectory_plot(positions, mercury_orbit_number):
    """
    Produce a plot of the trajectory of Mercury.

    Parameters
    ----------
    positions : TYPE
        DESCRIPTION.
    mercury_orbit_number : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(positions[:, 0] / AU, positions[:, 1] / AU, linewidth=2)
    ax.plot([0], [0], "yo", markersize=15)
    ax.set_xlabel("x [AU]")
    ax.set_ylabel("y [AU]")
    ax.grid(True)
    ax.set_title(mercury_orbit_number + " Mercury orbits")
    ax.set_aspect("equal")
    # fig.savefig(vis_dir+'Trajectory_'+str(int(number_of_orbits*MERCURY_YEAR/YEAR))+'.pdf')
    plt.show()


def periapsis_angle_plot(periapsides_angles):
    """
    Produce a plot of the angle of the periapsis.

    Parameters
    ----------
    times : TYPE
        DESCRIPTION.
    periapsides_angles : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    times = np.arange(0, len(periapsides_angles))
    rate = precession_rate(times, periapsides_angles)
    plt.plot(times, periapsides_angles, times, rate * times + periapsides_angles[0])
    plt.xlabel("Time (Earth YEARs)")
    plt.ylabel("Periapsis angle (arcseconds)")
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    plt.title("Precession rate: " + str(round(rate * 100, 0)) + " arcseconds/century")
    plt.tight_layout()
    plt.show()


def distance_plot(distances):
    """
    Produce a plot of the distance of Mercury over time.

    Parameters
    ----------
    times : TYPE
        DESCRIPTION.
    distances : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    plt.plot(np.around(distances / MERCURY_PEREHELION_DISTANCE, 10))  # *MERCURY_YEAR/YEAR
    plt.xlabel("Time (Mercury YEARs)")
    plt.ylabel("Distance from Sun (Mercury Perihelion distance)")
    plt.tight_layout()
    # plt.title(str(int(number_of_orbits))+" Earth YEARs")
    # plt.savefig(vis_dir+'Distance_'+str(int(number_of_orbits))+'.pdf')
    plt.show()


def grav_wave_freq_plot(grav_wave_freq_list):
    """
    Produce a plot of the frequencies of the gravitational waves emitted by Mercury.

    Parameters
    ----------
    grav_wave_freq_list : list of floats
        list of gravitational wave frequencies.

    Returns
    -------
    Plot of emitted GW frequencies.

    """
    plt.plot(grav_wave_freq_list.round(10) * 10**6)
    plt.xlabel("Orbits")
    plt.ylabel("f_GW [micro-Hz]")
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    # 'Precession of Mercury over '+str(int(number_of_orbits*MERCURY_YEAR/YEAR))+" Earth YEARs")
    # plt.title('Precession rate: '+str(round(rate*100, 0))+' arcseconds/century')
    plt.tight_layout()
    # plt.savefig(vis_dir+'Precession_Angle_'+str(int(number_of_orbits*MERCURY_YEAR/YEAR))+'.pdf')
    plt.show()
