"""PLOTTING FUNCTIONS."""
import os
import numpy as np
import matplotlib.pyplot as plt
from constants import AU, MERCURY_PEREHELION_DISTANCE, MERCURY_YEAR, YEAR
from functions import precession_rate

cwd_parent = os.path.dirname(os.getcwd())
vis_dir = cwd_parent + '\\Plots\\'


def trajectory_plot(positions):
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
    plt.plot(positions[:, 0] / AU, positions[:, 1] / AU, linewidth=2)
    plt.plot([0], [0], "yo", markersize=15)
    plt.xlim([-0.1, 0.35])
    plt.ylim([-0.1, 0.1])
    plt.xlabel("x [AU]")
    plt.ylabel("y [AU]")
    plt.grid(True)
    plt.title("Periapsides over " + str(len(positions)) + " Mercury orbits")
    plt.tight_layout()
    plt.savefig(vis_dir + 'Periapsides_' + str(int(len(positions) * MERCURY_YEAR / YEAR)) + '.pdf')
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
    plt.xlabel("Time [Orbits]")
    plt.ylabel("Periapsis angle (arcseconds)")
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    plt.title("Precession rate: " + str(round((rate * MERCURY_YEAR / YEAR) * 100, 0)) + " arcseconds/century")
    plt.tight_layout()
    plt.savefig(vis_dir + 'Periapsis_angle_' + str(int(len(times) * MERCURY_YEAR / YEAR)) + '.pdf')
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
    plt.plot(np.around(distances / MERCURY_PEREHELION_DISTANCE, 9))  # *MERCURY_YEAR/YEAR
    plt.xlabel("Time (Mercury YEARs)")
    plt.ylabel("% Periapsis distance")
    plt.tight_layout()
    plt.savefig(vis_dir + 'Periapsis_distance_' + str(int(len(distances) * MERCURY_YEAR / YEAR)) + '.pdf')
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
    plt.savefig(vis_dir + 'GW_freq_' + str(int(len(grav_wave_freq_list) * MERCURY_YEAR / YEAR)) + '.pdf')
    plt.show()
