"""Particle class."""
import numpy as np
from numpy import linalg as LA
import constants as c


class Particle:
    """Particle class."""

    def __init__(self, orbit_mass, mass, position, velocity):
        """
        Construct particle object.

        Parameters
        ----------
        central_mass : float
            Mass of the body around which the particle is orbiting.
        mass : float
            Particle mass.
        position : [float, float]
            Particle position.
        velocity : [float, float]
            Particle velocity.

        Returns
        -------
        None.

        """
        self.orbit_mass = orbit_mass
        self.mass = mass
        self.position = position
        self.velocity = velocity
        self.set_schwarz_radius(orbit_mass)
        self.set_acceleration()
        self.set_angular_momentum()
        self.set_reduced_mass()
        self.set_pop()
        self.eccentricity = (self.pop / LA.norm(position)) - 1
        self.set_semi_major_axis()
        self.set_period()
        print(self)

    def __str__(self):
        """
        Print string representation of object.

        Returns
        -------
        str
            String representation of object.

        """
        return (
            f"{self.orbit_mass} = orbit_mass\n"
            f"{self.mass} = mass\n"
            f"{self.position} = position\n"
            f"{self.velocity} = velocity\n"
            f"{self.schwarz_radius} = schwarz rad\n"
            f"{self.get_acceleration()} = acceleration\n"
            f"{self.get_angular_momentum()} = ang_mom\n"
            f"{self.get_reduced_mass()} = reduced mass\n"
            f"{self.get_pop()} = pop\n"
            f"{self.eccentricity} = eccentricity\n"
            f"{self.get_semi_major_axis()} = sma\n"
            f"{self.get_period()} = period\n"
        )

    def set_schwarz_radius(self, central_mass):
        """
        Schwarzchild radius for a the central mass.

        Parameters
        ----------
        mass : float
            Mass of the object around which the particle is being attracted.

        """
        self.schwarz_radius = 2 * c.NEWTON_G * central_mass / c.LIGHT_SPEED_2

# ==========VELOCITY=============
    def get_position(self):
        """
        Get particle position.

        Returns
        -------
        [float, float]
            Particle position.

        """
        return self.position

    def set_position(self, position):
        """
        Set the position of the particle.

        Parameters
        ----------
        position : [float, float]
            Position of the particle.

        Returns
        -------
        None.

        """
        self.position = position

# ==========VELOCITY=============
    def get_velocity(self):
        """
        Get particle velocity.

        Returns
        -------
        [float, float]
            Particle velocity.

        """
        return self.velocity

    def set_velocity(self, velocity):
        """
        Set the velocity of the particle.

        Parameters
        ----------
        velocity : [float, float]
            Velocity of the particle.

        Returns
        -------
        None.

        """
        self.velocity = velocity

# ==========R_SQUARE_ANG_MOM=============
    def set_r_square_ang_mom_(self):
        """
        Square-norm of angular momentum / (m c)^2.

        Parameters
        ----------
        position : [float, float]
            2D position.
        velocity : [float, float]
            2D velocity.

        """
        self.r_square_ang_mom = np.sum(np.cross(self.position, self.velocity) ** 2) / c.LIGHT_SPEED_2

# ==========ACCELERATION=============
    def get_acceleration(self):
        """
        Get particle acceleration.

        Returns
        -------
        [float, float]
            Particle acceleration.

        """
        return self.acceleration

    def set_acceleration(self):
        """
        Set the acceleration of the particle due to gravity.

        Returns
        -------
        None.

        """
        distance = LA.norm(self.position)  # This implicitly defines the sun as the origin
        self.set_r_square_ang_mom_()
        relativistic_correction = (1 + 3 * self.r_square_ang_mom / distance**2)
        self.acceleration = -0.5 * c.LIGHT_SPEED_2 * self.schwarz_radius * relativistic_correction * self.position / distance**3

# ==========ANGULAR MOMENTUM=============
    def get_angular_momentum(self):
        """
        Get the angular momentum of the particle.

        Returns
        -------
        [float, float]
            Angular momentum of the particle.

        """
        return self.angular_momentum

    def set_angular_momentum(self):
        """
        Set the angular momentum of the particle.

        Returns
        -------
        None.

        """
        # ASSUMING CENTRAL BODY IS FIXED
        self.angular_momentum = self.orbit_mass * np.cross(self.position, self.velocity)

# ==========REDUCED MASS=============
    def get_reduced_mass(self):
        """
        Get the reduced mass of the particle and its orbited body.

        Returns
        -------
        float
            Reduced mass.

        """
        return self.reduced_mass

    def set_reduced_mass(self):
        """
        Set the reduced mass between the particle and its orbited body.

        Returns
        -------
        None.

        """
        self.reduced_mass = self.orbit_mass * self.mass / (self.orbit_mass + self.mass)

# ==========PRINCIPAL ORBITAL PARAMETER=============
    def get_pop(self):
        """
        Get the principal orbital parameter of the particle and its orbited body.

        Returns
        -------
        float
            Principal orbital parameter.

        """
        return self.pop

    def set_pop(self):
        """
        Set the reduced mass between the particle and its orbited body.

        Returns
        -------
        None.

        """
        self.pop = self.angular_momentum ** 2 / (c.NEWTON_G * self.orbit_mass * self.mass * self.reduced_mass)

# ==========SEMI MAJOR AXIS=============
    def get_semi_major_axis(self):
        """
        Get the semi_major_axis of the particle and its orbited body.

        Returns
        -------
        float
            semi_major_axis.

        """
        return self.semi_major_axis

    def set_semi_major_axis(self):
        """
        Set the reduced mass between the particle and its orbited body.

        Returns
        -------
        None.

        """
        self.semi_major_axis = self.pop / (1 - self.eccentricity**2)

# ==========ORBITAL PERIOD=============
    def get_period(self):
        """
        Get the orbital period of the particle.

        Returns
        -------
        float
            Period.

        """
        return self.period

    def set_period(self):
        """
        Set the orbital period of the particle.

        Returns
        -------
        None.

        """
        self.period = np.sqrt(
            4 * np.pi**2 * self.semi_major_axis**3 * self.reduced_mass / (
                c.NEWTON_G * self.orbit_mass * self.mass
            )
        )


def particle_test():
    """
    Test particle class.

    Returns
    -------
    None.

    """
    p = Particle(c.SUN_MASS, c.MERCURY_MASS, np.array([c.MERCURY_PEREHELION_DISTANCE, 0]), np.array([0, c.MERCURY_PEREHELION_SPEED]))
    print("MASS:", "PASSED" if p.mass == c.MERCURY_MASS else "FAILED")
    print("POSITION:", "PASSED" if all(p.position == np.array([c.MERCURY_PEREHELION_DISTANCE, 0])) else "FAILED")
    print("VELOCITY:", "PASSED" if all(p.velocity ==  np.array([0, c.MERCURY_PEREHELION_SPEED])) else "FAILED")
    # print("ACCELERATION:", "PASSED" if all(p.acceleration)
