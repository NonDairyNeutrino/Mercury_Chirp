"""main."""
import numpy as np
# from numpy import linalg as LA
import matplotlib.pyplot as plt
# import matplotlib
# from mpl_toolkits.mplot3d import Axes3D
import multiprocessing as mp

from constants import *
from functions import *

modifier         =  "Newtonian"  # Can be GR or Newtonian
scale_factor     = 200
central_mass     = SUN_MASS
orbiting_mass    = MERCURY_MASS
initial_position = np.array([MERCURY_PEREHELION_DISTANCE, 0.])  # start at perihelion
initial_velocity = np.array([0., MERCURY_ORBITAL_VELOCITY])  # speed at perihelion
earth_years      = 1000
number_of_orbits = int(np.ceil(earth_years*(YEAR/MERCURY_YEAR))) # i.e. mercury_years
orbits           = np.arange(0, earth_years)

position, velocity, periapsides, periapsides_angle, eccentricity, grav_wave_freq_list = main(
    central_mass, orbiting_mass, initial_position, initial_velocity, number_of_orbits, modifier, scale_factor)

grav_wave_freq_plot(np.arange(1, number_of_orbits), grav_wave_freq_list[1:].round(6))

# SAVE JUST COMPUTED DATA
# np.savetxt(data_dir+test_modifier+str(int(number_of_orbits*mercury_year/year))+'positions.txt', test_position)
# np.savetxt(data_dir+test_modifier+str(int(number_of_orbits*mercury_year/year))+'periapsides.txt', test_periapsides)
# np.savetxt(data_dir+test_modifier+str(int(number_of_orbits*mercury_year/year))+'periapsides_angles.txt', test_periapsides_angle)

# np.savetxt(data_dir+test_modifier+'distances_'+str(int(number_of_orbits*mercury_year/year))+'.txt', test_distances_array)

# USE ALREADY COMPUTED DATA
# test_position = np.loadtxt(data_dir+test_modifier+str(int(number_of_orbits*mercury_year/year))+'positions.txt')
# test_periapsides = np.loadtxt(data_dir+test_modifier+str(int(number_of_orbits*mercury_year/year))+'periapsides.txt')
# test_periapsides_angle = np.loadtxt(data_dir+test_modifier+str(int(number_of_orbits*mercury_year/year))+'periapsides_angles.txt')

# test_distances_array = np.loadtxt(data_dir+'distances_'+str(int(number_of_orbits*mercury_year/year))+'.txt')

# orbits = np.arange(0, number_of_orbits)*mercury_year/year

trajectory_plot(position, str(number_of_orbits))
# periapsis_angle_plot(orbits, test_periapsides_angle)
# distance_plot(np.linspace(0., number_of_orbits, len(test_distances_array)), test_distances_array)

# TRAJECTORIES FOR DIFFERENT TIME STEPS
# initial_position = np.array([au, 0.]) #Earth distance from sun
# initial_velocity = np.array([0., 2*np.pi*au/year]) #Earth average orbital speed
# number_of_orbits = 10#int(np.ceil(1*(year/mercury_year)))
# test_modifier = "GR" #Can be GR or Newtonian
# scales = np.linspace(2, 20, 5)
# i=0
# fig = plt.figure()
# ax = fig.add_subplot(111)
#
# for time_scale in scales:
# 	test_position, test_velocity, test_periapsides, test_periapsides_angle, eccentricity = main(sun_mass, earth_mass, initial_position, initial_velocity, number_of_orbits, test_modifier, time_scale)
# 	np.savetxt(data_dir+'Earth'+test_modifier+str(number_of_orbits)+'_'+str(time_scale)+'positions.txt', test_position)
# 	np.savetxt(data_dir+'Earth'+test_modifier+str(number_of_orbits)+'_'+str(time_scale)+'periapsides.txt', test_periapsides)
# 	np.savetxt(data_dir+'Earth'+test_modifier+str(number_of_orbits)+'_'+str(time_scale)+'periapsides_angles.txt', test_periapsides_angle)
#
# 	ax.plot(test_position[:,0]/au, test_position[:,1]/au, linewidth = 1, label = str(time_scale))
#
# ax.plot([0], [0], 'yo', markersize = 15)
# ax.set_xlabel('x [au]')
# ax.set_ylabel('y [au]')
# ax.legend(loc = 'upper left')
# ax.grid(True)
# ax.set_title("Earth's orbit over "+str(number_of_orbits)+" Earth years")
# ax.set_aspect('equal')
# fig.savefig(vis_dir+'Earth_Trajectory_'+str(number_of_orbits)+'.pdf')
# plt.show()

# PERIAPSIDES ANGLES FOR DIFFERENT TIME STEPS
# initial_position = np.array([au, 0.]) #Earth distance from sun
# initial_velocity = np.array([0., 2*np.pi*au/year]) #Earth average orbital speed
# number_of_orbits = 10#int(np.ceil(1*(year/mercury_year)))
# test_modifier = "GR" #Can be GR or Newtonian
# scales = np.linspace(2, 20, 5)
# i=0
# fig = plt.figure()
# ax = fig.add_subplot(111)
#
# for time_scale in scales:
# 	test_position, test_velocity, test_periapsides, test_periapsides_angle, eccentricity = main(sun_mass, earth_mass, initial_position, initial_velocity, number_of_orbits, test_modifier, time_scale)
# 	np.savetxt(data_dir+'Earth'+test_modifier+str(number_of_orbits)+'_'+str(time_scale)+'positions.txt', test_position)
# 	np.savetxt(data_dir+'Earth'+test_modifier+str(number_of_orbits)+'_'+str(time_scale)+'periapsides.txt', test_periapsides)
# 	np.savetxt(data_dir+'Earth'+test_modifier+str(number_of_orbits)+'_'+str(time_scale)+'periapsides_angles.txt', test_periapsides_angle)
#
# 	ax.plot(test_position[:,0]/au, test_position[:,1]/au, linewidth = 1, label = str(time_scale))
#
# ax.plot([0], [0], 'yo', markersize = 15)
# ax.set_xlabel('x [au]')
# ax.set_ylabel('y [au]')
# ax.legend(loc = 'upper left')
# ax.grid(True)
# ax.set_title("Earth's orbit over "+str(number_of_orbits)+" Earth years")
# ax.set_aspect('equal')
# fig.savefig(vis_dir+'Earth_Trajectory_'+str(number_of_orbits)+'.pdf')
# plt.show()

# RATE COMPARISON FOR DIFFERENT RUN TIMES
# initial_position = np.array([46.*10.**9., 0.]) #start at perihelion
# initial_velocity = np.array([0., (59*10**3)]) #speed at perihelion
# test_modifier = "GR" #Can be GR or Newtonian
# times = np.linspace(1, 100, 20)
# rates = np.zeros(len(times))
# i=0
#
# for Earth_years in times:
# 	number_of_orbits = int(np.ceil(Earth_years*(year/mercury_year)))
# 	test_position, test_velocity, test_periapsides, test_periapsides_angle, eccentricity = main(sun_mass, mercury_mass, initial_position, initial_velocity, number_of_orbits, test_modifier)
# 	orbits = np.arange(0, number_of_orbits)*mercury_year/year
# 	rates[i] = precession_rate(orbits, test_periapsides_angle)
# 	i += 1
#
# plt.plot(times, rates)
# plt.xlabel('End Time (Earth years)')
# plt.ylabel('Angle Precession (arcseconds/century)')
# plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
# #plt.title('Precession rate: '+str(round(rate*100,0))+' arcseconds/century')#'Precession of Mercury over '+str(int(number_of_orbits*mercury_year/year))+" Earth years")
# plt.tight_layout()
# plt.savefig(vis_dir+'Precession_rate_Runtime_Comparison.pdf')
# plt.show()

# PRECESSION RATE VS ECCENTRICITY
# number_of_factors = 10
# initial_position = np.array([46.*10.**9., 0.])  # start at perihelion
# number_of_orbits = int(np.ceil(20*(YEAR/MERCURY_YEAR)))
# test_modifier = "GR"  # Can be GR or Newtonian
# orbits = np.arange(0, number_of_orbits)*MERCURY_YEAR/YEAR


# def foo(factor):
# 	initial_velocity = np.array([0., (59*10**3)*factor])  # speed at perihelion
# 	test_position, test_velocity, test_periapsides, test_periapsides_angle, eccentricity = main(
# 		SUN_MASS, MERCURY_MASS, initial_position, initial_velocity, number_of_orbits, test_modifier, 2000)
# 	rate = precession_rate(orbits, test_periapsides_angle)*100

# 	return eccentricity, rate


# if __name__ == '__main__':
# 	pool = mp.Pool(mp.cpu_count()-1)
# 	results = np.transpose(pool.map(foo, np.linspace(0.915, 1.2, number_of_factors)))  # 1.285 gives an eccentricity of about 0.99
# 	pool.close()
# 	plt.plot(results[0], results[1])
# 	plt.xlabel('Eccentricity')
# 	plt.ylabel('Precession Rate (arcseconds/century)')
# 	plt.tight_layout()
# 	plt.savefig(vis_dir+'rate_v_eccentricity_final.pdf')
# 	plt.show()

# i = 0
# number_of_factors = 2
# initial_position = np.array([46.*10.**9., 0.]) #start at perihelion
# number_of_orbits = int(np.ceil(100*(year/mercury_year)))
# test_modifier = "GR" #Can be GR or Newtonian
# orbits = np.arange(0, number_of_orbits)*mercury_year/year
#
# eccentricities = np.zeros(number_of_factors)
# rates = np.zeros(number_of_factors)
#
# for factor in [1.2]: #np.linspace(0.915, 1.285, number_of_factors): #Starting and stopping points are hard coded because I just manually found the factors that gave eccentricities close to 0 and 1
# 	initial_velocity = np.array([0., (59*10**3)*factor]) #speed at perihelion
# 	test_position, test_velocity, test_periapsides, test_periapsides_angle, eccentricity = main(sun_mass, mercury_mass, initial_position, initial_velocity, number_of_orbits, test_modifier)
# 	eccentricities[i] = eccentricity
# 	rates[i] = precession_rate(orbits, test_periapsides_angle)*100
# 	i += 1
#
# plt.plot(eccentricities, rates)
# plt.xlabel('Eccentricity')
# plt.ylabel('Precession Rate (arcseconds/century)')
# plt.show()

# COMPARISON BETWEEN NEWTONIAN GRAVITY AND GR FOR INTERESTING REGIMES
# gravity_comparison(0.3, 0.48, sun_mass, np.array([0, mercury_orbital_velocity]))

#3D PLOT OF ORBIT
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot(test_position[0]/au, test_position[1]/au, test_position[2]/au)
# plt.show()
