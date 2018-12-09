import numpy as np
import scipy
from matplotlib import pyplot as plt

def initialize_lattice(num_atoms, beta, lim):
    ''' parameters:
            num_atoms: number of atoms we're simulating with. Should be cubic number.
            beta: temperature parameter, hbar * omega / T
            lim: highest energy state occupied
        returns:
               lattice: (num_atoms**1/3, num_atoms**1/3, num_atoms**1/3, 3) array with energy levels distributed according to BE distribution
      NOTE: Assumes energy levels of each particle same as that of quantum harmonic oscillator hbar * omega (n + 1/2)         
      '''
    
    nx = int(round(num_atoms ** 1/3))
    lattice = np.ones([nx, nx, nx, 3], dtype=float)
    lat = np.random.rand(nx, nx, nx)
    prob = []
    for i in range(lim):
        compare = 1 / (np.exp(beta * (i+0.5))-1)
        prob.append(compare)
    prob = np.cumsum(prob/np.sum(prob))
    sel = lat < prob[0]
    lattice[sel] = 0
    for i in range(lim - 1):
        sel = (lat < prob[i+1]) & (lat > prob[i])
        lattice[sel] = i
    return lattice


def energy(lattice):
   ''' parameters:
   '''
    


###
# parameters:
#    lattice: the lattice we're working with
#    site: the site on the lattice being changed
#    e_coord: the energy coordinate being changed e.g. (0,0,1) for z, or (1,0,0) for x
#
# returns:
#    Difference in energy from changing specifiec site
###
def energy_diff(lattice, site, e_coord):
    old_en = np.sum(np.square(lattice[site]))
    new_en = np.sum(np.square(lattice[site]+e_coord))
    return old_en - new_en


if __name__ == '__main__':

    # Initial conditions in reduced units
    beta = 1.
    mass = 48.
    lim = 10
    num_atoms = 100
    nsweep = 1000
    lat_dat = 'lat.dat'
    scalar_dat = 'scalar.dat'
    nconf = 10
    ntherm = 1
    therm_fmt = '%8d %14.8f %14.8f %14.8f\n'
    therm_header = '# isweep  m2  mavg  pe\n'

    # Initialize the lattice at specified temperature
    lattice = initialize_lattice(num_atoms, beta, lim)
    
    restart = False
    if restart:
        lat = np.loadtxt('last_lat.dat')
        nx = lat.shape[0]
    else:
        with open(lat_dat, 'w') as flat:
            flat.write('# nx=%d\n' % nx)
        with open(scalar_dat, 'w') as fsca:
            fsca.write(therm_header)
    flat = open(lat_dat, 'a')
    fsca = open(scalar_dat, 'a')
    
    msg = 'starting %s sweeps with %d 
    if restart:
           msg = 're'+msg

