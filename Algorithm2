import numpy as np
import random
import scipy
from matplotlib import pyplot as plt

def neighbor_list(i, j, k, nx):
  """ Find all neighbors of site (i, j).

  Args:
    i (int): site index along x
    j (int): site index along y
    nx (int): number of sites along each dimension
  Return:
    list: a list of 2-tuples, [(i_left, j_left), (i_top, j_top),
     (i_right, j_right), (i_bottom, j_bottom)]
  """
  left_center   = (i-1, j, k)
  right_center  = (i+1, j, k)
  top_center    = (i, j+1, k)
  bottom_center = (i, j-1, k)
  left_up = (i, j, k + 1)
  left_down = (i, j, k -1)
  return np.mod([left_center, right_center, top_center, bottom_center, left_up, left_down], nx)


def initialize_lattice(beta, lim, n_atoms):
    lattice = np.ones([lim-1, lim-1, lim-1])
    for i in range(lim-1):
        for j in range(lim-1):
            for k in range(lim-1):
                occupation = 1 / (np.exp(beta * ((i+1)**2 + (k+1)**2 + (j+1)**2)) - 1)
                lattice[i][j][k] = occupation
    lattice = np.round(lattice/np.sum(lattice) * n_atoms)
    return lattice, np.sum(lattice)

def energy(lattice, lim):
    en = 0
    for i in range(lim-1):
        for j in range(lim-1):
            for k in range(lim-1):
                en = en + lattice[i,j,k] * ((i+1)**2 + (j+1)**2 + (k*i)**2)
    return en

def move_probability(lattice, lattice_temp, site, neighbor, beta, lim):
    energy1 = energy(lattice, lim)
    energy2 = energy(lattice_temp, lim)
    oc1 = lattice[site[0] , site[1], site[2]]
    oc2 = lattice[neighbor[0] , neighbor[1], neighbor[2]]
    return oc2 / oc1 * np.exp(-beta *(energy2 - energy1))

if __name__ == '__main__':

    # Initial conditions in reduced units
    betas = [5., 4., 3., 2., 1., 0.5, 0.25, 0.1]
    mass = 48.
    lim = 10                                    # allowed k states
    mu = 0.1                                     # chemical potential
    nsweep = 1000
    lat_dat = 'lat.dat'
    scalar_dat = 'scalar.dat'
    nconf = 10
    num_atoms = 1000
    ntherm = 1
    therm_fmt = '%8d %14.8f %14.8f\n'
    therm_header = '# isweep  e\n'
    for beta in betas:
        # Initialize the lattice at specified temperature
        lattice, num_atoms = initialize_lattice(beta, lim, num_atoms)
        restart = False
        if restart:
            lattice = np.loadtxt('last_lat.dat')
            num_atoms = np.sum(lattice)
        else:
            with open(lat_dat, 'w') as flat:
                flat.write('# n_atoms=%d\n' % num_atoms)
            with open(scalar_dat, 'w') as fsca:
                fsca.write(therm_header)
        flat = open(lat_dat, 'a')
        fsca = open(scalar_dat, 'a')

        msg = 'starting %s sweeps using temperature %s and number of particles %s'% (nsweep, 1/beta, num_atoms)
        if restart:
            msg = 're'+msg
        etot = energy(lattice, lim)
        from time import clock
        start = clock()
        print(msg)
        naccept = 0
        nattempt = 0
        for isweep in range(nsweep):
            if (isweep % ntherm == 0):
                etot = energy(lattice, lim)
                occu = lattice[0][0][0]
                fsca.write(therm_fmt % (isweep, etot, occu))
            site = np.random.randint(0, lim-1, (1, 3))[0]
            while lattice[site[0] , site[1], site[2]] < 1:
                site = np.random.randint(0, lim -1, (1, 3))[0]
            neighbors = neighbor_list(site[0], site[1], site[2], lim-1)
            neighbor = random.choice(neighbors)
            lattice_temp = np.copy(lattice)
            lattice_temp[site[0],site[1], site[2]] = lattice_temp[site[0],site[1], site[2]] - 1
            lattice_temp[neighbor[0], neighbor[1], neighbor[2]] = lattice_temp[neighbor[0], neighbor[1], neighbor[2]] + 1
            prob = move_probability(lattice, lattice_temp, site, neighbor, beta, lim)
            nattempt += 1
            rand = np.random.random()
            if rand < prob:
                naccept += 1
                lattice[site[0], site[1], site[2]] = lattice[site[0], site[1], site[2]] - 1
                lattice[neighbor[0], neighbor[1], neighbor[2]] = lattice[neighbor[0], neighbor[1], neighbor[2]] + 1
        print(lattice[0][0][0]/num_atoms, energy(lattice, lim))

