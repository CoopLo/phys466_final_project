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
    list_to_return = []
    if i > 0:
        left_center   = (i-1, j, k)
        list_to_return.append(left_center)
    if i<nx-1:
        right_center  = (i+1, j, k)
        list_to_return.append(right_center)
    if j>0:
        top_center    = (i, j+1, k)
        list_to_return.append(top_center)
    if j<nx-1:
        bottom_center = (i, j-1, k)
        list_to_return.append(bottom_center)
    if k>0:
        left_up = (i, j, k + 1)
        list_to_return.append(left_up)
    if j<nx-1:
        left_down = (i, j, k -1)
        list_to_return.append(left_down)

    return np.mod(list_to_return, nx)

def initialize_lattice(beta, lim, n_atoms,energy_table,dist):
    lattice = np.ones([lim-1, lim-1, lim-1])
    if dist == 'boltz':
        lattice[:]=np.exp(-energy_table[:]*beta)
    elif dist == 'rand':
        lattice[:]=np.random.random()
    else:
        lattice[:]=1 / (np.exp(beta * (energy_table[:])) - 1)
    # for i in range(lim-1):
    #     for j in range(lim-1):
    #         for k in range(lim-1):
    #             occupation = 1 / (np.exp(beta * (energy_table[i,j,k])) - 1)
    #             lattice[i][j][k] = occupation
    lattice = np.round(lattice/np.sum(lattice) * n_atoms)
    #lattice = np.zeros([lim-1,lim-1,lim-1])
    #lattice[0,0,0]=n_atoms
    return lattice, np.sum(lattice)

def energy(lattice, lim, energy_table):
    en = 0
    for i in range(lim-1):
        for j in range(lim-1):
            for k in range(lim-1):
                en = en + lattice[i,j,k] * energy_table[i,j,k]
    return en

def move_probability(lattice, lattice_temp, site, neighbor, beta, lim,energy_table):
    energy1 = energy(lattice, lim,energy_table)
    energy2 = energy(lattice_temp, lim,energy_table)
    oc1 = lattice[site[0] , site[1], site[2]]
    oc2 = lattice[neighbor[0] , neighbor[1], neighbor[2]]
    return oc2 / oc1 * np.exp(-beta *(energy2 - energy1))

if __name__ == '__main__':
    
    # Initial conditions in reduced units
    #betas = [5., 4., 3., 2., 1., 0.5, 0.25, 0.1]
    betas=[0.005]
    mass = 48.
    lim = 10                                # allowed k states
    mu = 0.1                                     # chemical potential
    nsweep = 30000
    dist='bose' #'boltz' for boltzmann distribution starting lattice, 'rand' for random, 'bose' for bose-einstein
    energy_model='harmonic' #'free', 'harmonic' defaults to free particle
    lat_dat = 'lat.dat'
    scalar_dat = 'scalar.dat'
    nconf = 10
    num_atoms = 333
    ntherm = 1
    therm_fmt = '%8d %14.8f %14.8f\n'
    therm_header = '# isweep  e\n'
    energy_table = np.zeros([lim-1,lim-1,lim-1])
    for i in range(lim-1):
        for j in range(lim-1):
            for k in range(lim-1):
                if energy_model =='free':
                    energy_table[i,j,k] = ((i+1)**2 + (j+1)**2 + (k+1)**2)
                elif energy_model =='harmonic':
                    energy_table[i,j,k] = (i + j + k +3/2)
                else:
                    energy_table[i,j,k] = ((i+1)**2 + (j+1)**2 + (k+1)**2)
    for beta in betas:
        # Initialize the lattice at specified temperature
        lattice, num_atoms = initialize_lattice(beta, lim, num_atoms,energy_table,dist)
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

        msg = 'starting %s sweeps using temperature %s and number of particles %s with a %s starting distribution'% (nsweep, 1/beta, num_atoms, dist)
        if restart:
            msg = 're'+msg
        etot = energy(lattice, lim,energy_table)
        from time import clock
        start = clock()
        print(msg)
        naccept = 0
        nattempt = 0
        indexes = np.ones([lim-1,lim-1,lim-1,3],dtype=int)
        for i in range(lim-1):
            for j in range(lim-1):
                for k in range(lim-1):
                    indexes[i,j,k]=[i,j,k]
        for isweep in range(nsweep):
            if (isweep % ntherm == 0):
                etot = energy(lattice, lim,energy_table)
                occu = lattice[0][0][0]
                fsca.write(therm_fmt % (isweep, etot, occu))
            # site = np.random.randint(0, lim-1, (1, 3))[0]
            # while lattice[site[0] , site[1], site[2]] < 1:
            #     site = np.random.randint(0, lim -1, (1, 3))[0]
            truth = lattice>=1
            choices = indexes[np.where(truth)]
            size_choice = np.shape(choices)[0]
            choice = np.random.choice(size_choice)
            site = choices[choice]
            
            neighbors = neighbor_list(site[0], site[1], site[2], lim-1)
            neighbor = random.choice(neighbors)
            lattice_temp = np.copy(lattice)
            if lattice_temp[neighbor[0], neighbor[1], neighbor[2]] + 1 >= lim: #prevents a move to an energy above lim
                continue
            lattice_temp[site[0],site[1], site[2]] = lattice_temp[site[0],site[1], site[2]] - 1
            lattice_temp[neighbor[0], neighbor[1], neighbor[2]] = lattice_temp[neighbor[0], neighbor[1], neighbor[2]] + 1
            prob = move_probability(lattice, lattice_temp, site, neighbor, beta, lim, energy_table)
            nattempt += 1
            rand = np.random.random()
            if rand < prob:
                naccept += 1
                lattice[site[0], site[1], site[2]] = lattice[site[0], site[1], site[2]] - 1
                lattice[neighbor[0], neighbor[1], neighbor[2]] = lattice[neighbor[0], neighbor[1], neighbor[2]] + 1
            if isweep%100 ==0:
                print(isweep)
        accept_ratio = naccept/nattempt
        print(lattice[0][0][0]/num_atoms, energy(lattice, lim,energy_table),accept_ratio)