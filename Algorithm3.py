import numpy as np
import random
import scipy
from matplotlib import pyplot as plt
from time import clock

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

def initialize_lattice(beta, lim, n_atoms,energy_table,dist):

    # Lattice all 1's to start (not sure if it matters)
    lattice = np.zeros([lim, lim, lim])

    # Fill according to selected distribution
    #if dist == 'boltz':
    #    lattice[:]=np.exp(-energy_table[:]*beta)

    #elif dist == 'rand':
    #    lattice[:]=np.random.random()

    #else: # bose-einstein distribution
    #    lattice[:]= 1 / (np.exp(beta * (energy_table[:])) - 1)

    # Maybe not the best thing to do here?
    #lattice = np.round(lattice/np.sum(lattice) * n_atoms)

    #if(n_atoms != np.sum(lattice)):
        #print(lattice)
        #exit(1)
    lattice[0][0][0] = n_atoms

    #if(np.sum(lattice) != n_atoms):
    #    print(lattice)
    #    exit(1)

    return lattice, np.sum(lattice)


def energy(lattice, lim, energy_table):
    en = 0
    for i in range(lim):
        for j in range(lim):
            for k in range(lim):
                en = en + (lattice[i,j,k] * energy_table[i,j,k])
    return en


def energy_difference(energy_table, s, n):
    return energy_table[n[0]][n[1]][n[2]] - energy_table[s[0]][s[1]][s[2]]


def move_probability(lattice, lattice_temp, site, neighbor, beta, lim,energy_table):
    energy1 = energy(lattice, lim,energy_table)
    energy2 = energy(lattice_temp, lim,energy_table)
    oc1 = lattice[site[0] , site[1], site[2]]
    oc2 = lattice[neighbor[0] , neighbor[1], neighbor[2]]
    return oc2 / oc1 * np.exp(-beta *(energy2 - energy1))

if __name__ == '__main__':

    # Initial conditions in reduced units
    betas = np.array([100, 10, 5., 4., 3., 2., 1., 0.5, .25, 0.1, 0.05, 0.01, 0.005, 0.001])
    #betas = np.array([1.0])
    mass = 48.
    lim = 50                                # allowed k states
    mu = 0.1                                     # chemical potential
    nsweep = 500
    num_atoms = 50

    # 'boltz' for boltzmann distribution starting lattice, 'rand' for random,
    # 'bose' for bose-einstein
    dist='bose' 
    energy_model='free' #'free', 'harmonic' defaults to free particle
    lat_dat = 'lat.dat'
    scalar_dat = 'scalar.dat'
    nconf = 10
    ntherm = 10
    therm_fmt = '%8d %14.8f %14.8f\n'
    therm_header = '# isweep  e\n'
    energy_table = np.zeros([lim,lim,lim])
    highest_val = -1

    # Initialize energy table
    for i in range(lim):
        for j in range(lim):
            for k in range(lim):
                if energy_model =='free':
                    energy_table[i,j,k] = ((i+1)**2 + (j+1)**2 + (k+1)**2)
                elif energy_model =='harmonic':
                    energy_table[i,j,k] = (i + j + k +3/2)
                else:
                    energy_table[i,j,k] = ((i+1)**2 + (j+1)**2 + (k+1)**2)


    # Loop over each beta
    av_occ = []
    av_en = []
    for beta in betas:

        # Yeet
        en = []
        occ_num = []

        # Initialize the lattice at specified temperature
        lattice, num_atoms = initialize_lattice(beta, lim, num_atoms, energy_table, dist)
        restart = False

        with open(lat_dat, 'w') as flat:
            flat.write('# n_atoms=%d\n' % num_atoms)
        with open(scalar_dat, 'w') as fsca:
            fsca.write(therm_header)
        flat = open(lat_dat, 'a')
        fsca = open(scalar_dat, 'a')

        msg = '%s sweeps, temperature %s, %s particles, %s distribution'% \
              (nsweep, 1/beta, num_atoms, dist)
        if restart:
            msg = 're'+msg
        etot = energy(lattice, lim,energy_table)
        start = clock()
        print(msg)
        naccept = 0
        nattempt = 0
        indexes = np.ones([lim,lim,lim,3],dtype=int)

        lm = ''
        for i in range(lim):
            for j in range(lim):
                for k in range(lim):
                    indexes[i,j,k]=[i,j,k]

        for isweep in range(nsweep):

            max_val = np.max(np.where(lattice))

            # Check highest state that gets occupied
            if(max_val > highest_val):
                highest_val = max_val

            if(isweep % (nsweep/20) == 0):
                left = '['
                right = ']'
                rm = ' '*(20 - len(lm))
                print("{}{}{}{} {}%".format(left,lm,rm,right, isweep/nsweep*100), end='\r')
                lm += '='

            for atom in range(int(num_atoms)):
                #if (isweep % ntherm == 0):
                #    etot = energy(lattice, lim,energy_table)
                #    occu = lattice[0][0][0]
                #    fsca.write(therm_fmt % (isweep, etot, occu))

                # Picks out energies that have at least one particle
                truth = lattice>=1

                # Index of non zero occupancies, need to weight by particle number
                choices = indexes[np.where(truth)]
                choice = np.ceil(num_atoms*np.random.random())

                #print("Choices: {}".format(choices))
                particles = 0
                for c in choices:
                    particles += lattice[c[0]][c[1]][c[2]]
                    if(choice <= particles):
                        site = c
                        break
                
                # Choose neighbor
                neighbors = neighbor_list(site[0], site[1], site[2], lim)

                # Get rid of unphysical neighbors
                neighbors[np.sum((neighbors - site), axis=1) > 1] = [0,0,0]

                # Get rid of neighbors that are the same as site
                neighbors = neighbors[np.sum(neighbors==[site[0],site[1],site[2]],axis=1)!=3]

                # Select neighbor
                neighbor = random.choice(neighbors)

                # Calculate energy difference
                energy_diff = energy_table[neighbor[0]][neighbor[1]][neighbor[2]] - \
                              energy_table[site[0]][site[1]][site[2]]

                # Get transition probability
                old_oc = lattice[site[0]][site[1]][site[2]] - 1
                new_oc = lattice[neighbor[0]][neighbor[1]][neighbor[2]] + 1
                prob =  np.exp(-beta*energy_diff)

                # Move if accepted
                acc = np.random.random()
                if(prob > acc):
                    lattice[site[0]][site[1]][site[2]] -= 1
                    lattice[neighbor[0]][neighbor[1]][neighbor[2]] += 1
                    naccept += 1

                nattempt += 1

            occ_num.append(lattice[0][0][0]/num_atoms)
            en.append(np.sum(lattice * energy_table))

        #print(len(en))
        accept_ratio = naccept/nattempt

        np.save("./run_dat/en_trace_{}_{}_{}".format(int(num_atoms), beta, energy_model), en)
        np.save("./run_dat/occ_trace_{}_{}_{}".format(int(num_atoms), beta, energy_model),
                                                      occ_num)
        np.save("./run_dat/acc_rate_{}_{}_{}".format(int(num_atoms), beta, energy_model),
                                                     np.array(accept_ratio))

        print("Occupation number: {}, Energy: {}, Acceptance Ratio: {}, Atoms: {}\n".format(
              lattice[0][0][0]/num_atoms, energy(lattice, lim, energy_table), accept_ratio,
              np.sum(lattice)))
        #exit(1)
        #occ_num.append(lattice[0][0][0]/num_atoms)
        av_occ.append(np.average(occ_num[int(nsweep/2):]))
        av_en.append(np.average(en[int(nsweep/2):]))
    
    fig, ax = plt.subplots(2)
    print(av_en)
    print(av_occ)
    ax[0].plot(1/betas, av_en)
    ax[1].plot(1/betas, av_occ)
    #plt.show()

