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

def initialize_lattice(lim, n_atoms):
    lattice = np.zeros((lim, lim, lim))
    lattice[0][0][0] = n_atoms
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


def de_broglie(mass, beta):
    return (6.626*10**(-34))/np.sqrt(2*mass*np.pi*(1./beta)) # ??? CHECK CONSTANT


def get_mu(mass, beta, volume, num_atoms):
    return 1./beta * np.log(de_broglie(mass, beta)**3) * num_atoms/volume


if __name__ == '__main__':

    # Initial conditions in reduced units
    betas = np.array([50, 10, 5., 4., 3., 2., 1., 0.75, 0.7, 0.6, 0.5, 0.4, 0.3, .25, 0.1, 
                      0.08, 0.05, 0.03, 0.01, 0.008, 0.005])
    betas = np.array([0.1])

    mass = 6.6465 * 10**(-27)  # Heluim
    length = 10*10**(-9)
    volume = (length)**(3)  # Because it just feels right
    lim = 10                                # allowed k states
    nsweep = 1000
    #n_atoms = [5, 10, 50, 100, 200, 500]
    n_atoms = [1000]

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
                    energy_table[i,j,k] = (1.054*10**(-34))**2 * np.pi**2/(2*mass*(length**2))* \
                                          ((i+1)**2 + (j+1)**2 + (k+1)**2)

                elif energy_model =='harmonic':
                    energy_table[i,j,k] = (i + j + k +3/2)

    # Loop over each beta
    av_occ = []
    av_en = []
    start = clock()

    for num_atoms in n_atoms:
        init_num_atoms = num_atoms
        for beta in betas:

            # Yeet
            en = []
            occ_num = []
            atom_trace = []

            # Initialize the lattice at specified temperature
            lattice, num_atoms = initialize_lattice(lim, init_num_atoms)
            restart = False

            with open(lat_dat, 'w') as flat:
                flat.write('# n_atoms=%d\n' % num_atoms)
            with open(scalar_dat, 'w') as fsca:
                fsca.write(therm_header)
            flat = open(lat_dat, 'a')
            fsca = open(scalar_dat, 'a')

            msg = '%s sweeps, temperature %s, %s particles'% \
                  (nsweep, 1/beta, num_atoms)
            if restart:
                msg = 're'+msg
            etot = energy(lattice, lim,energy_table)
            start = clock()
            print(msg)
            naccept = 0
            nattempt = 0
            indexes = np.ones([lim,lim,lim,3],dtype=int)

            for i in range(lim):
                for j in range(lim):
                    for k in range(lim):
                        indexes[i,j,k]=[i,j,k]

            lm = ''
            print("LATTICE SIZE: {}".format(np.sum(lattice)))
            for isweep in range(nsweep):

                #max_val = np.max(np.where(lattice))

                # Check highest state that gets occupied
                #if(max_val > highest_val):
                #    highest_val = max_val

                if(isweep % (nsweep/20) == 0):
                    left = '['
                    right = ']'
                    rm = ' '*(20 - len(lm))
                    print("{}{}{}{} {}%".format(left,lm,rm,right, isweep/nsweep*100), end='\r')
                    lm += '='

                for atom in range(int(init_num_atoms)):
                    #if (isweep % ntherm == 0):
                    #    etot = energy(lattice, lim,energy_table)
                    #    occu = lattice[0][0][0]
                    #    fsca.write(therm_fmt % (isweep, etot, occu))

                    # 50% chance of trying to exchange
                    if(np.random.random() < 0.5):

                        # Picks out energies that have at least one particle
                        truth = lattice>=1

                        # Index of non zero occupancies, need to weight by particle number
                        choices = indexes[np.where(truth)]
                        choice = np.ceil(num_atoms*np.random.random())

                        particles = 0
                        for c in choices:
                            particles += lattice[c[0]][c[1]][c[2]]
                            if(choice <= particles):
                                site = c
                                break
                        
                        # Choose neighbor
                        neighbors = neighbor_list(site[0], site[1], site[2], lim)

                        # Get rid of unphysical neighbors
                        neighbors[np.sum((neighbors - site), axis=1) > 1] = site

                        # Get rid of neighbors that are the same as site
                        neighbors = neighbors[np.sum(neighbors==[site[0],site[1],site[2]],
                                              axis=1)!=3]

                        # Select neighbor
                        neighbor = random.choice(neighbors)

                        # Calculate energy difference
                        energy_diff = energy_table[neighbor[0]][neighbor[1]][neighbor[2]] - \
                                      energy_table[site[0]][site[1]][site[2]]

                        # Get transition probability
                        old_oc = lattice[site[0]][site[1]][site[2]]
                        new_oc = lattice[neighbor[0]][neighbor[1]][neighbor[2]] + 1
                        
                        # Maybe different factor out front?
                        prob =  new_oc / old_oc * np.exp(-beta*energy_diff)

                        # Move if accepted
                        acc = np.random.random()
                        if(prob > acc):
                            lattice[site[0]][site[1]][site[2]] -= 1
                            lattice[neighbor[0]][neighbor[1]][neighbor[2]] += 1
                            naccept += 1

                        nattempt += 1

                    else: # Gonna particle exchange
                        #print(get_mu(mass, beta, volume, num_atoms), num_atoms)

                        if(np.random.random() < 0.5): # 50/50 to insert particle

                            # select site to insert, (random, then calculates energy)
                            site = np.random.randint(0, lim, (1,3))[0]
                            energy_diff = energy_table[site[0]][site[1]][site[2]]

                            # compute probability
                            #print(volume/de_broglie(mass, beta)**3 * (num_atoms+1))
                            #print(np.exp(beta * (get_mu(mass, beta, volume, num_atoms)+\
                            #                    energy_diff)))
                            prob = volume/(de_broglie(mass, beta)**3 * (num_atoms+1)) * \
                                   np.exp(beta * (get_mu(mass, beta, volume, num_atoms) + \
                                                  energy_diff))
                            #print("Insertion probability: {}".format(prob))

                            # insert or not
                            if(prob > np.random.random()):
                                lattice[site[0]][site[1]][site[2]] += 1
                                #print("INSERTED")
                                num_atoms = np.sum(lattice)

                        else: # delete particle

                            # select site to delete, (random)
                            # Index of non zero occupancies, need to weight by particle number
                            truth = lattice>=1
                            choices = indexes[np.where(truth)]
                            choice = np.ceil(num_atoms*np.random.random())

                            particles = 0
                            for c in choices:
                                particles += lattice[c[0]][c[1]][c[2]]
                                if(choice <= particles):
                                    site = c
                                    break

                            energy_diff = -energy_table[site[0]][site[1]][site[2]]

                            # compute probability
                            #print(de_broglie(mass, beta)**3 * num_atoms/volume)
                            #print(np.exp(-beta * (get_mu(mass, beta, volume, num_atoms) + 
                            #                      energy_diff)))
                            prob = de_broglie(mass, beta)**3 * num_atoms/volume * \
                                   np.exp(-beta * (get_mu(mass, beta, volume, num_atoms) - \
                                                  energy_diff))
                            #print(prob)
                            #print("Deletion probability: {}".format(prob))

                            #print("DELETION PROBABILITY: {}".format(prob))
                            # insert or not
                            if(prob > np.random.random()):
                                #print("DELETED")
                                lattice[site[0]][site[1]][site[2]] -= 1
                                num_atoms = np.sum(lattice)

                occ_num.append(lattice[0][0][0]/num_atoms)
                en.append(np.sum(lattice * energy_table))
                atom_trace.append(num_atoms)

                #try:
                #    if((en[isweep] == en[isweep-1]) and (en[isweep] == en[isweep-2])):
                #        print(lattice)
                #        exit(1)
                #except:
                #    print(isweep)

            #fig, ax = plt.subplots()
            #ax.plot(occ_num)
            #plt.show()
            #print()
            #print(naccept/nattempt)
            #print(np.average(occ_num[100:]))

            #print(len(en))
            try:
                accept_ratio = naccept/nattempt
            except:
                print("No moves attempted")

            np.save("./gc_run_dat/en_trace_{}_{}_{}".format(int(init_num_atoms), beta,
                                                   energy_model), en)
            np.save("./gc_run_dat/occ_trace_{}_{}_{}".format(int(init_num_atoms), beta,
                                                   energy_model), occ_num)
            np.save("./gc_run_dat/acc_rate_{}_{}_{}".format(int(init_num_atoms), beta,
                                                   energy_model), np.array(accept_ratio))
            np.save("./gc_run_dat/atom_trace_{}_{}_{}".format(int(init_num_atoms), beta,
                                                   energy_model), atom_trace)

            print("Occupation number: {}, Energy: {}, Acceptance Ratio: {}, Atoms: {}".format(
                  lattice[0][0][0]/num_atoms, energy(lattice, lim, energy_table), accept_ratio,
                  np.sum(lattice)))
            print("Number of atoms: {}\n".format(num_atoms))
            end = clock()
            #print("Total time: {}".format(end-start))
            #exit(1)
            #occ_num.append(lattice[0][0][0]/num_atoms)
            av_occ.append(np.average(occ_num[int(nsweep/2):]))
            av_en.append(np.average(en[int(nsweep/2):]))

            fig, ax = plt.subplots()
            ax.plot(atom_trace)
            plt.show()
            exit(1)
    
    fig, ax = plt.subplots(2)
    #print(av_en)
    #print(av_occ)
    ax[0].plot(1/betas, av_en)
    ax[1].plot(1/betas, av_occ)
    plt.show()


