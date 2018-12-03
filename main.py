import numpy as np
import scipy
from matplotlib import pyplot as plt


###
# parameters:
#    num_atoms: number of atoms we're simulating with. Should be cubic number
#    beta: temperature parameter, 1/T
#
# returns:
#    lattice with energy on sites distributed according to BE distribution
###
def initialize_lattice(num_atoms, beta):
    lattice = np.empty(num_atoms)

    # Initialize each site with energy such that average
    # is given temperature
    
    
    # Take into account energy levels and degeneracy for this?
    return lattice


# total energy of lattice
def energy(lattice):
    pass


###
# parameters:
#    lattice: the lattice we're working with
#    site: the site on the lattice being changed
#    e_coord: the energy coordinate being changed e.g. (0,0,1) for z, or (1,0,0) for x
#    pm: energy +- 1
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
    num_atoms = 10

    # Initialize the lattice at specified temperature
    lattice = initialize_lattice(num_atoms, beta)

    pass

