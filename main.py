import numpy as np
from matplotlib import pyplot as plt


def initialize_lattice(num_atoms, beta):
    lattice = np.empty(num_atoms)

    # Initialize each site with energy such that average
    # is given temperature
    
    # Take into account energy levels and degeneracy for this?
    return lattice

def energy(lattice):
    pass

if __name__ == '__main__':

    # Initial conditions in reduced units
    beta = 1.
    mass = 48.
    num_atoms = 10

    # Initialize the lattice at specified temperature
    lattice = initialize_lattice(num_atoms, beta)

    pass
