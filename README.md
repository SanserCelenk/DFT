# DFT

%%writefile LDA_Li_graphite_h.py


import numpy as np
from ase import Atom, Atoms
from ase.io import Trajectory
from ase.calculators.emt import EMT
from ase.lattice.hexagonal import Graphene
from ase.parallel import paropen
from gpaw import GPAW


xc='LDA' #xc-functional. Change to 'PBE' or 'BEEF-vdw' to run with other functionals

d = 1.56  # C-C distance. 
a = d * np.sqrt(3)  # Translate C-C distance to lattice parameter

traj = Trajectory('{}/LiC6_h.traj'.format(xc), 'w')

calc = GPAW(h=0.18,
         kpts=(5, 5, 6),
         txt='{}/result_LiC6_h.txt'.format(xc),
         xc=xc)

with paropen('{}/energies_LiC6_h.dat'.format(xc), 'w') as f:
    
    for h in np.linspace(3.5, 3.75, 8):
        
        
        system = Atoms('C6Li', 
                       positions=[[0, 0, 0],        # Carbon
                                  [0, d, 0],        # Carbon
                                  [a, 0, 0],        # Carbon
                                  [-a, 0, 0],       # Carbon
                                  [-a/2, -d/2, 0],  # Carbon
                                  [a/2, -d/2, 0],   # Carbon
                                  [0, -d, h/2]],   # Lithium
                       cell=([1.5*a, -1.5*d, 0], [1.5*a, 1.5*d, 0] , [0, 0, h]),
                       pbc=True) # Enable periodic boundary conditions
        system.wrap()
        
        system.set_calculator(calc) 
        
        energy = system.get_potential_energy() 
        
        traj.write(system)  
        
        print('Interlayer spacing: {:.3f}\tEnergy: {:.5f} eV'.format(h, energy), file=f)
        f.flush()  
