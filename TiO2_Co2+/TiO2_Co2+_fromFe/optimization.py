#This is a heavy calculation that requires the use of a supercomputer. It must be submitted using the run.sh script.
#Import all necessary modules
from ase import io
from ase.optimize import QuasiNewton #geometry optimization algorithm; QuasiNewton links to BFGS line search, which is the best general-purpose optimizer, but other options are available: https://wiki.fysik.dtu.dk/ase/ase/optimize.html
from espresso import espresso
import numpy as np

#initializing variables
slab_path = 'TiO2_Co2+_slab.traj'
output_path = 'results/'
#setup calculator
calcargs = dict(xc='BEEF-vdW',
        kpts=(4, 4, 1), #only need 1 kpt in z-direction
        pw=400.,
        dw=4000.,
        spinpol=True,
        beefensemble=True,
        printensemble=True,
        convergence={'energy':1e-6,
                    'mixing':0.05,
                    'maxsteps':1000,
                    'diag':'david'},
        startingwfc='atomic',
        smearing='fd', #fermi-dirac electron smearing
        sigma=0.1, #smearing width
        dipole={'status':True}, #dipole corrections True turns them on
        #parflags='-nk 2',
        outdir =output_path+'esp.log')
calc = espresso(**calcargs)
#atoms = io.read('POSCAR') #Read in the structure built by the other script
atoms = io.read(slab_path)
atoms.set_calculator(calc)
init_mag = np.zeros(len(atoms))
init_mag[-1] = 0.1
atoms.set_initial_magnetic_moments(magmoms=init_mag) #perturbs the magnetic moment of the system, this is required in spinpolarized systems but should not be done in spin paired systems
relax = QuasiNewton(atoms,logfile=output_path+'opt.log',trajectory=output_path+'opt.json',restart=output_path+'opt.pckl')
relax.run(fmax=0.05) #execute the relaxation algorithm. It will run until the maximum force on any atom is <0.05 eV/Angstrom.
energy = atoms.get_potential_energy() #this is the potential energy of the electrons as computed by DFT. It will be closely related to the enthalpy.
atoms.write(output_path+'converged_slab.traj')
f = open(output_path+'converged.log','w')
f.write(str(energy))
f.close()
