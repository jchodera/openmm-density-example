#!/usr/bin/env python
"""
Compute density of liquid chloroform using OpenMM

"""

################################################################################
# IMPORTS
################################################################################

from simtk import unit, openmm
from simtk.openmm import app

################################################################################
# PARAMETERS
################################################################################

temperature = 298.0 * unit.kelvin
pressure = 1.0 * unit.atmospheres
timestep = 2.0 * unit.femtoseconds
collision_rate = 5.0 / unit.picoseconds
cutoff = 11.0*unit.angstroms
constraints = app.HBonds
nsteps = 10000
nsteps_per_snapshot = 500
nsteps_per_density = 50
output_pdbfile = 'output.pdb'
density_outfile = 'statedata.out'

################################################################################
# MAIN
################################################################################

# Read initial snapshot
pdbfile = app.PDBFile('67-66-3-liq.pdb')

# Read gromacs input file and create system
print('Reading gromacs topology and creating system...')
top = app.GromacsTopFile('67-66-3-liq.top', periodicBoxVectors=pdbfile.topology.getPeriodicBoxVectors(), includeDir='.')
system = top.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=cutoff, constraints=constraints)

# Add a barostat
print('Adding barostat...')
barostat = openmm.MonteCarloBarostat(pressure, temperature)
system.addForce(barostat)

# Create simulation
print('Creating simulation...')
integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
simulation = app.Simulation(top.topology, system, integrator)
simulation.context.setPositions(pdbfile.positions)

# Minimize energy
print('Minimizing energy...')
simulation.minimizeEnergy()

# Append reporters
print('Creating reporters...')
simulation.reporters.append(app.PDBReporter(output_pdbfile, nsteps_per_snapshot))
simulation.reporters.append(app.StateDataReporter(density_outfile, nsteps_per_density, step=True, density=True))

# Run the simulation
print('Running the simulation...')
simulation.step(nsteps)
