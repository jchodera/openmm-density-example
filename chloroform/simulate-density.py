#!/usr/bin/env python
"""
Compute density of liquid chloroform using OpenMM

"""

################################################################################
# IMPORTS
################################################################################

from simtk import unit, openmm
from simtk.openmm import app
import numpy as np

################################################################################
# PARAMETERS
################################################################################

temperature = 298.15 * unit.kelvin
pressure = 1.0 * unit.atmospheres
timestep = 1.0 * unit.femtoseconds
collision_rate = 5.0 / unit.picoseconds
switch_width = 0.5*unit.angstroms
cutoff = 11.0*unit.angstroms
constraints = app.HBonds
simulation_length = 100 * unit.nanoseconds
nsteps = int(np.ceil(simulation_length / timestep))
nsteps_per_snapshot = 50000
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

# Ensure we are using a switching function
forces = { system.getForce(index).__class__.__name__ : system.getForce(index) for index in range(system.getNumForces()) }
forces['NonbondedForce'].setUseSwitchingFunction(True)
forces['NonbondedForce'].setSwitchingDistance(cutoff - switch_width)

# Ensure we are using mixed precision
print('Ensuring mixed precision...')
# Select platform automatically; use mixed precision
integrator = openmm.VerletIntegrator(timestep)
context = openmm.Context(system, integrator)
platform = context.getPlatform()
del context
platform.setPropertyDefaultValue('Precision', 'mixed')

# Create simulation
print('Creating simulation...')
integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
simulation = app.Simulation(top.topology, system, integrator, platform)
simulation.context.setPositions(pdbfile.positions)

# Minimize energy
print('Minimizing energy...')
simulation.minimizeEnergy()

# Append reporters
print('Creating reporters...')
simulation.reporters.append(app.PDBReporter(output_pdbfile, nsteps_per_snapshot))
simulation.reporters.append(app.StateDataReporter(density_outfile, nsteps_per_density, step=True, density=True))

# Run the simulation
print('Running the simulation for %d steps...' % nsteps)
simulation.step(nsteps)
