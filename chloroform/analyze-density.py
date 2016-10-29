#!/usr/local/env python

"""
Analyze density data.

"""

# Work with no DISPLAY set
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

# Imports
from simtk import unit
from pymbar import timeseries
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from simtk import unit

#"Step","Time (ps)","Potential Energy (kJ/mole)","Temperature (K)","Box Volume (nm^3)","Density (g/mL)","Speed (ns/day)"
data = np.genfromtxt('statedata.out', delimiter=',')

# Slice out density data
time_t = data[:,1] * unit.picoseconds 
density_t = data[:,5] * unit.grams/unit.milliliter
print('%.3f ns of simulation data read' % (time_t.max() / unit.nanoseconds))

# Compute expectation of density using <rho> = <M/V>
density_unit = (unit.grams/unit.milliliter)
print('Detecting equilibration...')
[t0, g, Neff] = timeseries.detectEquilibration(density_t / density_unit, nskip=10)
print('Discarding initial %d samples to equilibration (g = %.1f, Neff = %.1f)' % (t0, g, Neff))
density = (density_t / density_unit)[t0:].mean()
ddensity = (density_t / density_unit)[t0:].std() / np.sqrt(Neff)
print('density estimate from <M/V> = %.6f +- %.6f g/cm3' % (density, ddensity))

# Compute expectation of density using M/<V> = 1/<V/M> = 1/<1/rho>
print('Detecting equilibration...')
A_t = 1.0 / (density_t / density_unit)
[t0, g, Neff] = timeseries.detectEquilibration(A_t, nskip=10)
print('Discarding initial %d samples to equilibration (g = %.1f, Neff = %.1f)' % (t0, g, Neff))
EA = A_t[t0:].mean()
dEA = A_t[t0:].std() / np.sqrt(Neff)
density2 = 1.0 / EA
ddensity2 = dEA / EA
print('density estimate from M/<V> = %.6f +- %.6f g/cm3' % (density2, ddensity2))

# Plot
sns.set_style('white')
fig = plt.figure(figsize=[10,8])
plt.plot(time_t / unit.nanoseconds, density_t / (unit.grams/unit.milliliter), '.')
plt.hold(True)
plt.plot([0, (time_t / unit.nanoseconds).max()], [density, density], 'r-')
plt.plot([0, (time_t / unit.nanoseconds).max()], [density2, density2], 'b-')
plt.xlabel('time (ns)')
plt.ylabel('instantaneous density (g/mL)')
sns.despine()
fig.savefig('density.png')
