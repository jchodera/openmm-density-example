# Example illustrating the calculation of chloroform density using gromacs input files

## Prereqiusites

* Install [`miniconda`](http://conda.pydata.org/miniconda.html)
* Install OpenMM and other dependencies
```bash
conda config --add channels omnia
conda intall --yes openmm mdtraj numpy pymbar matplotlib seaborn
```
* Run the simulation (on a machine with a GPU)
```bash
python compute-density.py
```

## Manifest

* `trichloromethane.itp` - chloroform parameter file
* `67-66-3-liq.pdb` - initial liquid snapshot
* `grompp_LIQ.mdp` - gromacs MDP input file (unused)
