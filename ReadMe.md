# Greedy Conformer

Given a directory of xyz files which represent conformers of a given molecule,
this will look for duplicate conformers and prune them. This program employs
a simple greedy algorithm to do the pruning, and assess similarity in terms
of energy, RMSD, and maximum deviation of a given atom.

This code depends on the following packages: openbabel (with python bindings),
tqdm, and pyyaml. All of these can be installed with pip.

## Usage

First, create a configure file such as the included exaple config.yaml. Then
you can run the code as:
```
python main.py config.yaml
```
The `config.yaml` file has the following parameters:
* input: the path to the directory of xyz files to analyze.
* rmsd_cutoff: the cutoff value when comparing RMSD in angstrom.
* max_distance_cutoff: the cutoff value when comparing the maximum distance
  between atoms in angstrom.
* energy_cutoff: the cutoff value for the relative energy difference between 
  two structures (in %).
* aligned_directory: when a duplicate is detected, the structure is aligned
  to what it was matched to, and written as a file in this directory.

After the program runs, the following files are produced:
* failures.txt - a list of conformer pairs that were too similar, and the
  related metrics.
* successes.txt - 
* energies.txt - the energies of each system.

