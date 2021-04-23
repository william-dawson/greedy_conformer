# Greedy Conformer

Given a directory of xyz files which represent conformers of a given molecule,
this will look for duplicate conformers and prune them. This program employs
a simple greedy algorithm to do the pruning, and assess similarity in terms
of energy, RMSD, and maximum deviation of a given atom.

This code depends on the following packages: openbabel (with python bindings),
tqdm, numpy, and pyyaml. All of these can be installed with pip.

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
* successes.txt - the list of files that don't match anything else.
* energies.txt - the energies of each system.

## OpenBabel Memory Exception

When using large input files, it is possible to encounter the following error
```
==============================
*** Open Babel Error  in operator()
  memory limit exceeded...
```
When aligning molecules, openbabel will first produce a graph representation
of each molecule. It will then try to generate a mapping between the two graphs
and try aligning for all of those mappings. This is called the
"Graph Isomorphism" problem and belongs to the class NP. This is thus a very
hard problem to solve and I believe they rely on enumerating all possibilities.
Fearing that the number of conformers might grow too big, exhausting your
system's memory and leading to a crash, they have hard coded a memory limit.

To overcome this, first download openbabel from the 
[main repository](https://github.com/openbabel/openbabel). Then, open up
the file `include/openbabel/isomorphism.h`. On lines 140, 226, and 236, you
will see the value `30000000` corresponding to 300MB. You should increase
this, to say 4 or 8 GB. Then, you need to compile openbabel:
```
mkdir build
cd build
cmake .. -DPYTHON_BINDINGS=Yes -DRUN_SWIG=ON -DCMAKE_INSTALL_PREFIX=/path/to/install
make
make install
```
I recommend installing it somewhere locally so it doesn't conflict with a
more stable version of openbabel. Then, before you run this script, you should
set the environment variable:
```
export PYTHONPATH=/path/to/install/lib/python3.7/site-packages/$PYTHONPATH
```
The value `python3.7` will depend on your system. Just do ls in the 
install directory and you will see the right value. Then you can run
`main.py` as usual.
