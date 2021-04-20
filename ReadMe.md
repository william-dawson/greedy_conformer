# Greedy Conformer

This will greedily prune a set of conformers. To use it:
```
python main.py input 0.1
```
Inside `input` should be your list of xyz files. The second argument is
the RMSD value that you use as a tolerance. The output will be written
to failures.txt and successes.txt. Failures details which file matched
to something else in the set.

To use this you need to install the python bindings for openbabel
(either with pip or conda). The first time you run the script might
be slow as it JIT compiles the code

