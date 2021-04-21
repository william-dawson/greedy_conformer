"""
Greedily prune a set of conformers.

Usage:
    python main.py input 0.1
"""


def read_babel(f):
    from openbabel.openbabel import OBMol, OBConversion

    conv = OBConversion()
    conv.SetInFormat("xyz")

    with open(f) as ifile:
        sval = "".join([x for x in ifile])

    mol = OBMol()
    conv.ReadString(mol, sval)

    return mol


def filter(systems):
    from openbabel.openbabel import OBAlign

    failures = {}
    passing = {}

    passing[flist[0]] = systems[flist[0]]
    for f1 in tqdm(list(systems)):
        mol1 = systems[f1]
        found = False
        for f2, mol2 in passing.items():
            if f1 == f2:
                continue
            align = OBAlign(mol1, mol2, True, True)
            align.Align()
            score = align.GetRMSD()
            if score < cutoff:
                failures[f1 + " " + f2] = score
                found = True
                break
        if not found:
            passing[f1] = mol1

    return passing, failures


if __name__ == "__main__":
    from sys import argv
    from glob import glob
    from os.path import join
    from tqdm import tqdm

    # Get the input path
    if len(argv) < 3:
        raise Exception("Requires input path and cutoff")
    flist = glob(join(argv[1], "*"))
    cutoff = float(argv[2])

    # Read in the systems
    print("Reading in systems")
    systems = {}
    for i in tqdm(range(len(flist))):
        f = flist[i]
        systems[f] = read_babel(f)

    # Filter
    print("Filtering")
    passing, failures = filter(systems)

    print("Writing Output")
    # Write out the failures
    with open(join("failures.txt"), "w") as ofile:
        for k, rmsd in failures.items():
            ofile.write(k + "\t" + str(rmsd) + "\n")

    # Write out the successes
    with open(join("successes.txt"), "w") as ofile:
        ofile.write("\n".join(list(passing)))
