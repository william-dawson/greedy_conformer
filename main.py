"""
Greedily prune a set of conformers.

Usage:
    python main.py input 0.1 10
"""


def read_babel(f):
    """
    Read in an xyz file into an openbabel molecule.
    """
    from openbabel.openbabel import OBMol, OBConversion

    conv = OBConversion()
    conv.SetInFormat("xyz")

    with open(f) as ifile:
        sval = "".join([x for x in ifile])

    mol = OBMol()
    conv.ReadString(mol, sval)

    return mol


def energy(mol):
    """
    Compute the energy of a molecule using GAFF.
    """
    from openbabel.openbabel import OBForceField
    ff = OBForceField.FindForceField("gaff")
    ff.Setup(mol)

    return ff.Energy()


def filter(systems, energies, rmsd_cutoff, energy_cutoff):
    """
    Filter the systems based on an RMSD criteria.
    """
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
            ediff = 100 * abs(energies[f1] - energies[f2]) / abs(energies[f1])
            if ediff > energy_cutoff:
                continue
            align = OBAlign(mol1, mol2, True, True)
            align.Align()
            score = align.GetRMSD()
            if score < rmsd_cutoff:
                failures[f1 + "\t" + f2] = str(score)+"\t"+str(ediff)
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
    if len(argv) < 4:
        raise Exception("Requires input path, rmsd cutoff, energy cutoff")
    flist = glob(join(argv[1], "*"))
    rmsd_cutoff = float(argv[2])
    energy_cutoff = float(argv[3])

    # Read in the systems
    print("Reading in systems")
    systems = {}
    for i in tqdm(range(len(flist))):
        f = flist[i]
        systems[f] = read_babel(f)

    # Compute the energies
    print("Computing Energies")
    energies = {}
    for f in tqdm(list(systems)):
        energies[f] = energy(systems[f])

    # Filter
    print("Filtering")
    passing, failures = filter(systems, energies, rmsd_cutoff, energy_cutoff)

    print("Writing Output")
    # Write out the failures
    with open(join("failures.txt"), "w") as ofile:
        ofile.write("System 1\tSystem 2\tRMSD\tEnergy Difference (%)\n")
        for k, rmsd in failures.items():
            ofile.write(k + "\t" + str(rmsd) + "\n")

    # Write out the successes
    with open(join("successes.txt"), "w") as ofile:
        ofile.write("\n".join(list(passing)))

    # Write out the energies
    with open(join("energies.txt"), "w") as ofile:
        for k, ene in energies.items():
            ofile.write(k + "\t" + str(ene) + "\n")
