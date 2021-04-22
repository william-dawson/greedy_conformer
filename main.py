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


def write_babel(mol, fname):
    """
    Write an xyz file from an openbabel molecule
    """
    from openbabel.openbabel import OBConversion

    conv = OBConversion()
    conv.SetOutFormat("xyz")

    sval = conv.WriteString(mol)
    with open(fname, "w") as ofile:
        ofile.write(sval)


def energy(mol):
    """
    Compute the energy of a molecule using GAFF.
    """
    from openbabel.openbabel import OBForceField
    ff = OBForceField.FindForceField("gaff")
    ff.Setup(mol)

    return ff.Energy()


def max_deviation(mol1, mol2):
    """
    Return the furthest atom distance.
    """
    from numpy import array
    from numpy.linalg import norm
    plist = []
    numlist = []

    for i in range(mol2.NumAtoms()):
        at = mol2.GetAtom(i+1)
        plist.append([at.GetX(), at.GetY(), at.GetZ()])
        numlist.append(at.GetAtomicNum())

    minlist = []
    for i in range(mol1.NumAtoms()):
        at = mol1.GetAtom(i+1)
        pos = array([at.GetX(), at.GetY(), at.GetZ()])
        anum = at.GetAtomicNum()

        minlist.append(min([norm(pos - pos2)
                           for pos2, num2 in zip(plist, numlist)
                           if anum == num2]))

    return max(minlist)


def filter(systems, energies, rmsd_cutoff, max_cutoff, energy_cutoff,
           align_dir):
    """
    Filter the systems based on an RMSD criteria.
    """
    from openbabel.openbabel import OBAlign, OBMol

    failures = {}
    passing = {}

    align = OBAlign(True, True)

    for f1 in tqdm(list(systems)):
        mol1 = systems[f1]
        found = False
        align.SetRefMol(mol1)
        for f2, mol2 in passing.items():
            if f1 == f2:
                continue
            ediff = 100 * abs(energies[f1] - energies[f2]) / abs(energies[f1])
            if ediff > energy_cutoff:
                continue
            align.SetTargetMol(mol2)
            align.Align()
            rmsd = align.GetRMSD()
            if rmsd < rmsd_cutoff:
                align_mol = OBMol(mol2)
                align.UpdateCoords(align_mol)
                maxdev = max_deviation(mol1, align_mol)
                if maxdev < max_cutoff:
                    failures[(f1, f2)] = {"rmsd": rmsd, "maxdev": maxdev,
                                          "ediff": ediff}
                    write_babel(align_mol, join(align_dir, f2))
                    found = True
                    break
        if not found:
            passing[f1] = mol1

    return passing, failures


if __name__ == "__main__":
    from sys import argv
    from glob import glob
    from os.path import join, basename, exists
    from tqdm import tqdm
    from yaml import load, SafeLoader

    # Get the input path
    if len(argv) < 2:
        raise Exception("Please specify the configure file.")
    with open(argv[1]) as ifile:
        parameters = load(ifile, Loader=SafeLoader)

    if not exists(parameters["aligned_directory"]):
        raise Exception("directory", parameters["aligned_directory"],
                        "does not exist.")

    flist = glob(join(parameters["input"], "*"))

    # Read in the systems
    print("Reading in systems")
    systems = {}
    for i in tqdm(range(len(flist))):
        f = flist[i]
        systems[basename(f)] = read_babel(f)

    # Compute the energies
    print("Computing Energies")
    energies = {}
    for f in tqdm(list(systems)):
        energies[f] = energy(systems[f])

    # Filter
    print("Filtering")
    passing, failures = filter(systems, energies,
                               float(parameters["rmsd_cutoff"]),
                               float(parameters["max_distance_cutoff"]),
                               float(parameters["energy_cutoff"]),
                               parameters["aligned_directory"])

    print("Writing Output")
    # Write out the failures
    with open(join("failures.txt"), "w") as ofile:
        ofile.write("System 1\tSystem 2\t")
        ofile.write("RMSD\tMax Dist\tEnergy Difference (%)\n")
        for k, v in failures.items():
            ofile.write(k[0] + "\t" + k[1] + "\t")
            ofile.write(str(v["rmsd"]) + "\t")
            ofile.write(str(v["maxdev"]) + "\t")
            ofile.write(str(v["ediff"]) + "\n")

    # Write out the successes
    with open(join("successes.txt"), "w") as ofile:
        ofile.write("\n".join(list(passing)))

    # Write out the energies
    with open(join("energies.txt"), "w") as ofile:
        for k, ene in energies.items():
            ofile.write(k + "\t" + str(ene) + "\n")
