
# https://github.com/biofold/ddgun
DDGUN3D_SCRIPT="/ddgun/ddgun_3d.py"
DDGUN_ENV = "ddgun" # DDGun requires python 2.7

LEGACY_PDB_DUMP="../DataGeneration/data/LegacyPDB/"

def fetch_DDGun3D(sample):
    import os
    pdb_path = os.path.join(LEGACY_PDB_DUMP, "pdb" + sample["pdb_id"] + ".ent")

    import tempfile
    mutation_file = tempfile.NamedTemporaryFile().name

    mutation = sample['ref_residue'] + str(sample["resseq_position"]) + sample["mutated_residue"]
    with open(mutation_file, "w") as f:
        f.write(mutation)
    
    result_path = os.path.join("dump", sample["pdb_id"] + '_' + mutation)
    if not os.path.exists(result_path):
        cmd = "%s %s %c %s -o %s" % (DDGUN3D_SCRIPT, pdb_path, sample["chain_id"], mutation_file, result_path)
        os.system(". ~/miniconda3/etc/profile.d/conda.sh && conda activate %s && %s" % (DDGUN_ENV, cmd))
    with open(result_path, "r") as f:
        line = f.readlines()[1]
    ddG = float(line.split()[3])
    return { "ddG_DDGun3D" : ddG }

from DataFetch import DataFetch
fetcher = DataFetch()
fetcher.fetch("../DataGeneration/data/Stransitive/Stransitive.csv", "results/DDGun3D_Stransitive.csv", fetch_DDGun3D, web=False)

fetcher = DataFetch()
fetcher.fetch("../DataGeneration/data/S768/S768_SCONES.csv", "results/DDGun3D_S768.csv", fetch_DDGun3D, web=False)

fetcher = DataFetch()
fetcher.fetch("../DataGeneration/data/S350/S350_SCONES.csv", "results/DDGun3D_S350.csv", fetch_DDGun3D, web=False)