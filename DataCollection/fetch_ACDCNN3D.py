# https://github.com/compbiomed-unito/acdc-nn/
ACDCNN_PDBS="/acdc-nn/results_replication/pdbs/"
ACDCNN_PROFILES="/acdc-nn/results_replication/profiles/"

def fetch_ACDCNN3D(sample):
    import os
    pdb_path = os.path.join(ACDCNN_PDBS, sample["pdb_id"] + ".pdb.gz")
    profile_path = os.path.join(ACDCNN_PROFILES, sample["pdb_id"] + sample["chain_id"] + ".prof.gz")
    mutation = sample['ref_residue'] + str(sample["resseq_position"]) + sample["mutated_residue"]

    import subprocess
    p = subprocess.Popen("acdc-nn struct %s %s %s %c" % (mutation, profile_path, pdb_path, sample["chain_id"]), stdout=subprocess.PIPE, shell=True)
    stdout, stderr = p.communicate()
    print("stdout:", stdout)
    print("stderr:", stderr)

    ddG = float(stdout)
    print("ddG:", ddG)

    return { "ddG_ACDCNN3D" : ddG }

from DataFetch import DataFetch
fetcher = DataFetch()
fetcher.fetch("../DataGeneration/data/Stransitive/Stransitive.csv", "results/ACDCNN3D_Stransitive.csv", fetch_ACDCNN3D, web=False)