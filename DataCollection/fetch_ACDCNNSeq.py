# https://github.com/compbiomed-unito/acdc-nn/
ACDCNN_PROFILES="/acdc-nn/results_replication/profiles/"

def fetch_ACDCNNSeq(sample):
    import os
    profile_path = os.path.join(ACDCNN_PROFILES, sample["pdb_id"] + sample["chain_id"] + ".prof.gz")
    mutation = sample['ref_residue'] + str(sample["resseq_position"]) + sample["mutated_residue"]

    import subprocess
    p = subprocess.Popen("acdc-nn seq %s %s" % (mutation, profile_path), stdout=subprocess.PIPE, shell=True)
    stdout, stderr = p.communicate()
    print("stdout:", stdout)
    print("stderr:", stderr)

    ddG = float(stdout)
    print("ddG:", ddG)

    return { "ddG_ACDCNNSeq" : ddG }

from DataFetch import DataFetch
fetcher = DataFetch()
fetcher.fetch("../DataGeneration/data/Stransitive/Stransitive.csv", "results/ACDCNNSeq_Stransitive.csv", fetch_ACDCNNSeq, web=False)