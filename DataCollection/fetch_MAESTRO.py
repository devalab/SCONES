MAESTRO_DIR="/MAESTRO_linux_x64/"
LEGACY_PDB_DUMP="../DataGeneration/data/LegacyPDB/"

def fetch_MAESTRO(sample):
    import os
    pdb_path = os.path.join(LEGACY_PDB_DUMP, "pdb" + sample["pdb_id"] + ".ent")
    maestro_bin = os.path.join(MAESTRO_DIR, "maestro")
    maestro_config = os.path.join(MAESTRO_DIR, "config.xml")

    mutation = "%s%d.%s{%s}" % (sample['ref_residue'], sample["resseq_position"], sample["chain_id"], sample["mutated_residue"])    

    import numpy as np
    if 'pH' in sample and not np.isnan(sample["pH"]):
        cmd = "%s %s %s --evalmut=\"%s\" --bu --ph=%.2f" % (maestro_bin, maestro_config, pdb_path, mutation, float(sample["pH"]))
    else:
        cmd = "%s %s %s --evalmut=\"%s\" --bu" % (maestro_bin, maestro_config, pdb_path, mutation)

    import subprocess
    result = subprocess.check_output(cmd, shell=True)

    from io import BytesIO
    strio = BytesIO(result)

    import pandas as pd
    df = pd.read_csv(strio, sep="\t", header=0, index_col=False)
    df = df.set_index("mutation")
    print(df.to_markdown())
    row = df.iloc[1]
    assert(len(df) == 2)
    
    ddG = row.ddG
    ddG_confidence = row.ddG_confidence

    return { "ddG_MAESTRO" : ddG, "ddG_confidence" : ddG_confidence }

from DataFetch import DataFetch
fetcher = DataFetch()
fetcher.fetch("../DataGeneration/data/Stransitive/Stransitive.csv", "results/MAESTRO_Stransitive.csv", fetch_MAESTRO, web=False)

fetcher = DataFetch()
fetcher.fetch("../DataGeneration/data/S768/S768_SCONES.csv", "results/MAESTRO_S768.csv", fetch_MAESTRO, web=False)

fetcher = DataFetch()
fetcher.fetch("../DataGeneration/data/S350/S350_SCONES.csv", "results/MAESTRO_S350.csv", fetch_MAESTRO, web=False)