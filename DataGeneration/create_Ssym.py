import os
import pandas as pd
import numpy as np


SSYM_CSV_PATH = "data/Ssym/Ssym.csv"
PREMPS_SSYM_CSV_PATH = "data/Ssym/Ssym_PremPS.txt"
LEGACY_PDB_DUMP_PATH = "data/LegacyPDB"
PDB_DUMP_PATH = "data/PDB"
UNIPROT_DUMP_PATH = "data/UniProt"
RESSEQ_DUMP_PATH = "data/RESSEQ"

from datasets.SsymDatasetGenerator import SsymDatasetGenerator
generator = SsymDatasetGenerator()
generator.read_Ssym(SSYM_CSV_PATH)
generator.read_Ssym_PremPS(PREMPS_SSYM_CSV_PATH)
generator.download_pdb_files(PDB_DUMP_PATH)
generator.download_pdb_files(LEGACY_PDB_DUMP_PATH, file_format="pdb")
generator.download_uniprot(UNIPROT_DUMP_PATH)
generator.create_resseq(PDB_DUMP_PATH, RESSEQ_DUMP_PATH)
generator.calculate_seq_from_resseq(RESSEQ_DUMP_PATH)
generator.validate_data(RESSEQ_DUMP_PATH)
generator.export("data/Ssym/Ssym_SCONES.csv", "data/Ssym/Ssym_SCONES_idx.csv")