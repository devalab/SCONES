SDM_S350_CSV_PATH = "data/S350/S350_SDM2.csv"
PremPS_S350_CSV_PATH = "data/S350/S350_PremPS.txt"
PDB_DUMP_PATH = "data/PDB"
LEGACY_PDB_DUMP_PATH = "data/LegacyPDB"
UNIPROT_DUMP_PATH = "data/UniProt"
RESSEQ_DUMP_PATH = "data/RESSEQ"

from datasets.S350DatasetGenerator import S350DatasetGenerator
generator = S350DatasetGenerator()
generator.read_S350(SDM_S350_CSV_PATH)
generator.read_S350_PremPS(PremPS_S350_CSV_PATH)
generator.download_pdb_files(PDB_DUMP_PATH)
generator.download_pdb_files(LEGACY_PDB_DUMP_PATH, file_format="pdb")
generator.download_uniprot(UNIPROT_DUMP_PATH)
generator.create_resseq(PDB_DUMP_PATH, RESSEQ_DUMP_PATH)
generator.calculate_seq_from_resseq(RESSEQ_DUMP_PATH)
generator.validate_data(RESSEQ_DUMP_PATH)
generator.export("data/S350/S350_SCONES.csv")