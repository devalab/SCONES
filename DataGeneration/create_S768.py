S768_CSV_PATH = "data/S768/Data_Sheet_2_Prediction of Protein Mutational Free Energy_ Benchmark and Sampling Improvements Increase Classification Accuracy.CSV"
PDB_DUMP_PATH = "data/PDB"
LEGACY_PDB_DUMP_PATH = "data/LegacyPDB"
UNIPROT_DUMP_PATH = "data/UniProt"
RESSEQ_DUMP_PATH = "data/RESSEQ"

from datasets.S768DatasetGenerator import S768DatasetGenerator
generator = S768DatasetGenerator()
generator.read_S768(S768_CSV_PATH)
generator.download_pdb_files(PDB_DUMP_PATH)
generator.download_pdb_files(LEGACY_PDB_DUMP_PATH, file_format="pdb")
generator.download_uniprot(UNIPROT_DUMP_PATH)
generator.create_resseq(PDB_DUMP_PATH, RESSEQ_DUMP_PATH)
generator.calculate_seq_from_resseq(RESSEQ_DUMP_PATH)
generator.validate_data(RESSEQ_DUMP_PATH)
generator.export("data/S768/S768_SCONES.csv")