import os
import pandas as pd
import typing

class DatasetGenerator:
    """Base class to perform common tasks on all datasets such as downloading PDB structures and UniProt sequences.

    Each dataset will have a separate dataset generator class that derives from this class. The derived classes
    will handle dataset specific details and populate the `self.data` object using the 'add' method.
    The 'data' member object is a dictionary that is indexed by a key that can uniquely identify every sample in
    the dataset. Every sample is interally stored as a dictionary with predefined fieldnames for common attributes.
    The end-user can invoke the base class methods to download the mmCIF files, UniProt sequences, preprocess and validate them.
    """

    def __init__(self):
        self.data = {}

    def add(self, key,
        pdb_id: typing.Optional[str], chain_id: typing.Optional[str],
        ref_residue: typing.Optional[str], mutated_residue: typing.Optional[str],
        uniprot_id: typing.Optional[str]=None, sequence: typing.Optional[str]=None,
        seq_position: typing.Optional[int]=None, resseq_position: typing.Optional[int]=None,
        pH: typing.Optional[float]=None, T: typing.Optional[float]=None,
        ddG: typing.Optional[float]=None, dTm: typing.Optional[float]=None,
        extra: dict={}
    ) -> None:
        """Adds a sample to the dataset.

        Every sample is internally stored as a dictionary. All the positional arguments are stored with predefined
        attribute names. The user may optionally provide additional dataset specific attributes in `extra` argument. 

        Args:
            key: hashable key that can uniquely identify the sample
            pdb_id: PDB id (will be converted to lowercase internally)
            chain_id: PDB chain id (will be converted to uppercase internally)
            ref_residue: single-letter amino acid code for the reference residue
            mutated_residue: single-letter amino acid code for the mutated residue
            uniprot_id: UniProt sequence id (will be converted to uppercase internally) (Optional)
            sequence: string of single-letter amino acid codes representing the protein (Optional)
            seq_position: position of the mutating residue in the'sequence' argument (Optional)
            resseq_position: position of the mutating residue according to resseq id in the PDB file (Optional)
            pH: pH at which the sample was taken (Optional)
            T: temperature at which the sample was taken (Optional)
            ddG: target dG change (Optional)
            dTm: target Tm change (Optional)
            extra: dictionary of additional key-value pairs (Optional)

        Notes:
            - one of 'seq_position' or 'resseq_position' must always be given
            - 'pdb_id' must be provided if 'resseq_position' is given
            - 'sequence' must be given if 'seq_position' is given
        """
        
        assert(seq_position is not None or resseq_position is not None)
        assert(resseq_position is None or pdb_id is not None)
        assert(seq_position is None or sequence is not None)

        if key in self.data:
            raise ValueError("Duplicate Key: " + key)

        sample = { "id" : key }
        if pdb_id is not None:
            assert(chain_id is not None)
            sample["pdb_id"] = pdb_id.lower()
            sample["chain_id"] = chain_id.upper()
        if uniprot_id is not None:
            sample["uniprot_id"] = uniprot_id.upper()
        if sequence is not None:
            sample["sequence"] = sequence
        sample["ref_residue"] = ref_residue
        if seq_position is not None:
            sample["seq_position"] = int(seq_position)
        if resseq_position is not None:
            sample["resseq_position"] = int(resseq_position)
        sample["mutated_residue"] = mutated_residue

        if pH is not None:
            sample["pH"] = pH
        if T is not None:
            sample["T"] = T
        if ddG is not None:
            sample["ddG"] = ddG
        if dTm is not None:
            sample["dTm"] = dTm

        if extra:
            sample = {**sample, **extra}

        self.data[key] = sample

    def download_pdb_files(self, dump_path: str, file_format: str="mmCif") -> None:
        """Download the structure in mmCIF format using the 'pdb_id' field.
        
        Only considers samples which have 'pdb_id' field. The field may contain '-' to
        indicate no PDB id.

        Args:
            dump_path: directory where the mmCIF files must be downloaded to
        """

        pdb_ids = [sample["pdb_id"] for sample in self.data.values() if "pdb_id" in sample]
        missing_ids = pdb_ids.count('-')
        if missing_ids:
            print("[SKIP] PDB >> %d samples do not have pdb_id" % missing_ids)
        pdb_ids = list(filter(('-').__ne__, pdb_ids))
   
        from Bio.PDB import PDBList
        pdb_list = PDBList()
        pdb_list.download_pdb_files(pdb_ids, obsolete=False, pdir=dump_path, file_format=file_format)
        pdb_list.download_pdb_files(pdb_ids, obsolete=True, pdir=dump_path, file_format=file_format)

    def download_uniprot(self, dump_path: str) -> None:
        """Download UniProt sequences using 'uniprot_id' field and save them in FASTA format.
        
        Only considers samples which have 'uniprot_id' field. The field may contain '-' to
        indicate no UniProt id.

        Args:
            dump_path: directory where the sequences in FASTA file format must be downloaded to
        """

        uniprot_ids = [sample["uniprot_id"] for sample in self.data.values() if "uniprot_id" in sample]
        missing_ids = uniprot_ids.count('-')
        if missing_ids:
            print("[SKIP] UNIPROT >> %d samples do not have uniprot_id" % missing_ids)
        uniprot_ids = list(filter(('-').__ne__, uniprot_ids))

        for uniprot_id in set(uniprot_ids):
            fasta_path = os.path.join(dump_path, uniprot_id + ".fasta")
            if not os.path.exists(fasta_path):
                from utils.retrive_uniprot_sequence import retrive_uniprot_sequence
                status = retrive_uniprot_sequence(uniprot_id, fasta_path)
                if not status:
                    print("[SKIP] UNIPROT >> Failed to retrive sequence for", uniprot_id)

    def create_resseq(self, pdb_dir_path: str, resseq_dump_dir_path: str) -> None:
        """Creates FASTA and RESSEQ files for each chain.
        
           The FASTA file contains the polypetide sequence obtained from structural information.
           The RESSEQ file, also in FASTA format, contains a mapping of sequence positions to resseq ids.
        
        Args:
            pdb_dir_path: directory where the mmCIF files can be found
            resseq_dump_dir_path: directory where the FASTA and RESSEQ files must be stored
        """

        extended_pdb_ids = set([sample["pdb_id"] + sample["chain_id"] for sample in self.data.values() if "pdb_id" in sample])
        for extended_pdb_id in extended_pdb_ids:
            pdb_id, chain_id = extended_pdb_id[:4], extended_pdb_id[4]
            pdb_path = os.path.join(pdb_dir_path, pdb_id + ".cif")
            if not os.path.exists(pdb_path):
                print("[SKIP] PDB structure not found for", extended_pdb_id)
                continue

            fasta_path = os.path.join(resseq_dump_dir_path, extended_pdb_id +  ".fasta")
            resseq_path = os.path.join(resseq_dump_dir_path, extended_pdb_id + ".resseq")
            if not os.path.exists(fasta_path) or not os.path.exists(resseq_path):
                from utils.extract_resseq import extract_resseq
                extract_resseq(pdb_path, chain_id, fasta_path, resseq_path)

            if not os.path.exists(fasta_path) or not os.path.exists(resseq_path):
                print("[SKIP] could not create RESSEQ/FASTA for", extended_pdb_id)
                continue

    def calculate_seq_from_resseq(self, resseq_dump_dir_path: str) -> None:
        """Calculates `sequence` and `seq_position` fields from `resseq_position` field.
        
        This method stores the structure sequence in the `sequence` field and calculates
        the `seq_position` that corresponds to the `resseq_position`.
        
        Args:
            resseq_dump_dir_path: directory where the FASTA and RESSEQ files can be found
        """

        for sample in self.data.values():
            pdb_id, chain_id = sample["pdb_id"], sample["chain_id"]
            extended_pdb_id = pdb_id + chain_id

            fasta_path = os.path.join(resseq_dump_dir_path, extended_pdb_id +  ".fasta")
            resseq_path = os.path.join(resseq_dump_dir_path, extended_pdb_id + ".resseq")

            from Bio import SeqIO
            sequence = next(SeqIO.parse(fasta_path, "fasta"))
            sequence = str(sequence.seq)

            resseq = next(SeqIO.parse(resseq_path, "fasta"))
            resseq = str(resseq.seq)
            resseq = [int(i) for i in resseq.split(',')]
            resseq = {k : v for v, k in enumerate(resseq)}

            resseq_position = sample["resseq_position"]
            seq_position = resseq[resseq_position]

            sample["sequence"] = sequence
            sample["seq_position"] = seq_position

    def calculate_resseq_from_seq(self, resseq_dump_dir_path: str) -> None:
        """Calculates `resseq_position` from `seq_position` field.
        
        This method computes the `resseq_position` that corresponds to `seq_position`.
        Note that it is not possible to unambigiously calculate `resseq_position` in
        some cases. The first consistent position number is used instead.
        
        Args:
            resseq_dump_dir_path: directory where the FASTA and RESSEQ files can be found
        """

        for sample in self.data.values():
            pdb_id, chain_id = sample["pdb_id"], sample["chain_id"]
            extended_pdb_id = pdb_id + chain_id

            fasta_path = os.path.join(resseq_dump_dir_path, extended_pdb_id +  ".fasta")
            resseq_path = os.path.join(resseq_dump_dir_path, extended_pdb_id + ".resseq")

            from Bio import SeqIO
            resseq_sequence = next(SeqIO.parse(fasta_path, "fasta"))
            resseq_sequence = str(resseq_sequence.seq)

            seq_sequence = sample["sequence"]
            seq_position = sample["seq_position"]

            from utils.get_aligned_position import get_aligned_position
            score, pdb_position = get_aligned_position(seq_sequence, seq_position, resseq_sequence)
            if score == None:
                print("[FAIL] failed to calculate resseq position from sequence position for", sample["id"])
                continue

            resseq = next(SeqIO.parse(resseq_path, "fasta"))
            resseq = str(resseq.seq)
            resseq = [int(i) for i in resseq.split(',')]
            assert(len(resseq) == len(resseq_sequence))

            resseq_position = resseq[pdb_position]
            sample["resseq_position"] = resseq_position

    def validate_data(self, resseq_dump_dir_path: typing.Optional[str]) -> None:
        """Verify that the mutations and positions in the samples are valid for their pdb_id and chain_id.
        
        Args:
            resseq_dump_dir_path: directory where the FASTA and RESSEQ files can be found (Optional)

        `resseq_dump_dir_path` is only required if samples have `resseq_position` field.
        """

        for sample in self.data.values():
            ref_residue = sample["ref_residue"]
            mutated_residue = sample["mutated_residue"]
            assert(ref_residue != mutated_residue)

            if "resseq_position" in sample:
                pdb_id, chain_id = sample["pdb_id"], sample["chain_id"]
                extended_pdb_id = pdb_id + chain_id

                fasta_path = os.path.join(resseq_dump_dir_path, extended_pdb_id +  ".fasta")
                resseq_path = os.path.join(resseq_dump_dir_path, extended_pdb_id + ".resseq")

                from Bio import SeqIO
                sequence = next(SeqIO.parse(fasta_path, "fasta"))
                sequence = str(sequence.seq)

                resseq = next(SeqIO.parse(resseq_path, "fasta"))
                resseq = str(resseq.seq)
                resseq = [int(i) for i in resseq.split(',')]
                assert(len(resseq) == len(sequence))

                resseq = {k : v for v, k in enumerate(resseq)}

                resseq_position = sample["resseq_position"]
                assert(resseq_position in resseq)
                seq_position = resseq[resseq_position]

                seq_residue = sequence[seq_position]
                if seq_residue != ref_residue:
                    print("[SKIP] residue mismatch >> extended_pdb_id: %s, resseq_position: %d, seq_position: %d, expected_residue: %c, seq_residue: %c" % (extended_pdb_id, resseq_position, seq_position, ref_residue, seq_residue))
                    print("\t", sequence)

            if "seq_position" in sample:
                sequence = sample["sequence"]
                seq_position = sample["seq_position"]
                seq_residue = sequence[seq_position]
                if seq_residue != ref_residue:
                    print("[SKIP] residue mismatch >> seq_position: %d, expected_residue: %c, seq_residue: %c" % (seq_position, ref_residue, seq_residue))
                    print("\t", sequence)           

    def export(self, csv_path: str, columns: typing.Optional[list]=None) -> None:
        """Export data to CSV format.
        
        Args:
            csv_path: path to CSV file
            columns: field names (and their order) that have to be stored
        """

        records = self.data.values()
        df = pd.DataFrame.from_dict(records)
        if columns:
            for column in columns:
                assert(column in df.columns)
            df = df.reindex(columns, axis=1)
            assert(df.columns.size == len(columns))
        df.to_csv(csv_path, index=False)

