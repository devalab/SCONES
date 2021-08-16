import os

def extract_resseq(pdb_path: str, chain_id: str, fasta_path: str, resseq_path: str) -> None:
    from Bio.PDB.MMCIFParser import MMCIFParser
    pdb_id = os.path.basename(pdb_path).split('.')[0]
    struct = MMCIFParser(QUIET=True).get_structure(structure_id=pdb_id, filename=pdb_path)

    model = struct[0]
    chain = model[chain_id]

    from Bio.PDB import is_aa
    residues = [r for r in chain.get_residues() if is_aa(r)]
    resseq_ids = [r.get_id()[1] for r in residues]

    from Bio.SeqUtils import seq1
    sequence = "".join([seq1(r.get_resname()) for r in residues])

    from Bio.Seq import Seq
    sequence = Seq(sequence)

    seqid = pdb_id + ':' + chain_id

    from Bio import SeqRecord, SeqIO
    record = SeqRecord.SeqRecord(seq=sequence, id=seqid, description='')
    with open(fasta_path, 'w') as f:
        SeqIO.write(record, f, 'fasta')

    with open(resseq_path, 'w') as f:
        f.write('> ' + seqid + '\n')
        f.write(','.join(str(j) for j in resseq_ids))