def retrive_uniprot_sequence(uniprot_id: str, fasta_path: str) -> None:
    import requests as r
    response = r.post("http://www.uniprot.org/uniprot/" + uniprot_id + ".fasta")
    data = ''.join(response.text)

    if response.status_code != 200:
        return False
    
    from Bio import SeqIO
    from io import StringIO
    records = SeqIO.parse(StringIO(data), "fasta")
    SeqIO.write(records, fasta_path, "fasta")
    return True