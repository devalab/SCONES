def get_aligned_position(sequence1: str, position1: int, sequence2: str) -> tuple:
    """Finds position of a residue in sequence2 given the position in sequence1.

       The two sequences are expected to be similar. An arbitary alignment parameters
       are used to match the two sequences. The algorithm fails whenever the two sequences
       can be definitely said to be dissimilar: residue mismatch. Otherwise, it greedily
       finds the lowest consistent position of the residue in sequence2.

       Suppose position1 is 4 (zero based indexing) in sequence1. Position2 in sequence2
       can be 4 or 5 or 6 but the algorithm will choose the lowest consistent position.
       ALKTA--YLRL
       |||||  ||||
       ALKTAAAYLRL
    """

    from Bio import pairwise2
    alignments = pairwise2.align.globalms(
        sequence1, sequence2,
        1, -1, -3, -1,
        penalize_end_gaps=False,
        one_alignment_only=True
    )

    alignment = alignments[0]
    seqA = alignment.seqA
    seqB = alignment.seqB
    seqlen = len(seqA)
    assert(len(seqA) == len(seqB))
    assert(alignment.start == 0 and alignment.end == seqlen)

    mismatch = False
    posA, posB = 0, 0 # positions in sequence1 and sequence2 = how many residues have been seen so far
    aligned_position, position2 = 0, 0
    for i in range(len(seqA)):
        if seqA[i] != '-':
            posA += 1
        if seqB[i] != '-':
            posB += 1
        if not (seqA[i] == '-' or seqB[i] == '-') and seqA[i] != seqB[i]:
            mismatch=True
            break
        if posA == position1 + 1:
            position2 = posB - 1
            aligned_position = i

    if mismatch:
        return None, None

    if seqA[aligned_position] == '-' or seqB[aligned_position] == '-':
        return None, None

    if sequence1[position1] != sequence2[position2]:
        return None, None

    score = alignment.score / min(len(sequence1), len(sequence2))
    return score, position2
