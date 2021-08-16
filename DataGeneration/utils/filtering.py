def is_similar_easy(ref_residue1: str, ref_residue2: str, mutant_residue1: str, mutant_residue2: str, seqstrip1: str, seqstrip2: str) -> bool:
    is_same_fwd = (ref_residue1 == ref_residue2 and mutant_residue1 == mutant_residue2)
    is_same_bwd = (ref_residue1 == mutant_residue2 and ref_residue2 == mutant_residue1)

    if is_same_fwd or is_same_bwd:
        from Bio import pairwise2
        alignments = pairwise2.align.globalms(
            seqstrip1, seqstrip2,
            1, -1, 0, 0,
            one_alignment_only=True
        )

        alignment = alignments[0]
        seqA = alignment.seqA
        seqB = alignment.seqB
        seqlen = len(seqA)
        assert(len(seqA) == len(seqB))
        assert(alignment.start == 0 and alignment.end == seqlen)

        normalized_alignment_score = alignment.score / seqlen
        if normalized_alignment_score > 0.5:
            return True
    return False

def is_similar_hard(ref_residue1: str, ref_residue2: str, mutant_residue1: str, mutant_residue2: str, seqstrip1: str, seqstrip2: str) -> bool:
    from Bio import pairwise2
    alignments = pairwise2.align.globalms(
        seqstrip1, seqstrip2,
        1, -1, 0, 0,
        one_alignment_only=True
    )

    alignment = alignments[0]
    seqA = alignment.seqA
    seqB = alignment.seqB
    seqlen = len(seqA)
    assert(len(seqA) == len(seqB))
    assert(alignment.start == 0 and alignment.end == seqlen)

    normalized_alignment_score = alignment.score / len(seqstrip1)
    if normalized_alignment_score > 0.5:
        return True
    return False

def remove_similar_samples(train_dataset: dict, test_dataset: dict, similarity_function) -> bool:
    dataset = {}
    for idx, sample1 in train_dataset.items():
        primary1 = sample1["sequence"]
        position1 = sample1["seq_position"]
        ref_residue1 = sample1["ref_residue"]
        mutant_residue1 = sample1["mutated_residue"]
        primary1 = list(primary1)
        primary1[position1] = ref_residue1
        primary1 = "".join(primary1)

        start1 = max(0, position1 - 5)
        end1 = min(len(primary1), position1 + 5)
        seqstrip1 = primary1[start1 : end1]

        removed = False
        for sample2 in test_dataset.values():
            primary2 = sample2["sequence"]
            position2 = sample2["seq_position"]
            ref_residue2 = sample2["ref_residue"]
            mutant_residue2 = sample2["mutated_residue"]
            primary2 = list(primary2)
            primary2[position2] = ref_residue2
            primary2 = "".join(primary2)

            start2 = max(0, position2 - 5)
            end2 = min(len(primary2), position2 + 5)
            seqstrip2 = primary2[start2 : end2]

            if similarity_function(ref_residue1, ref_residue2, mutant_residue1, mutant_residue2, seqstrip1, seqstrip2):
                removed = True
                break

        if not removed:
            dataset[idx] = sample1
    return dataset
