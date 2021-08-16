from utils.get_aligned_position import get_aligned_position

PH_DELTA = 0.1
T_DELTA = 3
ALIGNMENT_SCORE_MIN = 0.98

alignment_cache = {}

class SampleKey:
    """Represents essential properties of a sample required for equality comparisions"""

    def __init__(self, primary_rp: str, ref_residue: str, position: int, mutated_residue: str, pH: float, T: float):
        self.primary_rp = primary_rp
        self.mutation = (ref_residue, position, mutated_residue)
        self.pH = pH
        self.T = T

    @staticmethod
    def is_pH_same(pH1: float, pH2: float) -> bool:
        return abs(pH1 - pH2) < PH_DELTA

    @staticmethod
    def is_T_same(T1: float, T2: float) -> bool:
        return abs(T1 - T2) < T_DELTA

    def __eq__(self, other):
        primary_rp1, (ref_residue1, position1, mutated_residue1) = self.primary_rp, self.mutation
        primary_rp2, (ref_residue2, position2, mutated_residue2) = other.primary_rp, other.mutation

        if ref_residue1 != ref_residue2 or mutated_residue1 != mutated_residue2:
            return False

        if not self.is_pH_same(self.pH, other.pH) or not self.is_T_same(self.T, other.T):
            return False

        if primary_rp1 == primary_rp2:
            if position1 != position2:
                return False
        else:
            cache_key = (primary_rp1, position1, primary_rp2)
            result_tuple = alignment_cache.get(cache_key, None)
            if result_tuple is None:
                result_tuple = get_aligned_position(primary_rp1, position1, primary_rp2)
                alignment_cache[cache_key] = result_tuple
            score, new_position2 = result_tuple
            
            if score is None or score < ALIGNMENT_SCORE_MIN:
                return False
            assert(primary_rp1[position1] == primary_rp2[position2])
            if new_position2 != position2:
                return False
        return True

    def __hash__(self):
        ref_residue, _, mutated_residue = self.mutation
        return int(ord(ref_residue) + ord(mutated_residue) * 100)

    def __repr__(self):
        ref_residue, position, mutated_residue = self.mutation
        repr = "Key:\n" +                                                           \
               "\tpH: %.2f, T: %.2f\n" % (self.pH, self.T) +                        \
               "\t%c -> %c @ %d\n" % (ref_residue, mutated_residue, position) +     \
               "\t%s\n" % self.primary_rp
        return repr
