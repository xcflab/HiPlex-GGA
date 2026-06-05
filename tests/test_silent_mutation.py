from oligos_gg_assembly.sequence_utils import translate
from oligos_gg_assembly.silent_mutation import mutate_restriction_sites


def test_mutate_restriction_sites_preserves_translation_and_removes_bsai():
    seq = "GCTGGTCTCGCT"
    changed = mutate_restriction_sites(seq, ("BsaI5", "BsaI3"))
    assert "GGTCTC" not in changed
    assert translate(changed) == translate(seq)

