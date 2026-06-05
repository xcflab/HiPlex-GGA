"""Shared DNA sequence utilities."""

CODON_TABLE = {
    "ATA": "I", "ATC": "I", "ATT": "I", "ATG": "M",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
    "AAC": "N", "AAT": "N", "AAA": "K", "AAG": "K",
    "AGC": "S", "AGT": "S", "AGA": "R", "AGG": "R",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CAC": "H", "CAT": "H", "CAA": "Q", "CAG": "Q",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
    "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
    "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
    "TTC": "F", "TTT": "F", "TTA": "L", "TTG": "L",
    "TAC": "Y", "TAT": "Y", "TAA": "_", "TAG": "_",
    "TGC": "C", "TGT": "C", "TGA": "_", "TGG": "W",
}

PREFERRED_CODONS = {
    "I": ["ATC", "ATT"],
    "M": ["ATG"],
    "T": ["ACA", "ACC", "ACG", "ACT"],
    "N": ["AAC", "AAT"],
    "K": ["AAA"],
    "S": ["AGC", "AGT", "TCA", "TCC", "TCG", "TCT"],
    "R": ["CGC", "CGT"],
    "L": ["CTC", "CTG", "CTT", "TTA", "TTG"],
    "P": ["CCA", "CCG", "CCT"],
    "H": ["CAC", "CAT"],
    "Q": ["CAA", "CAG"],
    "V": ["GTA", "GTC", "GTG", "GTT"],
    "A": ["GCA", "GCC", "GCG", "GCT"],
    "D": ["GAC", "GAT"],
    "E": ["GAA", "GAG"],
    "G": ["GGC", "GGG", "GGT"],
    "F": ["TTC", "TTT"],
    "Y": ["TAC", "TAT"],
    "_": ["TAA", "TAG", "TGA"],
    "W": ["TGG"],
}


def normalize_dna(seq: str) -> str:
    return seq.strip().upper()


def translate(seq: str) -> str:
    """Translate a DNA sequence in frame 0."""
    seq = normalize_dna(seq)
    if len(seq) % 3 != 0:
        raise ValueError(f"DNA sequence length must be divisible by 3: {len(seq)}")
    protein = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        try:
            protein.append(CODON_TABLE[codon])
        except KeyError as exc:
            raise ValueError(f"Unsupported codon: {codon}") from exc
    return "".join(protein)


def find_new_codon(aa: str, old_codon: str) -> str | None:
    """Return a synonymous codon using the historical E. coli-biased table."""
    old_codon = normalize_dna(old_codon)
    choices = list(PREFERRED_CODONS.get(aa.upper(), []))
    if not choices:
        return None
    if old_codon in choices and len(choices) > 1:
        choices.remove(old_codon)
    return choices[0]


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(table)[::-1]

