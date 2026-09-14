"""Split DNA sequences into Golden Gate oligos.

This module is based on the historical `oligo_split_251111.py` workflow.
"""

from __future__ import annotations

import argparse
import random
import sys
from dataclasses import dataclass
from importlib import resources
from pathlib import Path

from Bio.Seq import Seq
from tqdm import tqdm

from .add_adapters import add_adapter_list
from .io import (
    read_codon_table,
    read_designs,
    read_sequence_barcodes,
    read_subpool_barcodes,
    write_fragment_table,
    write_name_sequence_table,
    write_oligo_table,
)
from .sequence_utils import find_new_codon, translate

GOOD_OVERHANGS = [
    "TCAA", "ATAA", "TAGA", "GTTA", "AGGG", "CCTA", "AAGA", "TCCA",
    "AGAA", "AAAG", "ACAT", "GGGA", "GCAA", "TAAA", "TGAA", "GACA",
    "ACGA", "GGTA", "ATCC", "ATTG", "CTAC", "AAGC", "CATC", "ACTC",
    "CACA", "CTAA", "GAAA", "AGAC", "AGCA", "CGTA", "ACCG", "CAGA",
    "AACT", "AATA", "GCAC", "CCAG", "CAAG", "AAAT", "ATCA", "CAGG",
    "CATA", "GGAA", "AGGA", "ACGC", "ATAC", "CTCA", "GCCA", "CCGA",
    "ACAG", "AATC", "CAGC", "AAAA", "AGTG", "CGCA", "AACG", "GAGA",
    "ACTA", "TACA", "ATGA", "CGAC", "CGAA", "AGCC"
]

BAD_OVERHANGS = [
    "GGGG", "CCCC", "GGGC", "GCCC", "GGCG", "CGCC", "GCGG", "CCGC",
    "CGGG", "CCCG", "GGCC", "CCGG", "GCCG", "CGGC", "CGCG", "GCGC",
    "ATTA", "TAAT", "AATT", "ACGT", "AGCT", "ATAT", "CATG", "CTAG",
    "GATC", "GTAC", "TATA", "TCGA", "TGCA", "TTAA"
]


@dataclass
class SplitOutputs:
    final_oligos: dict[str, dict[str, str]]
    naked_fragments: dict[str, list[str]]
    whole_dna: dict[str, str]


RESTRICTION_SITES = {
    "NdeI": "CATATG",
    "XhoI": "CTCGAG",
    "HindIII": "AAGCTT",
    "EcoRI": "GAATTC",
    "NcoI": "CCATGG",
    "BamHI": "GGATCC",
    "BsaI5": "GGTCTC",
    "BsmBI5": "CGTCTC",
    "BsaI3": "GAGACC",
    "BsmBI3": "GAGACG",
    "NheI": "GCTAGC",
    "XbaI": "TCTAGA",
    "BbsI5": "GAAGAC",
    "BbsI3": "GTCTTC",
}


class Frag:
    """State for splitting one DNA sequence into fragments."""

    def __init__(
        self,
        frag_num: int,
        dna_seq: str,
        aa_seq: str | None,
        max_len_5: int,
        max_inner_len: int,
        max_len_3: int,
        min_oligo_len: int,
    ):
        self.frag_num = frag_num
        self.dna_seq = dna_seq
        self.max_len_5 = max_len_5
        self.max_inner_len = max_inner_len
        self.max_len_3 = max_len_3
        self.min_oligo_len = min_oligo_len
        self.aa_seq = aa_seq
        self.init_dna_seq = dna_seq
        self.overhang_len = 4
        self.end_pos = 0
        self.cut_order = 0
        self.frag_min_len: int | None = None

    def count_frag_lenth(self, cut_order: int, min_len: int) -> int:
        seq_len = len(self.dna_seq)
        frag_num_left = self.frag_num - cut_order
        frag_min_len = (
            seq_len
            - (self.max_inner_len - 4) * (frag_num_left - 1)
            - (self.max_len_3 - 4)
        )
        if self.frag_num == 8:
            if self.cut_order < 5:
                frag_min_len += 2
            if self.cut_order == 5:
                frag_min_len += 1

        if frag_min_len < min_len:
            frag_min_len = min_len

        if cut_order == 1 and frag_min_len > self.max_len_5:
            raise ValueError(
                "Each fragment should be shorter than the longest length of first "
                f"fragment. Minimum fragment length: {frag_min_len}; longest first "
                f"fragment length: {self.max_len_5}"
            )
        if cut_order != 1 and frag_min_len > self.max_inner_len:
            raise ValueError(
                "Each fragment should be shorter than the longest length of inner "
                f"fragments. Minimum fragment length: {frag_min_len}; longest inner "
                f"fragment length: {self.max_inner_len}"
            )
        return frag_min_len

    def update_dna_seq(self) -> None:
        new_dna_seq_begin = len(self.init_dna_seq) - len(self.dna_seq)
        self.dna_seq = self.init_dna_seq[new_dna_seq_begin:]


def default_codon_table_path() -> Path:
    return Path(str(resources.files("oligos_gg_assembly").joinpath("data/codontable.tab")))


def sort_by_length(input_dict: dict[str, str]) -> dict[str, str]:
    return {
        key: value
        for key, value in sorted(input_dict.items(), key=lambda item: len(item[1]), reverse=True)
    }


def shuffle_dict(input_dict: dict[str, str]) -> dict[str, str]:
    keys = list(input_dict.keys())
    random.shuffle(keys)
    return {key: input_dict[key] for key in keys}


def check_overhang_unique(end_list: list[str | Seq], to_check_end: str | Seq) -> bool:
    return to_check_end not in end_list


def split_with_pos(sequence: str, pos: int, overhang_len: int = 4) -> tuple[str, str, Seq, Seq]:
    seq_1 = sequence[0:pos]
    seq_2 = sequence[pos - overhang_len:]
    overhang = Seq(sequence[pos - overhang_len:pos])
    overhang_pair = overhang.reverse_complement()
    return seq_1, seq_2, overhang, overhang_pair


def get_frame_end(startpoint: int) -> int:
    frameend = startpoint
    if frameend % 3 == 0:
        frameend -= 1
    elif frameend % 3 == 1:
        frameend += 1
    return frameend


def change_the_site(seq: str, enzyme_types_seq: list[str], sites: dict[str, str]) -> str:
    selected = {
        enzyme_type: sites[enzyme_type]
        for enzyme_type in enzyme_types_seq
        if enzyme_type in sites
    }
    enzyme_site_free = False
    while not enzyme_site_free:
        for enzyme, pattern in selected.items():
            while pattern in seq:
                site = seq.find(pattern)
                remainder = site % 3
                if remainder == 0:
                    old_codon = seq[site:site + 3]
                elif remainder == 1:
                    old_codon = seq[site - 1:site + 2] if enzyme == "NcoI" else seq[site + 2:site + 5]
                else:
                    old_codon = seq[site + 1:site + 4]

                aa = translate(old_codon)
                new_codon = find_new_codon(aa, old_codon)
                if new_codon is None:
                    raise ValueError(f"Cannot find synonymous codon for {old_codon}")

                if remainder == 0:
                    seq = seq[:site] + new_codon + seq[site + 3:]
                elif remainder == 1:
                    if enzyme == "NcoI":
                        seq = seq[:site - 1] + new_codon + seq[site + 2:]
                    else:
                        seq = seq[:site + 2] + new_codon + seq[site + 5:]
                else:
                    seq = seq[:site + 1] + new_codon + seq[site + 4:]
        if all(site not in seq for site in sites.values()):
            enzyme_site_free = True
        else:
            selected = {key: value for key, value in sites.items() if value in seq}
    return seq


def replace_codons(s: str, frame_end: int, change_len: int, codons: dict[str, list[str]]) -> str:
    window_size = 6
    change_start = frame_end + 1
    change_end = get_frame_end(frame_end + 1 + change_len) + 2
    window_end = change_end + window_size
    window_start = change_start - window_size

    window_start_seq = s[window_start:change_start]
    window_end_seq = s[change_end + 2:window_end + 2]

    randomized_positions = [i for i in range(change_start // 3, change_end // 3)]
    random.shuffle(randomized_positions)
    poss_to_change = min(6, len(randomized_positions) - 1)
    changed = sorted(randomized_positions[:poss_to_change])

    sh = ""
    for p in range(change_start, change_end, 3):
        if p // 3 in changed:
            codon = s[p:p + 3]
            aa = str(Seq(codon).translate())
            r = random.randint(1, len(codons[aa]))
            sh += codons[aa][r - 1]
        else:
            sh += s[p:p + 3]

    new_seq = window_start_seq + sh + window_end_seq
    if any(site in new_seq for site in RESTRICTION_SITES.values()):
        matched_sites = [
            key for key, value in RESTRICTION_SITES.items() if value in new_seq
        ]
        new_seq = change_the_site(new_seq, matched_sites, RESTRICTION_SITES)
    return new_seq


def find_overhang_from_appropriate_range(
    gene: Frag,
    overhang_list: list[str | Seq],
    trial: int = -1,
) -> bool | tuple[str, str, Seq, Seq]:
    gene_max_len = gene.max_len_5 if gene.cut_order == 1 else gene.max_inner_len
    max_len = min(gene_max_len, gene.end_pos + 35)
    min_len = max(gene.frag_min_len or 0, gene.end_pos - 15)

    if max_len < gene.end_pos:
        return False

    pos_list = list(range(min_len, max_len + 1))
    random.shuffle(pos_list)
    last_result: tuple[str, str, Seq, Seq] | None = None
    for j in pos_list:
        frag, rest_seq, overhang, overhang_pair = split_with_pos(gene.dna_seq, j, gene.overhang_len)
        last_result = (frag, rest_seq, overhang, overhang_pair)
        overhang_str = str(overhang)
        if overhang_str in GOOD_OVERHANGS and check_overhang_unique(overhang_list, overhang):
            break
        if trial > 499 and overhang_str not in BAD_OVERHANGS and check_overhang_unique(overhang_list, overhang):
            break
    if last_result is None:
        return False
    return last_result


def find_overhang(
    gene: Frag,
    overhang_list: list[str | Seq],
    codons: dict[str, list[str]],
    overhang_len: int = 4,
    enable_replace_codons: bool = True,
) -> bool | tuple[str, str, Seq, Seq]:
    split_successful_tag = False
    if gene.max_len_5 < gene.max_inner_len and gene.cut_order == 1 and gene.end_pos >= gene.max_len_5:
        gene.end_pos = gene.frag_min_len or gene.end_pos
    if gene.frag_min_len is not None and gene.end_pos < gene.frag_min_len:
        gene.end_pos = gene.frag_min_len

    overhang: Seq | None = None
    frag = rest_seq = ""
    overhang_pair: Seq | None = None
    while not split_successful_tag:
        result = find_overhang_from_appropriate_range(gene, overhang_list, trial=-1)
        if isinstance(result, bool):
            break
        frag, rest_seq, overhang, overhang_pair = result
        overhang_str = str(overhang)

        if overhang_str in GOOD_OVERHANGS and check_overhang_unique(overhang_list, overhang):
            split_successful_tag = True
        elif not enable_replace_codons:
            break
        else:
            for trial in range(1000):
                overhang_str = str(overhang)
                if overhang_str in GOOD_OVERHANGS and check_overhang_unique(overhang_list, overhang):
                    split_successful_tag = True
                    break
                if trial > 499 and overhang_str not in BAD_OVERHANGS and check_overhang_unique(overhang_list, overhang):
                    split_successful_tag = True
                    break

                split_posi = len(gene.init_dna_seq) - len(gene.dna_seq) + gene.end_pos
                fiveprime = gene.init_dna_seq[:split_posi]
                frame_end = get_frame_end(len(fiveprime))
                gene_max_len = gene.max_len_5 if gene.cut_order == 1 else gene.max_inner_len
                max_len = min(gene_max_len, gene.end_pos + 35)
                min_len = max(gene.frag_min_len or 0, gene.end_pos - 15)
                change_len = max_len - min_len
                newseq = replace_codons(gene.init_dna_seq, frame_end, change_len, codons)
                old_seq = gene.init_dna_seq
                replace_start = frame_end + 1 - 6
                gene.init_dna_seq = (
                    gene.init_dna_seq[:replace_start]
                    + newseq
                    + gene.init_dna_seq[replace_start + len(newseq):]
                )
                if translate(gene.init_dna_seq) != translate(old_seq):
                    print("The amino acid sequence has changed!")
                gene.update_dna_seq()
                result = find_overhang_from_appropriate_range(gene, overhang_list, trial)
                if isinstance(result, bool):
                    break
                frag, rest_seq, overhang, overhang_pair = result
            if overhang is None or not check_overhang_unique(overhang_list, overhang):
                break

    if split_successful_tag and overhang is not None and overhang_pair is not None:
        return frag, rest_seq, overhang, overhang_pair
    return False


def split_sequences(
    designs: dict[str, str],
    protein_seqs: dict[str, str | None],
    codons: dict[str, list[str]],
    frag_num: int,
    subpool_barcodes: list[list[str]],
    subpool_index: int,
    adapter_f: str,
    adapter_r: str,
    sequence_barcodes: list[list[str]] | None,
    max_oligo_size: int = 300,
    min_oligo_len: int = 250,
    enable_replace_codons: bool = True,
) -> bool | tuple[SplitOutputs, list[list[str | Seq]]]:
    if len(designs) >= 96:
        raise ValueError("Each subpool should contain fewer than 96 designs")

    subp_barc_5_len = len(subpool_barcodes[0][subpool_index])
    subp_barc_3_len = len(subpool_barcodes[-1][subpool_index])
    subp_barc_in_len = len(subpool_barcodes[1][subpool_index])
    adapter_f_len = len(adapter_f)
    adapter_r_len = len(adapter_r)
    seq_barc_len = 0 if sequence_barcodes is None else len(sequence_barcodes[0][0]) + 6

    max_len_5 = max_oligo_size - (
        (adapter_f_len + (subp_barc_5_len + 6) + seq_barc_len)
        + ((subp_barc_in_len + 7) + adapter_r_len)
    )
    max_inner_len = max_oligo_size - (
        adapter_f_len + (subp_barc_in_len + 7) * 2 + adapter_r_len
    )
    max_len_3 = max_oligo_size - (
        (adapter_f_len + (subp_barc_in_len + 7))
        + (seq_barc_len + (6 + subp_barc_3_len) + adapter_r_len)
    )

    overhang_list: list[list[str | Seq]] = [[] for _ in range(1, frag_num)]
    split_result: dict[str, dict[str, str]] = {}
    naked_fragment_result: dict[str, list[str]] = {}
    whole_dna_result: dict[str, str] = {}

    seq_num = 0
    designs = sort_by_length(designs)
    for name, init_dna_seq in tqdm(designs.items()):
        gene = Frag(
            frag_num=frag_num,
            dna_seq=init_dna_seq,
            max_len_5=max_len_5,
            max_inner_len=max_inner_len,
            max_len_3=max_len_3,
            min_oligo_len=min_oligo_len,
            aa_seq=protein_seqs[name],
        )

        cut_successful_flag = False
        not_unique_overhang = False
        result: bool | tuple[str, str, Seq, Seq] = False
        try_times = 0
        while not cut_successful_flag:
            frag_seq: list[str] = []
            gene.dna_seq = gene.init_dna_seq
            for i in range(1, frag_num):
                gene.cut_order = i
                frag_num_left = frag_num - i
                gene.frag_min_len = gene.count_frag_lenth(cut_order=i, min_len=135)
                gene.end_pos = round((len(gene.dna_seq) + frag_num_left * 4) / (frag_num_left + 1))

                result = find_overhang(
                    gene,
                    overhang_list[i - 1],
                    codons,
                    overhang_len=4,
                    enable_replace_codons=enable_replace_codons,
                )
                if isinstance(result, bool):
                    try_times += 1
                    break

                frag, rest_seq, overhang, overhang_pair = result
                gene.dna_seq = rest_seq
                frag_seq.append(frag)
                if i == frag_num - 1:
                    frag_seq.append(rest_seq)
                overhang_list[i - 1].append(overhang)
                overhang_list[i - 1].append(overhang_pair)

            if isinstance(result, bool) and try_times < 5:
                continue
            if try_times == 5:
                not_unique_overhang = True
            else:
                cut_successful_flag = True

        if not_unique_overhang:
            print(f"Oligo {name}")
        else:
            naked_fragment_result[name] = frag_seq
            dna_seq = ""
            for i, fragment in enumerate(frag_seq):
                dna_seq += fragment if i == 0 else fragment[4:]
            whole_dna_result[name] = dna_seq

            seq_barcode: list[str] | str = "" if sequence_barcodes is None else sequence_barcodes[seq_num]
            frag_list = add_adapter_list(
                name,
                frag_seq,
                subpool_barcodes,
                subpool_index,
                frag_num,
                adapter_f,
                adapter_r,
                seq_barcode,
                enzyme_scheme="standard",
            )
            split_result[name] = frag_list
        seq_num += 1

    not_suitable_length_flag = any(
        len(value) < min_oligo_len or len(value) > max_oligo_size
        for frag_seq in split_result.values()
        for value in frag_seq.values()
    )
    if len(split_result) < len(designs) or not_suitable_length_flag:
        print("The split process failed because some genes are not splited. Try again!")
        return False
    return SplitOutputs(split_result, naked_fragment_result, whole_dna_result), overhang_list


def split_design_file(
    input_path: str | Path,
    subpool_barcode_path: str | Path,
    subpool_index: int,
    frag_num: int,
    adapter_f: str,
    adapter_r: str,
    sequence_barcode_path: str | Path | None,
    codon_table_path: str | Path | None,
    max_oligo_length: int,
    min_oligo_length: int,
    seed: int | None = None,
    enable_replace_codons: bool = True,
    return_intermediates: bool = False,
) -> dict[str, dict[str, str]] | SplitOutputs:
    if seed is not None:
        random.seed(seed)

    designs, protein_seqs = read_designs(input_path)
    codons = read_codon_table(codon_table_path or default_codon_table_path())
    subpool_barcodes = read_subpool_barcodes(subpool_barcode_path, frag_num)
    sequence_barcodes = (
        None
        if sequence_barcode_path is None
        else read_sequence_barcodes(sequence_barcode_path, len(designs))
    )

    result = split_sequences(
        designs,
        protein_seqs,
        codons,
        frag_num,
        subpool_barcodes,
        subpool_index - 1,
        adapter_f,
        adapter_r,
        sequence_barcodes,
        max_oligo_size=max_oligo_length,
        min_oligo_len=min_oligo_length,
        enable_replace_codons=enable_replace_codons,
    )
    split_attempts = 1
    max_split_attempts = 100
    while isinstance(result, bool):
        split_attempts += 1
        if split_attempts > max_split_attempts:
            raise ValueError(
                f"Unable to split all sequences after {max_split_attempts} attempts. "
                "Try changing fragment count, oligo length limits, barcode lengths, "
                "or re-enable codon replacement."
            )
        designs = shuffle_dict(designs)
        result = split_sequences(
            designs,
            protein_seqs,
            codons,
            frag_num,
            subpool_barcodes,
            subpool_index - 1,
            adapter_f,
            adapter_r,
            sequence_barcodes,
            max_oligo_size=max_oligo_length,
            min_oligo_len=min_oligo_length,
            enable_replace_codons=enable_replace_codons,
        )
    split_outputs, _overhang_list = result
    if return_intermediates:
        return split_outputs
    return split_outputs.final_oligos


def default_intermediate_output_paths(output_path: str | Path) -> tuple[Path, Path]:
    output_path = Path(output_path)
    output_prefix = output_path.with_suffix("") if output_path.suffix else output_path
    return (
        output_prefix.with_name(f"{output_prefix.name}_naked_fragments.tab"),
        output_prefix.with_name(f"{output_prefix.name}_whole_dna.tsv"),
    )


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--input", "--input-list", "--input_list", dest="input", required=True)
    parser.add_argument("--output", "-o", required=True)
    parser.add_argument(
        "--mute",
        action="store_true",
        help="Do not write naked fragments or reconstructed whole DNA outputs.",
    )
    parser.add_argument("--naked-output", "--naked_output", dest="naked_output")
    parser.add_argument("--whole-dna-output", "--whole_dna_output", dest="whole_dna_output")
    parser.add_argument("--subpool-barcodes", "--subpool_barcode_fname", dest="subpool_barcodes", required=True)
    parser.add_argument("--adapter-f", "--adapter_f", dest="adapter_f", default="F" * 20)
    parser.add_argument("--adapter-r", "--adapter_r", dest="adapter_r", default="R" * 20)
    parser.add_argument("--sequence-barcodes", "--sequence_barcode_fname", dest="sequence_barcodes")
    parser.add_argument("--subpool-index", "--subp_barc_index", dest="subpool_index", type=int, required=True)
    parser.add_argument("--max-oligo-length", "--max_oligo_length", dest="max_oligo_length", type=int, default=300)
    parser.add_argument("--min-oligo-length", "--min_oligo_length", dest="min_oligo_length", type=int, default=251)
    parser.add_argument("--frag-num", "--frag_num", dest="frag_num", type=int, default=3)
    parser.add_argument("--codon-table", "--codontable_fname", dest="codon_table")
    parser.add_argument("--seed", type=int, help="Random seed for reproducible examples/tests.")
    parser.add_argument(
        "--disable-replace-codons",
        "--no-codon-redesign",
        dest="enable_replace_codons",
        action="store_false",
        default=True,
        help="Disable synonymous codon redesign when searching for acceptable overhangs.",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Split genes into orthogonal pieces for HiPlex-GGA hierarchical multiplexed Golden Gate assembly."
    )
    add_arguments(parser)
    return parser


def run_from_args(args: argparse.Namespace) -> int:
    try:
        result = split_design_file(
            input_path=args.input,
            subpool_barcode_path=args.subpool_barcodes,
            subpool_index=args.subpool_index,
            frag_num=args.frag_num,
            adapter_f=args.adapter_f,
            adapter_r=args.adapter_r,
            sequence_barcode_path=args.sequence_barcodes,
            codon_table_path=args.codon_table,
            max_oligo_length=args.max_oligo_length,
            min_oligo_length=args.min_oligo_length,
            seed=args.seed,
            enable_replace_codons=args.enable_replace_codons,
            return_intermediates=True,
        )
    except ValueError as exc:
        print(exc)
        return 1
    write_oligo_table(result.final_oligos, args.output)

    naked_output = None
    whole_dna_output = None
    if not getattr(args, "mute", False):
        default_naked_output, default_whole_dna_output = default_intermediate_output_paths(
            args.output
        )
        naked_output = args.naked_output or default_naked_output
        whole_dna_output = args.whole_dna_output or default_whole_dna_output
    if naked_output:
        write_fragment_table(result.naked_fragments, naked_output)
    if whole_dna_output:
        write_name_sequence_table(result.whole_dna, whole_dna_output)
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return run_from_args(args)


if __name__ == "__main__":
    raise SystemExit(main())

