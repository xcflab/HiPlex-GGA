from pathlib import Path

from oligos_gg_assembly.cli import build_parser
from oligos_gg_assembly.split_oligos import default_intermediate_output_paths


def test_cli_has_split_subcommand():
    parser = build_parser()
    args = parser.parse_args([
        "split",
        "--input", "input.tsv",
        "--output", "out.tab",
        "--subpool-barcodes", "barcodes.tsv",
        "--subpool-index", "1",
    ])
    assert args.command == "split"
    assert args.frag_num == 3
    assert args.write_intermediate_outputs is False


def test_cli_can_disable_replace_codons():
    parser = build_parser()
    args = parser.parse_args([
        "split",
        "--input", "input.tsv",
        "--output", "out.tab",
        "--subpool-barcodes", "barcodes.tsv",
        "--subpool-index", "1",
        "--disable-replace-codons",
    ])
    assert args.enable_replace_codons is False


def test_cli_accepts_no_codon_redesign_alias():
    parser = build_parser()
    args = parser.parse_args([
        "split",
        "--input", "input.tsv",
        "--output", "out.tab",
        "--subpool-barcodes", "barcodes.tsv",
        "--subpool-index", "1",
        "--no-codon-redesign",
    ])
    assert args.enable_replace_codons is False


def test_cli_accepts_intermediate_split_outputs():
    parser = build_parser()
    args = parser.parse_args([
        "split",
        "--input", "input.tsv",
        "--output", "out.tab",
        "--write-intermediate-outputs",
        "--naked-output", "naked.tab",
        "--whole-dna-output", "whole.tsv",
        "--subpool-barcodes", "barcodes.tsv",
        "--subpool-index", "1",
    ])
    assert args.write_intermediate_outputs is True
    assert args.naked_output == "naked.tab"
    assert args.whole_dna_output == "whole.tsv"


def test_default_intermediate_output_paths_use_output_prefix():
    naked_output, whole_dna_output = default_intermediate_output_paths(
        "results/final_oligos.tab"
    )

    assert naked_output == Path("results/final_oligos_naked_fragments.tab")
    assert whole_dna_output == Path("results/final_oligos_whole_dna.tsv")
