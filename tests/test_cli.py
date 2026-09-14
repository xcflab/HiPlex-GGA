from pathlib import Path
from types import SimpleNamespace

from oligos_gg_assembly import split_oligos

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
    assert args.mute is False


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


def test_cli_accepts_mute():
    parser = build_parser()
    args = parser.parse_args([
        "split",
        "--input", "input.tsv",
        "--output", "out.tab",
        "--mute",
        "--subpool-barcodes", "barcodes.tsv",
        "--subpool-index", "1",
    ])
    assert args.mute is True


def test_cli_rejects_removed_intermediate_outputs_flag():
    parser = build_parser()
    try:
        parser.parse_args([
            "split",
            "--input", "input.tsv",
            "--output", "out.tab",
            "--write-intermediate-outputs",
            "--subpool-barcodes", "barcodes.tsv",
            "--subpool-index", "1",
        ])
    except SystemExit as exc:
        assert exc.code == 2
    else:
        raise AssertionError("removed --write-intermediate-outputs was accepted")


def test_default_intermediate_output_paths_use_output_prefix():
    naked_output, whole_dna_output = default_intermediate_output_paths(
        "results/final_oligos.tab"
    )

    assert naked_output == Path("results/final_oligos_naked_fragments.tab")
    assert whole_dna_output == Path("results/final_oligos_whole_dna.tsv")


def test_split_writes_intermediates_by_default(monkeypatch, tmp_path):
    parser = build_parser()
    args = parser.parse_args([
        "split",
        "--input", "input.tsv",
        "--output", str(tmp_path / "final.tab"),
        "--subpool-barcodes", "barcodes.tsv",
        "--subpool-index", "1",
    ])
    result = SimpleNamespace(final_oligos={}, naked_fragments={}, whole_dna={})
    written = []
    monkeypatch.setattr(split_oligos, "split_design_file", lambda **kwargs: result)
    monkeypatch.setattr(split_oligos, "write_oligo_table", lambda data, path: written.append(("final", Path(path))))
    monkeypatch.setattr(split_oligos, "write_fragment_table", lambda data, path: written.append(("naked", Path(path))))
    monkeypatch.setattr(split_oligos, "write_name_sequence_table", lambda data, path: written.append(("whole", Path(path))))

    assert split_oligos.run_from_args(args) == 0
    assert written == [
        ("final", tmp_path / "final.tab"),
        ("naked", tmp_path / "final_naked_fragments.tab"),
        ("whole", tmp_path / "final_whole_dna.tsv"),
    ]


def test_split_mute_skips_intermediate_outputs(monkeypatch, tmp_path):
    parser = build_parser()
    args = parser.parse_args([
        "split",
        "--input", "input.tsv",
        "--output", str(tmp_path / "final.tab"),
        "--subpool-barcodes", "barcodes.tsv",
        "--subpool-index", "1",
        "--mute",
    ])
    result = SimpleNamespace(final_oligos={}, naked_fragments={}, whole_dna={})
    written = []
    monkeypatch.setattr(split_oligos, "split_design_file", lambda **kwargs: result)
    monkeypatch.setattr(split_oligos, "write_oligo_table", lambda data, path: written.append(("final", Path(path))))
    monkeypatch.setattr(split_oligos, "write_fragment_table", lambda data, path: written.append(("naked", Path(path))))
    monkeypatch.setattr(split_oligos, "write_name_sequence_table", lambda data, path: written.append(("whole", Path(path))))

    assert split_oligos.run_from_args(args) == 0
    assert written == [("final", tmp_path / "final.tab")]
