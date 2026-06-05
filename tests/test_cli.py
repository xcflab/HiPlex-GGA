from oligos_gg_assembly.cli import build_parser


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

