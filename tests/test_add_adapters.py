from oligos_gg_assembly.add_adapters import add_adapters_to_fragments


def test_add_adapters_outputs_one_oligo_per_fragment():
    fragments = {"toy": ["ATG" + "GCT" * 30, "GCT" * 31, "GCT" * 31, "GCT" * 30 + "TAA"]}
    subpool_barcodes = [
        ["AACGTT"], ["TTCGAA"], ["CCGTAA"], ["TTACGG"],
        ["GGATTC"], ["GAATCC"], ["CTAGGA"], ["TCCTAG"],
    ]
    sequence_barcodes = [["ATGCGA", "TCGCAT"]]
    result = add_adapters_to_fragments(
        fragments,
        subpool_barcodes,
        subpool_index=0,
        frag_num=4,
        adapter_f="ccactccattcgtatcccacgtg",
        adapter_r="cggaatggctaggctgtacggat",
        sequence_barcodes=sequence_barcodes,
    )
    assert set(result["toy"]) == {"a_toy", "b_toy", "c_toy", "d_toy"}
    assert result["toy"]["a_toy"].startswith("ccactccattcgtatcccacgtgAACGTT")

