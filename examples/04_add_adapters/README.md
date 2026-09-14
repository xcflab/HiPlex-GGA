# Example 04: add adapters to naked fragments

```bash
hiplex-gga add-adapters \
  --input examples/04_add_adapters/naked_fragments.tab \
  --subpool-barcodes examples/04_add_adapters/subpool_barcode_list_250410.tsv \
  --sequence-barcodes examples/04_add_adapters/sequence_barcode_13bp_96pairs_list_250506.tsv \
  --subpool-index 1 \
  --frag-num 4 \
  --adapter-f ccactccattcgtatcccacgtg \
  --adapter-r cggaatggctaggctgtacggat \
  --enzyme-scheme standard \
  --output examples/04_add_adapters/oligos_with_adapters.tab
```

The input fragment list may be JSON-style as shown here. Historical Python-list
format is also accepted through `ast.literal_eval` for compatibility.

