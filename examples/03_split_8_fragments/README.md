# Example 03: split one design into eight Golden Gate oligos

```bash
oligos-gg split \
  --input examples/03_split_8_fragments/input.tsv \
  --subpool-barcodes examples/03_split_8_fragments/subpool_barcodes.tsv \
  --sequence-barcodes examples/03_split_8_fragments/sequence_barcodes.tsv \
  --subpool-index 1 \
  --frag-num 8 \
  --adapter-f ccactccattcgtatcccacgtg \
  --adapter-r cggaatggctaggctgtacggat \
  --seed 1 \
  --output examples/03_split_8_fragments/oligos.tab
```

This example uses synthetic DNA and synthetic barcode tables.

