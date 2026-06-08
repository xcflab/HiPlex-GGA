# Example 03: split one design into eight Golden Gate oligos

```bash
mkdir -p ./tmp

oligos-gg split \
  --input examples/03_split_8_fragments/input.tsv \
  --subpool-barcodes examples/03_split_8_fragments/subpool_barcodes.tsv \
  --sequence-barcodes examples/03_split_8_fragments/sequence_barcodes.tsv \
  --subpool-index 1 \
  --frag-num 8 \
  --adapter-f ccactccattcgtatcccacgtg \
  --adapter-r cggaatggctaggctgtacggat \
  --seed 1 \
  --output ./tmp/split8_s.tab \
  --write-intermediate-outputs
```

This example uses synthetic DNA and synthetic barcode tables. The generated
files are `./tmp/split8_s.tab`, `./tmp/split8_s_naked_fragments.tab`, and
`./tmp/split8_s_whole_dna.tsv`.

`replace_codons` is enabled by default during overhang search. Add
`--disable-replace-codons` or `--no-codon-redesign` to require splitting without
synonymous codon redesign.
