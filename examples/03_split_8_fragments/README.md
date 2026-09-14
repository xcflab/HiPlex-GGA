# Example 03: split one design into eight Golden Gate oligos

```bash
mkdir -p ./tmp

hiplex-gga split \
  --input examples/03_split_8_fragments/input.tsv \
  --subpool-barcodes examples/03_split_8_fragments/subpool_barcode_list_250410.tsv \
  --sequence-barcodes examples/03_split_8_fragments/sequence_barcode_13bp_96pairs_list_250506.tsv \
  --subpool-index 1 \
  --frag-num 8 \
  --adapter-f ccactccattcgtatcccacgtg \
  --adapter-r cggaatggctaggctgtacggat \
  --seed 1 \
  --output ./tmp/split8_s.tab
```

This example uses synthetic DNA and synthetic barcode tables. The generated
files are `./tmp/split8_s.tab`, `./tmp/split8_s_naked_fragments.tab`, and
`./tmp/split8_s_whole_dna.tsv`. Use `--mute` to skip both intermediate files.

`replace_codons` is enabled by default during overhang search. Add
`--disable-replace-codons` to require splitting without
synonymous codon redesign.
