# Example 02: split one design into four Golden Gate oligos

```bash
oligos-gg split \
  --input examples/02_split_4_fragments/input.tsv \
  --subpool-barcodes examples/02_split_4_fragments/subpool_barcodes.tsv \
  --sequence-barcodes examples/02_split_4_fragments/sequence_barcodes.tsv \
  --subpool-index 1 \
  --frag-num 4 \
  --adapter-f ccactccattcgtatcccacgtg \
  --adapter-r cggaatggctaggctgtacggat \
  --seed 1 \
  --output examples/02_split_4_fragments/oligos.tab \
  --write-intermediate-outputs
```

The main output contains one `fragment_name, oligo_sequence` line per designed
oligo. The intermediate outputs are
`examples/02_split_4_fragments/oligos_naked_fragments.tab` and
`examples/02_split_4_fragments/oligos_whole_dna.tsv`.

`replace_codons` is enabled by default during overhang search. Add
`--disable-replace-codons` or `--no-codon-redesign` to require splitting without
synonymous codon redesign.
