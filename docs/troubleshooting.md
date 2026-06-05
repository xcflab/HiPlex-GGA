# Troubleshooting

## Split process failed

The sequence may be too short or too long for the selected fragment number,
adapter length, barcode length, and oligo length limits. Try a different
fragment number or adjust `--max-oligo-length` / `--min-oligo-length`.

## Oligo shorter than the minimum

The workflow adds short filler sequences for very short fragments, but some
parameter combinations can still fail. Check adapter and barcode lengths first.

## Enzyme site remains after mutation

Confirm that the sequence is a coding sequence in frame 0 and contains only
standard DNA bases.

## Results differ between runs

Overhang selection and synonymous codon redesign use random choices. Use
`--seed` for reproducible examples and tests.

