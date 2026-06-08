# File formats

## Design input

Two supported whitespace-delimited formats:

```text
name DNA_sequence
name amino_acid_sequence DNA_sequence
```

DNA sequences should be coding sequences in frame 0 when silent mutation or
codon redesign is used.

## Subpool barcode table

The first column is a pool name. The remaining columns are paired fragment
barcodes. A four-fragment design needs eight barcode columns:

```text
pool A_F A_R B_F B_R C_F C_R D_F D_R
pool_1 AACGTT TTCGAA CCGTAA TTACGG GGATTC GAATCC CTAGGA TCCTAG
```

An eight-fragment design needs sixteen barcode columns.

## Sequence barcode table

```text
name seq5 seq3
toy_1 ATGCGA TCGCAT
```

CSV-style comma separators are also accepted.

## Oligo output

```text
a_design_name, oligo_sequence
b_design_name, oligo_sequence
```

Fragment prefixes follow the fragment order: `a_`, `b_`, `c_`, etc.

## Naked fragment output

When `--write-intermediate-outputs` or `--naked-output` is supplied, the
splitter writes one design per line. With `--output oligos.tab`, the default
path is `oligos_naked_fragments.tab`:

```text
name, ["fragment_1", "fragment_2", "fragment_3", "fragment_4"]
```

The fragment list is JSON-compatible and can be used as input for
`oligos-gg add-adapters`. Historical Python-list format is also accepted by
`add-adapters`.

## Whole DNA output

When `--write-intermediate-outputs` or `--whole-dna-output` is supplied, the
splitter writes reconstructed whole DNA after removing repeated four-base
overhangs from downstream fragments. With `--output oligos.tab`, the default
path is `oligos_whole_dna.tsv`:

```text
name<TAB>whole_DNA_sequence
```
