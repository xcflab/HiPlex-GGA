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

