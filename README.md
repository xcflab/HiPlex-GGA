# Oligos GG Assembly

Tools for designing oligo pools for Golden Gate assembly workflows.

The package supports:

- silent mutation of restriction enzyme sites while preserving protein sequence;
- splitting coding DNA sequences into Golden Gate fragments;
- adding pool adapters, subpool barcodes, sequence barcodes, and enzyme sites.

This public release focuses only on oligo / Golden Gate assembly design. It does
not include real experimental datasets, NGS/PacBio analysis, DNWorks batch
pipelines, or overlap PCR / two-oligo assembly workflows.

## Installation

Clone the repository and install it into a Python environment:

```bash
python -m pip install .
```

Minimum requirements:

- Python 3.10 or newer
- Biopython
- NumPy
- tqdm

## Command Line Usage

After installation, the package provides the `oligos-gg` command.

### Silent mutation

```bash
oligos-gg mutate-sites \
  --input examples/01_silent_mutation/input.tsv \
  --enzyme-type All \
  --output changed.tsv
```

### Split into Golden Gate oligos

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
  --output oligos.tab \
  --write-intermediate-outputs
```

With `--output oligos.tab`, intermediate outputs are written as
`oligos_naked_fragments.tab` and `oligos_whole_dna.tsv`. To override those
paths explicitly:

```bash
oligos-gg split ... \
  --output final_oligos.tab \
  --write-intermediate-outputs \
  --naked-output custom_naked_fragments.tab \
  --whole-dna-output custom_whole_dna.tsv
```

By default, `split` may use synonymous codon redesign through `replace_codons`
to search for acceptable and unique overhangs. To disable this behavior:

```bash
oligos-gg split ... --disable-replace-codons
```

The alias `--no-codon-redesign` is also accepted.

### Add adapters to naked fragments

```bash
oligos-gg add-adapters \
  --input examples/04_add_adapters/naked_fragments.tab \
  --subpool-barcodes examples/04_add_adapters/subpool_barcodes.tsv \
  --sequence-barcodes examples/04_add_adapters/sequence_barcodes.tsv \
  --subpool-index 1 \
  --frag-num 4 \
  --adapter-f ccactccattcgtatcccacgtg \
  --adapter-r cggaatggctaggctgtacggat \
  --enzyme-scheme standard \
  --output oligos_with_adapters.tab
```

## Input Formats

Design input is whitespace-delimited:

```text
name DNA_sequence
```

or:

```text
name amino_acid_sequence DNA_sequence
```

See [docs/file_formats.md](docs/file_formats.md) for barcode and output formats.

## Examples

- [Silent mutation](examples/01_silent_mutation/README.md)
- [Four-fragment split](examples/02_split_4_fragments/README.md)
- [Eight-fragment split](examples/03_split_8_fragments/README.md)
- [Add adapters](examples/04_add_adapters/README.md)

All examples use synthetic data.

## Notes

The public `split` implementation is based on the historical
`oligo_split_251111.py` script and is exposed as
`oligos_gg_assembly.split_oligos`.
