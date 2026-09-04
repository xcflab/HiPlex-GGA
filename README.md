# Oligos GG Assembly

Tools for designing oligo pools for Golden Gate assembly workflows.

The package supports:

- silent mutation of restriction enzyme sites while preserving protein sequence;
- splitting coding DNA sequences into Golden Gate fragments;
- adding pool adapters, subpool barcodes, sequence barcodes, and enzyme sites.

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

The `mutate-sites` command removes restriction enzyme recognition sites from
coding DNA by making synonymous codon changes, preserving the translated amino
acid sequence.

```bash
oligos-gg mutate-sites \
  --input examples/01_silent_mutation/input.tsv \
  --enzyme-type All \
  --output changed.tsv
```

Supported `--enzyme-type` values:

| Value | Restriction sites included |
| --- | --- |
| `BsaI` | BsaI forward and reverse-complement sites: `GGTCTC`, `GAGACC` |
| `BsmBI` | BsmBI forward and reverse-complement sites: `CGTCTC`, `GAGACG` |
| `TypeII` | BsaI and BsmBI sites: `GGTCTC`, `GAGACC`, `CGTCTC`, `GAGACG` |
| `TypeI` | Common cloning sites: NdeI `CATATG`, XhoI `CTCGAG`, HindIII `AAGCTT`, EcoRI `GAATTC`, NcoI `CCATGG`, BamHI `GGATCC`, NheI `GCTAGC`, XbaI `TCTAGA` |
| `All` | All supported sites in `TypeI` and `TypeII`, plus BbsI forward and reverse-complement sites: `GAAGAC`, `GTCTTC` |

Input sequences should be coding DNA in frame 0 so synonymous codon changes can
be checked by translation.

### Split into Golden Gate oligos

```bash
oligos-gg split \
  --input examples/02_split_4_fragments/input.tsv \
  --subpool-barcodes examples/02_split_4_fragments/subpool_barcode_list_250410.tsv \
  --sequence-barcodes examples/02_split_4_fragments/sequence_barcode_13bp_96pairs_list_250506.tsv \
  --subpool-index 1 \
  --frag-num 4 \
  --adapter-f ccactccattcgtatcccacgtg \
  --adapter-r cggaatggctaggctgtacggat \
  --seed 1 \
  --output oligos.tab \
  --write-intermediate-outputs
```

Split parameters:

| Parameter | Required | Description |
| --- | --- | --- |
| `--input` | Yes | Input design table. Each row contains a design name and coding DNA sequence, optionally with an amino acid sequence column. |
| `--subpool-barcodes` | Yes | Barcode table that provides the subpool barcode pairs used on each fragment end. The table must contain enough barcode columns for `--frag-num`. |
| `--subpool-index` | Yes | 1-based index of the subpool barcode set to use from the barcode table. Use `1` for the first subpool row. |
| `--frag-num` | No | Number of Golden Gate fragments to split each design into. Default: `3`. Examples use `4` or `8`. |
| `--adapter-f` | No | Forward outer adapter sequence added to every final oligo. |
| `--adapter-r` | No | Reverse outer adapter sequence added to every final oligo. |
| `--sequence-barcodes` | No | Optional per-design barcode table. When supplied, each design receives its own 5-prime and 3-prime sequence barcode. |
| `--output` | Yes | Main output file for final oligos with adapters, subpool barcodes, sequence barcodes, and enzyme sites. |
| `--write-intermediate-outputs` | No | Also write naked fragments and reconstructed whole DNA using the `--output` prefix. |
| `--naked-output` | No | Custom path for naked fragments before adapters are added. Overrides the default path from `--write-intermediate-outputs`. |
| `--whole-dna-output` | No | Custom path for reconstructed whole DNA. Overrides the default path from `--write-intermediate-outputs`. |
| `--max-oligo-length` | No | Maximum allowed final oligo length after adapters/barcodes/enzyme sites are added. Default: `300`. |
| `--min-oligo-length` | No | Minimum allowed final oligo length after adapters/barcodes/enzyme sites are added. Default: `251`. |
| `--codon-table` | No | Custom codon table used for synonymous codon redesign. If omitted, the bundled codon table is used. |
| `--seed` | No | Random seed for reproducible overhang search and example outputs. |
| `--disable-replace-codons`, `--no-codon-redesign` | No | Disable synonymous codon redesign during overhang search. |

Several historical option aliases are also accepted, including names with
underscores such as `--frag_num`, `--adapter_f`, and `--subpool_barcode_fname`.
The hyphenated option names shown above are recommended for new use.

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

### Add adapters to naked fragments

The `add-adapters` command converts pre-split naked DNA fragments into final
Golden Gate oligos. It adds the forward and reverse outer adapters, subpool
barcode pairs, optional per-design sequence barcodes, and the restriction enzyme
sites needed for Golden Gate assembly. This is useful when you already have a
naked fragment file, such as the `*_naked_fragments.tab` output from
`oligos-gg split --write-intermediate-outputs`, and want to regenerate final
oligos with a different barcode/adapter setup without re-running the splitting
step.

```bash
oligos-gg add-adapters \
  --input examples/04_add_adapters/naked_fragments.tab \
  --subpool-barcodes examples/04_add_adapters/subpool_barcode_list_250410.tsv \
  --sequence-barcodes examples/04_add_adapters/sequence_barcode_13bp_96pairs_list_250506.tsv \
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
