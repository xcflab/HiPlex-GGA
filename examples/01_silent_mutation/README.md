# Example 01: silent mutation

This example removes Type IIS enzyme sites while preserving the translated amino
acid sequence.

```bash
oligos-gg mutate-sites \
  --input examples/01_silent_mutation/input.tsv \
  --enzyme-type All \
  --output examples/01_silent_mutation/changed.tsv
```

Input format:

```text
name<TAB>DNA_sequence
```

The command writes one `name<TAB>mutated_DNA_sequence` record per input design.

