# Method overview

The public workflow has three steps.

1. Silently mutate restriction enzyme sites without changing the encoded amino
   acid sequence.
2. Split each coding DNA sequence into fragments with four-base Golden Gate
   overhangs.
3. Add adapters, subpool barcodes, optional sequence barcodes, and enzyme sites
   to produce orderable oligos.

The splitter searches for overhangs from a curated good-overhang list. If a
candidate overhang is not suitable or not unique within the same fragment
position, nearby codons can be synonymously redesigned and the search is retried.

