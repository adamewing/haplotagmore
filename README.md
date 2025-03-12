## haplotagmore

```
usage: haplotagmore.py [-h] -v VCF -b BAM -r REF [--minvotes MINVOTES] [--minmargin MINMARGIN] [--debug] [--all_snps] [--loose_indels] [--ignore_indels]
                       [--phased_only]

assign alignments to parent of origin based on snps and indels

optional arguments:
  -h, --help            show this help message and exit
  -v VCF, --vcf VCF     vcf
  -b BAM, --bam BAM     bam
  -r REF, --ref REF     indexed fasta
  --minvotes MINVOTES   default = 1
  --minmargin MINMARGIN
                        default = 1
  --debug
  --all_snps            ignore vcf filters (default: PASS only)
  --loose_indels        allow left end of deletion or right end of insertion matching (useful for nanopore indels)
  --ignore_indels       do not consider indels
  --phased_only         only output phased reads
  ```
