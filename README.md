# DNA Read Mapper

A short-read DNA mapper that aligns FASTQ reads to a reference genome using seed-based search and Levenshtein distance alignment. Written in pure Python with no external dependencies.

## How it works

Mapping each read against a 4.6 Mbp genome by brute-force edit distance would be too slow. This tool uses a two-stage approach:

**1. k-mer index (seed lookup)**

Before processing any reads, the reference genome is indexed: every 15-nucleotide substring (k-mer) is mapped to its position(s) in the genome. This gives O(1) candidate lookup.

**2. Seed-based candidate selection**

For each read, three seeds are extracted — from the prefix, middle, and suffix. The seed with the fewest hits in the index is chosen to minimize false candidates. This handles cases where one region of the read falls in a repetitive part of the genome.

**3. Band-limited Levenshtein alignment**

Each candidate position is verified by computing the Levenshtein (edit) distance between the read and the corresponding genomic window. The DP is band-limited to a maximum of 5 errors, which allows early exit and keeps the inner loop to O(k·n) instead of O(n²). Common prefix/suffix is trimmed before DP to further reduce work.

**4. Unique mapping filter**

A read is recorded only if it maps to exactly one position (unique mapping). Multi-mapped reads are counted in stats but excluded from coverage.

**5. Coverage calculation**

Mapped intervals are merged (overlapping intervals collapsed) to compute the total genome coverage without double-counting.

## Results

Tested on *E. coli* K-12 (4.6 Mbp genome, `GCF_000005845.2`) with Illumina reads `ERR022075_1.fastq`:

| Metric | Value |
|---|---|
| Reads processed | 22,720,100 (100%) |
| Mapping rate | 48.2% |
| Uniquely mapped | 97.0% |
| Multi-mapped | 3.0% |
| Avg. errors per read | 0.197 |
| Genome coverage | 98.48% |

## Usage

**Requirements:** Python 3.8+, no external packages.

**Input files** (place in the same directory):
- Reference genome in FASTA format: `GCF_000005845.2_ASM584v2_genomic.fna`
- Reads in FASTQ format: `ERR022075_1.fastq`

Both files can be downloaded from NCBI:
- Genome: [GCF_000005845.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005845.2/)
- Reads: [ERR022075](https://www.ebi.ac.uk/ena/browser/view/ERR022075) (ENA)

**Run:**
```bash
python match_reads_and_dna.py
```

Progress is printed every 100k reads.

## Parameters

| Parameter | Value | Location |
|---|---|---|
| k-mer seed size | 15 | `main()` |
| Max edit distance | 5 | `process_reads()` |
