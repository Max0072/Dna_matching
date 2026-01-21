from collections import defaultdict


# Levenshtein distance
def lev_dist(a, b, k):
    INF = k + 1
    # --- trim common prefix ---
    n, m, pref_c = len(a), len(b), 0
    lim = min(n, m)
    while pref_c < lim and a[pref_c] == b[pref_c]:
        pref_c += 1
    if pref_c: a, b = a[pref_c:], b[pref_c:]
    # --- trim common suffix ---
    n, m, suf_c = len(a), len(b), 0
    lim = min(n, m)
    while suf_c < lim and a[n - 1 - suf_c] == b[m - 1 - suf_c]:
        suf_c += 1
    if suf_c: a, b = a[:n - suf_c], b[:m - suf_c]
    # -----------------
    if b == "":
        return len(a) if len(a) <= k else k + 1
    # -----------------
    a = a.encode()
    b = b.encode()
    n, m = len(a), len(b)
    if m > n:
        a, b = b, a
        n, m = m, n
    # -----------------
    prev = [INF] * (2 * k + 3)
    cur  = [INF] * (2 * k + 3)
    prev_lo = 0
    prev_hi = min(m, k)
    prev[0] = INF
    for j in range(prev_hi + 1):
        prev[j + 1] = j
    prev[prev_hi + 2] = INF
    # -----------------
    for i in range(1, n + 1):
        if min(prev[1:prev_hi-prev_lo+3]) > k:
            return k + 1
        ai = a[i - 1]
        cur_lo = max(0, i - k)
        cur_hi = min(m, i + k)
        cur_len = cur_hi - cur_lo + 3
        cur[0] = INF
        cur[cur_len - 1] = INF
        base_prev = 1 - prev_lo
        base_cur  = 1 - cur_lo
        for j in range(cur_lo, cur_hi + 1):
            pt = j + base_prev
            t  = j + base_cur
            delete  = prev[pt] + 1
            insert  = cur[t - 1] + 1
            if j == 0:
                cur[t] = delete
            else:
                replace = prev[pt - 1] + (ai != b[j - 1])
                cur[t] = min(delete, insert, replace)
        prev, cur = cur, prev
        prev_lo, prev_hi = cur_lo, cur_hi
    return prev[m - prev_lo + 1]

# Read dna from the file
def read_fna(filename, max_size=0):
    dna = ""
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            if line[0] == ">":
                continue
            if max_size and len(dna) >= max_size:
                break
            dna += line.strip()
    print(f"File loaded, dna size is {len(dna)}\n")
    return dna

# Get reads
def read_fastq(filename):
    with open(filename) as f:
        while f.readline():
            seq = f.readline().strip()
            f.readline()
            f.readline()
            yield seq

# Build index for kmers
def build_kmer_index(dna, k):
    index = defaultdict(list)
    n = len(dna)
    for i in range(n - k + 1):
        index[dna[i:i + k]] += [i]
    return index

# Print statistics
def print_info(intervals, total_reads_processed, total_reads_mapped, total_unique_matchings, total_alignment_errors, dna_size):
    print(f"Total number of reads processed: {total_reads_processed} ({total_reads_processed / 22720100 * 100:.1f}%)")
    print(f"Mapping rates: {total_reads_mapped / total_reads_processed * 100 if total_reads_processed > 0 else 0:.1f}%")
    print(f"Percentage of unique reads: {total_unique_matchings / total_reads_mapped * 100 if total_reads_mapped > 0 else 0:.1f}%")
    print(f"Percentage of multi-mapped reads: {(total_reads_mapped - total_unique_matchings) / total_reads_mapped * 100 if total_reads_mapped > 0 else 0:.1f}%")
    print(f"Alignment quality: {total_alignment_errors/total_unique_matchings if total_unique_matchings > 0 else 0} errors per read")

# Find reads in dna
def process_reads(S, kmer_index, k, filename):
    intervals = []
    total_reads_processed = 0
    total_reads_mapped = 0
    total_unique_matchings = 0
    total_alignment_errors = 0
    dna_size = len(S)
    for i, r in enumerate(read_fastq(filename)):
        total_reads_processed += 1

        # print info
        if i % 100000 == 99999:
            print_info(intervals, total_reads_processed, total_reads_mapped, total_unique_matchings, total_alignment_errors, dna_size)
            coverage, intervals = count_total_matchings_and_update(intervals)
            print(f"Genome read coverage: {coverage / dna_size * 100:.4f}% ({coverage} units)")
            print("-" * 30)

        # see candidates by small seed
        r_size = len(r)
        seeds = ((0, r[:k]), #prefix
                 (r_size//2, r[r_size//2-k//2:r_size//2-k//2+k]), #middle
                 (r_size-k, r[r_size-k:r_size])) #suffix
        best = None
        for off, seed in seeds:
            cand = kmer_index.get(seed, [])
            if not cand:
                continue
            if best is None or len(cand) < len(best[2]):
                best = (off, seed, cand)
        if not best:
            continue

        # check full match
        off, seed, candidates = best
        place = None
        m_count = 0
        for c in candidates:
            max_dist = 5
            if c + r_size <= len(S) and ((distance:=lev_dist(S[c-off:c+r_size-off], r, max_dist)) <= max_dist):
                m_count += 1
                place = c
                if m_count == 1:
                    saved_dist = distance
                if m_count > 1:
                    break
        if m_count == 1:
            intervals.append((place - off, place + r_size - off))

        # collect stats
        if m_count == 1:
            total_unique_matchings += 1
            total_alignment_errors += saved_dist
        if m_count >= 1: total_reads_mapped += 1
    print_info(intervals, total_reads_processed, total_reads_mapped, total_unique_matchings, total_alignment_errors, dna_size)
    coverage, intervals = count_total_matchings_and_update(intervals)
    print(f"Genome read coverage: {coverage / dna_size * 100:.4f}% ({coverage} units)")
    print("-" * 30)
    return intervals, total_reads_processed, total_reads_mapped, total_unique_matchings

# count total amount of matchings in dna
def count_total_matchings_and_update(intervals):
    if not intervals:
        return 0
    new_intervals = []
    intervals.sort()
    count = 0
    cur_l, cur_r = intervals[0]
    for l, r in intervals[1:]:
        if l > cur_r:
            new_intervals.append((cur_l, cur_r))
            count += cur_r - cur_l
            cur_l, cur_r = l, r
        elif r > cur_r:
            cur_r = r
    new_intervals.append((cur_l, cur_r))
    count += cur_r - cur_l
    return count, new_intervals


def main():
    S = read_fna("./GCF_000005845.2_ASM584v2_genomic.fna")
    k = 15
    kmer_index = build_kmer_index(S, k)
    process_reads(S, kmer_index, k, "./ERR022075_1.fastq")


if __name__ == "__main__":
    main()













