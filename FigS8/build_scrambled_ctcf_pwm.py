#!/usr/bin/env python3
"""
Build a PWM from CTCF anchor coordinates and produce a position-shuffled control PWM.

Pipeline
--------
1. Read BED6 of single-bp CTCF motif anchors (chrom, start, end, name, score, strand).
2. Extend each anchor +/- FLANK bp around its center  ->  (2*FLANK + 1) bp window
   (default FLANK=9  ->  19 bp).
3. Extract the genomic sequence of every window from a FASTA (faidx-indexable),
   reverse-complementing windows on the '-' strand so all sequences share the
   motif orientation.
4. Tally a 4 x W base-count matrix and convert it to a position-probability matrix (PPM).
5. Shuffle: randomly PERMUTE THE W COLUMNS (positions). Each position's four-base
   frequency vector is swapped, as a block, to another position. This keeps the
   per-position composition but destroys the positional pattern  ->  a null control.
6. Write both the original and shuffled PWMs in MEME minimal-motif format (ready for FIMO),
   plus the permutation used and the raw counts for reproducibility.

Usage
-----
    python build_scrambled_ctcf_pwm.py \
        --bed   Hs_CTCF.motif \
        --genome /path/to/genome.fa \
        --flank 9 \
        --seed  0 \
        --out-prefix out/CTCF

Dependencies: numpy, and one of {pysam, pyfaidx}.
    pip install numpy pysam      # (or)  pip install numpy pyfaidx
The genome FASTA must have a .fai index (samtools faidx genome.fa). pysam/pyfaidx
build it automatically on first use if the file is writable.
"""

import argparse
import os
import sys

import numpy as np

BASES = "ACGT"
B2I = {b: i for i, b in enumerate(BASES)}
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


# ----------------------------------------------------------------------------- #
# Pure-Python indexed FASTA reader (no third-party deps). Works on a plain       #
# (uncompressed) FASTA; reads an existing .fai or builds one. Used as a fallback #
# so the script runs even where pysam/pyfaidx are unavailable.                   #
# ----------------------------------------------------------------------------- #
class FaidxLite:
    def __init__(self, path):
        self.path = path
        self.index = {}      # name -> (length, offset, linebases, linewidth)
        fai = path + ".fai"
        if os.path.exists(fai) and os.path.getmtime(fai) >= os.path.getmtime(path):
            self._load_fai(fai)
        else:
            self._build_fai(fai)
        self._fh = open(path, "rb")

    def _load_fai(self, fai):
        with open(fai) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                self.index[p[0]] = (int(p[1]), int(p[2]), int(p[3]), int(p[4]))

    def _build_fai(self, fai):
        rows = []
        with open(self.path, "rb") as f:
            name = None
            length = offset = linebases = linewidth = 0
            seen_short = False
            while True:
                pos = f.tell()
                raw = f.readline()
                if not raw:
                    break
                if raw[:1] == b">":
                    if name is not None:
                        rows.append((name, length, offset, linebases, linewidth))
                    name = raw[1:].split()[0].decode() if len(raw) > 1 else ""
                    length = 0
                    offset = f.tell()
                    linebases = linewidth = 0
                    seen_short = False
                    self.index[name] = None
                else:
                    nb = len(raw.rstrip(b"\r\n"))
                    if linebases == 0:
                        linebases, linewidth = nb, len(raw)
                    elif nb and not seen_short and nb != linebases:
                        # only the final (short) line of a record may differ
                        seen_short = True
                    length += nb
            if name is not None:
                rows.append((name, length, offset, linebases, linewidth))
        with open(fai, "w") as out:
            for name, length, offset, lb, lw in rows:
                self.index[name] = (length, offset, lb, lw)
                out.write(f"{name}\t{length}\t{offset}\t{lb}\t{lw}\n")

    def references(self):
        return self.index.keys()

    def length(self, chrom):
        return self.index[chrom][0]

    def fetch(self, chrom, start, end):
        length, offset, lb, lw = self.index[chrom]
        start = max(0, start)
        end = min(end, length)
        if end <= start:
            return ""
        b0 = offset + (start // lb) * lw + (start % lb)
        b1 = offset + (end // lb) * lw + (end % lb)
        self._fh.seek(b0)
        raw = self._fh.read(b1 - b0)
        return raw.replace(b"\n", b"").replace(b"\r", b"").decode().upper()


# ----------------------------------------------------------------------------- #
# Genome access: prefer pysam (handles bgzipped .fa.gz + .fai), then pyfaidx,    #
# then the built-in FaidxLite. Uniform .fetch(chrom, start, end) 0-based API.    #
# ----------------------------------------------------------------------------- #
class Genome:
    def __init__(self, path):
        self.backend = None
        try:
            import pysam
            self._fa = pysam.FastaFile(path)
            self._chroms = set(self._fa.references)
            self.backend = "pysam"
            return
        except Exception:
            pass
        try:
            from pyfaidx import Fasta
            self._fa = Fasta(path, sequence_always_upper=True,
                             rebuild=False, as_raw=True)
            self._chroms = set(self._fa.keys())
            self.backend = "pyfaidx"
            return
        except Exception:
            pass
        try:
            self._fa = FaidxLite(path)
            self._chroms = set(self._fa.references())
            self.backend = "builtin-faidx"
        except Exception as e:
            sys.exit(
                "ERROR: could not open genome with pysam, pyfaidx, or builtin reader.\n"
                "  (builtin reader needs an UNCOMPRESSED .fa; for .fa.gz use pysam)\n"
                f"  underlying error: {e}"
            )

    def has(self, chrom):
        return chrom in self._chroms

    def length(self, chrom):
        if self.backend == "pysam":
            return self._fa.get_reference_length(chrom)
        if self.backend == "builtin-faidx":
            return self._fa.length(chrom)
        return len(self._fa[chrom])

    def fetch(self, chrom, start, end):
        """0-based, half-open. Returns uppercase string."""
        if self.backend == "pysam":
            return self._fa.fetch(chrom, start, end).upper()
        if self.backend == "builtin-faidx":
            return self._fa.fetch(chrom, start, end)
        return str(self._fa[chrom][start:end]).upper()


# ----------------------------------------------------------------------------- #
def build_count_matrix(bed_path, genome, flank, site_cap=None):
    """Return (counts[4,W], n_used, n_skipped). W = 2*flank + 1."""
    width = 2 * flank + 1
    counts = np.zeros((4, width), dtype=np.int64)
    n_used = n_skip = 0
    warned_chroms = set()

    with open(bed_path) as fh:
        for ln, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            if len(f) < 3:
                continue
            chrom = f[0]
            try:
                start = int(f[1])  # 0-based start of the 1-bp anchor
            except ValueError:
                continue
            strand = f[5] if len(f) >= 6 and f[5] in ("+", "-") else "+"

            # center = the anchor base; window = [center-flank, center+flank+1)
            center = start
            ws, we = center - flank, center + flank + 1

            if not genome.has(chrom):
                if chrom not in warned_chroms:
                    warned_chroms.add(chrom)
                    sys.stderr.write(
                        f"[warn] chrom '{chrom}' not in genome FASTA "
                        f"(check chr-naming); skipping such records.\n")
                n_skip += 1
                continue
            if ws < 0 or we > genome.length(chrom):
                n_skip += 1
                continue

            seq = genome.fetch(chrom, ws, we)
            if len(seq) != width:
                n_skip += 1
                continue
            if strand == "-":
                seq = revcomp(seq)

            ok = True
            col_idx = np.empty(width, dtype=np.int64)
            for j, b in enumerate(seq):
                bi = B2I.get(b, -1)
                if bi < 0:          # N or other ambiguous base -> drop window
                    ok = False
                    break
                col_idx[j] = bi
            if not ok:
                n_skip += 1
                continue

            counts[col_idx, np.arange(width)] += 1
            n_used += 1
            if site_cap and n_used >= site_cap:
                break

    return counts, n_used, n_skip


def counts_to_ppm(counts, pseudocount):
    """Column-normalize counts (+pseudocount) -> probability matrix [4,W]."""
    c = counts.astype(np.float64) + pseudocount
    return c / c.sum(axis=0, keepdims=True)


def write_meme(path, ppm, motif_name, nsites, bg):
    """Write a single motif in MEME minimal-motif format (FIMO-ready)."""
    W = ppm.shape[1]
    with open(path, "w") as o:
        o.write("MEME version 4\n\n")
        o.write("ALPHABET= ACGT\n\n")
        o.write("strands: + -\n\n")
        o.write("Background letter frequencies\n")
        o.write("A %.5f C %.5f G %.5f T %.5f\n\n" % tuple(bg))
        o.write(f"MOTIF {motif_name}\n")
        o.write("letter-probability matrix: alength= 4 w= %d nsites= %d E= 0\n"
                % (W, nsites))
        for j in range(W):
            o.write("  %.6f  %.6f  %.6f  %.6f\n" % tuple(ppm[:, j]))
        o.write("\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bed", required=True, help="BED6 of 1-bp CTCF anchors")
    ap.add_argument("--genome", required=True, help="genome FASTA (faidx-indexable)")
    ap.add_argument("--flank", type=int, default=9, help="bp each side (default 9 -> 19 bp)")
    ap.add_argument("--seed", type=int, default=0, help="RNG seed for the column shuffle")
    ap.add_argument("--pseudocount", type=float, default=0.25,
                    help="added to each base count before normalizing (default 0.25)")
    ap.add_argument("--site-cap", type=int, default=None,
                    help="optional: only use the first N valid windows")
    ap.add_argument("--out-prefix", required=True, help="output path prefix")
    args = ap.parse_args()

    outdir = os.path.dirname(os.path.abspath(args.out_prefix))
    os.makedirs(outdir, exist_ok=True)

    sys.stderr.write("[1/4] opening genome ...\n")
    genome = Genome(args.genome)
    sys.stderr.write(f"      backend = {genome.backend}\n")

    sys.stderr.write("[2/4] extracting windows & tallying counts ...\n")
    counts, n_used, n_skip = build_count_matrix(
        args.bed, genome, args.flank, args.site_cap)
    width = 2 * args.flank + 1
    sys.stderr.write(f"      width = {width} bp | windows used = {n_used} | "
                     f"skipped (edge/N/chrom) = {n_skip}\n")
    if n_used == 0:
        sys.exit("ERROR: no usable windows. Check chrom naming and FASTA path.")

    # background = empirical mononucleotide freq across all kept windows
    bg = counts.sum(axis=1).astype(np.float64)
    bg = bg / bg.sum()

    sys.stderr.write("[3/4] building PPM and shuffling positions ...\n")
    ppm = counts_to_ppm(counts, args.pseudocount)

    rng = np.random.default_rng(args.seed)
    perm = rng.permutation(width)            # new order of the W columns
    ppm_shuf = ppm[:, perm]
    counts_shuf = counts[:, perm]

    sys.stderr.write("[4/4] writing outputs ...\n")
    pfx = args.out_prefix
    write_meme(pfx + "_original.meme", ppm, "CTCF_original", n_used, bg)
    write_meme(pfx + "_shuffled.meme", ppm_shuf, "CTCF_posShuffled", n_used, bg)

    # permutation record (new_position <- old_position)
    with open(pfx + "_perm.txt", "w") as o:
        o.write("# column shuffle: new_index\told_index (0-based)\n")
        o.write(f"# seed={args.seed}\n")
        for new_i, old_i in enumerate(perm):
            o.write(f"{new_i}\t{old_i}\n")

    # raw counts (rows = A,C,G,T; cols = positions), original orientation
    np.savetxt(pfx + "_counts.tsv", counts, fmt="%d", delimiter="\t",
               header="rows=A,C,G,T  cols=pos0..pos%d (motif orientation)" % (width - 1))
    np.savetxt(pfx + "_counts_shuffled.tsv", counts_shuf, fmt="%d", delimiter="\t",
               header="rows=A,C,G,T  cols=shuffled positions")

    sys.stderr.write(
        "done.\n"
        f"  original PWM : {pfx}_original.meme\n"
        f"  shuffled PWM : {pfx}_shuffled.meme   <-- scrambled control motif for FIMO or logo plotting\n"
        f"  permutation  : {pfx}_perm.txt\n"
        f"  bg (A,C,G,T) : {bg[0]:.4f} {bg[1]:.4f} {bg[2]:.4f} {bg[3]:.4f}\n")


if __name__ == "__main__":
    main()
