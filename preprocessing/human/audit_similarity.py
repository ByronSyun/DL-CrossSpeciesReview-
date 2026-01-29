import argparse
import json
import os
import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
import subprocess
from bisect import bisect_left, bisect_right

import pandas as pd


VARIANT_ID_PATTERNS = [
    # chr1:1234:C>A
    re.compile(r"^(chr)?(?P<chrom>[^:]+):(?P<pos>\d+):(?P<ref>[ACGT])>(?P<alt>[ACGT])$", re.IGNORECASE),
    # 1:1234_C/A
    re.compile(r"^(chr)?(?P<chrom>[^:]+):(?P<pos>\d+)[_/](?P<ref>[ACGT])[/>(](?P<alt>[ACGT])\)?$", re.IGNORECASE),
    # chr1-1234-C-A (SpliceVarDB hg38 column)
    re.compile(r"^(chr)?(?P<chrom>[^\-]+)-(?P<pos>\d+)-(?P<ref>[ACGT])-(?P<alt>[ACGT])$", re.IGNORECASE),
    # chr1_1234_A_T (Vex-seq processed tables)
    re.compile(r"^(chr)?(?P<chrom>[^_]+)_(?P<pos>\d+)_(?P<ref>[ACGT]+)_(?P<alt>[ACGT]+)$", re.IGNORECASE),
]


def parse_variant_id(s: str) -> Optional[Tuple[str, int, str, str]]:
    s = str(s).strip()
    for pat in VARIANT_ID_PATTERNS:
        m = pat.match(s)
        if m:
            chrom = m.group("chrom")
            if not chrom.startswith("chr"):
                chrom = f"chr{chrom}"
            return chrom, int(m.group("pos")), m.group("ref").upper(), m.group("alt").upper()
    return None


def normalize_row(row: Dict) -> Optional[Tuple[str, int, str, str]]:
    chrom = row.get("chrom")
    pos = row.get("pos")
    ref = row.get("ref")
    alt = row.get("alt")
    if all(k is not None for k in [chrom, pos, ref, alt]):
        chrom = str(chrom)
        if not chrom.startswith("chr"):
            chrom = f"chr{chrom}"
        return chrom, int(pos), str(ref).upper(), str(alt).upper()
    # common variant_id-like columns
    v = (row.get("variant_id") or row.get("variant") or row.get("id") or
         row.get("hg38") or row.get("hg19") or row.get("hg37"))
    if v is not None:
        return parse_variant_id(v)
    return None


def load_splicevardb(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    cols = {c.lower(): c for c in df.columns}
    out = []
    for _, r in df.iterrows():
        row = {k: r[v] for k, v in cols.items()}
        parsed = normalize_row(row)
        if parsed:
            out.append(parsed)
    return pd.DataFrame(out, columns=["chrom", "pos", "ref", "alt"]).drop_duplicates()


def load_dataset_generic(file: str, fmt: str, mapping: Dict[str, str]) -> pd.DataFrame:
    if fmt.lower() == "vcf":
        chroms = []
        poss = []
        refs = []
        alts = []
        with open(file, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 5:
                    continue
                chrom, pos, _id, ref, alt = parts[:5]
                # Skip multi-allelic beyond first alt for strict match
                alt_first = alt.split(",")[0]
                if not chrom.startswith("chr"):
                    chrom = f"chr{chrom}"
                chroms.append(chrom)
                poss.append(int(pos))
                refs.append(ref.upper())
                alts.append(alt_first.upper())
        return pd.DataFrame({"chrom": chroms, "pos": poss, "ref": refs, "alt": alts}).drop_duplicates()
    
    sep = "\t" if fmt.lower() == "tsv" else ","
    df = pd.read_csv(file, sep=sep, dtype=str)
    cols = {c.lower(): c for c in df.columns}
    out = []
    for _, r in df.iterrows():
        row = {k: r[v] for k, v in cols.items()}
        
        for k in ("chrom", "pos", "ref", "alt", "variant_id"):
            if mapping.get(k):
                row[k] = r[mapping[k]] if mapping[k] in df.columns else None
        
        # Auto-map common MFASS headers if not provided
        if row.get("chrom") is None and "chr" in cols:
            row["chrom"] = r[cols["chr"]]
        if row.get("pos") is None:
            for kpos in ("snp_position_hg38_1based", "snp_position_hg37_1based", "position", "pos"):
                if kpos in cols:
                    row["pos"] = r[cols[kpos]]
                    break
        if row.get("ref") is None:
            for kref in ("ref_allele", "refallele", "ref"):
                if kref in cols:
                    row["ref"] = r[cols[kref]]
                    break
        if row.get("alt") is None:
            for kalt in ("alt_allele", "altallele", "alt"):
                if kalt in cols:
                    row["alt"] = r[cols[kalt]]
                    break
        parsed = normalize_row(row)
        if parsed:
            out.append(parsed)
    return pd.DataFrame(out, columns=["chrom", "pos", "ref", "alt"]).drop_duplicates()


def exact_and_positional_overlap(a: pd.DataFrame, b: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    a_key = a.assign(key=a["chrom"] + ":" + a["pos"].astype(str) + ":" + a["ref"] + ":" + a["alt"])
    b_key = b.assign(key=b["chrom"] + ":" + b["pos"].astype(str) + ":" + b["ref"] + ":" + b["alt"])
    exact = a_key.merge(b_key[["key"]], on="key", how="inner").drop_duplicates()
    
    a_pos = a.assign(
        pkey=a["chrom"] + ":" + a["pos"].astype(str),
        key=a["chrom"] + ":" + a["pos"].astype(str) + ":" + a["ref"] + ":" + a["alt"]
    )
    b_pos = b.assign(pkey=b["chrom"] + ":" + b["pos"].astype(str))
    pos_overlap = a_pos.merge(b_pos[["pkey"]], on="pkey", how="inner")
    
    pos_overlap = pos_overlap.merge(exact[["key"]], on="key", how="left", indicator=True)
    pos_overlap = pos_overlap[pos_overlap["_merge"] == "left_only"].drop(columns=["_merge"]).drop_duplicates()
    return exact, pos_overlap


def try_import_minhash():
    try:
        from datasketch import MinHash
        return MinHash
    except Exception:
        return None


def build_kmers(seq: str, k: int = 21) -> List[str]:
    seq = seq.upper()
    return [seq[i:i+k] for i in range(0, max(0, len(seq)-k+1)) if "N" not in seq[i:i+k]]


def minhash_jaccard(kmers_a: List[str], kmers_b: List[str], num_perm: int = 256) -> float:
    MinHash = try_import_minhash()
    if MinHash is None:
        # fallback: exact Jaccard on sets (slower but ok for small batches)
        set_a, set_b = set(kmers_a), set(kmers_b)
        if not set_a or not set_b:
            return 0.0
        return len(set_a & set_b) / len(set_a | set_b)
    mh_a = MinHash(num_perm=num_perm)
    for kmer in kmers_a:
        mh_a.update(kmer.encode("utf8"))
    mh_b = MinHash(num_perm=num_perm)
    for kmer in kmers_b:
        mh_b.update(kmer.encode("utf8"))
    return mh_a.jaccard(mh_b)


def extract_window_sequences(df: pd.DataFrame, fasta_path: str, window: int) -> Dict[str, str]:
    try:
        from pyfaidx import Fasta
    except Exception as e:
        raise SystemExit(f"pyfaidx is required for sequence extraction: {e}")
    fa = Fasta(fasta_path, as_raw=True, sequence_always_upper=True)
    seqs = {}
    half = window
    for _, r in df.iterrows():
        chrom, pos = r["chrom"], int(r["pos"])  # 1-based
        start = max(1, pos - half)
        end = pos + half
        try:
            seq = str(fa[chrom][start-1:end])
            key = f"{chrom}:{pos}:{r['ref']}:{r['alt']}"
            seqs[key] = seq
        except KeyError:
            continue
    return seqs


def main():
    p = argparse.ArgumentParser(description="Audit overlap/similarity between SpliceVarDB and risk datasets")
    p.add_argument("--spdb", required=True, help="SpliceVarDB TSV path (filtered benchmark)")
    p.add_argument("--dataset", action="append", default=[],
                   help="Dataset spec: name=<n> file=<path> fmt=<vcf|tsv|csv> [chrom= col] [pos= col] [ref= col] [alt= col]")
    p.add_argument("--fasta", help="Reference FASTA for sequence extraction (optional)")
    p.add_argument("--window", type=int, default=1000, help="Half window for sequence extraction (default 1000bp)")
    p.add_argument("--jaccard", type=float, default=0.9, help="MinHash/Set Jaccard threshold for similarity candidates")
    p.add_argument("--blastn", action="store_true", help="Run BLASTN on extracted windows for high-confidence similarity")
    p.add_argument("--blastn_cmd", default="blastn", help="blastn binary (default: blastn)")
    p.add_argument("--makeblastdb_cmd", default="makeblastdb", help="makeblastdb binary (default: makeblastdb)")
    p.add_argument("--perc_identity", type=float, default=95.0, help="BLASTN percent identity threshold")
    p.add_argument("--qcov", type=float, default=80.0, help="BLASTN query coverage percent threshold")
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    sp = load_splicevardb(args.spdb)

    datasets = []
    for spec in args.dataset:
        parts = dict(kv.split("=", 1) for kv in spec.split() if "=" in kv)
        if not {"name", "file", "fmt"} <= parts.keys():
            raise SystemExit(f"Bad --dataset spec: {spec}")
        mapping = {k: parts.get(k) for k in ("chrom", "pos", "ref", "alt", "variant_id")}
        ds = load_dataset_generic(parts["file"], parts["fmt"], mapping)
        datasets.append((parts["name"], ds))

    summary = {}
    sp_seqs = None
    if args.fasta:
        sp_seqs = extract_window_sequences(sp, args.fasta, args.window)

    for name, ds in datasets:
        exact, pos_only = exact_and_positional_overlap(sp, ds)
        exact.to_csv(os.path.join(args.outdir, f"exact_overlap_{name}.csv"), index=False)
        pos_only.to_csv(os.path.join(args.outdir, f"positional_overlap_{name}.csv"), index=False)

        entry = {
            "splicevardb_n": int(len(sp)),
            "dataset_n": int(len(ds)),
            "exact_n": int(len(exact)),
            "positional_n": int(len(pos_only)),
        }

        sim_n = 0
        blast_n = 0
        if args.fasta:
            ds_seqs = extract_window_sequences(ds, args.fasta, args.window)
            # Index ds by chrom and sorted positions to prune comparisons
            ds_index: Dict[str, Dict[str, List]] = {}
            for _, r in ds.iterrows():
                chrom, pos = r["chrom"], int(r["pos"])
                key = f"{chrom}:{pos}:{r['ref']}:{r['alt']}"
                if key not in ds_seqs:
                    continue
                if chrom not in ds_index:
                    ds_index[chrom] = {"pos": [], "key": [], "kmers": []}
                ds_index[chrom]["pos"].append(pos)
                ds_index[chrom]["key"].append(key)
                ds_index[chrom]["kmers"].append(build_kmers(ds_seqs[key], k=21))
            
            for chrom in ds_index:
                order = sorted(range(len(ds_index[chrom]["pos"])), key=lambda i: ds_index[chrom]["pos"][i])
                ds_index[chrom] = {
                    "pos": [ds_index[chrom]["pos"][i] for i in order],
                    "key": [ds_index[chrom]["key"][i] for i in order],
                    "kmers": [ds_index[chrom]["kmers"][i] for i in order],
                }

            rows = []
            pos_window = args.window * 2
            for _, r in sp.iterrows():
                chrom, pos = r["chrom"], int(r["pos"])
                sp_key = f"{chrom}:{pos}:{r['ref']}:{r['alt']}"
                if sp_key not in sp_seqs:
                    continue
                km_sp = build_kmers(sp_seqs[sp_key], k=21)
                if not km_sp:
                    continue
                if chrom not in ds_index:
                    continue
                pos_list = ds_index[chrom]["pos"]
                lo = bisect_left(pos_list, pos - pos_window)
                hi = bisect_right(pos_list, pos + pos_window)
                for i in range(lo, hi):
                    jac = minhash_jaccard(km_sp, ds_index[chrom]["kmers"][i])
                    if jac >= args.jaccard:
                        rows.append({
                            "sp": sp_key,
                            "ds": ds_index[chrom]["key"][i],
                            "jaccard": jac
                        })
            if rows:
                out = pd.DataFrame(rows).sort_values("jaccard", ascending=False)
                out.to_csv(os.path.join(args.outdir, f"similarity_candidates_{name}.csv"), index=False)
                sim_n = int(len(out))

            if args.blastn and rows:
                workdir = os.path.join(args.outdir, f"blast_{name}")
                os.makedirs(workdir, exist_ok=True)
                sp_fa = os.path.join(workdir, "sp_windows.fa")
                ds_fa = os.path.join(workdir, "ds_windows.fa")
                
                cand_sp_keys = sorted(set([r["sp"] for r in rows]))
                cand_ds_keys = sorted(set([r["ds"] for r in rows]))
                with open(sp_fa, "w") as f:
                    for key_sp in cand_sp_keys:
                        seq = sp_seqs.get(key_sp)
                        if seq:
                            f.write(f">{key_sp}\n{seq}\n")
                with open(ds_fa, "w") as f:
                    for key_ds in cand_ds_keys:
                        seq = ds_seqs.get(key_ds)
                        if seq:
                            f.write(f">{key_ds}\n{seq}\n")
                
                subprocess.run([args.makeblastdb_cmd, "-in", ds_fa, "-dbtype", "nucl"], check=True)
                blast_out = os.path.join(workdir, "blast_out.tsv")
                cmd = [
                    args.blastn_cmd,
                    "-query", sp_fa,
                    "-db", ds_fa,
                    "-outfmt", "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore",
                    "-perc_identity", str(args.perc_identity),
                    "-max_hsps", "1",
                    "-max_target_seqs", "1",
                ]
                subprocess.run(cmd, check=True, stdout=open(blast_out, "w"))
                
                cols = ["qseqid","sseqid","pident","length","qlen","slen","qstart","qend","sstart","send","evalue","bitscore"]
                bo = pd.read_csv(blast_out, sep="\t", names=cols)
                if not bo.empty:
                    bo["qcov"] = 100.0 * bo["length"] / bo["qlen"].replace(0, pd.NA)
                    bo = bo[(bo["pident"] >= args.perc_identity) & (bo["qcov"] >= args.qcov)]
                    bo.to_csv(os.path.join(args.outdir, f"blast_hits_{name}.csv"), index=False)
                    blast_n = int(len(bo))
        entry["similarity_candidates_n"] = sim_n
        entry["blast_hits_n"] = blast_n
        summary[name] = entry

    with open(os.path.join(args.outdir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()


