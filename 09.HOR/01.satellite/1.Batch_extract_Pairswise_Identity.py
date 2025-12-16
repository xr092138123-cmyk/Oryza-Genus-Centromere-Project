#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse, re, sys, os, glob
from collections import defaultdict

def log(msg):
    print(msg, flush=True)

# ---------- 基础工具 ----------
HDR_COLS = ["query_name","query_start","query_end",
            "reference_name","reference_start","reference_end","perID_by_events"]

def parse_region(name):
    m = re.match(r"([^:]+):(\d+)-(\d+)$", name)
    if not m:
        raise ValueError(f"Bad region: {name}")
    return m.group(1), int(m.group(2)), int(m.group(3))

def bin_rel_start(rel_start, win):
    # 1-based 相对起点按窗口分箱：1, 1+win, 1+2*win, ...
    return ((rel_start - 1) // win) * win + 1

def make_id(chrom, abs_start, win):
    return f"{chrom}:{abs_start}-{abs_start + win - 1}"

def read_fasta_ids_and_win(fa_path):
    ids=set(); lens={}
    with open(fa_path) as f:
        h=None; L=0
        for line in f:
            if line.startswith(">"):
                if h is not None:
                    ids.add(h); lens[h]=L
                h=line[1:].strip().split()[0]
                L=0
            else:
                L += len(line.strip())
        if h is not None:
            ids.add(h); lens[h]=L
    win = next((v for v in lens.values() if v>0), None)
    if not win:
        raise ValueError(f"Cannot infer window size from {fa_path}")
    return ids, win

def parse_header_or_fixed(first_hash_line):
    idx = {k: None for k in HDR_COLS}
    if first_hash_line is not None:
        header = first_hash_line.lstrip("#").strip()
        cols = re.split(r"\s+|\t", header)
        for i,c in enumerate(cols):
            if c in idx and idx[c] is None:
                idx[c]=i
    if any(v is None for v in idx.values()):
        idx = {k:i for i,k in enumerate(HDR_COLS)}  # 回退固定顺序
    return idx

# ---------- 核心：处理单个 FASTA ----------
def process_one(fasta_path, moddot_root, out_dir, min_id=None):
    fasta_name = os.path.basename(fasta_path)
    m = re.match(r"^(?P<prefix>.+\.Chr\d+)\.w(?P<w>\d+)\.Top1\.fa$", fasta_name)
    if not m:
        raise ValueError(f"FASTA 名称不符合规则：{fasta_name}")
    prefix = m.group("prefix")           # 例如 AA_Ogla_hap1.Chr01
    w_str  = m.group("w")                # 例如 100
    w_int  = int(w_str)

    # 期望的 moddot 目录
    mod_dir = os.path.join(moddot_root, f"{prefix}.centromere.w{w_str}")
    if not os.path.isdir(mod_dir):
        raise FileNotFoundError(f"找不到 moddot 目录：{mod_dir}")

    # 目录内查找唯一 .bed
    bed_list = [p for p in glob.glob(os.path.join(mod_dir, "*.bed")) if os.path.isfile(p)]
    if len(bed_list) != 1:
        raise RuntimeError(f"在 {mod_dir} 中必须且仅有 1 个 .bed，当前找到 {len(bed_list)} 个：{bed_list}")
    bed_path = bed_list[0]

    # 读取 FASTA 的窗口 ID 集合 & 窗口大小（用于双保险）
    fa_ids, win_from_fa = read_fasta_ids_and_win(fasta_path)
    if win_from_fa != w_int:
        log(f"[WARN] {fasta_name}: FASTA 推断窗口={win_from_fa}, 文件名窗口={w_int}，将使用 {win_from_fa}")

    win = win_from_fa
    log(f"[RUN ] {fasta_name}  ->  {os.path.basename(bed_path)}   (win={win})")

    # 输出文件
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{prefix}.w{win}.pairwise.tsv")

    # 首行可能是带 # 的表头；我们先探测一下
    with open(bed_path) as f:
        first = f.readline()
    idx = parse_header_or_fixed(first if first.startswith("#") else None)

    # 主动再打开并逐行处理（大文件友好）
    pair2idty = defaultdict(float)
    total = 0; kept = 0
    with open(bed_path) as f:
        for ln in f:
            if not ln.strip(): continue
            if ln.startswith("#"):  # 跳过表头或注释
                continue
            parts = re.split(r"\s+|\t", ln.strip())
            # 坏行保护
            if len(parts) < max(idx.values())+1:
                continue

            try:
                qn = parts[idx["query_name"]]
                qs = int(parts[idx["query_start"]])
                rn = parts[idx["reference_name"]]
                rs = int(parts[idx["reference_start"]])
                idty = float(parts[idx["perID_by_events"]])
            except Exception:
                continue

            total += 1
            if (min_id is not None) and (idty < min_id):
                continue

            q_chr, q_s, q_e = parse_region(qn)
            r_chr, r_s, r_e = parse_region(rn)

            qs_bin = bin_rel_start(qs, win)
            rs_bin = bin_rel_start(rs, win)

            q_abs = q_s + (qs_bin - 1)
            r_abs = r_s + (rs_bin - 1)

            q_id = make_id(q_chr, q_abs, win)
            r_id = make_id(r_chr, r_abs, win)

            # 仅收 FASTA 集合内部的 pair
            if (q_id in fa_ids) and (r_id in fa_ids):
                a, b = (q_id, r_id) if q_id <= r_id else (r_id, q_id)
                if idty > pair2idty[(a,b)]:
                    pair2idty[(a,b)] = idty
                kept += 1

    # 写出
    with open(out_path, "w") as w:
        w.write("win_id_a\twin_id_b\tidentity\n")
        for (a,b),v in sorted(pair2idty.items()):
            w.write(f"{a}\t{b}\t{v}\n")

    log(f"[OK  ] {fasta_name}: 读入 {total:,} 行，命中 {kept:,} 行 → pairs={len(pair2idty):,} → {out_path}")
    return out_path

# ---------- 批量 ----------
def main():
    ap = argparse.ArgumentParser(description="批量把 moddotplot .bed 映射为 FASTA ID 对（仅保留 FASTA 内部两两）")
    ap.add_argument("--moddot-root", required=True, help="02.moddotplot_analysis 根目录")
    ap.add_argument("--fasta-dir",   required=True, help="06.repeat_top1_fa 目录")
    ap.add_argument("--out-dir",     required=True, help="输出目录（自动创建）")
    ap.add_argument("--min-id",      type=float, default=None, help="相似度下限（与输入单位一致；bed 若是百分比，填 86 即可）")
    ap.add_argument("--pattern",     default="*.Top1.fa", help="匹配 FASTA 文件名的通配符（默认 *.Top1.fa）")
    args = ap.parse_args()

    fas = sorted(glob.glob(os.path.join(args.fasta_dir, args.pattern)))
    if not fas:
        log(f"[ERR ] 在 {args.fasta_dir} 下找不到 {args.pattern}")
        sys.exit(1)

    ok, fail = 0, 0
    for fa in fas:
        try:
            process_one(fa, args.moddot_root, args.out_dir, args.min_id)
            ok += 1
        except Exception as e:
            log(f"[FAIL] {os.path.basename(fa)}: {e}")
            fail += 1
    log(f"[DONE] 成功 {ok} 个，失败 {fail} 个")

if __name__ == "__main__":
    main()
