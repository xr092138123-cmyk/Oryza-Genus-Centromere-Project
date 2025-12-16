#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
批量：按 blocks 文件的每一行 mono_start–mono_end 区间，
在对应 score 文件中截取 pos 落入区间的数据并作图：
- 合并图（所有区间叠加；PNG+PDF）
- 单独区间图（默认 PNG；可选导出一个多页 PDF）
- 统计汇总 TSV（每个区间的点数、范围与统计量）
- 缺失区间（若某区间没有命中点，记录在 blocks_no_points.tsv）

匹配原则：
- 在 blocks_dir 中递归寻找所有以 blocks_suffix 结尾的文件（默认 ".hor_pairs_blocks.tsv"）
- 公共前缀 = 文件名去掉 blocks_suffix（例如：AA_Ogla_hap1.Chr01.w100.Top1）
- 对应 score 文件 = 在 score_dir 中递归查找 (公共前缀 + score_suffix)（默认 ".hor_score.tsv"）

用法（示例）：
python plot_scores_by_blocks.py \
  --blocks_dir /share/org/YOUR/.../blocks/ \
  --score_dir  /share/org/YOUR/.../scores/ \
  --out_root   /share/org/YOUR/.../hor_plots_batch/ \
  --individual-pdf

依赖：pandas, matplotlib（可选：Pillow 用于个别转换，但本脚本不强制）
"""

import argparse
import os
from pathlib import Path
import sys
import glob
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")  # 后端适配服务器无显示环境
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def log(msg: str):
    print(msg, flush=True)


# ---------- 读取与清洗 ----------
def read_blocks(p: Path) -> pd.DataFrame:
    df = pd.read_csv(p, sep="\t", dtype=str)
    df.columns = [c.strip() for c in df.columns]

    # 候选列名（大小写无关）
    cand_start = [c for c in df.columns if c.lower() in {"mono_start", "start", "apos_start"}]
    cand_end   = [c for c in df.columns if c.lower() in {"mono_end", "end", "bpos_end"}]
    if not cand_start or not cand_end:
        raise ValueError(f"[{p}] 缺少 mono_start/mono_end 相关列；当前列：{df.columns.tolist()}")

    s_col, e_col = cand_start[0], cand_end[0]
    df = df[[s_col, e_col]].rename(columns={s_col: "mono_start", e_col: "mono_end"})

    # 转为整数 + 规范化区间
    df["mono_start"] = pd.to_numeric(df["mono_start"], errors="coerce")
    df["mono_end"]   = pd.to_numeric(df["mono_end"], errors="coerce")
    df = df.dropna(subset=["mono_start", "mono_end"]).astype({"mono_start": int, "mono_end": int})
    df["mono_start"], df["mono_end"] = np.minimum(df["mono_start"], df["mono_end"]), np.maximum(df["mono_start"], df["mono_end"])

    # 添加 block_id
    df = df.reset_index(drop=True)
    df.insert(0, "block_id", np.arange(1, len(df) + 1))
    return df


def read_score(p: Path) -> tuple[pd.DataFrame, str]:
    df = pd.read_csv(p, sep="\t", dtype=str)
    df.columns = [c.strip() for c in df.columns]

    # 定位 pos 列
    pos_col = None
    for c in df.columns:
        if c.lower() in {"pos", "position", "bp", "loc", "site", "idx"}:
            pos_col = c
            break
    if pos_col is None:
        raise ValueError(f"[{p}] 未找到 pos 列；当前列：{df.columns.tolist()}")

    # 自动选择 score 列：优先包含 'score' 的数值列；否则第一个数值列
    score_col = None
    for c in df.columns:
        if c == pos_col:
            continue
        if "score" in c.lower():
            tmp = pd.to_numeric(df[c], errors="coerce")
            if tmp.notna().any():
                score_col = c
                break

    if score_col is None:
        for c in df.columns:
            if c == pos_col:
                continue
            tmp = pd.to_numeric(df[c], errors="coerce")
            if tmp.notna().any():
                score_col = c
                break

    if score_col is None:
        raise ValueError(f"[{p}] 未找到可用数值列作为 score；列：{df.columns.tolist()}")

    out = pd.DataFrame({
        "pos": pd.to_numeric(df[pos_col], errors="coerce"),
        "score": pd.to_numeric(df[score_col], errors="coerce")
    }).dropna().sort_values("pos").reset_index(drop=True)

    # 同位点平均
    out = out.groupby("pos", as_index=False, sort=True).agg({"score": "mean"})
    return out, score_col


# ---------- 核心处理 ----------
def process_one(prefix: str,
                blocks_file: Path,
                score_file: Path,
                out_dir: Path,
                max_legend_items: int = 30,
                export_individual_pdf: bool = False) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    log(f"\n=== 处理：{prefix} ===")
    log(f"Blocks: {blocks_file}")
    log(f"Score : {score_file}")
    log(f"Out   : {out_dir}")

    blocks = read_blocks(blocks_file)
    score_df, used_score_col = read_score(score_file)

    # 为每个区间切片
    slices = []
    missing = []
    for _, r in blocks.iterrows():
        b_id = int(r["block_id"])
        a, b = int(r["mono_start"]), int(r["mono_end"])
        sub = score_df[(score_df["pos"] >= a) & (score_df["pos"] <= b)].copy().sort_values("pos")
        sub["block_id"] = b_id
        sub["mono_start"] = a
        sub["mono_end"] = b
        if sub.empty:
            missing.append({"block_id": b_id, "mono_start": a, "mono_end": b})
        slices.append(sub)

    sliced_df = pd.concat(slices, ignore_index=True) if slices else pd.DataFrame(columns=["pos", "score", "block_id", "mono_start", "mono_end"])

    # 统计表
    meta = (
        sliced_df.groupby(["block_id", "mono_start", "mono_end"], as_index=False)
        .agg(n_points=("pos", "size"),
             pos_min=("pos", "min"),
             pos_max=("pos", "max"),
             score_min=("score", "min"),
             score_max=("score", "max"),
             score_mean=("score", "mean"))
        .sort_values("block_id")
    )
    meta_path = out_dir / "blocks_slice_summary.tsv"
    meta.to_csv(meta_path, sep="\t", index=False)

    if missing:
        pd.DataFrame(missing).to_csv(out_dir / "blocks_no_points.tsv", sep="\t", index=False)

    # 合并图：PNG + PDF
    if not sliced_df.empty:
        fig = plt.figure(figsize=(12, 6))
        ax = plt.gca()
        for b_id, dfb in sliced_df.groupby("block_id"):
            label = f"block_{b_id} [{int(dfb['mono_start'].iloc[0])}-{int(dfb['mono_end'].iloc[0])}]"
            ax.plot(dfb["pos"].values, dfb["score"].values, label=label, linewidth=1.0)
        ax.set_xlabel("pos")
        ax.set_ylabel(f"score (from column: {used_score_col})")
        ax.set_title("Score across mono_start–mono_end intervals (combined)")

        # 图例限制以避免过多项目导致布局失败
        handles, labels = ax.get_legend_handles_labels()
        if len(labels) > max_legend_items:
            show_labels = labels[:max_legend_items] + ["..."]
            show_handles = handles[:max_legend_items]
            ax.legend(show_handles, show_labels, loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=8, ncol=1)
        else:
            ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=8, ncol=1)

        plt.tight_layout()
        plt.savefig(out_dir / "combined_blocks.png", dpi=200)
        plt.savefig(out_dir / "combined_blocks.pdf")
        plt.close(fig)

    # 单独图：默认只导出 PNG；可选导出单个多页 PDF
    pdf_pages = None
    if export_individual_pdf and not sliced_df.empty:
        pdf_pages = PdfPages(out_dir / "individual_blocks.pdf")

    count_png = 0
    for b_id, dfb in sliced_df.groupby("block_id"):
        a = int(dfb["mono_start"].iloc[0])
        b = int(dfb["mono_end"].iloc[0])
        fig = plt.figure(figsize=(10, 4))
        ax = plt.gca()
        ax.plot(dfb["pos"].values, dfb["score"].values, linewidth=1.2)
        ax.set_xlabel("pos")
        ax.set_ylabel(f"score (from column: {used_score_col})")
        ax.set_title(f"{prefix} | block_{b_id}: [{a}-{b}]")
        ax.set_xlim(a, b)
        plt.tight_layout()

        # PNG
        png_path = out_dir / f"block_{b_id:04d}_{a}-{b}.png"
        plt.savefig(png_path, dpi=200)
        count_png += 1

        # 可选：多页 PDF
        if pdf_pages is not None:
            pdf_pages.savefig(fig)

        plt.close(fig)

    if pdf_pages is not None:
        pdf_pages.close()

    log(f"[完成] 统计表: {meta_path.name}")
    if not sliced_df.empty:
        log(f"[完成] 合并图: combined_blocks.png / combined_blocks.pdf")
        log(f"[完成] 单独区间 PNG 数量: {count_png}")
        if export_individual_pdf:
            log(f"[完成] 单独区间多页 PDF: individual_blocks.pdf")
    if missing:
        log(f"[提示] 有 {len(missing)} 个区间在 score 中没有命中点，详见 blocks_no_points.tsv")


# ---------- 递归查找文件 ----------
def find_files_recursive(directory: Path, pattern: str):
    """递归查找匹配模式的文件"""
    return list(directory.rglob(pattern))


# ---------- 批量匹配 ----------
def main():
    ap = argparse.ArgumentParser(description="Batch plot score~pos by mono_start–mono_end intervals (blocks).")
    ap.add_argument("--blocks_dir", required=True, help="包含 *.hor_pairs_blocks.tsv 的目录（递归搜索子文件夹）")
    ap.add_argument("--score_dir",  required=True, help="包含 *.hor_score.tsv 的目录（递归搜索子文件夹）")
    ap.add_argument("--out_root",   required=True, help="输出的根目录（每一对会生成一个子目录）")
    ap.add_argument("--blocks_suffix", default=".hor_pairs_blocks.tsv", help="blocks 文件后缀（默认 .hor_pairs_blocks.tsv）")
    ap.add_argument("--score_suffix",  default=".hor_score.tsv",        help="score 文件后缀（默认 .hor_score.tsv）")
    ap.add_argument("--max_legend_items", type=int, default=30, help="合并图图例最多显示多少条（默认 30）")
    ap.add_argument("--individual-pdf", action="store_true", help="同时导出单个多页 PDF（individual_blocks.pdf）")
    args = ap.parse_args()

    blocks_dir = Path(args.blocks_dir)
    score_dir  = Path(args.score_dir)
    out_root   = Path(args.out_root)

    if not blocks_dir.is_dir():
        log(f"[错误] blocks_dir 不存在：{blocks_dir}")
        sys.exit(1)
    if not score_dir.is_dir():
        log(f"[错误] score_dir 不存在：{score_dir}")
        sys.exit(1)

    out_root.mkdir(parents=True, exist_ok=True)

    # 递归查找 blocks 文件
    blocks_files = find_files_recursive(blocks_dir, f"*{args.blocks_suffix}")
    if not blocks_files:
        log(f"[警告] 未在 {blocks_dir} 及其子文件夹下找到匹配 {args.blocks_suffix} 的文件")
        sys.exit(0)

    log(f"发现 {len(blocks_files)} 个 blocks 文件，开始匹配与处理……")

    n_ok, n_miss = 0, 0
    for bf in blocks_files:
        stem = bf.name
        if not stem.endswith(args.blocks_suffix):
            continue
        prefix = stem[: -len(args.blocks_suffix)]
        
        # 递归查找对应的 score 文件
        score_pattern = f"*{prefix}*{args.score_suffix}"
        score_candidates = find_files_recursive(score_dir, score_pattern)
        
        score_file = None
        if score_candidates:
            # 优先选择与 blocks 文件相同相对路径的 score 文件
            relative_path = bf.relative_to(blocks_dir)
            expected_score_path = score_dir / relative_path.with_name(f"{prefix}{args.score_suffix}")
            
            for candidate in score_candidates:
                if candidate == expected_score_path:
                    score_file = candidate
                    break
            
            # 如果没有找到相同相对路径的，选择第一个候选
            if score_file is None:
                score_file = score_candidates[0]
        
        if score_file is None or not score_file.exists():
            log(f"[跳过] 无法为 {bf} 找到对应 score：期望 {prefix}{args.score_suffix}")
            n_miss += 1
            continue

        # 创建输出目录，保留原始目录结构
        relative_path = bf.relative_to(blocks_dir)
        out_dir = out_root / relative_path.parent / prefix
        
        try:
            process_one(prefix=prefix,
                        blocks_file=bf,
                        score_file=score_file,
                        out_dir=out_dir,
                        max_legend_items=args.max_legend_items,
                        export_individual_pdf=args.individual_pdf)
            n_ok += 1
        except Exception as e:
            log(f"[失败] 处理 {prefix} 时出错：{e}")

    log(f"\n=== 全部完成 ===")
    log(f"成功：{n_ok} ；未匹配到 score：{n_miss}")
    log(f"输出根目录：{out_root}")


if __name__ == "__main__":
    main()