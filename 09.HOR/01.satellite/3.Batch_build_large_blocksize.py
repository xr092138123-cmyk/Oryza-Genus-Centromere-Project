#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
批量处理 HOR pair 文件，规则：
- 仅使用 shift==1 的行；
- 从遇到的第一个 block_size==3 开始；在"见过 >3"之后再次遇到下一个 3，才结算上一段；
- 段的区间 = 起点：该段第一个 3 行的 A_pos 第一个数字；终点：段内"达到最大 block_size 的行（若多次出现取最后一次）"的 B_pos 最后一个数字；
- 输出列：block_idx, mono_start, mono_end, mono_count, blocksize_max, win_bp（= mono_count × 文件名中的 w 数值）

修改版：处理指定目录下的所有子文件夹
"""

import argparse
import re
import pandas as pd
from pathlib import Path


def first_int(s):
    if pd.isna(s):
        return None
    m = re.search(r"\d+", str(s))
    return int(m.group()) if m else None


def last_int(s):
    if pd.isna(s):
        return None
    m = list(re.finditer(r"\d+", str(s)))
    return int(m[-1].group()) if m else None


def get_w_from_filename(name: str):
    """从文件名中提取 w 数值，例如 AA.Chr01.w100.xxx.tsv -> 100"""
    m = re.search(r"w(\d+)", name)
    return int(m.group(1)) if m else None


def process_file(infile: Path) -> pd.DataFrame:
    """按规则处理单个 pair TSV，返回结果 DataFrame（若无结果返回空表）"""
    df = pd.read_csv(infile, sep="\t", dtype=str, engine="python")

    need = {"HOR_ID", "block_size", "shift", "A_pos", "B_pos"}
    if not need.issubset(df.columns):
        print(f"[跳过] {infile.name} 缺少列: {sorted(need - set(df.columns))}")
        return pd.DataFrame()

    # 仅保留 shift==1，保持原始顺序
    df = df[df["shift"].astype(str) == "1"].reset_index(drop=True)
    if df.empty:
        print(f"[空] {infile.name}（无 shift==1 行）")
        return pd.DataFrame()

    # 解析数字
    df["block_size_i"] = pd.to_numeric(df["block_size"], errors="coerce")
    df["A_first"] = df["A_pos"].map(first_int)
    df["B_last"] = df["B_pos"].map(last_int)
    df = df.dropna(subset=["block_size_i", "A_first", "B_last"]).reset_index(drop=True)
    if df.empty:
        print(f"[空] {infile.name}（关键列解析为空）")
        return pd.DataFrame()

    # 规则状态变量
    blocks = []
    in_block = False
    seen_non3 = False            # 是否在当前段里见过 >3
    start_a = None
    last_b_of_start = None       # 整段都只有 3 时使用
    max_bs = None
    end_b_at_max = None          # 达到最大 block_size 的那行的 B_last（若多次出现取最后一次）

    for _, r in df.iterrows():
        bs = int(r["block_size_i"])
        a = int(r["A_first"])
        b = int(r["B_last"])

        if bs == 3:
            if not in_block:
                # 开启新段（遇到第一个 3）
                in_block = True
                seen_non3 = False
                start_a = a
                last_b_of_start = b
                max_bs = 3
                end_b_at_max = b
            else:
                # 段内再次遇到 3
                if seen_non3:
                    # 见过 >3 后再次见到 3 -> 结算上一段
                    final_end = end_b_at_max if end_b_at_max is not None else last_b_of_start
                    blocks.append({
                        "mono_start": start_a,
                        "mono_end": final_end,
                        "mono_count": final_end - start_a + 1,
                        "blocksize_max": max_bs if max_bs is not None else 3,
                    })
                    # 以当前 3 作为新段起点
                    in_block = True
                    seen_non3 = False
                    start_a = a
                    last_b_of_start = b
                    max_bs = 3
                    end_b_at_max = b
                else:
                    # 连续的多个 3 仍属于同一段
                    last_b_of_start = b
        else:
            if in_block:
                seen_non3 = True
                if max_bs is None or bs > max_bs:
                    max_bs = bs
                    end_b_at_max = b
                elif bs == max_bs:
                    # 最大值重复出现，取最后一次的 B_last
                    end_b_at_max = b
            else:
                # 段外遇到 >3，忽略，直到遇到第一个 3 才开始
                pass

    # 文件末尾收尾
    if in_block:
        final_end = end_b_at_max if (seen_non3 and end_b_at_max is not None) else last_b_of_start
        blocks.append({
            "mono_start": start_a,
            "mono_end": final_end,
            "mono_count": final_end - start_a + 1,
            "blocksize_max": max_bs if max_bs is not None else 3,
        })

    out = pd.DataFrame(blocks)
    if out.empty:
        return out

    out.insert(0, "block_idx", range(1, len(out) + 1))

    # win_bp 列：mono_count × w
    w = get_w_from_filename(infile.name)
    out["w"] = w
    out["win_bp"] = out["mono_count"] * out["w"] if w else None

    # 附带文件名便于汇总
    out.insert(0, "file", infile.name)
    # 添加文件夹信息
    out.insert(0, "folder", infile.parent.name)
    return out[
        ["folder", "file", "block_idx", "mono_start", "mono_end", "mono_count", "blocksize_max", "w", "win_bp"]
    ]


def main():
    ap = argparse.ArgumentParser(description="大 block 提取（处理所有子文件夹）")
    ap.add_argument("--indir", required=True, help="输入目录（包含多个子文件夹，每个子文件夹包含 pair tsv）")
    ap.add_argument("--outdir", required=True, help="输出目录")
    ap.add_argument("--glob", default="*.tsv", help="文件匹配模式，默认 *.tsv")
    args = ap.parse_args()

    indir = Path(args.indir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 获取所有子文件夹
    subfolders = [f for f in indir.iterdir() if f.is_dir()]
    if not subfolders:
        print(f"[警告] {indir} 下未找到任何子文件夹")
        return

    print(f"📁 找到 {len(subfolders)} 个子文件夹")

    all_merged = []
    
    for folder in subfolders:
        print(f"\n🔍 处理文件夹: {folder.name}")
        
        # 在子文件夹中查找文件
        files = sorted(folder.glob(args.glob))
        if not files:
            print(f"  [跳过] {folder.name} 中没有匹配 {args.glob} 的文件")
            continue
            
        print(f"  找到 {len(files)} 个文件")
        
        folder_merged = []
        for f in files:
            out = process_file(f)
            if out is None or out.empty:
                continue
                
            # 为每个子文件夹创建输出目录
            folder_outdir = outdir / folder.name
            folder_outdir.mkdir(parents=True, exist_ok=True)
            
            # 每个文件各自保存一份
            per_file_path = folder_outdir / f"{f.stem}_blocks.tsv"
            out.to_csv(per_file_path, sep="\t", index=False)
            print(f"  [完成] {f.name} -> {folder.name}/{per_file_path.name}（段数 {len(out)}）")
            folder_merged.append(out)
        
        # 为每个子文件夹生成汇总表
        if folder_merged:
            folder_df = pd.concat(folder_merged, ignore_index=True)
            folder_summary_path = folder_outdir / f"{folder.name}_summary_blocks.tsv"
            folder_df.to_csv(folder_summary_path, sep="\t", index=False)
            print(f"  [汇总] {folder.name}_summary_blocks.tsv 已生成（{len(folder_df)} 行）")
            
            all_merged.append(folder_df)

    # 合并所有文件夹的总表
    if all_merged:
        all_df = pd.concat(all_merged, ignore_index=True)
        all_path = outdir / "all_folders_blocks.tsv"
        all_df.to_csv(all_path, sep="\t", index=False)
        print(f"\n🎉 全局汇总: all_folders_blocks.tsv 已生成（{len(all_df)} 行，来自 {len(subfolders)} 个文件夹）")
    else:
        print("\n[提示] 没有任何有效输出。")


if __name__ == "__main__":
    main()