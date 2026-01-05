#!/usr/bin/env python3
"""
蛋白质-蛋白质复合物界面过滤脚本 (基于PyMOL)
过滤出受体和配体中彼此距离在指定范围内的残基
"""
import sys
import argparse
from pathlib import Path

import pymol
from pymol import cmd


def filter_interface_pymol(receptor_file, ligand_file, cutoff=6.0, output_dir=None):
    """
    使用PyMOL过滤蛋白质-蛋白质复合物界面残基

    Args:
        receptor_file: 受体PDB文件路径
        ligand_file: 配体PDB文件路径
        cutoff: 距离阈值(Å)
        output_dir: 输出目录，默认为输入文件所在目录

    Returns:
        tuple: (受体输出文件, 配体输出文件, 保留的受体原子数, 保留的配体原子数)
    """
    # 静默模式启动PyMOL
    pymol.finish_launching(["pymol", "-cq"])

    try:
        # 对象命名
        rec_name = "receptor_obj"
        lig_name = "ligand_obj"

        print(f"[-] 正在加载: {receptor_file} 和 {ligand_file}")

        # 加载结构
        cmd.load(receptor_file, rec_name)
        cmd.load(ligand_file, lig_name)

        print(f"[-] 正在筛选 {cutoff} Å 范围内的界面残基...")

        # 选择界面残基（按残基选择）
        cmd.select("rec_interface", f"byres ({rec_name} within {cutoff} of {lig_name})")
        cmd.select("lig_interface", f"byres ({lig_name} within {cutoff} of {rec_name})")

        # 统计原子数
        rec_atom_count = cmd.count_atoms("rec_interface")
        lig_atom_count = cmd.count_atoms("lig_interface")
        print(f"    受体保留原子数: {rec_atom_count}")
        print(f"    配体保留原子数: {lig_atom_count}")

        # 确定输出路径
        receptor_path = Path(receptor_file)
        ligand_path = Path(ligand_file)

        if output_dir is None:
            output_dir = receptor_path.parent
        else:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

        # 生成输出文件名
        rec_stem = receptor_path.stem
        lig_stem = ligand_path.stem
        suffix = f"_{int(cutoff)}A"

        out_rec = output_dir / f"{rec_stem}{suffix}.pdb"
        out_lig = output_dir / f"{lig_stem}{suffix}.pdb"

        # 保存文件
        if rec_atom_count > 0:
            cmd.save(str(out_rec), "rec_interface")
            print(f"[+] 受体文件已生成: {out_rec}")
        else:
            print("[!] 受体无残基在范围内，未生成文件")
            out_rec = None

        if lig_atom_count > 0:
            cmd.save(str(out_lig), "lig_interface")
            print(f"[+] 配体文件已生成: {out_lig}")
        else:
            print("[!] 配体无残基在范围内，未生成文件")
            out_lig = None

        return out_rec, out_lig, rec_atom_count, lig_atom_count

    except Exception as e:
        print(f"[!] 错误: {e}")
        return None, None, 0, 0
    finally:
        # 确保退出PyMOL
        try:
            cmd.quit()
        except Exception:
            pass


def main():
    parser = argparse.ArgumentParser(
        description="使用PyMOL过滤蛋白质-蛋白质复合物界面残基"
    )
    parser.add_argument(
        "--receptor", type=str, required=True, help="受体蛋白PDB文件路径"
    )
    parser.add_argument(
        "--ligand", type=str, required=True, help="配体蛋白PDB文件路径"
    )
    parser.add_argument(
        "--cutoff", type=float, default=6.0, help="距离阈值(Å)，默认6.0"
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default=None,
        help="输出目录，默认为受体文件所在目录",
    )

    args = parser.parse_args()

    # 验证输入文件
    if not Path(args.receptor).is_file():
        print(f"错误: 受体文件不存在: {args.receptor}")
        return 1

    if not Path(args.ligand).is_file():
        print(f"错误: 配体文件不存在: {args.ligand}")
        return 1

    print("处理蛋白质-蛋白质复合物:")
    print(f"  受体: {args.receptor}")
    print(f"  配体: {args.ligand}")
    print(f"  距离阈值: {args.cutoff} Å")
    print(f"  输出目录: {args.output_dir or '默认(受体同目录)'}")
    print()

    # 执行过滤
    out_rec, out_lig, rec_count, lig_count = filter_interface_pymol(
        args.receptor, args.ligand, args.cutoff, args.output_dir
    )

    if out_rec is None and out_lig is None:
        print("\n[!] 处理失败: 无残基在指定范围内")
        return 1

    print("\n[+] 处理完成!")
    print(f"    受体保留原子: {rec_count}")
    print(f"    配体保留原子: {lig_count}")

    return 0


if __name__ == "__main__":
    sys.exit(main())