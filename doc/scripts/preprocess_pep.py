#!/usr/bin/env python3
"""
蛋白序列预处理脚本

功能：
1. 只保留序列ID（类似 seqkit seq -i）
2. 去掉序列末尾的 *（终止密码子符号）
3. 删除序列中间有 * 的蛋白（内部终止密码子）
4. 将非20个标准氨基酸替换为模糊码 X
"""

import os
import sys
import argparse
from pathlib import Path

# 20个标准氨基酸
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def parse_fasta(filepath):
    """解析FASTA文件，返回(id, sequence)列表"""
    sequences = []
    current_id = None
    current_seq = []

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_id is not None:
                    sequences.append((current_id, ''.join(current_seq)))
                # 只保留ID（第一个空格前的部分）
                header = line[1:]  # 去掉 >
                current_id = header.split()[0]  # 只取第一个空格前的部分
                current_seq = []
            else:
                current_seq.append(line)

        # 处理最后一条序列
        if current_id is not None:
            sequences.append((current_id, ''.join(current_seq)))

    return sequences


def process_sequence(seq_id, sequence):
    """
    处理单条序列
    返回: (processed_id, processed_seq) 或 None（如果序列应被过滤）
    """
    # 去掉末尾的 *
    if sequence.endswith('*'):
        sequence = sequence[:-1]

    # 检查序列中间是否有 *（内部终止密码子）
    if '*' in sequence:
        return None

    # 将非标准氨基酸替换为 X
    processed_seq = []
    for aa in sequence.upper():
        if aa in STANDARD_AA:
            processed_seq.append(aa)
        else:
            processed_seq.append('X')

    return (seq_id, ''.join(processed_seq))


def write_fasta(sequences, output_path, line_width=80):
    """写出FASTA文件"""
    with open(output_path, 'w') as f:
        for seq_id, sequence in sequences:
            f.write(f'>{seq_id}\n')
            # 按行宽输出序列
            for i in range(0, len(sequence), line_width):
                f.write(sequence[i:i+line_width] + '\n')


def process_file(input_path, output_path):
    """处理单个FASTA文件"""
    sequences = parse_fasta(input_path)

    processed = []
    filtered_count = 0
    modified_count = 0

    for seq_id, sequence in sequences:
        result = process_sequence(seq_id, sequence)
        if result is None:
            filtered_count += 1
        else:
            processed.append(result)
            # 检查是否有非标准氨基酸被替换
            if 'X' in result[1]:
                modified_count += 1

    write_fasta(processed, output_path)

    return len(sequences), len(processed), filtered_count, modified_count


def main():
    parser = argparse.ArgumentParser(
        description='蛋白序列预处理：简化ID、去除终止密码子、过滤内部终止密码子、替换非标准氨基酸'
    )
    parser.add_argument('-i', '--input', required=True,
                        help='输入目录或单个FASTA文件')
    parser.add_argument('-o', '--output', required=True,
                        help='输出目录或单个FASTA文件')
    parser.add_argument('-s', '--suffix', default='.faa',
                        help='输入文件后缀 (默认: .faa)')

    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)

    if input_path.is_file():
        # 处理单个文件
        output_path.parent.mkdir(parents=True, exist_ok=True)
        total, kept, filtered, modified = process_file(input_path, output_path)
        print(f"{input_path.name}: 总计 {total} 条, 保留 {kept} 条, "
              f"过滤 {filtered} 条(内部终止密码子), {modified} 条含非标准AA")

    elif input_path.is_dir():
        # 处理目录
        output_path.mkdir(parents=True, exist_ok=True)

        files = list(input_path.glob(f'*{args.suffix}'))
        if not files:
            print(f"警告: 在 {input_path} 中未找到 {args.suffix} 文件")
            return

        total_seqs = 0
        total_kept = 0
        total_filtered = 0
        total_modified = 0

        for fasta_file in sorted(files):
            out_file = output_path / fasta_file.name
            total, kept, filtered, modified = process_file(fasta_file, out_file)
            total_seqs += total
            total_kept += kept
            total_filtered += filtered
            total_modified += modified
            print(f"{fasta_file.name}: 总计 {total} 条, 保留 {kept} 条, "
                  f"过滤 {filtered} 条, {modified} 条含非标准AA")

        print(f"\n总计: {total_seqs} 条序列, 保留 {total_kept} 条, "
              f"过滤 {total_filtered} 条, {total_modified} 条含非标准AA被修正")

    else:
        print(f"错误: {input_path} 不是有效的文件或目录")
        sys.exit(1)


if __name__ == '__main__':
    main()