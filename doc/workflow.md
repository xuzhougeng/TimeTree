# TimeTree Workflow 完整调用链

## 概述

```
OrthoFinder → MSA → Supermatrix → IQ-TREE → IQ2MC → MCMCtree → TimeTree
```

## 详细流程图

```
┌─────────────────────────────────────────────────────────────────────────────────┐
│                              TimeTree Workflow                                   │
│                    OrthoFinder → MSA → IQ-TREE → MCMCtree                       │
└─────────────────────────────────────────────────────────────────────────────────┘

                              ┌──────────────────┐
                              │   OrthoFinder    │  (外部运行)
                              │    Results/      │
                              └────────┬─────────┘
                                       │
                                       ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│  1. extract_sco (checkpoint)                           [orthogroups.smk]         │
│     ────────────────────────────────────────────────────────────────────         │
│     输入: Orthogroups_SingleCopyOrthologues.txt + Orthogroup_Sequences/          │
│     脚本: extract_sco_fastas.py                                                  │
│     输出: work/orthogroups_sco/OG*.faa (单拷贝直系同源基因)                       │
└──────────────────────────────────────────────────────────────────────────────────┘
                                       │
                            ┌──────────┴──────────┐
                            ▼                     ▼
              ┌─────────────────────┐  ┌─────────────────────┐
              │   OG0008709.faa     │  │   OG0008714.faa     │  ... (N个基因家族)
              └──────────┬──────────┘  └──────────┬──────────┘
                         │                        │
                         ▼                        ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│  2. mafft_align (per OG)                                    [msa.smk]            │
│     ────────────────────────────────────────────────────────────────────         │
│     命令: mafft --auto --thread 2 {og}.faa > {og}.raw.aln                        │
│     输出: work/msa/OG*.raw.aln                                                   │
└──────────────────────────────────────────────────────────────────────────────────┘
                                       │
                                       ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│  3. trim_alignment (per OG)                                 [msa.smk]            │
│     ────────────────────────────────────────────────────────────────────         │
│     命令: trimal -in {og}.raw.aln -out {og}.aln.faa -automated1                  │
│     输出: work/msa/OG*.aln.faa                                                   │
└──────────────────────────────────────────────────────────────────────────────────┘
                                       │
                            (收集所有702个trimmed alignments)
                                       │
                                       ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│  4. concat_supermatrix                                   [supermatrix.smk]       │
│     ────────────────────────────────────────────────────────────────────         │
│     脚本: concat_alignments.py                                                   │
│     输出: results/timetree/supermatrix.phy    (超级矩阵 PHYLIP)                  │
│           results/timetree/supermatrix.faa    (超级矩阵 FASTA)                   │
│           results/timetree/partitions.nex     (分区文件 NEXUS)                   │
└──────────────────────────────────────────────────────────────────────────────────┘
                                       │
                                       ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│  5. iqtree_ml (IQ2MC Step 1)                               [iqtree.smk]          │
│     ────────────────────────────────────────────────────────────────────         │
│     命令: iqtree3 -s supermatrix.phy [-p partitions.nex] -m MFP -B 1000          │
│                                       ↑                                          │
│                          use_partition: true/false 控制是否使用                  │
│     输出: results/timetree/iqtree/species.treefile (ML树)                        │
└──────────────────────────────────────────────────────────────────────────────────┘
                                       │
                                       ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│  6. root_tree                                              [iqtree.smk]          │
│     ────────────────────────────────────────────────────────────────────         │
│     脚本: root_tree.py                                                           │
│     方法: outgroup rooting 或 midpoint rooting                                   │
│     输出: results/timetree/species_tree.rooted.nwk                               │
└──────────────────────────────────────────────────────────────────────────────────┘
                                       │
                                       ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│  7. annotate_calibrations                            [calibrate_tree.smk]        │
│     ────────────────────────────────────────────────────────────────────         │
│     脚本: annotate_calibrations.py                                               │
│     输入: rooted tree + calibrations.tsv                                         │
│     输出: results/timetree/species_tree.calibrated.nwk                           │
│           results/timetree/calibration_mapping.log                               │
│     格式: 节点标注 'B(min,max)' 格式的校准约束                                   │
└──────────────────────────────────────────────────────────────────────────────────┘
                                       │
                                       ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│  8. iqtree_dating_mcmctree (IQ2MC Step 2)                  [iq2mc.smk]           │
│     ────────────────────────────────────────────────────────────────────         │
│     命令: iqtree3 -s supermatrix.phy [-p partitions.nex] -m MFP                  │
│            -te species_tree.calibrated.nwk --dating mcmctree                     │
│            --mcmc-iter 10000,100,20000 --mcmc-bds 1.0,1.0,0.5 --mcmc-clock IND   │
│     输出: results/timetree/iq2mc/species.mcmctree.hessian  (Hessian矩阵)         │
│           results/timetree/iq2mc/species.mcmctree.ctl      (控制文件)            │
│           results/timetree/iq2mc/species.rooted.nwk        (树文件)              │
│           results/timetree/iq2mc/species.dummy.phy         (伪比对文件)          │
└──────────────────────────────────────────────────────────────────────────────────┘
                                       │
                                       ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│  9. run_mcmctree (IQ2MC Step 3)                          [mcmctree.smk]          │
│     ────────────────────────────────────────────────────────────────────         │
│     前处理脚本:                                                                  │
│       1. localize_mcmctree_ctl.py  - 路径转为相对路径                            │
│       2. clean_mcmctree_tree.py    - 'B_x_y_' → 'B(x,y)' 格式转换                │
│       3. adjust_rootage.py         - 自动设置 RootAge = max_calibration * 1.5   │
│       4. ensure_mcmctree_print.py  - 确保 print=1                                │
│     命令: mcmctree species.mcmctree.ctl                                          │
│     输出: results/timetree/mcmctree/mcmc.txt       (MCMC采样)                    │
│           results/timetree/mcmctree/FigTree.tre    (定年树)                      │
└──────────────────────────────────────────────────────────────────────────────────┘
                                       │
                                       ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│  10. postprocess_timetree                                [mcmctree.smk]          │
│      ───────────────────────────────────────────────────────────────────         │
│      脚本: mcmctree_postprocess.py                                               │
│      输出: results/timetree/timetree.final.nwk    (最终时间树 Newick)            │
│            results/timetree/timetree.final.nex    (最终时间树 NEXUS)             │
└──────────────────────────────────────────────────────────────────────────────────┘
                                       │
                                       ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│  11. plot_timetree                                    [visualization.smk]        │
│      ───────────────────────────────────────────────────────────────────         │
│      脚本: plot_timetree.R (ggtree)                                              │
│      输出: results/timetree/timetree.pdf                                         │
│            results/timetree/timetree.png                                         │
│            results/timetree/timetree.svg                                         │
└──────────────────────────────────────────────────────────────────────────────────┘

                              ┌──────────────────┐
                              │   rule all       │
                              │  (最终目标)      │
                              └──────────────────┘
```

## 简化版流程图

```
OrthoFinder结果
      │
      ▼
 extract_sco ──────────────────► 702个 OG*.faa
      │
      ├──► mafft_align (×702) ──► OG*.raw.aln
      │
      ├──► trim_alignment (×702) ► OG*.aln.faa
      │
      ▼
 concat_supermatrix ───────────► supermatrix.phy + partitions.nex
      │
      ▼
 iqtree_ml ────────────────────► species.treefile (ML树)
      │
      ▼
 root_tree ────────────────────► species_tree.rooted.nwk
      │
      ▼
 annotate_calibrations ────────► species_tree.calibrated.nwk
      │
      ▼
 iqtree_dating_mcmctree ───────► .hessian + .ctl (IQ2MC)
      │
      ▼
 run_mcmctree ─────────────────► FigTree.tre (定年结果)
      │
      ▼
 postprocess_timetree ─────────► timetree.final.nwk
      │
      ▼
 plot_timetree ────────────────► timetree.pdf/png/svg
```

## 规则文件结构

| 文件 | 包含的规则 | 功能 |
|------|-----------|------|
| `orthogroups.smk` | extract_sco | 提取单拷贝直系同源基因 |
| `msa.smk` | mafft_align, trim_alignment | 多序列比对和修剪 |
| `supermatrix.smk` | concat_supermatrix | 合并为超级矩阵 |
| `iqtree.smk` | iqtree_ml, root_tree | ML树推断和根化 |
| `calibrate_tree.smk` | annotate_calibrations | 注入化石校准点 |
| `iq2mc.smk` | iqtree_dating_mcmctree | IQ2MC Step 2: 生成Hessian |
| `mcmctree.smk` | run_mcmctree, postprocess_timetree | MCMCtree定年 |
| `visualization.smk` | plot_timetree | ggtree可视化 |

## 关键配置项

| 配置项 | 影响的规则 | 作用 |
|-------|-----------|------|
| `use_partition: true/false` | iqtree_ml, iqtree_dating_mcmctree | 是否使用分区模型 |
| `outgroup_taxa` | iqtree_ml, root_tree | 外群指定/根化方式 |
| `calibrations.tsv` | annotate_calibrations | 化石校准点 |
| `clock_model: IND/COR` | iqtree_dating_mcmctree | 分子钟模型 |
| `model: MFP` | iqtree_ml, iqtree_dating_mcmctree | 替代模型选择 |
| `bootstrap` | iqtree_ml | Bootstrap重复次数 |

## 输出文件说明

### 中间文件 (`work/`)
- `work/orthogroups_sco/` - 单拷贝基因FASTA
- `work/msa/` - 比对文件
- `work/logs/` - 各步骤日志

### 结果文件 (`results/timetree/`)
- `supermatrix.phy` - 超级矩阵
- `partitions.nex` - 分区定义
- `iqtree/species.treefile` - ML树
- `species_tree.rooted.nwk` - 根化树
- `species_tree.calibrated.nwk` - 校准树
- `iq2mc/species.mcmctree.hessian` - Hessian矩阵
- `mcmctree/FigTree.tre` - MCMCtree输出
- `timetree.final.nwk` - **最终定年树**
- `timetree.pdf/png/svg` - 可视化图片

## 常见问题排查

### 1. Hessian矩阵溢出
**症状**: 分支长度为0，lnL值为天文数字
**原因**: 分区数量过多（如702个）
**解决**: 设置 `use_partition: false`

### 2. 校准点未生效
**检查**: `calibration_mapping.log` 中是否有 FAIL
**原因**: 物种名不匹配
**解决**: 确保 calibrations.tsv 中的物种名与树中一致

### 3. MCMCtree未收敛
**检查**: 用Tracer查看 `mcmc.txt` 的ESS值
**解决**: 增加 burnin 和 nsample
