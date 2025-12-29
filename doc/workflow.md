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
│     目的: 从OrthoFinder结果中筛选单拷贝直系同源基因(SCO)，确保每个物种只有      │
│           一个拷贝，避免旁系同源基因干扰系统发育分析                             │
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
│     目的: 对每个基因家族进行多序列比对，识别同源位点，为后续系统发育分析         │
│           建立位置对应关系                                                       │
│     输入: work/orthogroups_sco/OG*.faa                                           │
│     命令: mafft --auto --thread 2 {og}.faa > {og}.raw.aln                        │
│     输出: work/msa/OG*.raw.aln                                                   │
└──────────────────────────────────────────────────────────────────────────────────┘
                                       │
                                       ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│  3. trim_alignment (per OG)                                 [msa.smk]            │
│     ────────────────────────────────────────────────────────────────────         │
│     目的: 移除比对中的低质量区域(高gap、高变异位点)，保留可靠的系统发育          │
│           信息位点，减少比对误差对树推断的影响                                   │
│     输入: work/msa/OG*.raw.aln                                                   │
│     命令: trimal -in {og}.raw.aln -out {og}.aln.faa -automated1                  │
│     输出: work/msa/OG*.aln.faa                                                   │
└──────────────────────────────────────────────────────────────────────────────────┘
                                       │
                            (收集所有trimmed alignments)
                                       │
                                       ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│  4. concat_supermatrix                                   [supermatrix.smk]       │
│     ────────────────────────────────────────────────────────────────────         │
│     目的: 将多个单基因比对横向拼接成超级矩阵，整合多基因信息提高系统发育         │
│           分辨率；同时生成分区文件记录每个基因在矩阵中的位置                     │
│     输入: work/msa/OG*.aln.faa (所有trimmed比对)                                 │
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
│     目的: 使用最大似然法推断物种树拓扑结构，为后续分歧时间估计提供树骨架；       │
│           可选择是否使用分区模型以考虑不同基因的进化速率异质性                   │
│     输入: results/timetree/supermatrix.phy                                       │
│           results/timetree/partitions.nex (可选)                                 │
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
│     目的: 为无根树指定根的位置，确定进化方向；可使用外群法(已知最早分化的        │
│           类群)或中点法(假设进化速率恒定)                                        │
│     输入: results/timetree/iqtree/species.treefile                               │
│     脚本: root_tree.py                                                           │
│     方法: outgroup rooting 或 midpoint rooting                                   │
│     输出: results/timetree/species_tree.rooted.nwk                               │
└──────────────────────────────────────────────────────────────────────────────────┘
                                       │
                                       ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│  7. annotate_calibrations                            [calibrate_tree.smk]        │
│     ────────────────────────────────────────────────────────────────────         │
│     目的: 将化石校准点(已知分歧时间约束)注入到树的对应节点上，为贝叶斯           │
│           定年提供时间先验信息，锚定绝对时间尺度                                 │
│     输入: results/timetree/species_tree.rooted.nwk                               │
│           calibrations.tsv (化石校准表)                                          │
│     脚本: annotate_calibrations.py                                               │
│     输出: results/timetree/species_tree.calibrated.nwk                           │
│           results/timetree/calibration_mapping.log                               │
│     格式: 节点标注 'B(min,max)' 格式的校准约束                                   │
└──────────────────────────────────────────────────────────────────────────────────┘
                                       │
                                       ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│  8. iqtree_dating_mcmctree (IQ2MC Step 2)                  [iq2mc.smk]           │
│     ────────────────────────────────────────────────────────────────────         │
│     目的: 使用IQ-TREE计算似然函数的梯度和Hessian矩阵，用于MCMCtree的近似         │
│           似然计算，大幅加速贝叶斯定年过程(IQ2MC方法核心步骤)                    │
│     输入: results/timetree/supermatrix.phy                                       │
│           results/timetree/partitions.nex (可选)                                 │
│           results/timetree/species_tree.calibrated.nwk                           │
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
│     目的: 使用贝叶斯MCMC方法估计各节点的分歧时间，结合分子数据、化石校准         │
│           和松弛分子钟模型，给出时间估计及置信区间                               │
│     输入: results/timetree/iq2mc/species.mcmctree.ctl                            │
│           results/timetree/iq2mc/species.mcmctree.hessian                        │
│           results/timetree/iq2mc/species.rooted.nwk                              │
│           results/timetree/iq2mc/species.dummy.phy                               │
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
│      目的: 将MCMCtree输出转换为标准格式，提取节点年龄和置信区间，生成             │
│            可用于下游分析和可视化的时间树文件                                    │
│      输入: results/timetree/mcmctree/FigTree.tre                                 │
│            results/timetree/mcmctree/mcmc.txt                                    │
│      脚本: mcmctree_postprocess.py                                               │
│      输出: results/timetree/timetree.final.nwk    (最终时间树 Newick)            │
│            results/timetree/timetree.final.nex    (最终时间树 NEXUS)             │
└──────────────────────────────────────────────────────────────────────────────────┘
                                       │
                                       ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│  11. plot_timetree                                    [visualization.smk]        │
│      ───────────────────────────────────────────────────────────────────         │
│      目的: 生成发表级别的时间树可视化图，展示物种关系、分歧时间和置信区间，       │
│            包含地质年代标尺                                                      │
│      输入: results/timetree/mcmctree/FigTree.tre                                 │
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
 extract_sco ──────────────────► N个 OG*.faa
      │                          (筛选单拷贝直系同源基因)
      │
      ├──► mafft_align (×N) ───► OG*.raw.aln
      │                          (多序列比对)
      │
      ├──► trim_alignment (×N) ► OG*.aln.faa
      │                          (去除低质量区域)
      ▼
 concat_supermatrix ───────────► supermatrix.phy + partitions.nex
      │                          (合并为超级矩阵)
      ▼
 iqtree_ml ────────────────────► species.treefile (ML树)
      │                          (最大似然树推断)
      ▼
 root_tree ────────────────────► species_tree.rooted.nwk
      │                          (树根化)
      ▼
 annotate_calibrations ────────► species_tree.calibrated.nwk
      │                          (注入化石校准点)
      ▼
 iqtree_dating_mcmctree ───────► .hessian + .ctl
      │                          (计算Hessian矩阵)
      ▼
 run_mcmctree ─────────────────► FigTree.tre
      │                          (贝叶斯分歧时间估计)
      ▼
 postprocess_timetree ─────────► timetree.final.nwk
      │                          (格式转换)
      ▼
 plot_timetree ────────────────► timetree.pdf/png/svg
                                 (可视化)
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
