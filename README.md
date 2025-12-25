
# TimeTree 工作流程

使用 Snakemake 工作流程，通过 IQ-TREE 和 MCMCtree（IQ2MC 方法）从 OrthoFinder 结果推断时间校准的系统发育树。

## 概述

该流程：
1. 从 OrthoFinder 输出中提取单拷贝直系同源群（SCO）
2. 使用 MAFFT 对序列进行比对，使用 trimAl 进行修剪
3. 将比对结果连接成带有分区的超矩阵
4. 使用 IQ-TREE 推断最大似然（ML）树
5. 从表中注入化石/节点校准
6. 通过 `iqtree --dating mcmctree` 生成 Hessian/控制文件
7. 运行 MCMCtree 进行贝叶斯时间估计
8. 以 Newick 和 NEXUS 格式导出最终时间树

参考文献：[IQ-TREE Phylogenetic Dating](https://iqtree.github.io/doc/Dating)

## 要求

- Snakemake >= 7.0
- Conda/Mamba（用于环境管理）
- OrthoFinder 结果
- 包含化石/节点年龄约束的校准表

## 环境设置

### 方法 1：使用 Conda 自动管理环境（推荐）

Snakemake 会自动为每个规则创建和激活所需的 conda 环境。只需确保安装了 conda/mamba：

```bash
# 安装 Mamba（推荐，比 conda 快）
conda install -n base -c conda-forge mamba

# 或者使用 Miniforge（已包含 mamba）
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
```

**重要**：即使使用自动环境，MCMCtree 仍需要手动安装 IQ2MC 修改版（见下方"方法 2"中的安装步骤）。Snakemake 创建 `mcmctree` 环境后，需要在该环境中手动编译并安装修改版 MCMCtree。

### 方法 2：手动创建统一环境

如果您希望手动管理环境而不是让 Snakemake 自动创建：

```bash
# 创建主环境
conda create -n timetree python=3.12
conda activate timetree

# 安装 Snakemake
pip install snakemake

# 安装 Python 依赖
pip install  biopython ete3 pandas dendropy

# 安装系统发育工具
conda install  -c conda-forge -c bioconda \
    mafft trimal clipkit seqkit iqtree

# 安装 IQ2MC 的修改版 MCMCTree（必需）
# IQ2MC 工作流需要来自 https://github.com/iqtree/paml 的修改版
git clone https://github.com/iqtree/paml
cd paml/src
make -j
mv baseml basemlg chi2 codeml evolver infinitesites pamp yn00 mcmctree $CONDA_PREFIX/bin
cd ../..

# 运行工作流程（不使用 --use-conda）
snakemake --cores 8
```

**注意**：如果你的 `iqtree3` 不在 PATH，请在 `config/config.yaml` 指定绝对路径：
```yaml
binaries:
  iqtree_dating: \"/ABS/PATH/TO/iqtree3\"
```

### 验证环境

```bash
# 检查关键工具是否可用
mafft --version
trimal --version
iqtree --version
python -c "import Bio; import ete3; import pandas; print('Python packages OK')"
which mcmctree
```

## 快速开始

### 使用自动 Conda 环境（推荐）

```bash
# 1. 编辑 config/config.yaml 中的路径
vim config/config.yaml

# 2. 编辑 config/calibrations.tsv 中的校准点
vim config/calibrations.tsv

# 3. 运行工作流程（Snakemake 自动管理所有环境）
snakemake --use-conda --cores 8

# 干运行查看将执行的步骤
snakemake --use-conda --cores 8 -n

# 查看 DAG 图（需要 graphviz）
snakemake --dag | dot -Tpdf > workflow_dag.pdf
```

### 使用手动环境

```bash
# 先激活手动创建的环境
conda activate timetree

# 运行工作流程（不使用 --use-conda）
snakemake --cores 8
```

## 输入文件

### OrthoFinder 结果

将 `config/config.yaml` 中的 `results_dir` 指向您的 OrthoFinder 结果文件夹，例如：

```yaml
results_dir: "data/orthofinder/Results_May01"
```

所需文件：
- `Orthogroups/Orthogroups_SingleCopyOrthologues.txt`
- `Orthogroup_Sequences/OG*.fa`（或提供 `proteomes_dir` 作为备用）

### 校准表

创建 `config/calibrations.tsv`，包含以下列：

| 列 | 描述 |
|--------|-------------|
| `taxa_csv` | 定义 MRCA 节点的逗号分隔分类单元 |
| `min_age` | 最小年龄边界（例如：Myr） |
| `max_age` | 最大年龄边界 |
| `label` | （可选）此校准的名称 |
| `prior` | （可选）先验类型：`B`（软边界）、`L`（下界）、`U`（上界） |

示例：
```tsv
taxa_csv	min_age	max_age	label	prior
Arabidopsis,Oryza	120	150	Monocot-Dicot	B
Zea,Sorghum	10	20	Grass_Crown	B
```

### 外群分类单元

为可靠地确定根，请在 `config.yaml` 中指定外群分类单元：

```yaml
outgroup_taxa: ["outgroup_species1", "outgroup_species2"]
```

## 配置

`config/config.yaml` 中的关键选项：

```yaml
# IQ-TREE 设置
iqtree:
  model: "MFP"        # ModelFinder 或特定模型
  use_partition: true # 按基因分区
  bootstrap: 1000     # UFBoot 重复次数（0 表示禁用）

# MCMCtree 设置  
mcmctree:
  clock_model: "IND"  # IND、CORR 或 EQUAL
  burnin: 20000
  nsample: 50000
```

## 输出文件

结果写入 `results/timetree/`：

| 文件 | 描述 |
|------|-------------|
| `supermatrix.phy` | 连接的比对（PHYLIP 格式） |
| `partitions.nex` | 分区定义（NEXUS 格式） |
| `species_tree.rooted.nwk` | 有根的 ML 树 |
| `species_tree.calibrated.nwk` | 带校准注释的 ML 树 |
| `iq2mc/species.mcmctree.hessian` | MCMCtree 的 Hessian |
| `iq2mc/species.mcmctree.ctl` | MCMCtree 控制文件 |
| `mcmctree/FigTree.tre` | MCMCtree 输出（FigTree 格式） |
| `timetree.final.nwk` | 最终时间树（Newick 格式） |
| `timetree.final.nex` | 最终时间树（NEXUS 格式，用于 FigTree） |

## 可视化结果

在 FigTree 中打开 `results/timetree/timetree.final.nex`：
1. 在左侧面板中点击 "Node Labels"
2. 将 "Display" 设置为 "date" 以显示估计的年龄

## 注意事项

### MCMCtree 二进制文件

IQ2MC 工作流程**必须**使用来自 [https://github.com/iqtree/paml](https://github.com/iqtree/paml) 的修改版 MCMCtree。标准 conda 的 `paml` 包可能不兼容。

**手动安装方法**（见上方"方法 2：手动创建统一环境"）：
```bash
git clone https://github.com/iqtree/paml
cd paml/src
make -j
mv baseml basemlg chi2 codeml evolver infinitesites pamp yn00 mcmctree $CONDA_PREFIX/bin
```

**验证安装**：
```bash
# 应显示用法信息而不是错误
mcmctree 2>&1 | head -10
which mcmctree  # 确认在 PATH 中
```

### MCMC 收敛性

通过以下方式检查收敛性：
1. 在 Tracer 中加载 `mcmctree/mcmc.txt`
2. 确保所有参数的 ESS > 200
3. 如果未收敛，请在配置中增加 `burnin` 和 `nsample`

## 故障排除

**未提取到 SCO**：检查 OrthoFinder 输出和蛋白质组文件之间的物种名称是否匹配。

**校准失败**：确保 `calibrations.tsv` 中的分类单元与树中的末端标签完全匹配。

**MCMCtree 错误**：验证校准树具有有效的 MCMCtree 语法。检查 `calibration_mapping.log`。

## 许可证

MIT