
[IQ2MC](./IQ2MC.md)核心在于获取 MCMCTree运行所需的 hessian 文件，里面包括

- 梯度向量（gradient vector）—— 一阶导数
- Hessian 矩阵 —— 二阶导数

> 两者都是泰勒展开近似似然函数所必需的

核心参数就是`--dating mcmctree`, 在此基础上

- 选择不同的进化模式，核苷酸序列对应GTR+G4, "MIX{GTR,HKY}+G4", 氨基酸序列LG+G4+C60， 或者直接用MFP都测试一遍
- 添加基因的分区信息，-Q  ( Like -p but edge-unlinked partition model)
- 额外的mcmc控制参数，如--mcmc-iter 20000,200,50000 --mcmc-bds 1,1,0.5 --mcmc-clock IND


```bash
# NT序列，无分区信息，GTR+GT4模型
iqtree3 -s example.phy     -te example_tree.nwk    --dating mcmctree -m GTR+G4  --prefix example

# NT序列，无分区信息，GTR+GT4模型, 额外控制参数
iqtree3 -s example.phy     -te example_tree.nwk    --dating mcmctree -m GTR+G4  --mcmc-iter 20000,200,50000 --mcmc-bds 1,1,0.5 --mcmc-clock IND

# NT序列,  有分区信息，GTR+GT4模型
iqtree3 -s example.phy     -te example_tree.nwk    --dating mcmctree -m GTR+G4 -Q example.nex

# NT序列，无分区信息， "MIX{GTR,HKY}+G4"模型
iqtree3 -s example.phy     -te example_tree.nwk    –-dating mcmctree -m "MIX{GTR,HKY}+G4"

# AA序列，无分区信息, LG+G4+C60模型
iqtree3 -s example_aa.phy  -te example_aa_tree.nwk –-dating mcmctree -m LG+G4+C60
```