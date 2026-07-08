# V_plot_pipeline (ATAC 为例)

从 ATAC 片段 bed 一路到"区分真 TF 足迹 vs 酶切偏好"的二维阈值 + 散点图。

## 流程

| 步 | 脚本 | 输入 → 输出 |
|----|------|------------|
| 0 | `0_prepare_fragments.sh` | 原始片段(3列 chr,start,end) → 中点+长度(6列) |
| 1 | `1_locate_bias.sh` | bias k-mer 列表 + 基因组fa → 各 k-mer 坐标 bed (seqkit) |
| 2 | `run_TF_scatter.sh` | TF motif 目录 + 片段 → 每个 TF 的 V-plot 距离分布 |
| 3 | `run_bias_scatter.sh` | bias 位点 + 片段 → V-plot 距离分布 (抽样/排除指定TF) |
| 4 | `fit_vplot_apex.py` | 距离分布 → 生成模型拟合: apex、V-channel 宽度、V内富集、LRT |
| 5 | `5_scatter_cutoffs.py` | TF/bias apex → Wilcoxon rank-sum 检验；显著时再输出 Youden 二维阈值 + LOO + 散点图 |

## 依赖
- `bedtools` (closestBed / intersectBed)、`seqkit`
- `python3` + numpy / pandas / scipy / matplotlib

## 用法

参数可**直接命令行传入**，也可把 `config.example.sh` 复制为 `config.sh` 后填写路径；
命令行参数覆盖 config，二者皆可省略其一。

**A. 命令行 (推荐)**
```bash
bash run_pipeline.sh \
  --genome hg38.fa \
  --fragment my_ATAC.bed \
  --tf-dir   motif_dir \
  --bias-kmers my_bias.txt \
  --exclude-tf "motif_dir/CTCF.bed motif_dir/NFIC.bed" \
  --mode ATAC \
  --out  results_myATAC
```
建议挂后台：`nohup bash run_pipeline.sh ... > run.log 2>&1 &`

**B. 用 config.sh**：`cp config.example.sh config.sh`，编辑同名变量后直接 `bash run_pipeline.sh`。

常用流程控制 (命令行)：
```bash
bash run_pipeline.sh ... -f 2        # 断点续跑: 从步骤2开始(复用已完成的步骤0/1)
bash run_pipeline.sh ... -f 5        # 只重画图(读已有 apex tsv)
bash run_pipeline.sh ... -k 0        # 跑完删除中间文件
bash run_pipeline.sh ... --shuf-n 0          # bias 不抽样, 用全部位点 (默认抽样20万)
bash run_pipeline.sh ... --shuf-seed none   # bias 抽样用真随机(默认固定种子42, 结果可复现)
bash run_pipeline.sh ... --rank-alpha 0.05   # step 5: 两个特征均需 Wilcoxon rank-sum p < alpha 才划阈值
bash run_pipeline.sh -h              # 全部可用参数
```

拟合参数(步骤4)用 `--mode` 选预设：
- `loMNase` / `DNase`：apex_x∈[-10,10], apex_y∈[0,60], frag 20–100, x-window 150, max-n 0, 不置换；
- `ATAC`：同上但 frag 20–**150**；
- `custom`：用命令行(或 config)的 `--apex-x-lo/hi`、`--apex-y-lo/hi`、`--frag-min/max`、
  `--x-window`、`--max-n`、`--permutations`、`--perm-n`；**任一项不给则回退到 ATAC 预设**。

- 输出目录 `OUTDIR` 没定义会自建；**中间文件全部放在 `OUTDIR/intermediate/`**
  (片段中点 bed、bias 坐标、各 V-plot 距离分布等)。`-k 0` 在跑完后删除该子目录，
  最终结果 (apex tsv / 阈值 csv / 散点图) 保留在 `OUTDIR/`。
- bias 与 motif 均带链信息，排除真 TF 位点时**按同链进行**。

## 输入格式
- **TF motif** (`TF_MOTIF_DIR/`)：每个文件须为 6 列 `chr start end name score strand`(strand=+/-/.)；
  不限后缀，非 6 列格式会直接报错。
- **bias k-mer 列表**：每行一条序列 (如 `GCC`)。
- **片段 bed**：3 列 `chr start end`。

## 输出 (`results/`)
- `<ASSAY>_TF_apex.tsv` / `<ASSAY>_bias_apex.tsv`：每个 motif 的 apex 打分。
- `<ASSAY>_methodB_mw_test.csv`：两特征的 Wilcoxon rank-sum 检验、AUC 和是否通过 `alpha`；文件名中的 `mw` 为兼容旧版本保留。
- `<ASSAY>_methodB_conclusion.txt`：是否通过 rank-sum gate；未通过时说明不划分阈值。
- `<ASSAY>_methodB_cutoffs.csv`：若两个特征均显著，输出两特征阈值 + margin + AUC + LOO，及二维联合规则的 LOO；未通过时写入 `no_cutoff_MW_gate_failed`，不含阈值。
- `<ASSAY>_methodB_scatter.png/.pdf`：width × E 散点；通过 rank-sum gate 时显示 cut-off 虚线和右上达标区。

## 判别规则
先对 V-channel 宽度和 V内富集 log2(V-in/V-out) 分别做双侧 Wilcoxon rank-sum 检验。只有两个特征都满足
`p < RANK_ALPHA` 时才继续划分阈值。通过后，判为 **TF-footprint-like** ⇔ V-channel 宽度 ≥ 阈值<sub>w</sub>
**且** V内富集 ≥ 阈值<sub>E</sub> (散点图右上阴影区)。阈值取 Youden threshold；当 bias 与 TF 完全可分时，
它等于两类分布间 gap 的中点。LOO 留一交叉验证给出 sens/spec/acc。
