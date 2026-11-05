# HCF · P3 工具验证 — 黑盒信效度复现包

本仓库伴随论文 **《Black-Box Reliability and Validity of a Proprietary Computational Assessment System in Two Cohorts》** 发布，用于复现公开范围内的黑盒信效度与下游统计结果。

> **公开边界**：仓库提供九个后台直出读数（四个 α、四个 carrier、一个 H）、外部变量、公开量表结果及下游统计。书籍已公开的 H_trend = [α1]² 仅作为 α1 的确定性下游变换复现，不是第十个后台读数，也不表示实测纵向变化。内部题目作答、内部中间量、H 的后台生成公式、映射规则、校准参数与专有计算组件均不公开；公开代码不会从原始作答或其他读数重构内部计分。

## 1. 版本标识

| 对象 | 对外名称 |
|---|---|
| 第一期问卷 | 科研-第一套题 |
| 第二期新题 | 科研-第二套题 |
| 两期算法版本 | 算法v1.1 |

版本字段仅用于标识数据口径，不包含内部规则。详见 `docs/版本三口径说明.md`。

## 2. 数据发布层级

| 层 | 内容 | 获取方式 |
|---|---|---|
| **L1 公开代码** | 本仓库中的 R、Python、文档与聚合结果 | 公开 |
| **L2 去标识分析数据** | 九个后台直出读数、外部健康结局、公开量表结果、人口学与版本标识 | 本仓库公开 |
| **L3 受控输入** | 受控原始作答、配对桥接与题目级核验材料 | 仅限受控审稿核验 |
| **L4 内部组件** | 专有计算实现与内部规则 | 不公开 |

公开 CSV 不含内部逐题字段、派生分类或内部派生计数。`gad7_item_count`、`ibs_symptom_count`、`sleep_symptom_count` 均为外部变量，不属于内部计分派生。

## 3. 目录结构

```text
03_GitHub代码包/
├── code/
│   ├── R/                         # 双队列 R 下游统计
│   ├── python/                    # Python 下游统计、E2/E3 与公开数据导出
│   └── data/
│       ├── derived/               # 公开 L2 基底
│       └── semantic/              # 编号相似度矩阵与聚合统计
├── data_output/
│   ├── analysis_results/          # 公开 Python 聚合结果
│   ├── E2_results/                # E2 聚合结果
│   ├── E3b_results/               # E3a/E3b 聚合结果与公开重复基底
│   ├── a42_semantic_confound/     # 语义聚合结果
│   └── llm_e1/                    # LLM 聚合结果
└── docs/                          # 公开边界、版本与偏离说明
```

## 4. 复现指南

### 4.1 公开数据核验与清理

```bash
python code/python/verify_phase2.py
python code/python/sanitize_public_data.py
```

`sanitize_public_data.py` 以白名单原地整理本仓库的六份公开派生基底及其公开重复副本，不读取受控原始数据。六份基底为 `p3_base_prolific.csv`、`p3_base_jianshu.csv`、`p3_base_gold.csv`、`p3_base_phase2_long.csv`、`p3_base_phase2_wide.csv` 和 `p3_base_e2_retest.csv`。`export_public_l2.py` 用于从受控权威基底重新导出这些公开文件；只有明确列入白名单的字段能够输出。E2 文件仅含164个匿名配对、两波五个直接读数及重测间隔，不含原始作答、平台标识或配对桥接。第二期受控 `h3_present` 仅改名为公开 `H_score`，原始平台标识只在内存中转换为匿名 token，宽表与长表保持同一参与者对应关系：

```bash
python code/python/export_public_l2.py \
  --source-root /path/to/controlled/project \
  --out-dir /path/to/public/code/data/derived
```

导出器拒绝把受控 `代码/data/derived` 或其子目录作为输出目录；固定校验匿名 token、第二期长宽表和 Gold 与全样本的全部共享字段。所有检查及临时文件回读通过后才原子替换目标 CSV，写入失败会回滚。公开清理脚本采用相同的分阶段写出机制。

### 4.2 Python 下游统计

```bash
python code/python/run_all.py
python code/python/e03b_phase2.py
```

`run_all.py` 使用固定的公开模块白名单，显式纳入 `e03b_phase2.py` 与 `e2_icc_dual_20261031.py`，不通过目录自动发现脚本；因此本机保留的受控或历史模块不会误入公开回归。白名单文件缺失会立即报错。运行清单区分 `success`、`partial_public_l2`、`skipped_controlled_l3_required` 与 `failed`：部分完成只表示列明的公开子分析已运行，不能解释为该模块的全部历史/受控子分析均已完成。

公开 Python 模块共 27 个：`a01`、`a02`、`a03`、`a04`、`a05`、`a06`、`a07`、`a09`、`a10`、`a11`、`a12`、`a15`、`a16`、`a18`、`a19`、`a24`、`a25`、`a27`、`a32`、`a36`、`a37`、`a38`、`a39`、`a41`、`a42`、`e03b`、`e2_icc_dual_20261031`。

隔离运行时，将 `P3_RESULT_DIR` 设为独立目录；模块结果写入该目录，运行清单位于其上一级。E3 使用 `P3_E3B_OUT`，E2 使用 `P3_E2_OUT` 指定隔离输出目录。E2 默认读取公开匿名配对基底并复算 ICC(3,1)、ICC(2,1)、Pearson、Spearman 和均值差；只有重新执行受控数据血缘核验时才需要 `P3_CONTROLLED_ROOT`。

### 4.3 R 下游统计

```bash
Rscript code/R/10_analysis_run_all.R
```

使用 `P3_OUT_DIR` 指定隔离输出目录。公开 R 总入口包含 28 个分析模块：`01`、`02`、`03`、`04`、`05`、`06`、`07`、`09`、`11`、`12`、`14`、`15`、`16`、`17`、`18`、`20`、`21`、`22`、`23`、`24`、`25`、`26`、`27`、`32`、`34`、`35`、`36`、`43`。需要受控题目级或数值锚点输入的整模块或子分支在公开 L2 下显式跳过，不生成占位结果；清单同样区分完整成功、部分完成、受控跳过和失败。

### 4.4 语义与 LLM 聚合结果

- `a42_semantic_confound.py` 只读取 `code/data/semantic/` 中的编号矩阵与聚合统计。
- `llm_e1_deepseek_redo.py` 仅执行公开输出级下游分析：读取 `code/data/derived/p3_llm_predictions.csv` 中的匿名 token 与八项模型数值预测，联结公开 L2 后复算五张 E1 聚合表。公开脚本不包含个案文本构造、原始作答读取、平台 ID 桥接、API 调用或原始模型返回；这些受控组件不进入公开仓库。
- 运行 `python code/python/llm_e1_deepseek_redo.py` 可复现 E1 五表。`p3_llm_predictions.csv` 不含题干、选项、个案文本、原始返回、response ID、输入哈希或平台 ID；它不是原始作答数据，也不用于重构后台读数。
- `verify_scoring.py` 在公开 L2 下只核验总分范围、公开衍生字段与二值字段，并将原始题项计分重算明确标为受控跳过；只有设置 `P3_CONTROLLED_ROOT` 且受控 L3 输入完整时才逐题重算 PHQ-9/GAD-7 与 Gold 外部量表。配置列出或总分范围通过不等于公开包已独立验证官方计分。`verify_phase0_statistics.py` 可直接执行公开模式：独立复算 C1/E3a、a41五表、a04与a25，E2原始桥接和a27逐题核验明确记为受控跳过；设置 `P3_CONTROLLED_ROOT` 后才执行全部受控步骤，失败时同样保留机读报告。受控核验可另设 `P0_REPORT_DIR` 与 `P3_RESULT_DIR`，将含受控条目级结果的报告写到公开仓库之外。

## 5. 关键结果入口

| 内容 | 脚本 | 公开结果 |
|---|---|---|
| E2 重测 ICC | `e2_icc_dual_20261031.py`（默认读取 `code/data/derived/p3_base_e2_retest.csv` 的164个匿名配对；直接分析最终读数） | `data_output/E2_results/E2_icc_dual_20261031.csv` |
| E3a 配对与顺序效应 | `e03b_phase2.py` | `data_output/E3b_results/` |
| E3b 注册 C1/C2/C3 | `e03b_phase2.py` / R `10_analysis_43_e3b_registered.R` | `C1_预注册标准化主结果.csv`、`C2_方向复现.csv`、`C3_选项级分布.csv`；C1 为第2批240对原中文队列504，abs(Cliff's δ)<.147 |
| H_trend 输出级内部一致性 | R `10_analysis_25_replication_extension.R` / Python `a25_replication_extension.py` | `data_output/analysis_results/a25_replication_extension/T25_htrend_replication.csv` |
| 收敛效度 | R `10_analysis_12_convergent_semantic.R` / Python `a12_convergent.py` | 聚合统计表 |
| 筛查与增量效度 | 保留的 R/Python 下游统计模块 | 聚合统计表 |

E3b 的 C1 使用 Brunner–Munzel 中秩推断，α3、α4 通过，α1、α2 未通过；C2 仅作 α3–IBS 与 α3–PHQ-9 的方向描述；C3 公开的是不含题目文本、选项文本或映射的100格聚合计数。第1批对第2批均值尺度分析明确标为事后分析，不属于注册 C1/C2。

## 6. 依赖与验收

- R：tidyverse、psych、lme4、boot、MASS 等，具体见各脚本。
- Python：pandas、numpy、scipy、statsmodels、scikit-learn。
- 清理后的验收以隔离目录中的语法检查、公开回归日志与机读清单为准；历史全模块验收数字不代表当前公开集合。

## 7. 伦理与预注册

- 伦理：南京师范大学生物医学伦理委员会，批件号 **NNU202605041**。
- 预注册：OSF **d4e7p**（V2.0，2026-05-04）；E3b 另有独立预注册。
- 偏离摘要见 `docs/偏离表摘要.md`，伦理与公开边界见 `docs/伦理声明.md`。
