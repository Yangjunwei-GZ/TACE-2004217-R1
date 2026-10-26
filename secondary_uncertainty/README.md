# 新队列次要性能评价入口

本目录只适配本轮新冻结队列的次要统计计算，不含旧 118 例的计算结果，不训练模型，不修改医学处置、主分析输入或已发布 v1.0.0。

`secondary_metrics.py` 与 2026-10-22 次要统计程序逐字节相同；`ORIGINAL_SECONDARY_METHODS_REFERENCE.json` 仅保存原统计定义的历史来源，其 118 例输入哈希、支持范围及结果不作为本轮证据。新方案由 `freeze` 命令根据实际新矩阵生成，写入真实日期和新输入指纹。

## 统计口径

- A：`clinical_diameter_mean_hu`，固定临床变量＋连续最大 3D 直径＋平均 HU。
- B：A＋每个外层训练折内部完成筛选、标准化和调参的 radiomics score。
- 12、24、36 月累计/动态 IPCW AUC 与生存 Brier；训练集 reverse-KM 删失权重；概率来自当前训练拟合的 Cox 基线，不使用风险分数冒充概率。
- 固定 IBS：每月 1–36 月，梯形积分除以 35。A/B、每折、每次重复使用相同窗口；不因新结果调整窗口。
- 点估计使用 point run 0、种子 20261014 的一次分组 5 折：AUC/Brier/IBS 为五折算术均值；校准为该次患者折外预测合并后五分位分组的 KM 风险误差。
- 原 10 次 point 重复只作为分折敏感性均值及 SD。500 次全流程患者 bootstrap 的原种子 20261015–20261514 保留；区间与 point 0 的单次 5 折统计量匹配，不给 10 次均值套同一区间。
- 各指标 A、B 与 B−A 使用共同可评估 bootstrap 集合的 2.5/97.5 百分位。点估计可评估且共同有效次数至少 450/500 才给 CI；其余保留点估计及 NA 原因，不改种子补齐。
- 分组 KM 的 Greenwood log-log CI 反映条件于分组的观察风险不确定性，不代表全模型开发过程；五组绝对误差的非加权与人数加权均值分别在 bootstrap 中计算。
- ΔBrier/ΔIBS 小于 0 表示误差减少。指标 CI 跨 0 不支持宣布改善，也不证明完全无效。

## 当前 checkpoint 要求

执行目录必须已含完整 `04_validation/VALIDATION_INTEGRITY.json`，状态为 `COMPLETE`，以及 10 point＋500 bootstrap 的真实新 checkpoint。适配器再次逐份运行当前 runtime 的完整 checkpoint 检查，核对矩阵/配置/源代码 SHA、原抽样和患者分组折次、所有适用模型及主比较。

新 checkpoint 必须保存 `fitted_baseline`、`monthly_survival`、训练结局及抽样行号。本程序优先直接使用这些实际拟合对象，并在固定系数下用精确训练样本独立重建 Breslow 作数值核验。该核验不重新拟合回归系数。Coxnet 的 `prediction_offset` 单独处理；第一观察时间前 H0=0，最后训练随访之后的预测输出保留缺失。

缺对象的旧 checkpoint 会明确报错，不自动取旧 118 例数据，也不隐式重训。所有程序异常写入 JSONL，返回非零；统计不可评估另存 NA 与原因。输入不完整、重复执行或计算失败不会生成 `STATISTICS_COMPLETE.json`。

## 执行命令（PowerShell）

使用项目锁定的 CPython 3.10.11。以下路径均可通过参数更换，输出必须位于主科学 revision 目录之外。

```powershell
$py = 'python' # Use the locked CPython3.10.11 executable
$entry = 'secondary_uncertainty/code/run_new_secondary.py'
$revision = 'new scientific run/revision_R1_20261014'
$runtime = 'code' # Exact directory used to generate those checkpoints
$out = 'new secondary statistics'

# 新矩阵和本轮 runtime 最终冻结之后，计算次要指标之前执行。
& $py -B -X utf8 $entry freeze --revision $revision --runtime $runtime --output $out

# 仅在真实主验证完整完成后执行。默认执行全体 510 个既定 checkpoint。
& $py -B -X utf8 $entry compute --revision $revision --runtime $runtime --output $out

# 对既有本轮输出重新作只读来源核对与独立数值重算。
& $py -B -X utf8 $entry audit --revision $revision --runtime $runtime --output $out
```

首次 freeze 后不修改本目录 Python 代码；任何代码、矩阵或配置变化均需新输出目录及新冻结记录。compute 保留已有结果，不覆盖运行台账；若中断或失败，先保留日志，另选输出目录重新执行本次次要计算。

## 结果文件

- `results/uncertainty_summary.csv`：13 个主比较次要指标/时间组合，含 A/B/差值、CI、共同有效和 NA 次数，及独立列示的 10 次敏感性统计。
- `results/paired_bootstrap_metrics.csv`、`run_metrics.csv`、`fold_metrics.csv.gz`：重复/折级追溯和配对重算。
- `results/point0_calibration_with_KM_CI.csv`、`calibration_group_details.csv.gz`：主校准图、各组人数/事件/支持范围及观察风险区间。
- `results/patient_OOF_fixed_time_predictions.csv.gz`、`original_train_validation_partitions.csv.gz`、`reconstructed_predictions/`：患者、抽样行、折次和完整月度概率。
- `results/NON_EVALUABLE_RUNS.csv`、`NON_EVALUABLE_REASON_COUNTS.csv`：真实 NA 分母和原因。
- `results/task_ledger.csv`、`logs/TASK_EVENTS.jsonl`、`PROGRAM_FAILURES.jsonl`、`ENTRY_FAILURES.jsonl`：计划/实际尝试/成功/失败证据。
- `audit/NEW_RUN_FOLLOWUP_SUPPORT.csv`、`source_binding/CURRENT_SCIENTIFIC_BINDINGS.json`、`SECONDARY_INDEPENDENT_RECOMPUTATION.json`：新队列支持审计、原 checkpoint 字节哈希及独立摘要复算。
- `results/COMPUTATION_STATUS.json` 与 `STATISTICS_COMPLETE.json`：前者为计算状态，后者只在全部执行完整且自动数值审计通过后生成。该标记不等于临床增量价值成立，也不替代稿件独立复核。

## 测试

```powershell
& $py -B -X utf8 -m unittest discover -s (Split-Path $entry) -p 'test_*.py' -v
```

合成测试与真实科学计算分开。测试覆盖 IPCW 与 scikit-survival 的支持域内一致性、KM 与 lifelines 的一致性、概率方向/单调性、Coxnet offset、保存基线与概率、原抽样/分组折次、共同可评估百分位区间、输入指纹漂移、499/500 失败和缺失/重复任务门禁。
