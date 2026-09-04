# MIAO 设计说明

本文记录当前 beta 版本中与 intORF 分析有关的有意设计取舍。它们不是待清理的测试功能。

## 1. 同一 stop 区段只保留最早起点

当前流程没有 TI-seq 等起始位点证据，因此不能可靠区分同一 reading frame、同一终止密码子上游的多个嵌套起点。扫描器在每个 stop-delimited segment 中保留最早起点，即最长的 start-to-stop ORF，减少高度重叠候选和多重检验负担。更下游起点可能才是真实起点，这是当前 beta 版本的定位限制。

## 2. 与注释 CDS 共享 stop 时使用完整 CDS

当扫描 ORF 与注释 CDS 共享终止位置时，当前版本使用完整注释 CDS 范围并标记为 `CDS`。CDS 在本项目中主要用于定义 intORF 的相对位置、重叠范围和 reading frame，并非主要研究对象。

## 3. 负链处理

规范的负链转录本已经在坐标映射和序列读取两处配套处理：外显子顺序按转录方向组织，序列按需要反向互补。因此当前 GENCODE `+/-` 注释下的负链结果是正常的。

## 4. gORF 验证范围

`validate_gorf_outputs()` 验证的是下游 intORF 分析所依赖的核心 tORF-gORF 映射，包括 ID 唯一性、每个 tORF 只属于一个 gORF、membership 映射一致以及 tORF 集合一致。它不是覆盖 TSV、FASTA 和所有汇总字段的完整一致性审计器。

## 5. intORF_altframe 与其他注释 CDS 共享 stop

候选首先按照所在转录本 CDS 判定 F1/F2。仅对已经判定为 `intORF_altframe` 的 ORF，再将完整三核苷酸 genomic stop signature 与 GENCODE 同基因及同链注释 CDS stop 比较。原始 F1/F2 分类不被覆盖，但同 stop 的 annotation-confounded 标记会进入 tORF/gORF 输出；DM caller 默认排除这类候选。

## 6. FASTA 句柄复用

扫描器按进程复用 `pyfaidx.Fasta` 句柄：单进程只打开一次，多进程由 pool initializer 为每个 worker 打开一次。该优化不改变序列、坐标、ORF 分类或输出内容。

## 7. 正式 DM 推断

正式 caller 只使用 host-only 与 host-plus-target 的 Dirichlet-multinomial mixture 检验产生统计 p 值。此前探索的 gene-shrunk noise、global-noise、联合 uniform/stop p 和 calibrated-LRT 没有进入正式模型，相关试验脚本已经删除。

最终可信判定在 BH 后同时应用预先固定的诊断门槛，包括：

- active-core codon 数量和比例；
- target residual breadth；
- lambda 下限与两方向 lambda 一致性；
- mixture segment distance；
- template separation。

这些 gate 不替代 DM p 值，也不重新定义 FDR family。

尾部概率计算提供两个正式模式：

- `accurate`：固定 9,999 次 importance sample，默认模式，适合正式结果与排序；
- `fast`：pilot 只决定独立 confirmation sample 的规模，适合常规筛查和初步运行。

两个模式使用完全相同的生物模型、候选集合、LRT、FDR 和 gates，只改变 Monte Carlo 资源分配。
