# Local and Historical Resources

## 摘要描述

本文件登记重要的本地代码、历史note、半成品文章和工具指南。具体读取范围和授权状态按条目记录；T004 只扩展了 `C01` 的 W02 实现范围和 `H01` 的数值演化附录。历史资源均不保证正确，不能直接作为可信物理依据。

## 具体内容

## C01

### 2026-09-08 · T009 授权补充

Owner 在 `TASK.md` 任务2明确授权主 agent 参考 `Codes/errorAnalyse` 与
`Codes/mma-expansion` 框架，同时禁止 coding 子 agent 阅读参考代码文件夹。
本次主 agent 已读相关分析说明、轨道频率定义、angular/complex-MST 接口片段；
仅用于形成新计划，没有执行历史代码或将其结果迁入本次证据。
下方 T004 记录为历史范围，不覆盖这项新的限定授权。

- 标题：Direct_Wave_in_Higher_order/Codes
- 类型：code collection
- 路径：`/Users/indigo/codes/Direct_Wave_in_Higher_order/Codes`
- 用途：既有direct-wave研究代码的历史入口
- 访问条件：未经仓库所有者明确允许，不得列目录、搜索、打开或读取其细节；T004 已明确授权读取 `NewPointParticlesm2Evolution` 及其直接引用的 plunge 初值实现
- 状态：T004 已审计上述限定范围内的当前主线推导、solver、tests、初值生成和保留结果，并区分目标仓库自身的 legacy outputs；结论是以当前主线修正迁移，内容尚未获本仓库 owner 批准；其余代码仍未读取

## H01

- 标题：Direct_Wave_in_Higher_order/Paper/main.tex
- 类型：manuscript draft
- 路径：`/Users/indigo/codes/Direct_Wave_in_Higher_order/Paper/main.tex`
- 用途：历史文章结构与论述记录
- 访问条件：仅在任务需要追溯历史写作时读取；T004 已授权审计数值演化附录
- 状态：T004 已读取相关附录；坐标块基本可用，但有源方程的 tetrad/signature 衔接、$\Psi_4$--strain 表述、horizon boundary-term 理由和若干记号存在缺口；半成品，未批准

## H02

- 标题：ExciseOfDW.md
- 类型：note
- 路径：`/Users/indigo/codes/Direct_Wave_in_Higher_order/Paper/ExciseOfDW.md`
- 用途：历史想法或推导记录
- 访问条件：仅在任务明确需要相关历史背景时读取
- 状态：本次未读取；不保证正确

## H03

- 标题：direct_wave_quadratic.tex
- 类型：manuscript draft
- 路径：`/Users/indigo/codes/Nonlinear_direct_wave/paper/manuscript/latex_brand/direct_wave_quadratic.tex`
- 当前任务参考PDF：`/Users/indigo/codes/Nonlinear_direct_wave/paper/manuscript/latex_brand/direct_wave_quadratic.pdf`
- 用途：非线性direct-wave工作的历史文章记录
- 访问条件：仅在任务需要追溯该方向历史内容时读取
- 状态：T010 经owner提供PDF核对第4页式(29)/(30)、(39)/(40)，后两式明确
  Lu的p推广模板；未读取tex或验证整篇标量非线性推导。具体版本hash及使用范围
  见 [W08任务3计划](../workspace/plan/W08_end_to_end_error/task3_plan.md)。

## H04

- 标题：Groupmenting824.md
- 类型：meeting note
- 路径：`/Users/indigo/codes/Direct_Wave_in_Higher_order/Paper/notes/Groupmenting824.md`
- 用途：组会讨论的历史记录
- 访问条件：仅在任务需要追溯相关讨论时读取
- 状态：本次未读取；不保证正确

## N01

- 标题：Mathematica MCP usage guide for this repository
- 类型：note
- 路径：`/Users/indigo/codes/Direct_Wave_in_Higher_order/Reference/mathematica_mcp_ai_guide.md`
- 用途：MCP自包含调用、紧凑输出、session状态、超时拆分、失败诊断和验证记录
- 访问条件：首次构造MMA MCP检查，或MCP连接、Kernel、环境、超时和状态问题再次出现时
- 状态：本次已核对全文；指南本身声明未实际调用MCP，其中配置、版本和环境事实需在使用前复查
- 适用边界：可复用工程经验；原仓库专用物理假设、固定参数和强制路由不能无条件移植

## 总结及下一步预告

真正可信的内容必须来自后续正文核对、明确推导或实际执行检查。历史材料中的结论在完成这些工作前均视为未验证。
