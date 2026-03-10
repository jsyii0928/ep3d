# Gpulse与EPOCH3D融合重构：当前进度与问题总结

## 1. 当前工作进度 (Progress Status)

截至目前，已基本完成 **"阶段 1：程序的精简、移植与编译适配 (Porting)"** 以及部分**"阶段 2"**。具体工作包括：

- **文件移植与拆分**：将 `fpulse` 目录下的核心逻辑移动到 `epoch3d/src/` 中。为适配 Fortran 的编译流与模块依赖，将合在一起的代码成功拆分为两个文件：
  - `epoch3d/src/gpulse_constants.f90` (独立存放物理常数和辅助计算结构)
  - `epoch3d/src/gpulse.f90` (包含主电磁场偏微分和积分求解器)
- **构建系统 (Makefile) 改造**：在 `epoch3d/Makefile` 的 `SOURCES` 中登记了上述两个新文件，并指定了正确的编译依赖链：
  - `gpulse_constants.o` $\rightarrow$ `gpulse.o` $\rightarrow$ `laser.o`
- **处理命名空间冲突**：Fpulse 原本的 `constants` 模块名称与 EPOCH 内部自带的 `constants` (通过 `shared_data` 或底层自带) 发生同名冲突。已安全地将新增的模块重命名为 `gpulse_constants` 以隔离命名空间。
- **Fortran 语言标准向后兼容降级**：修复了 Fpulse 原本使用了 Fortran 2008 标准的动态单元号 `open(newunit=iu...)` 语法特性。由于 EPOCH 整体强制要求使用 `-std=f2003` 进行编译，导致编译器报错。目前已将其降级处理为手动指定固定的静态通道号（`iu = 283`以及 `unit=283`）。

---

## 2. 遇到的主要问题与解决手段 (Issues Encountered & Solutions)

### 问题 A：模块同名（Namespace Clash）导致编译挂起
- **现象**：Fpulse 在自身的体系内定义了 `module constants`。当引入到 EPOCH 并用 `make` 编译时，发现 `USE constants` 引发了严重的模块读取解析错误，破坏了后续如 `shared_data` 甚至原版 `timer.f90` 的编译。
- **解决**：全局替换新增的文件模块名，从 `constants` 变为 `gpulse_constants`，并在 `gpulse.f90` 内做了相应的 `USE gpulse_constants` 追踪替换。

### 问题 B：Fortran 标准降级导致的语法不支持（F2008 vs F2003）
- **现象**：`gpulse.f90` 使用了 `OPEN(newunit=iu, ...)` 动态分配文件通道的功能，但 EPOCH 的 `Makefile` 里默认写死了 `-std=f2003`，编译器提示：`Fortran 2008: NEWUNIT specifier at (1)`。
- **解决**：回退到了传统的方法，声明了一个显式的整数并由我们自行控制可用单元口进行关联打开：`iu = 283`，`OPEN(unit=283, ...)`。

### 问题 C：变量声明与多次定义
- **现象**：在通过脚本修改 `iu=283` 时，一度出现 `Symbol 'iu' already has basic type of INTEGER`。
- **解决**：这是因为原本已在头部声明过 `integer :: iu`，需将变量声明与赋值分离，修改已完成。

---

## 3. 下一步迭代计划 (Next Steps)

一旦上一步的修复完全通过编译，接下来将马上进行后续的核心步骤：

1. **彻底应用双精度类型对齐 (精度安全)**：
   取代由于拆分模块未覆盖彻底的宏变量定义。必须确认 `gpulse_constants` 中的精度基准（`dp` 以及 `dpc`）使用 EPOCH 官方内置的 `num`（REAL(8)）进行接管。
2. **挂载生命周期 (阶段 2/3)**：
   在 `laser.f90` 的 setup 阶段调用 `gpulse_init()` 初始化，并绑定卸载函数 `gpulse_free()`。检查其内存及缓存（ `.gpulse_cache.bin` ）在集群 MPI 多进程架构下的防冲突逻辑。
3. **完成核心网格点的三维场注入替换 (阶段 5)**：
   进入 `outflow_bcs_x_min` 极化处理源码。通过 `CALL gpulse(...)`获取结果代替原本的高斯和平面包络解。映射正确的交错网络偏移坐标（`0.5*dx` 等步长错位）。
