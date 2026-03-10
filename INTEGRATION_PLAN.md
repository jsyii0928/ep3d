# Gpulse 与 EPOCH3D (laser.f90) 深度融合与重构计划

本文档记录了将用于计算复杂入射激光电磁场的 `gpulse` 模块，深度融合并替换 EPOCH 3D PIC 模拟程序中 `laser.f90` 内部边界条件的详细实施架构。主要目标是实现具有时空拓扑特性（如纵向场、相干复合）的三维大孔径精确光场的注入。

---

## 阶段 1：程序的精简、移植与编译适配 (Porting)

1. **核心文件剥离与转移**
   - 将 `fpulse/gpulse_mod.f90`（包含电磁场偏微分、积分与插值求解的核心模块）以及相关包含物理常数的模块直接**拷贝或移动**到 `epoch3d/src/` 文件夹下，重命名为标准的 EPOCH 约定命名（如 `gpulse.f90`）。
   - 剔除 `fpulse` 文件夹中原本用于独立测试、可视化的辅助代码（如 Bash 构建脚本、Python 画图脚本等），保持 EPOCH 工程干净。
2. **系统构建 Makefile 改造**
   - 在 `epoch3d/Makefile` 的 `SOURCES`（源文件列表）中添加新加入的 `gpulse.f90`。
   - 增加编译依赖关系：确保 `gpulse.o` 先于 `laser.o` 编译（即在 Makefile 下方添加 `laser.o: gpulse.o` 的依赖说明），同时在 `laser.f90` 开头引入 `USE gpulse_mod`。

---

## 阶段 2：数据类型对齐与内存生命周期管理

1. **统一双精度浮点类型 (Precision Match)**
   - EPOCH 使用预定义的 `num` 作为所有实数的精度基准（在 `shared_data` 模块中通常等于 `REAL(KIND=8)`）。
   - 全局替换 `gpulse` 中的自有双精度定义（如 `dp`, `gp_dp` 等）直接对齐到 EPOCH 的 `num`。避免在频繁的源项内外传递中引发类型强制转换甚至精度丢失。
2. **内存分配与销毁的同步挂载**
   - 将 `gpulse` 的启动（缓存读取、网格初始化等 `gpulse_init`）挂载到 EPOCH 的激光初始化过程（如 `setup_laser_phases` 或者模拟最前期的 setup 中）。
   - 将释放过程 (`gpulse_free`) 挂载到 `laser.f90` 的 `deallocate_lasers` 子程序中。
   - 需要特别注意 MPI 并行下的 IO：确保如果 `gpulse` 需要拉取大型缓存文件，仅通过 Master 节点读取然后广播（Broadcast），或者让各 MPI 节点独立读取，避免并发 IO 冲突。

---

## 阶段 3：物理单位制的自洽统一 (Unit System)

1. **维持 EPOCH 的 SI 国际单位制设定**
   - EPOCH 宏观使用标准的国际单位制（SI）：电场 $E$ 为 `V/m`，磁场 $B$ 为 纯 `Tesla (T)`，距离为 `m`，时间为 `s`。
   - 核对 `gpulse` 内部和输出层：确保从 `gpulse` 中得出的 `Ex, Ey, Ez, Bx, By, Bz` 严格且绝对为国际单位制下真实的电磁场强度值。如果有使用归一化单位（如高斯制或无量纲的 $a_0 = eE/m\omega c$），必须在传递回 EPOCH `laser.f90` 前完成换算因子的乘入。

---

## 阶段 4：坐标系统体系与 Yee-Grid 的精确映射

由于 EPOCH 基于 MPI 并行以及交错网格（Yee-Cell），边界处的 `(i,j,k)` 并不是空间绝对坐标。必须向 `gpulse` 传入**全局的精确物理坐标**。

1. **绝对宏观坐标的演算**
   - 使用定义在 `shared_data` 或 `mpi_routines` 中的局部框起点 `x_min, y_min, z_min`，以及网格步长 `dx, dy, dz`。
   - 对于网格索引 `i, j, k`，其实际坐标应为：`phys_y = y_min + REAL(i, num) * dy`，以此类推。
2. **兼顾交错网格的半步错位 (Yee-Grid Staggering)**
   - 磁场网格面心和电场网格边缘在时间、空间上天然相差半格。求值时要对宏观坐标进行偏移修正。
   - 例如要获取 S方向对齐的 $E_y$，其评估位置可能在 `(x, y+0.5*dy, z)`；
   - 而针对同时刻的纵向磁场 $B_x$，评估位置通常在面心 `(x, y+0.5*dy, z+0.5*dz)`。
   - 对于边界条件中时间半步错开的特征，传入 `gpulse` 的时间参数需要评估是否应该使用 `time + 0.5 * dt`。

---

## 阶段 5：3D电磁场全分量(纵波+横波)注入实施

原 EPOCH 激光注入仅简化出两个横向包络：`source1` 和 `source2`。但对于强聚焦的 `gpulse` 理论预言的光场，带有显著的纵向分量，漏掉它们会导致 $\nabla \cdot \mathbf{E} \neq 0$ 以及 $\nabla \cdot \mathbf{B} \neq 0$ 从而激发假震荡。

1. **替换原生高斯/平顶横向场生成**
   - 在 `outflow_bcs_x_min` 等处理边界条件循环 `DO j... DO i...` 中。彻底废弃其基于 `t_env`，`profile` 函数的求解方法。
   - 使用宏观物理量 `(phys_x, phys_y, phys_z, time)` 作为入参调用 `CALL gpulse(...)` 获得当地时空的数值结果解析解 `(gp_Ex, gp_Ey, gp_Ez, gp_Bx, gp_By, gp_Bz)`。
   - 以 $x\_min$ 边界为例：将 `source1` 的赋值改写为 `gp_Ey` 关联的值，`source2` 赋为 `gp_Ez` 关联的值（注意匹配原代码中的 S / P 偏振方向定义）。
2. **平行传播方向的纵向磁场叠加**
   - 原程序对顺传播方向场是假定为零的。
   - 以 $x\_min$ 为例：注入时增加一行关于平行磁场同步的修正代码：
     ```fortran
     ! 比如原来是 bx(...) = bx_x_min(...)
     ! 修改为叠加 gpulse 输出的对应面上纵向磁场
     bx(laserpos-1, 0:ny, 0:nz) = bx_x_min(0:ny, 0:nz) + b_pulse_x_matrix
     ```
   - 同理在 $y$ 和 $z$ 入射边界也要完成类似的 `by(...) = by_y_min(...) + b_y_pulse` 缝合。 

---

## 最终执行范例模式框架 (Pseudocode Example)

以替换 `outflow_bcs_x_min` 为核心参照点：

```fortran
! 在 laser.f90 : outflow_bcs_x_min 内部

IF (add_laser(n)) THEN
  ! ... existing variables ...
  REAL(num) :: p_x, p_y, p_z, p_t
  REAL(num) :: gp_Ex, gp_Ey, gp_Ez, gp_Bx, gp_By, gp_Bz

  current => lasers
  DO WHILE(ASSOCIATED(current))
    IF (current%boundary == c_bd_x_min) THEN
      IF (time >= current%t_start .AND. time <= current%t_end) THEN
        
        p_x = x_min ! 针对当前节点的实际 x_min 值
        p_t = time

        DO j = 0, nz
          ! 注意此处理论上应当加上对应场位置的偏移量，如 + 0.5*dz
          p_z = z_min + REAL(j, num) * dz 
          DO i = 0, ny
            p_y = y_min + REAL(i, num) * dy

            ! 从预计算好的库或者理论表达式中抽出真实解
            CALL gpulse(p_t, p_x, p_y, p_z, gp_Ex, gp_Ey, gp_Ez, gp_Bx, gp_By, gp_Bz)

            ! 替代横波源场 (S、P极化映射)
            source1(i,j) = gp_Ey
            source2(i,j) = gp_Ez
            
            ! 保存算出的纵波，如 Bx_pulse_matrix(i,j) = gp_Bx
          END DO
        END DO
      END IF
    END IF
    current => current%next
  END DO
END IF

! 将由于非平面波特性多出来的散度场注入边界磁场中
! bx(laserpos-1, 0:ny, 0:nz) = bx_x_min(0:ny, 0:nz) + Bx_pulse_matrix 
```
