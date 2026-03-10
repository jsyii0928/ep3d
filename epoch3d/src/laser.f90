! Copyright (C) 2009-2019 University of Warwick
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE laser

  USE custom_laser
  USE evaluator

  IMPLICIT NONE

CONTAINS

  ! 初始化单个激光对象的属性
  ! boundary: 入射边界索引
  ! laser: 待初始化的激光块结构体
  SUBROUTINE init_laser(boundary, laser)

    INTEGER, INTENT(IN) :: boundary
    TYPE(laser_block), INTENT(INOUT) :: laser

    laser%boundary = boundary
    laser%id = -1
    laser%use_time_function = .FALSE.
    laser%use_phase_function = .FALSE.
    laser%use_profile_function = .FALSE.
    laser%use_omega_function = .FALSE.
    laser%amp = -1.0_num
    laser%omega = -1.0_num
    laser%pol_angle = 0.0_num
    laser%t_start = 0.0_num
    laser%t_end = t_end
    ! 记录激光随时间累积的相位积分 (用于处理变频激光)
    laser%current_integral_phase = 0.0_num
    NULLIFY(laser%profile)
    NULLIFY(laser%phase)
    NULLIFY(laser%next)

    ! 根据边界类型分配空间分布和相位分布的数组空间
    CALL allocate_with_boundary(laser%profile, boundary)
    CALL allocate_with_boundary(laser%phase, boundary)
    ! 默认分布为均匀分布 (profile=1.0)，初始相位为 0
    laser%profile = 1.0_num
    laser%phase = 0.0_num

  END SUBROUTINE init_laser



  ! 设置所有激光的初始相位积分
  ! 通常在从重启点恢复模拟时使用
  SUBROUTINE setup_laser_phases(phases)

    REAL(num), DIMENSION(:), INTENT(IN) :: phases
    TYPE(laser_block), POINTER :: laser
    INTEGER :: ilas

    ilas = 1
    laser => lasers
    DO WHILE(ASSOCIATED(laser))
      laser%current_integral_phase = phases(ilas)
      ilas = ilas + 1
      laser => laser%next
    END DO

  END SUBROUTINE setup_laser_phases



  ! 释放单个激光对象的内存
  SUBROUTINE deallocate_laser(laser)

    TYPE(laser_block), POINTER :: laser

    IF (ASSOCIATED(laser%profile)) DEALLOCATE(laser%profile)
    IF (ASSOCIATED(laser%phase)) DEALLOCATE(laser%phase)
    IF (laser%use_profile_function) &
        CALL deallocate_stack(laser%profile_function)
    IF (laser%use_phase_function) &
        CALL deallocate_stack(laser%phase_function)
    IF (laser%use_time_function) &
        CALL deallocate_stack(laser%time_function)
    IF (laser%use_omega_function) &
        CALL deallocate_stack(laser%omega_function)
    DEALLOCATE(laser)

  END SUBROUTINE deallocate_laser



  ! 释放所有已定义的激光对象的内存链表
  SUBROUTINE deallocate_lasers

    TYPE(laser_block), POINTER :: current, next

    current => lasers
    DO WHILE(ASSOCIATED(current))
      next => current%next
      CALL deallocate_laser(current)
      current => next
    END DO

  END SUBROUTINE deallocate_lasers



  ! 将创建的激光对象挂载到全局激光链表的末尾
  SUBROUTINE attach_laser(laser)

    TYPE(laser_block), POINTER :: laser
    TYPE(laser_block), POINTER :: current
    INTEGER :: boundary

    boundary = laser%boundary

    ! 增加该边界上的激光计数
    n_lasers(boundary) = n_lasers(boundary) + 1

    IF (ASSOCIATED(lasers)) THEN
      current => lasers
      DO WHILE(ASSOCIATED(current%next))
        current => current%next
      END DO
      current%next => laser
    ELSE
      lasers => laser
    END IF

  END SUBROUTINE attach_laser



  ! 根据激光所在的边界，设置参数包中的坐标索引
  ! parameters: 用于表达式求值的参数包
  ! 这个函数的作用是为数学表达式及其求值器 (evaluator) 提供空间上下文。
  ! 因为激光位于特定的模拟框边界上，所以该边界法线方向的坐标是固定的。
  !
  ! 例如：如果激光位于 x_min 边界，那么所有激光点的 ix 坐标都固定为 0。
  ! 只有平行于边界平面的坐标 (例如 y 和 z) 会在后续的循环中变化。
  !
  ! 输入:
  !   laser: 指向激光对象的指针，包含了边界信息 (laser%boundary)。
  ! 输出:
  !   parameters: 一个参数包结构体，将被传递给求值器来计算数学表达式。
  SUBROUTINE populate_pack_from_laser(laser, parameters)

    TYPE(laser_block), POINTER :: laser
    TYPE(parameter_pack), INTENT(INOUT) :: parameters

    ! 初始化所有方向的网格索引为 0
    parameters%pack_ix = 0
    parameters%pack_iy = 0
    parameters%pack_iz = 0

    ! 根据激光具体位于哪个边界，设置对应的固定坐标索引
    SELECT CASE(laser%boundary)
      CASE(c_bd_x_min)
        ! 激光在 X 最小边界，固定 ix = 0
        parameters%pack_ix = 0
      CASE(c_bd_x_max)
        ! 激光在 X 最大边界，固定 ix = nx (X方向网格最大值)
        parameters%pack_ix = nx
      CASE(c_bd_y_min)
        ! 激光在 Y 最小边界，固定 iy = 0
        parameters%pack_iy = 0
      CASE(c_bd_y_max)
        ! 激光在 Y 最大边界，固定 iy = ny
        parameters%pack_iy = ny
      CASE(c_bd_z_min)
        ! 激光在 Z 最小边界，固定 iz = 0
        parameters%pack_iz = 0
      CASE(c_bd_z_max)
        ! 激光在 Z 最大边界，固定 iz = nz
        parameters%pack_iz = nz
    END SELECT

  END SUBROUTINE populate_pack_from_laser



  ! 计算激光的时间分布强度
  ! 如果定义了 t_profile 函数，则评估该函数
  ! 否则使用 custom_laser.f90 中的硬编码配置
  FUNCTION laser_time_profile(laser)

    TYPE(laser_block), POINTER :: laser
    REAL(num) :: laser_time_profile
    INTEGER :: err
    TYPE(parameter_pack) :: parameters

    err = 0
    CALL populate_pack_from_laser(laser, parameters)
    IF (laser%use_time_function) THEN
      laser_time_profile = evaluate_with_parameters(laser%time_function, &
          parameters, err)
      RETURN
    END IF

    laser_time_profile = custom_laser_time_profile(laser)

  END FUNCTION laser_time_profile



  ! 根据 user_deck 中定义的 phase 函数，更新激光在边界上的相位分布数组
  SUBROUTINE laser_update_phase(laser)

    TYPE(laser_block), POINTER :: laser
    INTEGER :: i, j, err
    TYPE(parameter_pack) :: parameters

    err = 0
    CALL populate_pack_from_laser(laser, parameters)
    SELECT CASE(laser%boundary)
      CASE(c_bd_x_min, c_bd_x_max)
        ! 在 Y-Z 平面上循环更新相位
        DO j = 0,nz
          parameters%pack_iz = j
          DO i = 0,ny
            parameters%pack_iy = i
            laser%phase(i,j) = &
                evaluate_with_parameters(laser%phase_function, parameters, err)
          END DO
        END DO
      CASE(c_bd_y_min, c_bd_y_max)
        ! 在 X-Z 平面上循环更新相位
        DO j = 0,nz
          parameters%pack_iz = j
          DO i = 0,nx
            parameters%pack_ix = i
            laser%phase(i,j) = &
                evaluate_with_parameters(laser%phase_function, parameters, err)
          END DO
        END DO
      CASE(c_bd_z_min, c_bd_z_max)
        ! 在 X-Y 平面上循环更新相位
        DO j = 0,ny
          parameters%pack_iy = j
          DO i = 0,nx
            parameters%pack_ix = i
            laser%phase(i,j) = &
                evaluate_with_parameters(laser%phase_function, parameters, err)
          END DO
        END DO
    END SELECT

  END SUBROUTINE laser_update_phase



  ! 根据 user_deck 中定义的 profile 函数，更新激光在边界上的空间强度分布数组
  SUBROUTINE laser_update_profile(laser)

    TYPE(laser_block), POINTER :: laser
    INTEGER :: i, j, err
    TYPE(parameter_pack) :: parameters

    err = 0
    CALL populate_pack_from_laser(laser, parameters)
    SELECT CASE(laser%boundary)
      CASE(c_bd_x_min, c_bd_x_max)
        DO j = 0,nz
          parameters%pack_iz = j
          DO i = 0,ny
            parameters%pack_iy = i
            laser%profile(i,j) = &
                evaluate_with_parameters(laser%profile_function, parameters, &
                    err)
          END DO
        END DO
      CASE(c_bd_y_min, c_bd_y_max)
        DO j = 0,nz
          parameters%pack_iz = j
          DO i = 0,nx
            parameters%pack_ix = i
            laser%profile(i,j) = &
                evaluate_with_parameters(laser%profile_function, parameters, &
                    err)
          END DO
        END DO
      CASE(c_bd_z_min, c_bd_z_max)
        DO j = 0,ny
          parameters%pack_iy = j
          DO i = 0,nx
            parameters%pack_ix = i
            laser%profile(i,j) = &
                evaluate_with_parameters(laser%profile_function, parameters, &
                    err)
          END DO
        END DO
    END SELECT

  END SUBROUTINE laser_update_profile



  ! 更新激光的角频率 omega
  ! 支持从频率 (Hz) 或波长 (m) 转换
  SUBROUTINE laser_update_omega(laser)

    TYPE(laser_block), POINTER :: laser
    INTEGER :: err
    TYPE(parameter_pack) :: parameters

    err = 0
    CALL populate_pack_from_laser(laser, parameters)
    laser%omega = &
        evaluate_with_parameters(laser%omega_function, parameters, err)
    ! 转换频率 (f -> 2*pi*f)
    IF (laser%omega_func_type == c_of_freq) &
        laser%omega = 2.0_num * pi * laser%omega
    ! 转换波长 (lambda -> 2*pi*c/lambda)
    IF (laser%omega_func_type == c_of_lambda) &
        laser%omega = 2.0_num * pi * c / laser%omega

  END SUBROUTINE laser_update_omega



  ! 在每个时间步更新所有启用变频的激光的频率和累积相位
  SUBROUTINE update_laser_omegas

    TYPE(laser_block), POINTER :: current

    current => lasers
    DO WHILE(ASSOCIATED(current))
      IF (current%use_omega_function) THEN
        CALL laser_update_omega(current)
        ! 累加相位： phi = phi + omega * dt
        current%current_integral_phase = current%current_integral_phase &
            + current%omega * dt
      END IF
      current => current%next
    END DO

  END SUBROUTINE update_laser_omegas



  ! 为特定的边界分配二维数组空间（用于 profile 和 phase）
  SUBROUTINE allocate_with_boundary(array, boundary)

    REAL(num), DIMENSION(:,:), POINTER :: array
    INTEGER, INTENT(IN) :: boundary

    IF (boundary == c_bd_x_min .OR. boundary == c_bd_x_max) THEN
      ! 如果是 X 边界，分配 Y-Z 平面
      ALLOCATE(array(1-ng:ny+ng, 1-ng:nz+ng))
    ELSE IF (boundary == c_bd_y_min .OR. boundary == c_bd_y_max) THEN
      ! 如果是 Y 边界，分配 X-Z 平面
      ALLOCATE(array(1-ng:nx+ng, 1-ng:nz+ng))
    ELSE IF (boundary == c_bd_z_min .OR. boundary == c_bd_z_max) THEN
      ! 如果是 Z 边界，分配 X-Y 平面
      ALLOCATE(array(1-ng:nx+ng, 1-ng:ny+ng))
    END IF

  END SUBROUTINE allocate_with_boundary



  ! 根据激光频率计算激光模拟所需的时间步长限制
  SUBROUTINE set_laser_dt

    REAL(num) :: dt_local
    TYPE(laser_block), POINTER :: current

    dt_laser = HUGE(1.0_num)

    current => lasers
    DO WHILE(ASSOCIATED(current))
      ! 计算激光周期 T = 2*pi / omega
      dt_local = 2.0_num * pi / current%omega
      dt_laser = MIN(dt_laser, dt_local)
      current => current%next
    END DO

    ! 根据奈奎斯特采样定理，每个周期至少需要 2 个时间步
    ! (Nyquist)
    dt_laser = dt_laser / 2.0_num

  END SUBROUTINE set_laser_dt



  ! 处理 X 最小边界 (x_min) 的出射边界条件和激光注入
  SUBROUTINE outflow_bcs_x_min

    REAL(num) :: t_env
    REAL(num) :: dtc2, lx, ly, lz, sum, diff, dt_eps, base
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: source1, source2
    INTEGER :: laserpos, n, i, j
    TYPE(laser_block), POINTER :: current

    n = c_bd_x_min

    laserpos = 1
    ! 如果使用了 CPML 吸收边界叠加激光
    IF (bc_field(n) == c_bc_cpml_laser) THEN
      laserpos = cpml_x_min_laser_idx
    END IF
    ! 预计算系数以优化边界条件方程
    dtc2 = dt * c**2
    lx = dtc2 / dx
    ly = dtc2 / dy
    lz = dtc2 / dz
    sum = 1.0_num / (lx + c)
    diff = lx - c
    dt_eps = dt / epsilon0

    ALLOCATE(source1(0:ny,0:nz))
    ALLOCATE(source2(0:ny,0:nz))
    source1 = 0.0_num
    source2 = 0.0_num

    ! 同步边界上的磁场 bx
    bx(laserpos-1, 0:ny, 0:nz) = bx_x_min(0:ny, 0:nz)

    ! 如果该边界上有激光需要注入
    IF (add_laser(n)) THEN
      current => lasers
      DO WHILE(ASSOCIATED(current))
        IF (current%boundary == c_bd_x_min) THEN
          ! 评估激光随时间的演化：只在 [t_start, t_end] 范围内有效
          IF (time >= current%t_start .AND. time <= current%t_end) THEN
            ! 如果设置了动态变化的相位或剖面函数，则更新它们
            IF (current%use_phase_function) CALL laser_update_phase(current)
            IF (current%use_profile_function) CALL laser_update_profile(current)
            ! 计算时间包络强度
            t_env = laser_time_profile(current) * current%amp
            DO j = 0,nz
            DO i = 0,ny
              ! 计算激光电场矢量： base = amp * profile * sin(phase_integral + phase_offset)
              base = t_env * current%profile(i,j) &
                * SIN(current%current_integral_phase + current%phase(i,j))
              ! 根据偏振角分解到两个分量 (S 与 P 偏振)
              source1(i,j) = source1(i,j) + base * COS(current%pol_angle)
              source2(i,j) = source2(i,j) + base * SIN(current%pol_angle)
            END DO
            END DO
          END IF
        END IF
        current => current%next
      END DO
    END IF

    ! 将激光源项代入 Mur 出射边界条件方程，计算边界处的磁场分量 bz 和 by
    bz(laserpos-1, 0:ny, 0:nz) = sum * ( 4.0_num * source1 &
        + 2.0_num * (ey_x_min(0:ny, 0:nz) + c * bz_x_min(0:ny, 0:nz)) &
        - 2.0_num * ey(laserpos, 0:ny, 0:nz) &
        - lz * (bx(laserpos, 0:ny, 0:nz) - bx(laserpos, 0:ny, -1:nz-1)) &
        + dt_eps * jy(laserpos, 0:ny, 0:nz) &
        + diff * bz(laserpos, 0:ny, 0:nz))

    by(laserpos-1, 0:ny, 0:nz) = sum * (-4.0_num * source2 &
        - 2.0_num * (ez_x_min(0:ny, 0:nz) - c * by_x_min(0:ny, 0:nz)) &
        + 2.0_num * ez(laserpos, 0:ny, 0:nz) &
        - ly * (bx(laserpos, 0:ny, 0:nz) - bx(laserpos, -1:ny-1, 0:nz)) &
        - dt_eps * jz(laserpos, 0:ny, 0:nz) &
        + diff * by(laserpos, 0:ny, 0:nz))

    DEALLOCATE(source1, source2)

    ! 如果启用了能量吸收计算
    IF (dump_absorption) THEN
      IF (add_laser(n)) THEN
        CALL calc_absorption(c_bd_x_min, lasers=lasers)
      ELSE
        CALL calc_absorption(c_bd_x_min)
      END IF
    END IF

  END SUBROUTINE outflow_bcs_x_min



  SUBROUTINE outflow_bcs_x_max

    REAL(num) :: t_env
    REAL(num) :: dtc2, lx, ly, lz, sum, diff, dt_eps, base
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: source1, source2
    INTEGER :: laserpos, n, i, j
    TYPE(laser_block), POINTER :: current

    n = c_bd_x_max

    laserpos = nx
    IF (bc_field(n) == c_bc_cpml_laser) THEN
      laserpos = cpml_x_max_laser_idx
    END IF
    dtc2 = dt * c**2
    lx = dtc2 / dx
    ly = dtc2 / dy
    lz = dtc2 / dz
    sum = 1.0_num / (lx + c)
    diff = lx - c
    dt_eps = dt / epsilon0

    ALLOCATE(source1(0:ny,0:nz))
    ALLOCATE(source2(0:ny,0:nz))
    source1 = 0.0_num
    source2 = 0.0_num

    ! 同步边界上的磁场 bx
    bx(laserpos+1, 0:ny, 0:nz) = bx_x_max(0:ny, 0:nz)

    ! 如果该边界上有激光需要注入
    IF (add_laser(n)) THEN
      current => lasers
      DO WHILE(ASSOCIATED(current))
        IF (current%boundary == c_bd_x_max) THEN
          ! 评估激光随时间的演化：只在 [t_start, t_end] 范围内有效
          IF (time >= current%t_start .AND. time <= current%t_end) THEN
            ! 如果设置了动态变化的相位或剖面函数，则更新它们
            IF (current%use_phase_function) CALL laser_update_phase(current)
            IF (current%use_profile_function) CALL laser_update_profile(current)
            ! 计算时间包络强度
            t_env = laser_time_profile(current) * current%amp
            DO j = 0,nz
            DO i = 0,ny
              ! 计算激光电场矢量： base = amp * profile * sin(phase_integral + phase_offset)
              base = t_env * current%profile(i,j) &
                * SIN(current%current_integral_phase + current%phase(i,j))
              ! 根据偏振角分解到两个分量 (S 与 P 偏振)
              source1(i,j) = source1(i,j) + base * COS(current%pol_angle)
              source2(i,j) = source2(i,j) + base * SIN(current%pol_angle)
            END DO
            END DO
          END IF
        END IF
        current => current%next
      END DO
    END IF

    ! 将激光源项代入 Mur 出射边界条件方程，计算边界处的磁场分量 bz 和 by
    bz(laserpos, 0:ny, 0:nz) = sum * (-4.0_num * source1 &
        - 2.0_num * (ey_x_max(0:ny, 0:nz) - c * bz_x_max(0:ny, 0:nz)) &
        + 2.0_num * ey(laserpos, 0:ny, 0:nz) &
        + lz * (bx(laserpos, 0:ny, 0:nz) - bx(laserpos, 0:ny, -1:nz-1)) &
        - dt_eps * jy(laserpos, 0:ny, 0:nz) &
        + diff * bz(laserpos-1, 0:ny, 0:nz))

    by(laserpos, 0:ny, 0:nz) = sum * ( 4.0_num * source2 &
        + 2.0_num * (ez_x_max(0:ny, 0:nz) + c * by_x_max(0:ny, 0:nz)) &
        - 2.0_num * ez(laserpos, 0:ny, 0:nz) &
        + ly * (bx(laserpos, 0:ny, 0:nz) - bx(laserpos, -1:ny-1, 0:nz)) &
        + dt_eps * jz(laserpos, 0:ny, 0:nz) &
        + diff * by(laserpos-1, 0:ny, 0:nz))

    DEALLOCATE(source1, source2)

    ! 如果启用了能量吸收计算
    IF (dump_absorption) THEN
      IF (add_laser(n)) THEN
        CALL calc_absorption(c_bd_x_max, lasers=lasers)
      ELSE
        CALL calc_absorption(c_bd_x_max)
      END IF
    END IF

  END SUBROUTINE outflow_bcs_x_max



  SUBROUTINE outflow_bcs_y_min

    REAL(num) :: t_env
    REAL(num) :: dtc2, lx, ly, lz, sum, diff, dt_eps, base
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: source1, source2
    INTEGER :: laserpos, n, i, j
    TYPE(laser_block), POINTER :: current

    n = c_bd_y_min

    laserpos = 1
    IF (bc_field(n) == c_bc_cpml_laser) THEN
      laserpos = cpml_y_min_laser_idx
    END IF
    dtc2 = dt * c**2
    lx = dtc2 / dx
    ly = dtc2 / dy
    lz = dtc2 / dz
    sum = 1.0_num / (ly + c)
    diff = ly - c
    dt_eps = dt / epsilon0

    ALLOCATE(source1(0:nx,0:nz))
    ALLOCATE(source2(0:nx,0:nz))
    source1 = 0.0_num
    source2 = 0.0_num

    ! 同步边界上的磁场 by
    by(0:nx, laserpos-1, 0:nz) = by_y_min(0:nx, 0:nz)

    ! 如果该边界上有激光需要注入
    IF (add_laser(n)) THEN
      current => lasers
      DO WHILE(ASSOCIATED(current))
        IF (current%boundary == c_bd_y_min) THEN
          ! 评估激光随时间的演化：只在 [t_start, t_end] 范围内有效
          IF (time >= current%t_start .AND. time <= current%t_end) THEN
            ! 如果设置了动态变化的相位或剖面函数，则更新它们
            IF (current%use_phase_function) CALL laser_update_phase(current)
            IF (current%use_profile_function) CALL laser_update_profile(current)
            ! 计算时间包络强度
            t_env = laser_time_profile(current) * current%amp
            DO j = 0,nz
            DO i = 0,nx
              ! 计算激光电场矢量： base = amp * profile * sin(phase_integral + phase_offset)
              base = t_env * current%profile(i,j) &
                * SIN(current%current_integral_phase + current%phase(i,j))
              ! 根据偏振角分解到两个分量 (S 与 P 偏振)
              source1(i,j) = source1(i,j) + base * COS(current%pol_angle)
              source2(i,j) = source2(i,j) + base * SIN(current%pol_angle)
            END DO
            END DO
          END IF
        END IF
        current => current%next
      END DO
    END IF

    ! 将激光源项代入 Mur 出射边界条件方程，计算边界处的磁场分量 bz 和 by
    bx(0:nx, laserpos-1, 0:nz) = sum * ( 4.0_num * source1 &
        + 2.0_num * (ez_y_min(0:nx, 0:nz) + c * bx_y_min(0:nx, 0:nz)) &
        - 2.0_num * ez(0:nx, laserpos, 0:nz) &
        - lx * (by(0:nx, laserpos, 0:nz) - by(-1:nx-1, laserpos, 0:nz)) &
        + dt_eps * jz(0:nx, laserpos, 0:nz) &
        + diff * bx(0:nx, laserpos, 0:nz))

    bz(0:nx, laserpos-1, 0:nz) = sum * (-4.0_num * source2 &
        - 2.0_num * (ex_y_min(0:nx, 0:nz) - c * bz_y_min(0:nx, 0:nz)) &
        + 2.0_num * ex(0:nx, laserpos, 0:nz) &
        - lz * (by(0:nx, laserpos, 0:nz) - by(0:nx, laserpos, -1:nz-1)) &
        - dt_eps * jx(0:nx, laserpos, 0:nz) &
        + diff * bz(0:nx, laserpos, 0:nz))

    DEALLOCATE(source1, source2)

    ! 如果启用了能量吸收计算
    IF (dump_absorption) THEN
      IF (add_laser(n)) THEN
        CALL calc_absorption(c_bd_y_min, lasers=lasers)
      ELSE
        CALL calc_absorption(c_bd_y_min)
      END IF
    END IF

  END SUBROUTINE outflow_bcs_y_min



  SUBROUTINE outflow_bcs_y_max

    REAL(num) :: t_env
    REAL(num) :: dtc2, lx, ly, lz, sum, diff, dt_eps, base
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: source1, source2
    INTEGER :: laserpos, n, i, j
    TYPE(laser_block), POINTER :: current

    n = c_bd_y_max

    laserpos = ny
    IF (bc_field(n) == c_bc_cpml_laser) THEN
      laserpos = cpml_y_max_laser_idx
    END IF
    dtc2 = dt * c**2
    lx = dtc2 / dx
    ly = dtc2 / dy
    lz = dtc2 / dz
    sum = 1.0_num / (ly + c)
    diff = ly - c
    dt_eps = dt / epsilon0

    ALLOCATE(source1(0:nx,0:nz))
    ALLOCATE(source2(0:nx,0:nz))
    source1 = 0.0_num
    source2 = 0.0_num

    ! 同步边界上的磁场 by
    by(0:nx, laserpos+1, 0:nz) = by_y_max(0:nx, 0:nz)

    ! 如果该边界上有激光需要注入
    IF (add_laser(n)) THEN
      current => lasers
      DO WHILE(ASSOCIATED(current))
        IF (current%boundary == c_bd_y_max) THEN
          ! 评估激光随时间的演化：只在 [t_start, t_end] 范围内有效
          IF (time >= current%t_start .AND. time <= current%t_end) THEN
            ! 如果设置了动态变化的相位或剖面函数，则更新它们
            IF (current%use_phase_function) CALL laser_update_phase(current)
            IF (current%use_profile_function) CALL laser_update_profile(current)
            ! 计算时间包络强度
            t_env = laser_time_profile(current) * current%amp
            DO j = 0,nz
            DO i = 0,nx
              ! 计算激光电场矢量： base = amp * profile * sin(phase_integral + phase_offset)
              base = t_env * current%profile(i,j) &
                * SIN(current%current_integral_phase + current%phase(i,j))
              ! 根据偏振角分解到两个分量 (S 与 P 偏振)
              source1(i,j) = source1(i,j) + base * COS(current%pol_angle)
              source2(i,j) = source2(i,j) + base * SIN(current%pol_angle)
            END DO
            END DO
          END IF
        END IF
        current => current%next
      END DO
    END IF

    ! 将激光源项代入 Mur 出射边界条件方程，计算边界处的磁场分量 bz 和 by
    bz(laserpos, 0:ny, 0:nz) = sum * (-4.0_num * source1 &
        - 2.0_num * (ey_x_max(0:ny, 0:nz) - c * bz_x_max(0:ny, 0:nz)) &
        + 2.0_num * ey(laserpos, 0:ny, 0:nz) &
        + lz * (bx(laserpos, 0:ny, 0:nz) - bx(laserpos, 0:ny, -1:nz-1)) &
        - dt_eps * jy(laserpos, 0:ny, 0:nz) &
        + diff * bz(laserpos-1, 0:ny, 0:nz))

    by(laserpos, 0:ny, 0:nz) = sum * ( 4.0_num * source2 &
        + 2.0_num * (ez_x_max(0:ny, 0:nz) + c * by_x_max(0:ny, 0:nz)) &
        - 2.0_num * ez(laserpos, 0:ny, 0:nz) &
        + ly * (bx(laserpos, 0:ny, 0:nz) - bx(laserpos, -1:ny-1, 0:nz)) &
        + dt_eps * jz(laserpos, 0:ny, 0:nz) &
        + diff * by(laserpos-1, 0:ny, 0:nz))

    DEALLOCATE(source1, source2)

    ! 如果启用了能量吸收计算
    IF (dump_absorption) THEN
      IF (add_laser(n)) THEN
        CALL calc_absorption(c_bd_y_max, lasers=lasers)
      ELSE
        CALL calc_absorption(c_bd_y_max)
      END IF
    END IF

  END SUBROUTINE outflow_bcs_y_max



  SUBROUTINE outflow_bcs_z_min

    REAL(num) :: t_env
    REAL(num) :: dtc2, lx, ly, lz, sum, diff, dt_eps, base
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: source1, source2
    INTEGER :: laserpos, n, i, j
    TYPE(laser_block), POINTER :: current

    n = c_bd_z_min

    laserpos = 1
    IF (bc_field(n) == c_bc_cpml_laser) THEN
      laserpos = cpml_z_min_laser_idx
    END IF
    dtc2 = dt * c**2
    lx = dtc2 / dx
    ly = dtc2 / dy
    lz = dtc2 / dz
    sum = 1.0_num / (lz + c)
    diff = lz - c
    dt_eps = dt / epsilon0

    ALLOCATE(source1(0:nx,0:ny))
    ALLOCATE(source2(0:nx,0:ny))
    source1 = 0.0_num
    source2 = 0.0_num

    ! 同步边界上的磁场 bz
    bz(0:nx, 0:ny, laserpos-1) = bz_z_min(0:nx, 0:ny)

    ! 如果该边界上有激光需要注入
    IF (add_laser(n)) THEN
      current => lasers
      DO WHILE(ASSOCIATED(current))
        IF (current%boundary == c_bd_z_min) THEN
          ! 评估激光随时间的演化：只在 [t_start, t_end] 范围内有效
          IF (time >= current%t_start .AND. time <= current%t_end) THEN
            ! 如果设置了动态变化的相位或剖面函数，则更新它们
            IF (current%use_phase_function) CALL laser_update_phase(current)
            IF (current%use_profile_function) CALL laser_update_profile(current)
            ! 计算时间包络强度
            t_env = laser_time_profile(current) * current%amp
            DO j = 0,ny
            DO i = 0,nx
              ! 计算激光电场矢量： base = amp * profile * sin(phase_integral + phase_offset)
              base = t_env * current%profile(i,j) &
                * SIN(current%current_integral_phase + current%phase(i,j))
              ! 根据偏振角分解到两个分量 (S 与 P 偏振)
              source1(i,j) = source1(i,j) + base * COS(current%pol_angle)
              source2(i,j) = source2(i,j) + base * SIN(current%pol_angle)
            END DO
            END DO
          END IF
        END IF
        current => current%next
      END DO
    END IF

    ! 将激光源项代入 Mur 出射边界条件方程，计算边界处的磁场分量 bz 和 by
    by(0:nx, 0:ny, laserpos-1) = sum * ( 4.0_num * source1 &
        + 2.0_num * (ex_z_min(0:nx, 0:ny) + c * by_z_min(0:nx, 0:ny)) &
        - 2.0_num * ex(0:nx, 0:ny, laserpos) &
        - ly * (bz(0:nx, 0:ny, laserpos) - bz(0:nx, -1:ny-1, laserpos)) &
        + dt_eps * jx(0:nx, 0:ny, laserpos) &
        + diff * by(0:nx, 0:ny, laserpos))

    bx(0:nx, 0:ny, laserpos-1) = sum * (-4.0_num * source2 &
        - 2.0_num * (ey_z_min(0:nx, 0:ny) - c * bx_z_min(0:nx, 0:ny)) &
        + 2.0_num * ey(0:nx, 0:ny, laserpos) &
        - lx * (bz(0:nx, 0:ny, laserpos) - bz(-1:nx-1, 0:ny, laserpos)) &
        - dt_eps * jy(0:nx, 0:ny, laserpos) &
        + diff * bx(0:nx, 0:ny, laserpos))

    DEALLOCATE(source1, source2)

    ! 如果启用了能量吸收计算
    IF (dump_absorption) THEN
      IF (add_laser(n)) THEN
        CALL calc_absorption(c_bd_z_min, lasers=lasers)
      ELSE
        CALL calc_absorption(c_bd_z_min)
      END IF
    END IF

  END SUBROUTINE outflow_bcs_z_min



  SUBROUTINE outflow_bcs_z_max

    REAL(num) :: t_env
    REAL(num) :: dtc2, lx, ly, lz, sum, diff, dt_eps, base
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: source1, source2
    INTEGER :: laserpos, n, i, j
    TYPE(laser_block), POINTER :: current

    n = c_bd_z_max

    laserpos = nz
    IF (bc_field(n) == c_bc_cpml_laser) THEN
      laserpos = cpml_z_max_laser_idx
    END IF
    dtc2 = dt * c**2
    lx = dtc2 / dx
    ly = dtc2 / dy
    lz = dtc2 / dz
    sum = 1.0_num / (lz + c)
    diff = lz - c
    dt_eps = dt / epsilon0

    ALLOCATE(source1(0:nx,0:ny))
    ALLOCATE(source2(0:nx,0:ny))
    source1 = 0.0_num
    source2 = 0.0_num

    ! 同步边界上的磁场 bz
    bz(0:nx, 0:ny, laserpos+1) = bz_z_max(0:nx, 0:ny)

    ! 如果该边界上有激光需要注入
    IF (add_laser(n)) THEN
      current => lasers
      DO WHILE(ASSOCIATED(current))
        IF (current%boundary == c_bd_z_max) THEN
          ! 评估激光随时间的演化：只在 [t_start, t_end] 范围内有效
          IF (time >= current%t_start .AND. time <= current%t_end) THEN
            ! 如果设置了动态变化的相位或剖面函数，则更新它们
            IF (current%use_phase_function) CALL laser_update_phase(current)
            IF (current%use_profile_function) CALL laser_update_profile(current)
            ! 计算时间包络强度
            t_env = laser_time_profile(current) * current%amp
            DO j = 0,ny
            DO i = 0,nx
              ! 计算激光电场矢量： base = amp * profile * sin(phase_integral + phase_offset)
              base = t_env * current%profile(i,j) &
                * SIN(current%current_integral_phase + current%phase(i,j))
              ! 根据偏振角分解到两个分量 (S 与 P 偏振)
              source1(i,j) = source1(i,j) + base * COS(current%pol_angle)
              source2(i,j) = source2(i,j) + base * SIN(current%pol_angle)
            END DO
            END DO
          END IF
        END IF
        current => current%next
      END DO
    END IF

    ! 将激光源项代入 Mur 出射边界条件方程，计算边界处的磁场分量 bz 和 by
    by(0:nx, 0:ny, laserpos) = sum * (-4.0_num * source1 &
        - 2.0_num * (ex_z_max(0:nx, 0:ny) - c * by_z_max(0:nx, 0:ny)) &
        + 2.0_num * ex(0:nx, 0:ny, laserpos) &
        + ly * (bz(0:nx, 0:ny, laserpos) - bz(0:nx, -1:ny-1, laserpos)) &
        - dt_eps * jx(0:nx, 0:ny, laserpos) &
        + diff * by(0:nx, 0:ny, laserpos-1))

    bx(0:nx, 0:ny, laserpos) = sum * ( 4.0_num * source2 &
        + 2.0_num * (ey_z_max(0:nx, 0:ny) + c * bx_z_max(0:nx, 0:ny)) &
        - 2.0_num * ey(0:nx, 0:ny, laserpos) &
        + lx * (bz(0:nx, 0:ny, laserpos) - bz(-1:nx-1, 0:ny, laserpos)) &
        + dt_eps * jy(0:nx, 0:ny, laserpos) &
        + diff * bx(0:nx, 0:ny, laserpos-1))

    DEALLOCATE(source1, source2)

    ! 如果启用了能量吸收计算
    IF (dump_absorption) THEN
      IF (add_laser(n)) THEN
        CALL calc_absorption(c_bd_z_max, lasers=lasers)
      ELSE
        CALL calc_absorption(c_bd_z_max)
      END IF
    END IF

  END SUBROUTINE outflow_bcs_z_max



  SUBROUTINE calc_absorption(bd, lasers)

    TYPE(laser_block), POINTER, OPTIONAL :: lasers
    INTEGER, INTENT(IN) :: bd
    TYPE(laser_block), POINTER :: current
    REAL(num) :: t_env, dir, dd, factor, lfactor, laser_inject_sum
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: e1, e2, b1, b2
    INTEGER :: mm, nn, ibc, icell, jcell

    ! Note: ideally e1, e2, b1, b2 should be face-centred. However, this is not
    ! possible with 'open' boundaries since E-fields are not defined in the
    ! ghost cell, so we use the cell-centred quantities in the first cell.

    dir = 1.0_num
    mm = 1
    nn = 1

    SELECT CASE(bd)
      CASE(c_bd_x_min, c_bd_x_max)
        dd = dy * dz
        mm = ny
        nn = nz
        ALLOCATE(e1(mm, nn), e2(mm, nn), b1(mm, nn), b2(mm, nn))

        ibc = 1
        IF (bd == c_bd_x_max) THEN
          dir = -1.0_num
          ibc = nx
        END IF

        e1 = 0.5_num  * (ey(ibc  , 0:ny-1, 1:nz  ) + ey(ibc, 1:ny  , 1:nz  ))
        e2 = 0.5_num  * (ez(ibc  , 1:ny  , 0:nz-1) + ez(ibc, 1:ny  , 1:nz  ))
        b1 = 0.25_num * (bz(ibc-1, 0:ny-1, 1:nz  ) + bz(ibc, 0:ny-1, 1:nz  ) &
                       + bz(ibc-1, 1:ny  , 1:nz  ) + bz(ibc, 1:ny  , 1:nz  ))
        b2 = 0.25_num * (by(ibc-1, 1:ny  , 0:nz-1) + by(ibc, 1:ny  , 0:nz-1) &
                       + by(ibc-1, 1:ny  , 1:nz  ) + by(ibc, 1:ny  , 1:nz  ))

      CASE(c_bd_y_min, c_bd_y_max)
        dd = dx * dz
        mm = nx
        nn = nz
        ALLOCATE(e1(mm, nn), e2(mm, nn), b1(mm, nn), b2(mm, nn))

        ibc = 1
        IF (bd == c_bd_y_max) THEN
          dir = -1.0_num
          ibc = ny
        END IF

        e1 = 0.5_num  * (ez(1:nx  , ibc  , 0:nz-1) + ez(1:nx  , ibc, 1:nz  ))
        e2 = 0.5_num  * (ex(0:nx-1, ibc  , 1:nz  ) + ex(1:nx  , ibc, 1:nz  ))
        b1 = 0.25_num * (bx(1:nx  , ibc-1, 0:nz-1) + bx(1:nx  , ibc, 0:nz-1) &
                       + bx(1:nx  , ibc-1, 1:nz  ) + bx(1:nx  , ibc, 1:nz  ))
        b2 = 0.25_num * (bz(0:nx-1, ibc-1, 1:nz  ) + bz(0:nx-1, ibc, 1:nz  ) &
                       + bz(1:nx  , ibc-1, 1:nz  ) + bz(1:nx  , ibc, 1:nz  ))

      CASE(c_bd_z_min, c_bd_z_max)
        dd = dx * dy
        mm = nx
        nn = ny
        ALLOCATE(e1(mm, nn), e2(mm, nn), b1(mm, nn), b2(mm, nn))

        ibc = 1
        IF (bd == c_bd_z_max) THEN
          dir = -1.0_num
          ibc = nz
        END IF

        e1 = 0.5_num  * (ex(0:nx-1, 1:ny  , ibc  ) + ex(1:nx  , 1:ny  , ibc))
        e2 = 0.5_num  * (ey(1:nx  , 0:ny-1, ibc  ) + ey(1:nx  , 1:ny  , ibc))
        b1 = 0.25_num * (by(0:nx-1, 1:ny  , ibc-1) + by(0:nx-1, 1:ny  , ibc) &
                       + by(1:nx  , 1:ny  , ibc-1) + by(1:nx  , 1:ny  , ibc))
        b2 = 0.25_num * (bx(1:nx  , 0:ny-1, ibc-1) + bx(1:nx  , 0:ny-1, ibc) &
                       + bx(1:nx  , 1:ny  , ibc-1) + bx(1:nx  , 1:ny  , ibc))

      CASE DEFAULT
        dd = 0.0_num
        ALLOCATE(e1(mm, nn), e2(mm, nn), b1(mm, nn), b2(mm, nn))

        e1 = 0.0_num
        e2 = 0.0_num
        b1 = 0.0_num
        b2 = 0.0_num
    END SELECT

    factor = dt * dd * dir
    laser_absorb_local = laser_absorb_local &
        + (factor / mu0) * SUM(e1 * b1 - e2 * b2)

    IF (PRESENT(lasers)) THEN
      current => lasers
      DO WHILE(ASSOCIATED(current))
        IF (current%boundary == bd) THEN
          laser_inject_sum = 0.0_num
          DO jcell = 1, nn
          DO icell = 1, mm
            laser_inject_sum = laser_inject_sum &
                + current%profile(icell, jcell)**2
          END DO
          END DO
          t_env = laser_time_profile(current)
          lfactor = 0.5_num * epsilon0 * c * factor * (t_env * current%amp)**2
          laser_inject_local = laser_inject_local + lfactor * laser_inject_sum
        END IF
        current => current%next
      END DO
    END IF

    DEALLOCATE(e1, e2, b1, b2)

  END SUBROUTINE calc_absorption

END MODULE laser
