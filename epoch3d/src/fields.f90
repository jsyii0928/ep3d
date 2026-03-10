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

MODULE fields

  USE boundary

  IMPLICIT NONE

  ! =========================
  ! 本模块用途（Fields 模块）
  ! =========================
  ! 负责推进电磁场（E/B），实现 Maxwell 方程在 Yee 网格上的离散更新。
  ! 主要提供：
  !   - set_field_order: 设置空间差分阶数（2/4/6）以及相应 CFL 系数
  !   - set_maxwell_solver: 选择并配置不同的 Maxwell 求解器（Yee / Lehe / Cowan / Pukhov / 自定义）
  !   - update_e_field / update_b_field: 根据离散旋度更新 E/B（并考虑电流项与 CPML）
  !   - update_eb_fields_half / update_eb_fields_final: 时间推进（半步/终步），并调用边界条件
  !
  ! 说明：本文件只添加注释，不改变任何数值实现。

  INTEGER :: field_order
  ! field_order: 场更新的空间差分阶数（典型值 2/4/6）；越高阶数值色散越小但 stencil 更宽
  REAL(num) :: hdt, fac
  REAL(num) :: hdtx, hdty, hdtz
  REAL(num) :: cnx, cny, cnz
  ! hdt  : 半时间步长 dt/2
  ! hdtx : (dt/2)/dx，用于 B 场更新中 curl(E) 的系数
  ! cnx  : c^2*(dt/2)/dx，用于 E 场更新中 curl(B) 的系数
  ! fac  : (dt/2)/epsilon0，用于 E 场更新中的电流项 -J/epsilon0

  ! 以下系数用于改进的 Maxwell 求解器（通过修改 curl 算子来抑制数值色散/数值切伦科夫辐射等）
  ! 对 Yee 求解器，这些系数通常不参与计算（或等效为标准值）。
  REAL(num) :: alphax = 1.0_num, alphay = 1.0_num, alphaz = 1.0_num
  REAL(num) :: betaxy = 0.0_num, betayx = 0.0_num, betazx = 0.0_num
  REAL(num) :: betaxz = 0.0_num, betazy = 0.0_num, betayz = 0.0_num
  REAL(num) :: deltax = 0.0_num, deltay = 0.0_num, deltaz = 0.0_num
  REAL(num) :: gammax = 0.0_num, gammay = 0.0_num, gammaz = 0.0_num

CONTAINS

  SUBROUTINE set_field_order(order)

    INTEGER, INTENT(IN) :: order

  ! 设置空间差分阶数，并据此确定守恒变量在边界处需要的 ghost cell 数量等（fng 在其它模块中使用）
  field_order = order
    fng = field_order / 2

  ! CFL 系数：用于限制 dt 以满足稳定性条件；高阶差分的 CFL 上限会更小
    IF (field_order == 2) THEN
      cfl = 1.0_num
    ELSE IF (field_order == 4) THEN
      cfl = 6.0_num / 7.0_num
    ELSE
      cfl = 120.0_num / 149.0_num
    END IF

  END SUBROUTINE set_field_order



  SUBROUTINE set_maxwell_solver

    REAL(num) :: delta, dx_cdt
    REAL(num) :: c1, c2, c3, cx1, cx2

  ! 根据输入 deck 设置 Maxwell 求解器的系数。
  ! 这些系数会在 update_b_field（非 Yee 分支）中用于构造“改进的 curl”算子。
  ! 参考文献见各分支注释。
  IF (maxwell_solver == c_maxwell_solver_custom) THEN
      alphax = 1.0_num - 2.0_num * betaxy - 2.0_num * betaxz &
                       - 4.0_num * gammax - 3.0_num * deltax
      alphay = 1.0_num - 2.0_num * betayx - 2.0_num * betayz &
                       - 4.0_num * gammay - 3.0_num * deltay
      alphaz = 1.0_num - 2.0_num * betazx - 2.0_num * betazy &
                       - 4.0_num * gammaz - 3.0_num * deltaz

    ELSE IF (maxwell_solver == c_maxwell_solver_lehe_x) THEN
      ! R. Lehe et al., Phys. Rev. ST Accel. Beams 16, 021301 (2013)
  ! Lehe 求解器：在一个首选传播方向（这里是 x）上优化离散色散关系，以减弱
  ! relativistic beam/等离子体中的数值切伦科夫不稳定性。
      dx_cdt = dx / (c * dt)
      betaxy = 0.125_num * (dx / dy)**2
      betaxz = 0.125_num * (dx / dz)**2
      betayx = 0.125_num
      betazx = 0.125_num
      betazy = 0.0_num
      betayz = 0.0_num
      deltax = 0.25_num * (1.0_num - dx_cdt**2 * SIN(0.5_num * pi / dx_cdt)**2)
      deltay = 0.0_num
      deltaz = 0.0_num
      gammax = 0.0_num
      gammay = 0.0_num
      gammaz = 0.0_num
      alphax = 1.0_num - 2.0_num * betaxy - 2.0_num * betaxz - 3.0_num * deltax
      alphay = 1.0_num - 2.0_num * betayx
      alphaz = 1.0_num - 2.0_num * betazx

    ELSE IF (maxwell_solver == c_maxwell_solver_lehe_y) THEN
  ! 同上，但首选方向为 y
      dx_cdt = dy / (c * dt)
      betayx = 0.125_num * (dy / dx)**2
      betayz = 0.125_num * (dy / dz)**2
      betaxy = 0.125_num
      betazy = 0.125_num
      betaxz = 0.0_num
      betazx = 0.0_num
      deltax = 0.0_num
      deltay = 0.25_num * (1.0_num - dx_cdt**2 * SIN(0.5_num * pi / dx_cdt)**2)
      deltaz = 0.0_num
      gammax = 0.0_num
      gammay = 0.0_num
      gammaz = 0.0_num
      alphax = 1.0_num - 2.0_num * betaxy
      alphay = 1.0_num - 2.0_num * betayx - 2.0_num * betayz - 3.0_num * deltay
      alphaz = 1.0_num - 2.0_num * betazx

    ELSE IF (maxwell_solver == c_maxwell_solver_lehe_z) THEN
  ! 同上，但首选方向为 z
      dx_cdt = dz / (c * dt)
      betazx = 0.125_num * (dz / dx)**2
      betazy = 0.125_num * (dz / dy)**2
      betaxz = 0.125_num
      betayz = 0.125_num
      betayx = 0.0_num
      betaxy = 0.0_num
      deltax = 0.0_num
      deltay = 0.0_num
      deltaz = 0.25_num * (1.0_num - dx_cdt**2 * SIN(0.5_num * pi / dx_cdt)**2)
      gammax = 0.0_num
      gammay = 0.0_num
      gammaz = 0.0_num
      alphax = 1.0_num - 2.0_num * betaxy
      alphay = 1.0_num - 2.0_num * betayx
      alphaz = 1.0_num - 2.0_num * betazx - 2.0_num * betazy - 3.0_num * deltaz

    ELSE IF (maxwell_solver == c_maxwell_solver_cowan) THEN
      ! Cowan et al., Phys. Rev. ST Accel. Beams 16, 041303 (2013)
  ! Cowan 求解器：构造各向同性的改进 curl 算子，兼顾色散与稳定性。
      delta = MIN(dx, dy, dz)
      c1 = (delta / dx)**2
      c2 = (delta / dy)**2
      c3 = (delta / dz)**2
      cx1 = 1.0_num / (c1 * c2 + c2 * c3 + c1 * c3)
      cx2 = 1.0_num - c1 * c2 * c3 * cx1

      betayx = 0.125_num * c1 * cx2
      betaxy = 0.125_num * c2 * cx2
      betaxz = 0.125_num * c3 * cx2
      betazx = betayx
      betazy = betaxy
      betayz = betaxz
      deltax = 0.0_num
      deltay = 0.0_num
      deltaz = 0.0_num
      gammax = c2 * c3 * (0.0625_num - 0.125_num * c2 * c3 * cx1)
      gammay = c1 * c3 * (0.0625_num - 0.125_num * c1 * c3 * cx1)
      gammaz = c1 * c2 * (0.0625_num - 0.125_num * c1 * c2 * cx1)
      alphax = 1.0_num - 2.0_num * betaxy - 2.0_num * betaxz - 4.0_num * gammax
      alphay = 1.0_num - 2.0_num * betayx - 2.0_num * betayz - 4.0_num * gammay
      alphaz = 1.0_num - 2.0_num * betazx - 2.0_num * betazy - 4.0_num * gammaz

    ELSE IF (maxwell_solver == c_maxwell_solver_pukhov) THEN
      ! A. Pukhov, Journal of Plasma Physics 61, 425-433 (1999)
  ! Pukhov 求解器：使用最小网格间距 delta 构造在三维上相对均衡的高精度 curl。
      delta = MIN(dx, dy, dz)

      betayx = 0.125_num * (delta / dx)**2
      betaxy = 0.125_num * (delta / dy)**2
      betaxz = 0.125_num * (delta / dz)**2
      betazx = betayx
      betazy = betaxy
      betayz = betaxz
      deltax = 0.0_num
      deltay = 0.0_num
      deltaz = 0.0_num
      gammax = 0.0_num
      gammay = 0.0_num
      gammaz = 0.0_num
      alphax = 1.0_num - 2.0_num * betaxy - 2.0_num * betaxz
      alphay = 1.0_num - 2.0_num * betayx - 2.0_num * betayz
      alphaz = 1.0_num - 2.0_num * betazx - 2.0_num * betazy
    END IF

  ! 仅在非 Yee 求解器下打印参数，便于检查 deck 设置是否生效
  IF (rank == 0 .AND. maxwell_solver /= c_maxwell_solver_yee) THEN
      PRINT*, 'Maxwell solver set to the following parameters:'
      PRINT'(A9, 3F14.9)', 'alpha =', alphax, alphay, alphaz
      PRINT'(A9, 2F14.9)', 'betax =', betaxy, betaxz
      PRINT'(A9, 2F14.9)', 'betay =', betayx, betayz
      PRINT'(A9, 2F14.9)', 'betaz =', betazx, betazy
      PRINT'(A9, 3F14.9)', 'gamma =', gammax, gammay, gammaz
      PRINT'(A9, 3F14.9)', 'delta =', deltax, deltay, deltaz
      PRINT'(A9, 1F14.9)', 'c*dt/dx = ', dt * c / dx
      PRINT*
    END IF

  END SUBROUTINE set_maxwell_solver



  SUBROUTINE update_e_field

    INTEGER :: ix, iy, iz
    REAL(num) :: cpml_x, cpml_y, cpml_z
    REAL(num) :: c1, c2, c3
    REAL(num) :: cx1, cx2, cx3
    REAL(num) :: cy1, cy2, cy3
    REAL(num) :: cz1, cz2, cz3

  ! 推进电场：
  !   E^{n+1/2} = E^{n-1/2} + c^2*(dt/2)*curl(B^n) - (dt/2)/epsilon0 * J^n
  ! 其中 curl(B) 的离散依赖 field_order（2/4/6 阶）。
  ! 若启用 CPML，则 curl 系数要除以对应方向的 kappa；并在末尾推进 CPML 辅助电流。

    IF (cpml_boundaries) THEN
      IF (field_order == 2) THEN
        DO iz = 0, nz
          cz1 = cnz / cpml_kappa_ez(iz)
          DO iy = 0, ny
            cy1 = cny / cpml_kappa_ey(iy)
            DO ix = 0, nx
              cx1 = cnx / cpml_kappa_ex(ix)

              ex(ix, iy, iz) = ex(ix, iy, iz) &
                  + cy1 * (bz(ix  , iy  , iz  ) - bz(ix  , iy-1, iz  )) &
                  - cz1 * (by(ix  , iy  , iz  ) - by(ix  , iy  , iz-1)) &
                  - fac * jx(ix, iy, iz)

              ey(ix, iy, iz) = ey(ix, iy, iz) &
                  + cz1 * (bx(ix  , iy  , iz  ) - bx(ix  , iy  , iz-1)) &
                  - cx1 * (bz(ix  , iy  , iz  ) - bz(ix-1, iy  , iz  )) &
                  - fac * jy(ix, iy, iz)

              ez(ix, iy, iz) = ez(ix, iy, iz) &
                  + cx1 * (by(ix  , iy  , iz  ) - by(ix-1, iy  , iz  )) &
                  - cy1 * (bx(ix  , iy  , iz  ) - bx(ix  , iy-1, iz  )) &
                  - fac * jz(ix, iy, iz)
            END DO
          END DO
        END DO
      ELSE IF (field_order == 4) THEN
        c1 = 9.0_num / 8.0_num
        c2 = -1.0_num / 24.0_num

        DO iz = 0, nz
          cpml_z = cnz / cpml_kappa_ez(iz)
          cz1 = c1 * cpml_z
          cz2 = c2 * cpml_z
          DO iy = 0, ny
            cpml_y = cny / cpml_kappa_ey(iy)
            cy1 = c1 * cpml_y
            cy2 = c2 * cpml_y
            DO ix = 0, nx
              cpml_x = cnx / cpml_kappa_ex(ix)
              cx1 = c1 * cpml_x
              cx2 = c2 * cpml_x

              ex(ix, iy, iz) = ex(ix, iy, iz) &
                  + cy1 * (bz(ix  , iy  , iz  ) - bz(ix  , iy-1, iz  )) &
                  + cy2 * (bz(ix  , iy+1, iz  ) - bz(ix  , iy-2, iz  )) &
                  - cz1 * (by(ix  , iy  , iz  ) - by(ix  , iy  , iz-1)) &
                  - cz2 * (by(ix  , iy  , iz+1) - by(ix  , iy  , iz-2)) &
                  - fac * jx(ix, iy, iz)

              ey(ix, iy, iz) = ey(ix, iy, iz) &
                  + cz1 * (bx(ix  , iy  , iz  ) - bx(ix  , iy  , iz-1)) &
                  + cz2 * (bx(ix  , iy  , iz+1) - bx(ix  , iy  , iz-2)) &
                  - cx1 * (bz(ix  , iy  , iz  ) - bz(ix-1, iy  , iz  )) &
                  - cx2 * (bz(ix+1, iy  , iz  ) - bz(ix-2, iy  , iz  )) &
                  - fac * jy(ix, iy, iz)

              ez(ix, iy, iz) = ez(ix, iy, iz) &
                  + cx1 * (by(ix  , iy  , iz  ) - by(ix-1, iy  , iz  )) &
                  + cx2 * (by(ix+1, iy  , iz  ) - by(ix-2, iy  , iz  )) &
                  - cy1 * (bx(ix  , iy  , iz  ) - bx(ix  , iy-1, iz  )) &
                  - cy2 * (bx(ix  , iy+1, iz  ) - bx(ix  , iy-2, iz  )) &
                  - fac * jz(ix, iy, iz)
            END DO
          END DO
        END DO
      ELSE
        c1 = 75.0_num / 64.0_num
        c2 = -25.0_num / 384.0_num
        c3 = 3.0_num / 640.0_num

        DO iz = 0, nz
          cpml_z = cnz / cpml_kappa_ez(iz)
          cz1 = c1 * cpml_z
          cz2 = c2 * cpml_z
          cz3 = c3 * cpml_z
          DO iy = 0, ny
            cpml_y = cny / cpml_kappa_ey(iy)
            cy1 = c1 * cpml_y
            cy2 = c2 * cpml_y
            cy3 = c3 * cpml_y
            DO ix = 0, nx
              cpml_x = cnx / cpml_kappa_ex(ix)
              cx1 = c1 * cpml_x
              cx2 = c2 * cpml_x
              cx3 = c3 * cpml_x

              ex(ix, iy, iz) = ex(ix, iy, iz) &
                  + cy1 * (bz(ix  , iy  , iz  ) - bz(ix  , iy-1, iz  )) &
                  + cy2 * (bz(ix  , iy+1, iz  ) - bz(ix  , iy-2, iz  )) &
                  + cy3 * (bz(ix  , iy+2, iz  ) - bz(ix  , iy-3, iz  )) &
                  - cz1 * (by(ix  , iy  , iz  ) - by(ix  , iy  , iz-1)) &
                  - cz2 * (by(ix  , iy  , iz+1) - by(ix  , iy  , iz-2)) &
                  - cz3 * (by(ix  , iy  , iz+2) - by(ix  , iy  , iz-3)) &
                  - fac * jx(ix, iy, iz)

              ey(ix, iy, iz) = ey(ix, iy, iz) &
                  + cz1 * (bx(ix  , iy  , iz  ) - bx(ix  , iy  , iz-1)) &
                  + cz2 * (bx(ix  , iy  , iz+1) - bx(ix  , iy  , iz-2)) &
                  + cz3 * (bx(ix  , iy  , iz+2) - bx(ix  , iy  , iz-3)) &
                  - cx1 * (bz(ix  , iy  , iz  ) - bz(ix-1, iy  , iz  )) &
                  - cx2 * (bz(ix+1, iy  , iz  ) - bz(ix-2, iy  , iz  )) &
                  - cx3 * (bz(ix+2, iy  , iz  ) - bz(ix-3, iy  , iz  )) &
                  - fac * jy(ix, iy, iz)

              ez(ix, iy, iz) = ez(ix, iy, iz) &
                  + cx1 * (by(ix  , iy  , iz  ) - by(ix-1, iy  , iz  )) &
                  + cx2 * (by(ix+1, iy  , iz  ) - by(ix-2, iy  , iz  )) &
                  + cx3 * (by(ix+2, iy  , iz  ) - by(ix-3, iy  , iz  )) &
                  - cy1 * (bx(ix  , iy  , iz  ) - bx(ix  , iy-1, iz  )) &
                  - cy2 * (bx(ix  , iy+1, iz  ) - bx(ix  , iy-2, iz  )) &
                  - cy3 * (bx(ix  , iy+2, iz  ) - bx(ix  , iy-3, iz  )) &
                  - fac * jz(ix, iy, iz)
            END DO
          END DO
        END DO
      END IF

  ! CPML 更新：推进用于吸收边界的辅助变量/电流（E 相关）
      CALL cpml_advance_e_currents(hdt)
    ELSE
      IF (field_order == 2) THEN
        cx1 = cnx
        cy1 = cny
        cz1 = cnz

        DO iz = 0, nz
          DO iy = 0, ny
            DO ix = 0, nx
              ex(ix, iy, iz) = ex(ix, iy, iz) &
                  + cy1 * (bz(ix  , iy  , iz  ) - bz(ix  , iy-1, iz  )) &
                  - cz1 * (by(ix  , iy  , iz  ) - by(ix  , iy  , iz-1)) &
                  - fac * jx(ix, iy, iz)

              ey(ix, iy, iz) = ey(ix, iy, iz) &
                  + cz1 * (bx(ix  , iy  , iz  ) - bx(ix  , iy  , iz-1)) &
                  - cx1 * (bz(ix  , iy  , iz  ) - bz(ix-1, iy  , iz  )) &
                  - fac * jy(ix, iy, iz)

              ez(ix, iy, iz) = ez(ix, iy, iz) &
                  + cx1 * (by(ix  , iy  , iz  ) - by(ix-1, iy  , iz  )) &
                  - cy1 * (bx(ix  , iy  , iz  ) - bx(ix  , iy-1, iz  )) &
                  - fac * jz(ix, iy, iz)
            END DO
          END DO
        END DO
      ELSE IF (field_order == 4) THEN
        c1 = 9.0_num / 8.0_num
        c2 = -1.0_num / 24.0_num

        cx1 = c1 * cnx
        cx2 = c2 * cnx
        cy1 = c1 * cny
        cy2 = c2 * cny
        cz1 = c1 * cnz
        cz2 = c2 * cnz

        DO iz = 0, nz
          DO iy = 0, ny
            DO ix = 0, nx
              ex(ix, iy, iz) = ex(ix, iy, iz) &
                  + cy1 * (bz(ix  , iy  , iz  ) - bz(ix  , iy-1, iz  )) &
                  + cy2 * (bz(ix  , iy+1, iz  ) - bz(ix  , iy-2, iz  )) &
                  - cz1 * (by(ix  , iy  , iz  ) - by(ix  , iy  , iz-1)) &
                  - cz2 * (by(ix  , iy  , iz+1) - by(ix  , iy  , iz-2)) &
                  - fac * jx(ix, iy, iz)

              ey(ix, iy, iz) = ey(ix, iy, iz) &
                  + cz1 * (bx(ix  , iy  , iz  ) - bx(ix  , iy  , iz-1)) &
                  + cz2 * (bx(ix  , iy  , iz+1) - bx(ix  , iy  , iz-2)) &
                  - cx1 * (bz(ix  , iy  , iz  ) - bz(ix-1, iy  , iz  )) &
                  - cx2 * (bz(ix+1, iy  , iz  ) - bz(ix-2, iy  , iz  )) &
                  - fac * jy(ix, iy, iz)

              ez(ix, iy, iz) = ez(ix, iy, iz) &
                  + cx1 * (by(ix  , iy  , iz  ) - by(ix-1, iy  , iz  )) &
                  + cx2 * (by(ix+1, iy  , iz  ) - by(ix-2, iy  , iz  )) &
                  - cy1 * (bx(ix  , iy  , iz  ) - bx(ix  , iy-1, iz  )) &
                  - cy2 * (bx(ix  , iy+1, iz  ) - bx(ix  , iy-2, iz  )) &
                  - fac * jz(ix, iy, iz)
            END DO
          END DO
        END DO
      ELSE
        c1 = 75.0_num / 64.0_num
        c2 = -25.0_num / 384.0_num
        c3 = 3.0_num / 640.0_num

        cx1 = c1 * cnx
        cx2 = c2 * cnx
        cx3 = c3 * cnx
        cy1 = c1 * cny
        cy2 = c2 * cny
        cy3 = c3 * cny
        cz1 = c1 * cnz
        cz2 = c2 * cnz
        cz3 = c3 * cnz

        DO iz = 0, nz
          DO iy = 0, ny
            DO ix = 0, nx
              ex(ix, iy, iz) = ex(ix, iy, iz) &
                  + cy1 * (bz(ix  , iy  , iz  ) - bz(ix  , iy-1, iz  )) &
                  + cy2 * (bz(ix  , iy+1, iz  ) - bz(ix  , iy-2, iz  )) &
                  + cy3 * (bz(ix  , iy+2, iz  ) - bz(ix  , iy-3, iz  )) &
                  - cz1 * (by(ix  , iy  , iz  ) - by(ix  , iy  , iz-1)) &
                  - cz2 * (by(ix  , iy  , iz+1) - by(ix  , iy  , iz-2)) &
                  - cz3 * (by(ix  , iy  , iz+2) - by(ix  , iy  , iz-3)) &
                  - fac * jx(ix, iy, iz)

              ey(ix, iy, iz) = ey(ix, iy, iz) &
                  + cz1 * (bx(ix  , iy  , iz  ) - bx(ix  , iy  , iz-1)) &
                  + cz2 * (bx(ix  , iy  , iz+1) - bx(ix  , iy  , iz-2)) &
                  + cz3 * (bx(ix  , iy  , iz+2) - bx(ix  , iy  , iz-3)) &
                  - cx1 * (bz(ix  , iy  , iz  ) - bz(ix-1, iy  , iz  )) &
                  - cx2 * (bz(ix+1, iy  , iz  ) - bz(ix-2, iy  , iz  )) &
                  - cx3 * (bz(ix+2, iy  , iz  ) - bz(ix-3, iy  , iz  )) &
                  - fac * jy(ix, iy, iz)

              ez(ix, iy, iz) = ez(ix, iy, iz) &
                  + cx1 * (by(ix  , iy  , iz  ) - by(ix-1, iy  , iz  )) &
                  + cx2 * (by(ix+1, iy  , iz  ) - by(ix-2, iy  , iz  )) &
                  + cx3 * (by(ix+2, iy  , iz  ) - by(ix-3, iy  , iz  )) &
                  - cy1 * (bx(ix  , iy  , iz  ) - bx(ix  , iy-1, iz  )) &
                  - cy2 * (bx(ix  , iy+1, iz  ) - bx(ix  , iy-2, iz  )) &
                  - cy3 * (bx(ix  , iy+2, iz  ) - bx(ix  , iy-3, iz  )) &
                  - fac * jz(ix, iy, iz)
            END DO
          END DO
        END DO
      END IF
    END IF

  END SUBROUTINE update_e_field



  SUBROUTINE update_b_field

    INTEGER :: ix, iy, iz
    REAL(num) :: cpml_x, cpml_y, cpml_z
    REAL(num) :: c1, c2, c3
    REAL(num) :: cx1, cx2, cx3
    REAL(num) :: cy1, cy2, cy3
    REAL(num) :: cz1, cz2, cz3

  ! 推进磁场：
  !   B^{n+1} = B^{n} - (dt/2)*curl(E^{n+1/2})
  ! field_order 控制 curl 离散阶数。
  ! 当 maxwell_solver==Yee 时使用标准 stencil；否则使用 set_maxwell_solver
  ! 预先计算的 (alpha/beta/gamma/delta) 组合构造改进的 curl。

    IF (cpml_boundaries) THEN
      IF (field_order == 2) THEN
        IF (maxwell_solver == c_maxwell_solver_yee) THEN
          DO iz = 0, nz
            cz1 = hdtz / cpml_kappa_bz(iz)
            DO iy = 0, ny
              cy1 = hdty / cpml_kappa_by(iy)
              DO ix = 0, nx
                cx1 = hdtx / cpml_kappa_bx(ix)

                bx(ix, iy, iz) = bx(ix, iy, iz) &
                    - cy1 * (ez(ix  , iy+1, iz  ) - ez(ix  , iy  , iz  )) &
                    + cz1 * (ey(ix  , iy  , iz+1) - ey(ix  , iy  , iz  ))

                by(ix, iy, iz) = by(ix, iy, iz) &
                    - cz1 * (ex(ix  , iy  , iz+1) - ex(ix  , iy  , iz  )) &
                    + cx1 * (ez(ix+1, iy  , iz  ) - ez(ix  , iy  , iz  ))

                bz(ix, iy, iz) = bz(ix, iy, iz) &
                    - cx1 * (ey(ix+1, iy  , iz  ) - ey(ix  , iy  , iz  )) &
                    + cy1 * (ex(ix  , iy+1, iz  ) - ex(ix  , iy  , iz  ))
              END DO
            END DO
          END DO
        ELSE
          DO iz = 0, nz
            cz1 = hdtz / cpml_kappa_bz(iz)
            DO iy = 0, ny
              cy1 = hdty / cpml_kappa_by(iy)
              DO ix = 0, nx
                cx1 = hdtx / cpml_kappa_bx(ix)

                bx(ix, iy, iz) = bx(ix, iy, iz) &
                    - cy1 * ( &
                       alphay * (ez(ix  , iy+1, iz  ) - ez(ix  , iy  , iz  ))  &
                     + betayx * (ez(ix+1, iy+1, iz  ) - ez(ix+1, iy  , iz  )   &
                               + ez(ix-1, iy+1, iz  ) - ez(ix-1, iy  , iz  ))  &
                     + betayz * (ez(ix  , iy+1, iz+1) - ez(ix  , iy  , iz+1)   &
                               + ez(ix  , iy+1, iz-1) - ez(ix  , iy  , iz-1))  &
                     + gammay * (ez(ix+1, iy+1, iz-1) - ez(ix+1, iy  , iz-1)   &
                               + ez(ix-1, iy+1, iz-1) - ez(ix-1, iy  , iz-1)   &
                               + ez(ix+1, iy+1, iz+1) - ez(ix+1, iy  , iz+1)   &
                               + ez(ix-1, iy+1, iz+1) - ez(ix-1, iy  , iz+1))  &
                     + deltay * (ez(ix  , iy+2, iz  ) - ez(ix  , iy-1, iz  ))) &
                    + cz1 * ( &
                       alphaz * (ey(ix  , iy  , iz+1) - ey(ix  , iy  , iz  ))  &
                     + betazx * (ey(ix+1, iy  , iz+1) - ey(ix+1, iy  , iz  )   &
                               + ey(ix-1, iy  , iz+1) - ey(ix-1, iy  , iz  ))  &
                     + betazy * (ey(ix  , iy+1, iz+1) - ey(ix  , iy+1, iz  )   &
                               + ey(ix  , iy-1, iz+1) - ey(ix  , iy-1, iz  ))  &
                     + gammaz * (ey(ix+1, iy-1, iz+1) - ey(ix+1, iy-1, iz  )   &
                               + ey(ix-1, iy-1, iz+1) - ey(ix-1, iy-1, iz  )   &
                               + ey(ix+1, iy+1, iz+1) - ey(ix+1, iy+1, iz  )   &
                               + ey(ix-1, iy+1, iz+1) - ey(ix-1, iy+1, iz  ))  &
                     + deltaz * (ey(ix  , iy  , iz+2) - ey(ix  , iy  , iz-1)))

                by(ix, iy, iz) = by(ix, iy, iz) &
                    - cz1 * ( &
                       alphaz * (ex(ix  , iy  , iz+1) - ex(ix  , iy  , iz  ))  &
                     + betazx * (ex(ix+1, iy  , iz+1) - ex(ix+1, iy  , iz  )   &
                               + ex(ix-1, iy  , iz+1) - ex(ix-1, iy  , iz  ))  &
                     + betazy * (ex(ix  , iy+1, iz+1) - ex(ix  , iy+1, iz  )   &
                               + ex(ix  , iy-1, iz+1) - ex(ix  , iy-1, iz  ))  &
                     + gammaz * (ex(ix+1, iy-1, iz+1) - ex(ix+1, iy-1, iz  )   &
                               + ex(ix-1, iy-1, iz+1) - ex(ix-1, iy-1, iz  )   &
                               + ex(ix+1, iy+1, iz+1) - ex(ix+1, iy+1, iz  )   &
                               + ex(ix-1, iy+1, iz+1) - ex(ix-1, iy+1, iz  ))  &
                     + deltaz * (ex(ix  , iy  , iz+2) - ex(ix  , iy  , iz-1))) &
                    + cx1 * ( &
                       alphax * (ez(ix+1, iy  , iz  ) - ez(ix  , iy  , iz  ))  &
                     + betaxy * (ez(ix+1, iy+1, iz  ) - ez(ix  , iy+1, iz  )   &
                               + ez(ix+1, iy-1, iz  ) - ez(ix  , iy-1, iz  ))  &
                     + betaxz * (ez(ix+1, iy  , iz+1) - ez(ix  , iy  , iz+1)   &
                               + ez(ix+1, iy  , iz-1) - ez(ix  , iy  , iz-1))  &
                     + gammax * (ez(ix+1, iy+1, iz-1) - ez(ix  , iy+1, iz-1)   &
                               + ez(ix+1, iy-1, iz-1) - ez(ix  , iy-1, iz-1)   &
                               + ez(ix+1, iy+1, iz+1) - ez(ix  , iy+1, iz+1)   &
                               + ez(ix+1, iy-1, iz+1) - ez(ix  , iy-1, iz+1))  &
                     + deltax * (ez(ix+2, iy  , iz  ) - ez(ix-1, iy  , iz  )))

                bz(ix, iy, iz) = bz(ix, iy, iz) &
                    - cx1 * ( &
                       alphax * (ey(ix+1, iy  , iz  ) - ey(ix  , iy  , iz  ))  &
                     + betaxy * (ey(ix+1, iy+1, iz  ) - ey(ix  , iy+1, iz  )   &
                               + ey(ix+1, iy-1, iz  ) - ey(ix  , iy-1, iz  ))  &
                     + betaxz * (ey(ix+1, iy  , iz+1) - ey(ix  , iy  , iz+1)   &
                               + ey(ix+1, iy  , iz-1) - ey(ix  , iy  , iz-1))  &
                     + gammax * (ey(ix+1, iy+1, iz-1) - ey(ix  , iy+1, iz-1)   &
                               + ey(ix+1, iy-1, iz-1) - ey(ix  , iy-1, iz-1)   &
                               + ey(ix+1, iy+1, iz+1) - ey(ix  , iy+1, iz+1)   &
                               + ey(ix+1, iy-1, iz+1) - ey(ix  , iy-1, iz+1))  &
                     + deltax * (ey(ix+2, iy  , iz  ) - ey(ix-1, iy  , iz  ))) &
                    + cy1 * ( &
                       alphay * (ex(ix  , iy+1, iz  ) - ex(ix  , iy  , iz  ))  &
                     + betayx * (ex(ix+1, iy+1, iz  ) - ex(ix+1, iy  , iz  )   &
                               + ex(ix-1, iy+1, iz  ) - ex(ix-1, iy  , iz  ))  &
                     + betayz * (ex(ix  , iy+1, iz+1) - ex(ix  , iy  , iz+1)   &
                               + ex(ix  , iy+1, iz-1) - ex(ix  , iy  , iz-1))  &
                     + gammay * (ex(ix+1, iy+1, iz-1) - ex(ix+1, iy  , iz-1)   &
                               + ex(ix-1, iy+1, iz-1) - ex(ix-1, iy  , iz-1)   &
                               + ex(ix+1, iy+1, iz+1) - ex(ix+1, iy  , iz+1)   &
                               + ex(ix-1, iy+1, iz+1) - ex(ix-1, iy  , iz+1))  &
                     + deltay * (ex(ix  , iy+2, iz  ) - ex(ix  , iy-1, iz  )))
              END DO
            END DO
          END DO
        END IF
      ELSE IF (field_order == 4) THEN
        c1 = 9.0_num / 8.0_num
        c2 = -1.0_num / 24.0_num

        DO iz = 0, nz
          cpml_z = hdtz / cpml_kappa_bz(iz)
          cz1 = c1 * cpml_z
          cz2 = c2 * cpml_z
          DO iy = 0, ny
            cpml_y = hdty / cpml_kappa_by(iy)
            cy1 = c1 * cpml_y
            cy2 = c2 * cpml_y
            DO ix = 0, nx
              cpml_x = hdtx / cpml_kappa_bx(ix)
              cx1 = c1 * cpml_x
              cx2 = c2 * cpml_x

              bx(ix, iy, iz) = bx(ix, iy, iz) &
                  - cy1 * (ez(ix  , iy+1, iz  ) - ez(ix  , iy  , iz  )) &
                  - cy2 * (ez(ix  , iy+2, iz  ) - ez(ix  , iy-1, iz  )) &
                  + cz1 * (ey(ix  , iy  , iz+1) - ey(ix  , iy  , iz  )) &
                  + cz2 * (ey(ix  , iy  , iz+2) - ey(ix  , iy  , iz-1))

              by(ix, iy, iz) = by(ix, iy, iz) &
                  - cz1 * (ex(ix  , iy  , iz+1) - ex(ix  , iy  , iz  )) &
                  - cz2 * (ex(ix  , iy  , iz+2) - ex(ix  , iy  , iz-1)) &
                  + cx1 * (ez(ix+1, iy  , iz  ) - ez(ix  , iy  , iz  )) &
                  + cx2 * (ez(ix+2, iy  , iz  ) - ez(ix-1, iy  , iz  ))

              bz(ix, iy, iz) = bz(ix, iy, iz) &
                  - cx1 * (ey(ix+1, iy  , iz  ) - ey(ix  , iy  , iz  )) &
                  - cx2 * (ey(ix+2, iy  , iz  ) - ey(ix-1, iy  , iz  )) &
                  + cy1 * (ex(ix  , iy+1, iz  ) - ex(ix  , iy  , iz  )) &
                  + cy2 * (ex(ix  , iy+2, iz  ) - ex(ix  , iy-1, iz  ))
            END DO
          END DO
        END DO
      ELSE
        c1 = 75.0_num / 64.0_num
        c2 = -25.0_num / 384.0_num
        c3 = 3.0_num / 640.0_num

        DO iz = 0, nz
          cpml_z = hdtz / cpml_kappa_bz(iz)
          cz1 = c1 * cpml_z
          cz2 = c2 * cpml_z
          cz3 = c3 * cpml_z
          DO iy = 0, ny
            cpml_y = hdty / cpml_kappa_by(iy)
            cy1 = c1 * cpml_y
            cy2 = c2 * cpml_y
            cy3 = c3 * cpml_y
            DO ix = 0, nx
              cpml_x = hdtx / cpml_kappa_bx(ix)
              cx1 = c1 * cpml_x
              cx2 = c2 * cpml_x
              cx3 = c3 * cpml_x

              bx(ix, iy, iz) = bx(ix, iy, iz) &
                  - cy1 * (ez(ix  , iy+1, iz  ) - ez(ix  , iy  , iz  )) &
                  - cy2 * (ez(ix  , iy+2, iz  ) - ez(ix  , iy-1, iz  )) &
                  - cy3 * (ez(ix  , iy+3, iz  ) - ez(ix  , iy-2, iz  )) &
                  + cz1 * (ey(ix  , iy  , iz+1) - ey(ix  , iy  , iz  )) &
                  + cz2 * (ey(ix  , iy  , iz+2) - ey(ix  , iy  , iz-1)) &
                  + cz3 * (ey(ix  , iy  , iz+3) - ey(ix  , iy  , iz-2))

              by(ix, iy, iz) = by(ix, iy, iz) &
                  - cz1 * (ex(ix  , iy  , iz+1) - ex(ix  , iy  , iz  )) &
                  - cz2 * (ex(ix  , iy  , iz+2) - ex(ix  , iy  , iz-1)) &
                  - cz3 * (ex(ix  , iy  , iz+3) - ex(ix  , iy  , iz-2)) &
                  + cx1 * (ez(ix+1, iy  , iz  ) - ez(ix  , iy  , iz  )) &
                  + cx2 * (ez(ix+2, iy  , iz  ) - ez(ix-1, iy  , iz  )) &
                  + cx3 * (ez(ix+3, iy  , iz  ) - ez(ix-2, iy  , iz  ))

              bz(ix, iy, iz) = bz(ix, iy, iz) &
                  - cx1 * (ey(ix+1, iy  , iz  ) - ey(ix  , iy  , iz  )) &
                  - cx2 * (ey(ix+2, iy  , iz  ) - ey(ix-1, iy  , iz  )) &
                  - cx3 * (ey(ix+3, iy  , iz  ) - ey(ix-2, iy  , iz  )) &
                  + cy1 * (ex(ix  , iy+1, iz  ) - ex(ix  , iy  , iz  )) &
                  + cy2 * (ex(ix  , iy+2, iz  ) - ex(ix  , iy-1, iz  )) &
                  + cy3 * (ex(ix  , iy+3, iz  ) - ex(ix  , iy-2, iz  ))
            END DO
          END DO
        END DO
      END IF

  ! CPML 更新：推进用于吸收边界的辅助变量/电流（B 相关）
      CALL cpml_advance_b_currents(hdt)
    ELSE
      IF (field_order == 2) THEN
        cx1 = hdtx
        cy1 = hdty
        cz1 = hdtz

        IF (maxwell_solver == c_maxwell_solver_yee) THEN
          DO iz = 0, nz
            DO iy = 0, ny
              DO ix = 0, nx
                bx(ix, iy, iz) = bx(ix, iy, iz) &
                    - cy1 * (ez(ix  , iy+1, iz  ) - ez(ix  , iy  , iz  )) &
                    + cz1 * (ey(ix  , iy  , iz+1) - ey(ix  , iy  , iz  ))

                by(ix, iy, iz) = by(ix, iy, iz) &
                    - cz1 * (ex(ix  , iy  , iz+1) - ex(ix  , iy  , iz  )) &
                    + cx1 * (ez(ix+1, iy  , iz  ) - ez(ix  , iy  , iz  ))

                bz(ix, iy, iz) = bz(ix, iy, iz) &
                    - cx1 * (ey(ix+1, iy  , iz  ) - ey(ix  , iy  , iz  )) &
                    + cy1 * (ex(ix  , iy+1, iz  ) - ex(ix  , iy  , iz  ))
              END DO
            END DO
          END DO
        ELSE
          DO iz = 0, nz
            DO iy = 0, ny
              DO ix = 0, nx
                bx(ix, iy, iz) = bx(ix, iy, iz) &
                    - cy1 * ( &
                       alphay * (ez(ix  , iy+1, iz  ) - ez(ix  , iy  , iz  ))  &
                     + betayx * (ez(ix+1, iy+1, iz  ) - ez(ix+1, iy  , iz  )   &
                               + ez(ix-1, iy+1, iz  ) - ez(ix-1, iy  , iz  ))  &
                     + betayz * (ez(ix  , iy+1, iz+1) - ez(ix  , iy  , iz+1)   &
                               + ez(ix  , iy+1, iz-1) - ez(ix  , iy  , iz-1))  &
                     + gammay * (ez(ix+1, iy+1, iz-1) - ez(ix+1, iy  , iz-1)   &
                               + ez(ix-1, iy+1, iz-1) - ez(ix-1, iy  , iz-1)   &
                               + ez(ix+1, iy+1, iz+1) - ez(ix+1, iy  , iz+1)   &
                               + ez(ix-1, iy+1, iz+1) - ez(ix-1, iy  , iz+1))  &
                     + deltay * (ez(ix  , iy+2, iz  ) - ez(ix  , iy-1, iz  ))) &
                    + cz1 * ( &
                       alphaz * (ey(ix  , iy  , iz+1) - ey(ix  , iy  , iz  ))  &
                     + betazx * (ey(ix+1, iy  , iz+1) - ey(ix+1, iy  , iz  )   &
                               + ey(ix-1, iy  , iz+1) - ey(ix-1, iy  , iz  ))  &
                     + betazy * (ey(ix  , iy+1, iz+1) - ey(ix  , iy+1, iz  )   &
                               + ey(ix  , iy-1, iz+1) - ey(ix  , iy-1, iz  ))  &
                     + gammaz * (ey(ix+1, iy-1, iz+1) - ey(ix+1, iy-1, iz  )   &
                               + ey(ix-1, iy-1, iz+1) - ey(ix-1, iy-1, iz  )   &
                               + ey(ix+1, iy+1, iz+1) - ey(ix+1, iy+1, iz  )   &
                               + ey(ix-1, iy+1, iz+1) - ey(ix-1, iy+1, iz  ))  &
                     + deltaz * (ey(ix  , iy  , iz+2) - ey(ix  , iy  , iz-1)))

                by(ix, iy, iz) = by(ix, iy, iz) &
                    - cz1 * ( &
                       alphaz * (ex(ix  , iy  , iz+1) - ex(ix  , iy  , iz  ))  &
                     + betazx * (ex(ix+1, iy  , iz+1) - ex(ix+1, iy  , iz  )   &
                               + ex(ix-1, iy  , iz+1) - ex(ix-1, iy  , iz  ))  &
                     + betazy * (ex(ix  , iy+1, iz+1) - ex(ix  , iy+1, iz  )   &
                               + ex(ix  , iy-1, iz+1) - ex(ix  , iy-1, iz  ))  &
                     + gammaz * (ex(ix+1, iy-1, iz+1) - ex(ix+1, iy-1, iz  )   &
                               + ex(ix-1, iy-1, iz+1) - ex(ix-1, iy-1, iz  )   &
                               + ex(ix+1, iy+1, iz+1) - ex(ix+1, iy+1, iz  )   &
                               + ex(ix-1, iy+1, iz+1) - ex(ix-1, iy+1, iz  ))  &
                     + deltaz * (ex(ix  , iy  , iz+2) - ex(ix  , iy  , iz-1))) &
                    + cx1 * ( &
                       alphax * (ez(ix+1, iy  , iz  ) - ez(ix  , iy  , iz  ))  &
                     + betaxy * (ez(ix+1, iy+1, iz  ) - ez(ix  , iy+1, iz  )   &
                               + ez(ix+1, iy-1, iz  ) - ez(ix  , iy-1, iz  ))  &
                     + betaxz * (ez(ix+1, iy  , iz+1) - ez(ix  , iy  , iz+1)   &
                               + ez(ix+1, iy  , iz-1) - ez(ix  , iy  , iz-1))  &
                     + gammax * (ez(ix+1, iy+1, iz-1) - ez(ix  , iy+1, iz-1)   &
                               + ez(ix+1, iy-1, iz-1) - ez(ix  , iy-1, iz-1)   &
                               + ez(ix+1, iy+1, iz+1) - ez(ix  , iy+1, iz+1)   &
                               + ez(ix+1, iy-1, iz+1) - ez(ix  , iy-1, iz+1))  &
                     + deltax * (ez(ix+2, iy  , iz  ) - ez(ix-1, iy  , iz  )))

                bz(ix, iy, iz) = bz(ix, iy, iz) &
                    - cx1 * ( &
                       alphax * (ey(ix+1, iy  , iz  ) - ey(ix  , iy  , iz  ))  &
                     + betaxy * (ey(ix+1, iy+1, iz  ) - ey(ix  , iy+1, iz  )   &
                               + ey(ix+1, iy-1, iz  ) - ey(ix  , iy-1, iz  ))  &
                     + betaxz * (ey(ix+1, iy  , iz+1) - ey(ix  , iy  , iz+1)   &
                               + ey(ix+1, iy  , iz-1) - ey(ix  , iy  , iz-1))  &
                     + gammax * (ey(ix+1, iy+1, iz-1) - ey(ix  , iy+1, iz-1)   &
                               + ey(ix+1, iy-1, iz-1) - ey(ix  , iy-1, iz-1)   &
                               + ey(ix+1, iy+1, iz+1) - ey(ix  , iy+1, iz+1)   &
                               + ey(ix+1, iy-1, iz+1) - ey(ix  , iy-1, iz+1))  &
                     + deltax * (ey(ix+2, iy  , iz  ) - ey(ix-1, iy  , iz  ))) &
                    + cy1 * ( &
                       alphay * (ex(ix  , iy+1, iz  ) - ex(ix  , iy  , iz  ))  &
                     + betayx * (ex(ix+1, iy+1, iz  ) - ex(ix+1, iy  , iz  )   &
                               + ex(ix-1, iy+1, iz  ) - ex(ix-1, iy  , iz  ))  &
                     + betayz * (ex(ix  , iy+1, iz+1) - ex(ix  , iy  , iz+1)   &
                               + ex(ix  , iy+1, iz-1) - ex(ix  , iy  , iz-1))  &
                     + gammay * (ex(ix+1, iy+1, iz-1) - ex(ix+1, iy  , iz-1)   &
                               + ex(ix-1, iy+1, iz-1) - ex(ix-1, iy  , iz-1)   &
                               + ex(ix+1, iy+1, iz+1) - ex(ix+1, iy  , iz+1)   &
                               + ex(ix-1, iy+1, iz+1) - ex(ix-1, iy  , iz+1))  &
                     + deltay * (ex(ix  , iy+2, iz  ) - ex(ix  , iy-1, iz  )))
              END DO
            END DO
          END DO
        END IF
      ELSE IF (field_order == 4) THEN
        c1 = 9.0_num / 8.0_num
        c2 = -1.0_num / 24.0_num

        cx1 = c1 * hdtx
        cx2 = c2 * hdtx
        cy1 = c1 * hdty
        cy2 = c2 * hdty
        cz1 = c1 * hdtz
        cz2 = c2 * hdtz

        DO iz = 0, nz
          DO iy = 0, ny
            DO ix = 0, nx
              bx(ix, iy, iz) = bx(ix, iy, iz) &
                  - cy1 * (ez(ix  , iy+1, iz  ) - ez(ix  , iy  , iz  )) &
                  - cy2 * (ez(ix  , iy+2, iz  ) - ez(ix  , iy-1, iz  )) &
                  + cz1 * (ey(ix  , iy  , iz+1) - ey(ix  , iy  , iz  )) &
                  + cz2 * (ey(ix  , iy  , iz+2) - ey(ix  , iy  , iz-1))

              by(ix, iy, iz) = by(ix, iy, iz) &
                  - cz1 * (ex(ix  , iy  , iz+1) - ex(ix  , iy  , iz  )) &
                  - cz2 * (ex(ix  , iy  , iz+2) - ex(ix  , iy  , iz-1)) &
                  + cx1 * (ez(ix+1, iy  , iz  ) - ez(ix  , iy  , iz  )) &
                  + cx2 * (ez(ix+2, iy  , iz  ) - ez(ix-1, iy  , iz  ))

              bz(ix, iy, iz) = bz(ix, iy, iz) &
                  - cx1 * (ey(ix+1, iy  , iz  ) - ey(ix  , iy  , iz  )) &
                  - cx2 * (ey(ix+2, iy  , iz  ) - ey(ix-1, iy  , iz  )) &
                  + cy1 * (ex(ix  , iy+1, iz  ) - ex(ix  , iy  , iz  )) &
                  + cy2 * (ex(ix  , iy+2, iz  ) - ex(ix  , iy-1, iz  ))
            END DO
          END DO
        END DO
      ELSE
        c1 = 75.0_num / 64.0_num
        c2 = -25.0_num / 384.0_num
        c3 = 3.0_num / 640.0_num

        cx1 = c1 * hdtx
        cx2 = c2 * hdtx
        cx3 = c3 * hdtx
        cy1 = c1 * hdty
        cy2 = c2 * hdty
        cy3 = c3 * hdty
        cz1 = c1 * hdtz
        cz2 = c2 * hdtz
        cz3 = c3 * hdtz

        DO iz = 0, nz
          DO iy = 0, ny
            DO ix = 0, nx
              bx(ix, iy, iz) = bx(ix, iy, iz) &
                  - cy1 * (ez(ix  , iy+1, iz  ) - ez(ix  , iy  , iz  )) &
                  - cy2 * (ez(ix  , iy+2, iz  ) - ez(ix  , iy-1, iz  )) &
                  - cy3 * (ez(ix  , iy+3, iz  ) - ez(ix  , iy-2, iz  )) &
                  + cz1 * (ey(ix  , iy  , iz+1) - ey(ix  , iy  , iz  )) &
                  + cz2 * (ey(ix  , iy  , iz+2) - ey(ix  , iy  , iz-1)) &
                  + cz3 * (ey(ix  , iy  , iz+3) - ey(ix  , iy  , iz-2))

              by(ix, iy, iz) = by(ix, iy, iz) &
                  - cz1 * (ex(ix  , iy  , iz+1) - ex(ix  , iy  , iz  )) &
                  - cz2 * (ex(ix  , iy  , iz+2) - ex(ix  , iy  , iz-1)) &
                  - cz3 * (ex(ix  , iy  , iz+3) - ex(ix  , iy  , iz-2)) &
                  + cx1 * (ez(ix+1, iy  , iz  ) - ez(ix  , iy  , iz  )) &
                  + cx2 * (ez(ix+2, iy  , iz  ) - ez(ix-1, iy  , iz  )) &
                  + cx3 * (ez(ix+3, iy  , iz  ) - ez(ix-2, iy  , iz  ))

              bz(ix, iy, iz) = bz(ix, iy, iz) &
                  - cx1 * (ey(ix+1, iy  , iz  ) - ey(ix  , iy  , iz  )) &
                  - cx2 * (ey(ix+2, iy  , iz  ) - ey(ix-1, iy  , iz  )) &
                  - cx3 * (ey(ix+3, iy  , iz  ) - ey(ix-2, iy  , iz  )) &
                  + cy1 * (ex(ix  , iy+1, iz  ) - ex(ix  , iy  , iz  )) &
                  + cy2 * (ex(ix  , iy+2, iz  ) - ex(ix  , iy-1, iz  )) &
                  + cy3 * (ex(ix  , iy+3, iz  ) - ex(ix  , iy-2, iz  ))
            END DO
          END DO
        END DO
      END IF
    END IF

  END SUBROUTINE update_b_field



  SUBROUTINE update_eb_fields_half

  ! 半步推进（leapfrog / Boris 推进中的“场半步”）：
  !  1) E: n -> n+1/2（包含 J）
  !  2) 对 E 做边界条件
  !  3) B: n -> n+1/2（使用 E^(n+1/2)）
  !  4) 对 B 做边界条件
  !  得到 t = n+1/2 的 E/B 后，粒子推进器可以使用它们。

    hdt  = 0.5_num * dt
    hdtx = hdt / dx
    hdty = hdt / dy
    hdtz = hdt / dz

    cnx = hdtx * c**2
    cny = hdty * c**2
    cnz = hdtz * c**2

    fac = hdt / epsilon0

  ! 更新 E 到 t+dt/2
    CALL update_e_field

  ! 已得到 E(t+dt/2)，施加 E 的边界条件
    CALL efield_bcs

  ! 使用 E(t+dt/2) 更新 B 到 t+dt/2
    CALL update_b_field

  ! 已得到 B(t+dt/2)，施加 B 的边界条件
    CALL bfield_bcs(.TRUE.)

    ! Now have E&B fields at t = t+dt/2
    ! Move to particle pusher

  END SUBROUTINE update_eb_fields_half



  SUBROUTINE update_eb_fields_final

  ! 终步推进：粒子完成推进并沉积电流后，补完该时间步的场更新。
  ! 顺序为：
  !   1) B: n+1/2 -> n+1（使用新的 E^(n+1/2)）并做“final”边界
  !   2) E: n+1/2 -> n+1（使用新的 J^(n+1)）并做边界

    hdt  = 0.5_num * dt
    hdtx = hdt / dx
    hdty = hdt / dy
    hdtz = hdt / dz

    cnx = hdtx * c**2
    cny = hdty * c**2
    cnz = hdtz * c**2

    fac = hdt / epsilon0

    CALL update_b_field

    CALL bfield_final_bcs

    CALL update_e_field

    CALL efield_bcs

  END SUBROUTINE update_eb_fields_final

END MODULE fields
