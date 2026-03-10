module gpulse_mod
   use gpulse_constants
   use omp_lib
   implicit none
   public

   type :: Cell
      complex(dpc) :: Ex, Ey, Ez
      complex(dpc) :: Bx, By, Bz
   end type Cell

   type :: SpectralField
      type(Cell), allocatable :: data(:)
      real(dp), allocatable :: omegas(:)
      integer :: N_pixels_x, N_pixels_y, N_freq
      real(dp) :: min_x, max_x, min_y, max_y, dx, dy
   end type SpectralField

   type(SpectralField), save, target :: g_field
   logical, save :: is_initialized = .false.
   character(*), parameter :: CACHE_FILENAME = '.gpulse_cache.bin'

contains

   ! wrapToPi: 将角度映射到 (-PI, PI]，与 C++ run_solver 完全匹配
   function wrapToPi(angle) result(res)
      real(dp), intent(in) :: angle
      real(dp) :: res
      res = mod(angle + PI, 2.0_dp * PI)
      if (res < 0.0_dp) res = res + 2.0_dp * PI
      res = res - PI
   end function wrapToPi

   subroutine gpulse_init()
      type(Geometry) :: geom
      type(Spectrum) :: spec
      integer :: n_x, n_y, total_src, total_pixels, n, s, k
      integer :: ix, iy, idx
      real(dp) :: theta_rot, aa, psi
      real(dp) :: range_val, t0, wn, lam_n, kn, E0_curr
      real(dp) :: factor, freq_weight
      real(dp) :: r, nx, ny, nz, obs_x, obs_y, obs_z
      real(dp) :: key_Rx, key_Ry, key_Rz, R_sq, R_norm
      complex(dpc) :: time_comp_phase, common_phase, gaussian_phase, phase_factor
      complex(dpc) :: E_scalar, Ex, Ey, Ez, Bx, By, Bz
      complex(dpc) :: G, dG, nxe, sEx, sEy, sEz, sBx, sBy, sBz

      real(dp), allocatable :: xf(:), yf(:), cos_2psi(:), sin_2psi(:)
      complex(dpc), allocatable :: vec_NxB_x(:), vec_NxB_y(:), vec_NxB_z(:), vec_N_dot_E(:)

      if (is_initialized) return

      ! 尝试从缓存加载
      if (load_cache()) then
         is_initialized = .true.
         return
      end if

      n_x = num_x
      n_y = num_y
      range_val = RANGE

      g_field%N_pixels_x = n_x
      g_field%N_pixels_y = n_y
      g_field%N_freq = N_FREQ
      allocate(g_field%data(n_x * n_y * N_FREQ))
      allocate(g_field%omegas(N_FREQ))

      g_field%min_x = -range_val
      g_field%max_x = range_val
      g_field%min_y = -range_val
      g_field%max_y = range_val
      g_field%dx = (2.0_dp * range_val) / real(max(1, n_x - 1), dp)
      g_field%dy = (2.0_dp * range_val) / real(max(1, n_y - 1), dp)

      print *, "[Gpulse] Initializing Field Solver (Fortran)..."

      geom = build_geometry()
      spec = calculate_spectrum()
      total_src = geom%total_points

      g_field%omegas = spec%omega_list

      allocate(xf(n_x * n_y), yf(n_x * n_y))
      do ix = 1, n_x
         do iy = 1, n_y
            idx = (ix - 1) * n_y + iy
            xf(idx) = -range_val + real(ix - 1, dp) * g_field%dx
            yf(idx) = -range_val + real(iy - 1, dp) * g_field%dy
         end do
      end do

      ! 偏振映射：与 C++ run_solver 完全匹配的 wrapToPi + 分支逻辑
      allocate(cos_2psi(total_src), sin_2psi(total_src))
      aa = PI / 4.0_dp
      do s = 1, total_src
         if (RADIAL == 1) then
            theta_rot = wrapToPi(geom%theta(s) + PI / 4.0_dp)
            if (theta_rot >= 0.0_dp .and. theta_rot < PI / 2.0_dp) then
               psi = aa           ! PI/4
            else if (theta_rot >= PI / 2.0_dp .and. theta_rot <= PI) then
               psi = 2.0_dp * aa  ! PI/2
            else if (theta_rot >= -PI .and. theta_rot < -PI / 2.0_dp) then
               psi = 3.0_dp * aa  ! 3*PI/4
            else
               psi = 0.0_dp       ! [-PI/2, 0)
            end if
         else
            psi = 0.0_dp
         end if
         cos_2psi(s) = cos(2.0_dp * psi)
         sin_2psi(s) = sin(2.0_dp * psi)
      end do

      allocate(vec_NxB_x(total_src), vec_NxB_y(total_src), vec_NxB_z(total_src))
      allocate(vec_N_dot_E(total_src))

      t0 = 2.0_dp * F0 / C_CONST
      total_pixels = n_x * n_y
      factor = 1.0_dp / (2.0_dp * PI)

      do n = 1, N_FREQ
         wn = spec%omega_list(n)
         lam_n = 2.0_dp * PI * C_CONST / wn
         kn = 2.0_dp * PI * 1.0_dp / lam_n
         E0_curr = spec%E_amp(n)

         time_comp_phase = exp(-I_COMP * wn * t0)
         common_phase = exp(-I_COMP * PI / 2.0_dp) * time_comp_phase

         !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(s, r, gaussian_phase, phase_factor, &
         !$OMP E_scalar, Ex, Ey, Ez, Bx, By, Bz, nx, ny, nz, freq_weight)
         do s = 1, total_src
            r = geom%r(s)
            gaussian_phase = exp(- (r / W0)**2)
            phase_factor = exp(-I_COMP * kn * geom%z(s)) * common_phase
            E_scalar = E0_curr * gaussian_phase * phase_factor

            nx = geom%Nx(s)
            ny = geom%Ny(s)
            nz = geom%Nz(s)

            if (RADIAL == 1) then
               Ex = E_scalar * sin_2psi(s)
               Ey = -E_scalar * cos_2psi(s)
               Ez = (0.0_dp, 0.0_dp)
               Bx = (E_scalar / C_CONST) * cos_2psi(s)
               By = (E_scalar / C_CONST) * sin_2psi(s)
               Bz = (0.0_dp, 0.0_dp)
            else
               Ex = (0.0_dp, 0.0_dp)
               Ey = E_scalar
               Ez = (0.0_dp, 0.0_dp)
               Bx = E_scalar / C_CONST
               By = (0.0_dp, 0.0_dp)
               Bz = (0.0_dp, 0.0_dp)
            end if

            ! Pre-multiply by dS(s) and other constants once
            freq_weight = spec%weights(n) * factor
            vec_NxB_x(s) = (ny * Bz - nz * By) * geom%dS(s) * freq_weight
            vec_NxB_y(s) = (nz * Bx - nx * Bz) * geom%dS(s) * freq_weight
            vec_NxB_z(s) = (nx * By - ny * Bx) * geom%dS(s) * freq_weight
            vec_N_dot_E(s) = (nx * Ex + ny * Ey + nz * Ez) * geom%dS(s) * freq_weight
         end do
         !$OMP END PARALLEL DO

         freq_weight = spec%weights(n)

         if (n == 1) then
            print *, "DEBUG F90: freq_weight=", freq_weight, " factor=", factor
         end if

         !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k, obs_x, obs_y, obs_z, &
         !$OMP sEx, sEy, sEz, sBx, sBy, sBz, s, key_Rx, key_Ry, key_Rz, &
         !$OMP R_sq, R_norm, G, dG, nxe, idx)
         do k = 1, total_pixels
            obs_x = xf(k)
            obs_y = yf(k)
            obs_z = 0.0_dp

            sEx = (0.0_dp, 0.0_dp); sEy = (0.0_dp, 0.0_dp); sEz = (0.0_dp, 0.0_dp)
            sBx = (0.0_dp, 0.0_dp); sBy = (0.0_dp, 0.0_dp); sBz = (0.0_dp, 0.0_dp)

            do s = 1, total_src
               key_Rx = obs_x - geom%x(s)
               key_Ry = obs_y - geom%y(s)
               key_Rz = obs_z - geom%z(s)
               R_sq = key_Rx * key_Rx + key_Ry * key_Ry + key_Rz * key_Rz
               R_norm = sqrt(R_sq)

               G = exp(I_COMP * kn * R_norm) / R_norm
               dG = (I_COMP * kn / R_norm - 1.0_dp / R_sq) * G
               nxe = vec_N_dot_E(s)

               sEx = sEx + (I_COMP * wn * vec_NxB_x(s) * G + nxe * dG * key_Rx)
               sEy = sEy + (I_COMP * wn * vec_NxB_y(s) * G + nxe * dG * key_Ry)
               sEz = sEz + (I_COMP * wn * vec_NxB_z(s) * G + nxe * dG * key_Rz)

               sBx = sBx + (vec_NxB_y(s) * (dG * key_Rz) - vec_NxB_z(s) * (dG * key_Ry))
               sBy = sBy + (vec_NxB_z(s) * (dG * key_Rx) - vec_NxB_x(s) * (dG * key_Rz))
               sBz = sBz + (vec_NxB_x(s) * (dG * key_Ry) - vec_NxB_y(s) * (dG * key_Rx))
            end do

            idx = (k - 1) * N_FREQ + n
            g_field%data(idx)%Ex = sEx
            g_field%data(idx)%Ey = sEy
            g_field%data(idx)%Ez = sEz
            g_field%data(idx)%Bx = sBx
            g_field%data(idx)%By = sBy
            g_field%data(idx)%Bz = sBz
         end do
         ! print progress
         if (mod(n, 10) == 0 .or. n == 1 .or. n == N_FREQ) then
            print *, "  Calculated freq ", n, "/", N_FREQ
         end if
      end do

      is_initialized = .true.
      call save_cache()
      print *, "[Gpulse] Initialization Complete."
   end subroutine gpulse_init

   subroutine gpulse_free()
      if (allocated(g_field%data)) deallocate(g_field%data)
      if (allocated(g_field%omegas)) deallocate(g_field%omegas)
      is_initialized = .false.
   end subroutine gpulse_free

   ! ==========================================
   ! 缓��� I/O (Cache I/O)
   ! ==========================================

   subroutine save_cache()
      USE shared_data, ONLY: rank
      integer :: iu, i, total_cells
      integer :: h_n_x, h_n_y, h_n_freq, h_order, h_radial
      real(dp) :: h_range, h_lambda0, h_w0, h_p_peak, h_width, h_f0, h_d, h_hole

      iu = 283
      if (.not. allocated(g_field%data)) return
      if (rank /= 0) return

      open(unit=283, file=CACHE_FILENAME, access='stream', &
         form='unformatted', status='replace', action='write')

      ! 写入 header 参数（用于后续验证）
      h_n_x = num_x;  h_n_y = num_y;  h_n_freq = N_FREQ
      h_range = RANGE; h_lambda0 = LAMBDA0; h_w0 = W0
      h_p_peak = P_PEAK; h_width = WIDTH; h_order = ORDER
      h_f0 = F0; h_d = D; h_hole = HOLE; h_radial = RADIAL

      write(iu) h_n_x, h_n_y, h_n_freq
      write(iu) h_range, h_lambda0, h_w0, h_p_peak, h_width
      write(iu) h_order
      write(iu) h_f0, h_d, h_hole
      write(iu) h_radial

      ! 写入 omegas 数组
      write(iu) g_field%omegas

      ! 逐个写入 Cell 数据（避免派生类型 padding 问题）
      total_cells = size(g_field%data)
      do i = 1, total_cells
         write(iu) g_field%data(i)%Ex, g_field%data(i)%Ey, g_field%data(i)%Ez, &
            g_field%data(i)%Bx, g_field%data(i)%By, g_field%data(i)%Bz
      end do

      close(iu)
      print *, '[Gpulse] Field data saved to cache: ', CACHE_FILENAME
   end subroutine save_cache

   function load_cache() result(success)
      logical :: success
      integer :: iu, i, total_cells, ios
      logical :: file_exists
      integer :: h_n_x, h_n_y, h_n_freq, h_order, h_radial
      real(dp) :: h_range, h_lambda0, h_w0, h_p_peak, h_width, h_f0, h_d, h_hole
      real(dp) :: range_val

      success = .false.
      iu = 283

      inquire(file=CACHE_FILENAME, exist=file_exists)
      if (.not. file_exists) return

      open(unit=283, file=CACHE_FILENAME, access='stream', &
         form='unformatted', status='old', action='read', iostat=ios)
      if (ios /= 0) return

      ! 读取 header
      read(iu, iostat=ios) h_n_x, h_n_y, h_n_freq
      if (ios /= 0) then; close(iu); return; end if
      read(iu, iostat=ios) h_range, h_lambda0, h_w0, h_p_peak, h_width
      if (ios /= 0) then; close(iu); return; end if
      read(iu, iostat=ios) h_order
      if (ios /= 0) then; close(iu); return; end if
      read(iu, iostat=ios) h_f0, h_d, h_hole
      if (ios /= 0) then; close(iu); return; end if
      read(iu, iostat=ios) h_radial
      if (ios /= 0) then; close(iu); return; end if

      ! 验证参数是否匹配
      if (h_n_x /= num_x .or. h_n_y /= num_y .or. h_n_freq /= N_FREQ .or. &
         abs(h_range - RANGE) > 1.0e-12_dp .or. &
         abs(h_lambda0 - LAMBDA0) > 1.0e-12_dp .or. &
         abs(h_w0 - W0) > 1.0e-12_dp .or. &
         abs(h_p_peak - P_PEAK) > 1.0e-12_dp .or. &
         abs(h_width - WIDTH) > 1.0e-12_dp .or. &
         h_order /= ORDER .or. &
         abs(h_f0 - F0) > 1.0e-12_dp .or. &
         abs(h_d - D) > 1.0e-12_dp .or. &
         abs(h_hole - HOLE) > 1.0e-12_dp .or. &
         h_radial /= RADIAL) then
         print *, '[Gpulse] Cache parameters mismatch. Recalculating...'
         close(iu)
         return
      end if

      ! 参数匹配，填充 g_field 元数据
      range_val = h_range
      g_field%N_pixels_x = h_n_x
      g_field%N_pixels_y = h_n_y
      g_field%N_freq = h_n_freq
      g_field%min_x = -range_val
      g_field%max_x = range_val
      g_field%min_y = -range_val
      g_field%max_y = range_val
      g_field%dx = (2.0_dp * range_val) / real(max(1, h_n_x - 1), dp)
      g_field%dy = (2.0_dp * range_val) / real(max(1, h_n_y - 1), dp)

      allocate(g_field%omegas(h_n_freq))
      allocate(g_field%data(h_n_x * h_n_y * h_n_freq))

      ! 读取 omegas
      read(iu, iostat=ios) g_field%omegas
      if (ios /= 0) then
         close(iu)
         call gpulse_free()
         return
      end if

      ! 逐个读取 Cell 数据
      total_cells = h_n_x * h_n_y * h_n_freq
      do i = 1, total_cells
         read(iu, iostat=ios) g_field%data(i)%Ex, g_field%data(i)%Ey, g_field%data(i)%Ez, &
            g_field%data(i)%Bx, g_field%data(i)%By, g_field%data(i)%Bz
         if (ios /= 0) then
            close(iu)
            call gpulse_free()
            return
         end if
      end do

      close(iu)
      success = .true.
      print *, '[Gpulse] Loaded field data from cache.'
   end function load_cache

   function bilinear_interp(c00, c10, c01, c11, tx, ty) result(res)
      complex(dpc), intent(in) :: c00, c10, c01, c11
      real(dp), intent(in) :: tx, ty
      complex(dpc) :: res, c0, c1

      c0 = c00 * (1.0_dp - tx) + c10 * tx
      c1 = c01 * (1.0_dp - tx) + c11 * tx
      res = c0 * (1.0_dp - ty) + c1 * ty
   end function bilinear_interp

   subroutine gpulse_complex(t, x, y, z, Ex, Ey, Ez, Bx, By, Bz)
      real(dp), intent(in) :: t, x, y, z
      complex(dpc), intent(out) :: Ex, Ey, Ez, Bx, By, Bz
      real(dp) :: u, v, tx, ty, theta_base
      integer :: ix, iy, ix0, ix1, iy0, iy1, k, n_y, n_f
      complex(dpc) :: phase, s(6, 0:1, 0:1)

      if (.not. is_initialized) call gpulse_init()

      if (x < g_field%min_x .or. x > g_field%max_x .or. &
         y < g_field%min_y .or. y > g_field%max_y) then
         Ex = (0.0_dp, 0.0_dp); Ey = (0.0_dp, 0.0_dp); Ez = (0.0_dp, 0.0_dp)
         Bx = (0.0_dp, 0.0_dp); By = (0.0_dp, 0.0_dp); Bz = (0.0_dp, 0.0_dp)
         return
      end if

      u = (x - g_field%min_x) / g_field%dx
      v = (y - g_field%min_y) / g_field%dy
      ix = floor(u); iy = floor(v)

      ix0 = max(0, min(ix, g_field%N_pixels_x - 1))
      ix1 = max(0, min(ix + 1, g_field%N_pixels_x - 1))
      iy0 = max(0, min(iy, g_field%N_pixels_y - 1))
      iy1 = max(0, min(iy + 1, g_field%N_pixels_y - 1))

      tx = max(0.0_dp, min(1.0_dp, u - real(ix, dp)))
      ty = max(0.0_dp, min(1.0_dp, v - real(iy, dp)))

      s = (0.0_dp, 0.0_dp)
      n_y = g_field%N_pixels_y; n_f = g_field%N_freq
      theta_base = t + z / C_CONST

      do k = 1, n_f
         phase = exp(-I_COMP * g_field%omegas(k) * theta_base)
         s(:,0,0) = s(:,0,0) + cell_to_vec(g_field%data((ix0 * n_y + iy0) * n_f + k)) * phase
         s(:,1,0) = s(:,1,0) + cell_to_vec(g_field%data((ix1 * n_y + iy0) * n_f + k)) * phase
         s(:,0,1) = s(:,0,1) + cell_to_vec(g_field%data((ix0 * n_y + iy1) * n_f + k)) * phase
         s(:,1,1) = s(:,1,1) + cell_to_vec(g_field%data((ix1 * n_y + iy1) * n_f + k)) * phase
      end do

      Ex = bilinear_interp(s(1,0,0), s(1,1,0), s(1,0,1), s(1,1,1), tx, ty)
      Ey = bilinear_interp(s(2,0,0), s(2,1,0), s(2,0,1), s(2,1,1), tx, ty)
      Ez = bilinear_interp(s(3,0,0), s(3,1,0), s(3,0,1), s(3,1,1), tx, ty)
      Bx = bilinear_interp(s(4,0,0), s(4,1,0), s(4,0,1), s(4,1,1), tx, ty)
      By = bilinear_interp(s(5,0,0), s(5,1,0), s(5,0,1), s(5,1,1), tx, ty)
      Bz = bilinear_interp(s(6,0,0), s(6,1,0), s(6,0,1), s(6,1,1), tx, ty)

   contains
      function cell_to_vec(c) result(vec_out)
         type(Cell), intent(in) :: c
         complex(dpc) :: vec_out(6)
         vec_out = [c%Ex, c%Ey, c%Ez, c%Bx, c%By, c%Bz]
      end function cell_to_vec
   end subroutine gpulse_complex

   subroutine gpulse(t, x, y, z, Ex, Ey, Ez, Bx, By, Bz)
      real(dp), intent(in) :: t, x, y, z
      real(dp), intent(out) :: Ex, Ey, Ez, Bx, By, Bz
      complex(dpc) :: cEx, cEy, cEz, cBx, cBy, cBz

      call gpulse_complex(t, x, y, z, cEx, cEy, cEz, cBx, cBy, cBz)

      if (USE_ABS_OUTPUT) then
         Ex = abs(cEx)
         Ey = abs(cEy)
         Ez = abs(cEz)
         Bx = abs(cBx)
         By = abs(cBy)
         Bz = abs(cBz)
      else
         Ex = real(cEx, dp)
         Ey = real(cEy, dp)
         Ez = real(cEz, dp)
         Bx = real(cBx, dp)
         By = real(cBy, dp)
         Bz = real(cBz, dp)
      end if
   end subroutine gpulse

end module gpulse_mod
