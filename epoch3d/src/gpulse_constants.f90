
module gpulse_constants
   use shared_data, only: num
   implicit none
   public
   ! Kinds
   integer, parameter :: dp = num
   ! complex kind
   integer, parameter :: dpc = num

   ! --- Basic Physical Constants ---
   real(dp), parameter :: PI = 3.14159265358979323846_dp
   complex(dpc), parameter :: I_COMP = (0.0_dp, 1.0_dp)

   ! Unit System Toggle
   ! .true.: Normalized Units, .false.: SI Units
   logical, parameter :: USE_NORMALIZED_UNITS = .false.

   ! --- SI Reference Values ---
   real(dp), parameter :: C_SI = 299792458.0_dp          ! (m/s)
   real(dp), parameter :: EPS0_SI = 8.854187817e-12_dp   ! (F/m)
   real(dp), parameter :: Me_SI = 9.10938356e-31_dp      ! (kg)
   real(dp), parameter :: Qe_SI = 1.60217663e-19_dp      ! (C)

   ! --- Input Parameters in SI ---
   real(dp), parameter :: LAMBDA0_SI = 1.8e-6_dp   ! (1.8 um)
   real(dp), parameter :: W0_SI = 8.9e-3_dp        ! (8.9 mm)
   real(dp), parameter :: P_PEAK_SI = 0.3e12_dp    ! (0.3 TW)
   real(dp), parameter :: WIDTH_SI = 700e-9_dp     ! (700 nm)
   real(dp), parameter :: F0_SI = 0.00635_dp       ! (6.35 mm)
   real(dp), parameter :: D_SI = 0.0254_dp         ! (25.4 mm)
   real(dp), parameter :: HOLE_SI = 0.005_dp       ! (5 mm)
   real(dp), parameter :: RANGE_SI = 4e-6_dp       ! (4 um)

   ! --- Reference Scales Calculation ---
   real(dp), parameter :: omega_r = 2.0_dp * PI * C_SI / LAMBDA0_SI

   real(dp), parameter :: T_r = 1.0_dp / omega_r
   real(dp), parameter :: L_r = C_SI / omega_r
   real(dp), parameter :: V_r = C_SI
   real(dp), parameter :: E_r = Me_SI * C_SI * omega_r / Qe_SI
   real(dp), parameter :: B_r = Me_SI * omega_r / Qe_SI

   real(dp), parameter :: I_r = 0.5_dp * C_SI * EPS0_SI * E_r * E_r
   real(dp), parameter :: Power_r = I_r * L_r * L_r

   ! --- Active Constants ---
   real(dp), parameter :: C_CONST = merge(1.0_dp, C_SI, USE_NORMALIZED_UNITS)
   real(dp), parameter :: EPS0 = merge(2.0_dp, EPS0_SI, USE_NORMALIZED_UNITS)

   real(dp), parameter :: LAMBDA0 = merge(LAMBDA0_SI / L_r, LAMBDA0_SI, USE_NORMALIZED_UNITS)
   real(dp), parameter :: W0      = merge(W0_SI / L_r,      W0_SI,      USE_NORMALIZED_UNITS)
   real(dp), parameter :: P_PEAK  = merge(P_PEAK_SI / Power_r,P_PEAK_SI, USE_NORMALIZED_UNITS)
   real(dp), parameter :: WIDTH   = merge(WIDTH_SI / L_r,   WIDTH_SI,   USE_NORMALIZED_UNITS)
   real(dp), parameter :: F0      = merge(F0_SI / L_r,      F0_SI,      USE_NORMALIZED_UNITS)
   real(dp), parameter :: D       = merge(D_SI / L_r,       D_SI,       USE_NORMALIZED_UNITS)
   real(dp), parameter :: HOLE    = merge(HOLE_SI / L_r,    HOLE_SI,    USE_NORMALIZED_UNITS)
   real(dp), parameter :: RANGE   = merge(RANGE_SI / L_r,   RANGE_SI,   USE_NORMALIZED_UNITS)

   integer, parameter :: ORDER = 7
   integer, parameter :: RADIAL = 1

   ! 输出模式开关：.true. 输出模值(abs)，.false. 输出实部(real)
   logical, parameter :: USE_ABS_OUTPUT = .false.

   integer, parameter :: NI = 8
   integer, parameter :: N_FREQ = 101

   integer, parameter :: num_x = 100
   integer, parameter :: num_y = 100

   ! --- Data Structures ---
   type :: Geometry
      real(dp), allocatable :: x(:), y(:), z(:)
      real(dp), allocatable :: Nx(:), Ny(:), Nz(:)
      real(dp), allocatable :: dS(:)
      real(dp), allocatable :: r(:), theta(:)
      integer :: total_points
   end type Geometry

   type :: Spectrum
      real(dp), allocatable :: omega_list(:)
      real(dp), allocatable :: weights(:)
      real(dp), allocatable :: E_amp(:)
      real(dp) :: d_omega
   end type Spectrum

contains

   ! --- Helper Functions ---

   function linspace(start_val, end_val, num) result(res)
      real(dp), intent(in) :: start_val, end_val
      integer, intent(in) :: num
      real(dp), allocatable :: res(:)
      integer :: i

      allocate(res(num))
      if (num > 1) then
         do i = 1, num
            res(i) = start_val + (end_val - start_val) * real(i - 1, dp) / real(num - 1, dp)
         end do
      else if (num == 1) then
         res(1) = start_val
      end if
   end function linspace

   function build_geometry() result(geom)
      type(Geometry) :: geom
      integer :: Nr, Nth, i, j, idx
      real(dp), allocatable :: r_vec(:), th_vec(:)
      real(dp), allocatable :: wr(:), wth(:)
      real(dp) :: dr, dth, r, th, dzdx, dzdy, Norm

      Nr = int(2**NI) + 1
      Nth = int(2**(NI + 1)) + 1
      geom%total_points = Nr * Nth

      allocate(geom%x(geom%total_points), geom%y(geom%total_points), geom%z(geom%total_points))
      allocate(geom%Nx(geom%total_points), geom%Ny(geom%total_points), geom%Nz(geom%total_points))
      allocate(geom%dS(geom%total_points))
      allocate(geom%r(geom%total_points), geom%theta(geom%total_points))

      r_vec = linspace(HOLE / 2.0_dp, D / 2.0_dp, Nr)
      th_vec = linspace(-PI, PI, Nth)
      dr = r_vec(2) - r_vec(1)
      dth = th_vec(2) - th_vec(1)

      allocate(wr(Nr), wth(Nth))
      wr = 1.0_dp
      do i = 2, Nr - 1, 2
         wr(i) = 4.0_dp
      end do
      do i = 3, Nr - 2, 2
         wr(i) = 2.0_dp
      end do

      wth = 1.0_dp
      do j = 2, Nth - 1, 2
         wth(j) = 4.0_dp
      end do
      do j = 3, Nth - 2, 2
         wth(j) = 2.0_dp
      end do

      idx = 1
      do j = 1, Nth
         do i = 1, Nr
            r = r_vec(i)
            th = th_vec(j)
            geom%r(idx) = r
            geom%theta(idx) = th

            geom%x(idx) = r * cos(th)
            geom%y(idx) = r * sin(th)
            geom%z(idx) = r * r / (4.0_dp * F0) - F0

            dzdx = geom%x(idx) / (2.0_dp * F0)
            dzdy = geom%y(idx) / (2.0_dp * F0)
            Norm = sqrt(dzdx * dzdx + dzdy * dzdy + 1.0_dp)

            geom%Nx(idx) = -dzdx / Norm
            geom%Ny(idx) = -dzdy / Norm
            geom%Nz(idx) = 1.0_dp / Norm

            geom%dS(idx) = (r * dr * dth) * (wth(j) * wr(i) / 9.0_dp) * Norm
            idx = idx + 1
         end do
      end do
   end function build_geometry

   function calculate_spectrum() result(spec)
      type(Spectrum) :: spec
      real(dp) :: omega0, dw_Range, m_val, width_p, I_peak, E_tgt, sum_num, val, ratio, w_simp
      integer :: i

      omega0 = 2.0_dp * PI * C_CONST / LAMBDA0
      dw_Range = WIDTH * omega0 * omega0 / (2.0_dp * PI * C_CONST)
      m_val = real(ORDER, dp)

      width_p = (dw_Range / 2.0_dp) / (log(2.0_dp)**(1.0_dp / m_val))

      allocate(spec%omega_list(N_FREQ))
      spec%omega_list = linspace(omega0 - dw_Range, omega0 + dw_Range, N_FREQ)
      spec%d_omega = spec%omega_list(2) - spec%omega_list(1)

      allocate(spec%weights(N_FREQ))
      do i = 1, N_FREQ
         if (i == 1 .or. i == N_FREQ) then
            w_simp = 1.0_dp
         else if (mod(i-1, 2) == 1) then
            w_simp = 4.0_dp
         else
            w_simp = 2.0_dp
         end if
         spec%weights(i) = (w_simp / 3.0_dp) * spec%d_omega
      end do

      I_peak = 2.0_dp * P_PEAK / (PI * W0**2)
      E_tgt = sqrt(2.0_dp * I_peak / (C_CONST * EPS0))

      allocate(spec%E_amp(N_FREQ))
      sum_num = 0.0_dp

      do i = 1, N_FREQ
         val = exp(-0.5_dp * (abs((spec%omega_list(i) - omega0) / width_p))**m_val)
         spec%E_amp(i) = val
         sum_num = sum_num + val * spec%weights(i)
      end do

      ratio = E_tgt / sum_num
      spec%E_amp = spec%E_amp * ratio

   end function calculate_spectrum

end module gpulse_constants
