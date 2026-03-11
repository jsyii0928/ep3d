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

MODULE deck_laser_block

   USE strings_advanced
   USE laser
   USE utilities

   IMPLICIT NONE
   SAVE

   PRIVATE
   PUBLIC :: laser_deck_initialise, laser_deck_finalise
   PUBLIC :: laser_block_start, laser_block_end
   PUBLIC :: laser_block_handle_element, laser_block_check

   TYPE(laser_block), POINTER :: working_laser
   LOGICAL :: boundary_set = .FALSE.
   INTEGER :: boundary

CONTAINS

   SUBROUTINE laser_deck_initialise

      n_lasers(:) = 0

   END SUBROUTINE laser_deck_initialise



   SUBROUTINE laser_deck_finalise

   END SUBROUTINE laser_deck_finalise



   SUBROUTINE laser_block_start

      IF (deck_state == c_ds_first) RETURN

      ! Every new laser uses the internal time function
      ! 当遇到一个新的 laser block 时，分配一个新的激光对象
      ! 并初始化默认行为：
      ! use_time_function = .FALSE. (默认不使用时间函数，除非指定 t_profile)
      ! use_phase_function = .TRUE. (默认使用相位函数)
      ! use_profile_function = .TRUE. (默认使用空间分布函数)
      ! use_omega_function = .FALSE. (默认不使用变频)
      ALLOCATE(working_laser)
      working_laser%use_time_function = .FALSE.
      working_laser%use_phase_function = .TRUE.
      working_laser%use_profile_function = .TRUE.
      working_laser%use_omega_function = .FALSE.

   END SUBROUTINE laser_block_start



   SUBROUTINE laser_block_end

      IF (deck_state == c_ds_first) RETURN

      ! laser block 结束时，将配置好的激光对象挂载到全局链表中
      CALL attach_laser(working_laser)
      boundary_set = .FALSE.

   END SUBROUTINE laser_block_end



   FUNCTION laser_block_handle_element(element, value) RESULT(errcode)
      ! 处理 laser block 中的配置元素。
      ! 这个函数会被多次调用，每次处理 input deck 中 laser block 的一行参数设置。
      ! element: 参数名 (key)，例如 "boundary", "amp", "lambda", "t_profile"
      ! value: 参数值 (value)，例如 "x", "1.0", "1e-6", "gauss(time, 10e-15, 20e-15)"
      ! errcode: 返回错误码

      CHARACTER(*), INTENT(IN) :: element, value
      INTEGER :: errcode
      REAL(num) :: dummy
      INTEGER :: io, iu

      errcode = c_err_none
      IF (deck_state == c_ds_first) RETURN
      IF (element == blank .OR. value == blank) RETURN

      ! 设置激光入射的边界 (boundary)
      ! 例如: boundary = x_min
      IF (str_cmp(element, 'boundary') .OR. str_cmp(element, 'direction')) THEN
         IF (rank == 0 .AND. str_cmp(element, 'direction')) THEN
            DO iu = 1, nio_units ! Print to stdout and to file (输出错误信息到标准输出和文件)
               io = io_units(iu)
               WRITE(io,*) '*** WARNING ***'
               WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
               WRITE(io,*) 'Element "direction" in the block "laser" is deprecated.'
               WRITE(io,*) 'Please use the element name "boundary" instead.'
            END DO
         END IF
         ! If the boundary has already been set, simply ignore further calls to it (防止重复设置边界)
         IF (boundary_set) RETURN
         ! 将字符串值转换为边界 ID (如 x_min, x_max 等)
         boundary = as_boundary_print(value, element, errcode)
         boundary_set = .TRUE.
         ! 初始化激光对象，将其绑定到指定边界
         CALL init_laser(boundary, working_laser)
         RETURN
      END IF

      ! 在设置其他属性之前必须先设置 boundary，否则报错
      IF (.NOT. boundary_set) THEN
         IF (rank == 0) THEN
            DO iu = 1, nio_units ! Print to stdout and to file
               io = io_units(iu)
               WRITE(io,*) '*** ERROR ***'
               WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
               WRITE(io,*) 'Cannot set laser properties before boundary is set'
            END DO
            CALL abort_code(c_err_required_element_not_set)
         END IF
         extended_error_string = 'boundary'
         errcode = c_err_required_element_not_set
         RETURN
      END IF

      ! 设置激光的归一化振幅 (amp)
      ! 例如: amp = 10.0
      IF (str_cmp(element, 'amp')) THEN
         working_laser%amp = as_real_print(value, element, errcode)
         RETURN
      END IF

      ! 设置激光强度/辐照度 (irradiance 或 intensity), 单位 W/m^2
      ! 自动转换为振幅存储
      ! SI (W/m^2)
      IF (str_cmp(element, 'irradiance') .OR. str_cmp(element, 'intensity')) THEN
         working_laser%amp = SQRT(as_real_print(value, element, errcode) &
            / (c*epsilon0/2.0_num))
         RETURN
      END IF

      ! 设置激光强度/辐照度 (irradiance_w_cm2 或 intensity_w_cm2), 单位 W/cm^2
      ! 自动转换为振幅存储
      IF (str_cmp(element, 'irradiance_w_cm2') &
         .OR. str_cmp(element, 'intensity_w_cm2')) THEN
         working_laser%amp = SQRT(as_real_print(value, element, errcode) &
            / (c*epsilon0/2.0_num)) * 100_num
         RETURN
      END IF

      ! 设置激光角频率 (omega)
      ! 也可以是一个随时间变化的函数
      IF (str_cmp(element, 'omega') .OR. str_cmp(element, 'freq')) THEN
         IF (rank == 0 .AND. str_cmp(element, 'freq')) THEN
            DO iu = 1, nio_units ! Print to stdout and to file
               io = io_units(iu)
               WRITE(io,*) '*** WARNING ***'
               WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
               WRITE(io,*) 'Element "freq" in the block "laser" is deprecated.'
               WRITE(io,*) 'Please use the element name "omega" instead.'
            END DO
         END IF
         ! 解析数学表达式
         CALL initialise_stack(working_laser%omega_function)
         CALL tokenize(value, working_laser%omega_function, errcode)
         working_laser%omega = 0.0_num
         working_laser%omega_func_type = c_of_omega
         ! 更新激光频率相关参数
         CALL laser_update_omega(working_laser)
         ! 标记频率是否随时变
         IF (working_laser%omega_function%is_time_varying) THEN
            working_laser%use_omega_function = .TRUE.
         ELSE
            CALL deallocate_stack(working_laser%omega_function)
         END IF
         RETURN
      END IF

      ! 设置频率 (frequency)，以 Hz 为单位
      IF (str_cmp(element, 'frequency')) THEN
         CALL initialise_stack(working_laser%omega_function)
         CALL tokenize(value, working_laser%omega_function, errcode)
         working_laser%omega = 0.0_num
         working_laser%omega_func_type = c_of_freq
         CALL laser_update_omega(working_laser)
         IF (working_laser%omega_function%is_time_varying) THEN
            working_laser%use_omega_function = .TRUE.
         ELSE
            CALL deallocate_stack(working_laser%omega_function)
         END IF
         RETURN
      END IF

      ! 设置波长 (lambda), 单位米
      IF (str_cmp(element, 'lambda')) THEN
         CALL initialise_stack(working_laser%omega_function)
         CALL tokenize(value, working_laser%omega_function, errcode)
         working_laser%omega = 0.0_num
         working_laser%omega_func_type = c_of_lambda
         CALL laser_update_omega(working_laser)
         IF (working_laser%omega_function%is_time_varying) THEN
            working_laser%use_omega_function = .TRUE.
         ELSE
            CALL deallocate_stack(working_laser%omega_function)
         END IF
         RETURN
      END IF

      ! 设置激光的空间分布 (profile)
      ! 这是一个空间坐标的函数，例如 y 或 z
      IF (str_cmp(element, 'profile')) THEN
         CALL initialise_stack(working_laser%profile_function)
         CALL tokenize(value, working_laser%profile_function, errcode)
         working_laser%profile = 0.0_num
         CALL laser_update_profile(working_laser)
         IF (working_laser%profile_function%is_time_varying) THEN
            working_laser%use_profile_function = .TRUE.
         ELSE
            CALL deallocate_stack(working_laser%profile_function)
         END IF
         RETURN
      END IF

      ! 设置激光的相位 (phase)
      IF (str_cmp(element, 'phase')) THEN
         CALL initialise_stack(working_laser%phase_function)
         CALL tokenize(value, working_laser%phase_function, errcode)
         working_laser%phase = 0.0_num
         CALL laser_update_phase(working_laser)
         IF (working_laser%phase_function%is_time_varying) THEN
            working_laser%use_phase_function = .TRUE.
         ELSE
            CALL deallocate_stack(working_laser%phase_function)
         END IF
         RETURN
      END IF

      ! 设置激光开启时间 (t_start)
      IF (str_cmp(element, 't_start')) THEN
         working_laser%t_start = as_time_print(value, element, errcode)
         RETURN
      END IF

      ! 设置相位偏移 (phase_offset)
      IF (str_cmp(element, 'phase_offset')) THEN
         working_laser%phase_offset = as_real_print(value, element, errcode)
         RETURN
      END IF

      ! 设置激光结束时间 (t_end)
      IF (str_cmp(element, 't_end')) THEN
         working_laser%t_end = as_time_print(value, element, errcode)
         RETURN
      END IF

      ! 设置激光的时间分布函数 (t_profile)
      ! 例如 t_profile = gauss(time, 0, 10e-15)
      IF (str_cmp(element, 't_profile')) THEN
         working_laser%use_time_function = .TRUE.
         CALL initialise_stack(working_laser%time_function)
         CALL tokenize(value, working_laser%time_function, errcode)
         ! evaluate it once to check that it's a valid block (尝试进行一次评估以检查语法)
         dummy = evaluate(working_laser%time_function, errcode)
         RETURN
      END IF

      ! 设置偏振角 (pol_angle 或 polarisation_angle)，单位是度
      ! 相对于模拟平面的偏振方向
      IF (str_cmp(element, 'pol_angle') &
         .OR. str_cmp(element, 'polarisation_angle')) THEN
         working_laser%pol_angle = as_real_print(value, element, errcode)
         RETURN
      END IF

      ! 设置偏振角 (pol 或 polarisation)，单位是度
      ! 内部转换为弧度
      IF (str_cmp(element, 'pol') &
         .OR. str_cmp(element, 'polarisation')) THEN
         ! Convert from degrees to radians
         working_laser%pol_angle = &
            pi * as_real_print(value, element, errcode) / 180.0_num
         RETURN
      END IF

      ! 设置激光的唯一标识符 ID (id)
      IF (str_cmp(element, 'id')) THEN
         working_laser%id = as_integer_print(value, element, errcode)
         RETURN
      END IF

      errcode = c_err_unknown_element

   END FUNCTION laser_block_handle_element



   FUNCTION laser_block_check() RESULT(errcode)

      INTEGER :: errcode
      TYPE(laser_block), POINTER :: current
      INTEGER :: error, io, iu

      errcode = c_err_none

      error = 0
      current => lasers
      DO WHILE(ASSOCIATED(current))
         ! 检查频率参数是否合法
         ! 如果频率、波长或角频率全为非法值，则报错
         IF (current%omega < 0.0_num) error = IOR(error, 1)
         ! 检查振幅参数是否合法
         ! 如果强度或振幅全为非法值，则报错
         IF (current%amp < 0.0_num) error = IOR(error, 2)
         current => current%next
      END DO

      IF (IAND(error, 1) /= 0) THEN
         IF (rank == 0) THEN
            DO iu = 1, nio_units ! Print to stdout and to file
               io = io_units(iu)
               WRITE(io,*) '*** ERROR ***'
               WRITE(io,*) 'Must define a "lambda" or "omega" for every laser.'
            END DO
         END IF
         errcode = c_err_missing_elements
      END IF

      IF (IAND(error, 2) /= 0) THEN
         IF (rank == 0) THEN
            DO iu = 1, nio_units ! Print to stdout and to file
               io = io_units(iu)
               WRITE(io,*) '*** ERROR ***'
               WRITE(io,*) 'Must define an "amp" or "irradiance" for every laser.'
            END DO
         END IF
         errcode = c_err_missing_elements
      END IF

   END FUNCTION laser_block_check

END MODULE deck_laser_block
