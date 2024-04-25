module mod_path_sampler
    use iso_fortran_env, only: &
        stdin => input_unit, &
        stdout => output_unit, &
        stderr => error_unit
!
    use mod_brownian_motion
    use mod_stochastic_logistic_model
    implicit none
contains
    subroutine get_sample_path(&
            &brownian_path, &
            &n_operative, &
            &h_op, &
            &x_0, &
            &alpha, &
            &m, &
            &sigma, &
            &sample_path, &
            &debug, &
            &filename&
    &)
!setup
        integer(kind=4), intent(in) :: n_operative
        real(kind=8), intent(in) :: brownian_path(n_operative + 1), h_op
        real(kind=8), intent(in) :: x_0, alpha, m, sigma
        real(kind=8) :: brownian_increment, current_state, &
            &current_drift, current_diffusion, &
            &diffusion_derivative, x_next, x_euler
        real(kind=8), allocatable, intent(out) :: sample_path(:)
        real(kind=8) :: euler_sample_path(n_operative + 1)
        logical, intent(in), optional :: debug
        character(len=100), intent(in), optional :: filename
        integer(kind=4) :: i, j, fileunit
        character(len=100) :: header
        character(len=100) :: default_filename
        character(len=1)   :: sep = ','
        character(len=:), allocatable :: csv_fmt, header_fmt
        sample_path = [(0.0, i=0, n_operative)]
        euler_sample_path = [(0.0, i=0, n_operative)]
        sample_path(1) = x_0
        euler_sample_path(1) = x_0
        x_next = x_0
        x_euler = x_0
        i = 0
        current_state = x_0
        header = 'i' // sep // 't_i' // sep // &
            &'|Milstein(X(t_{i}) - Euler(X(t_{i}))|'
        header_fmt ='(a4, a, a12, a, a12)'
        csv_fmt = '(i4, a, f12.8, a, f12.8)'
        default_filename = '../data/gen_log_sde_sample_path.csv'
!>
        if(present(debug) .AND. debug) then
            if(present(filename)) then
                open(newunit=fileunit, file=trim(filename))
                print *, filename
            else
                open(newunit=fileunit, file=trim(default_filename))
            end if
            if(debug) then
                write(stdout, fmt=header_fmt) 'i', sep,'t_i', sep, 'X(t_i)'
                write(stdout, fmt=csv_fmt) i, sep, i * h_op, sep, sample_path(0)
            end if
!        else
!            write(fileunit, fmt=header_fmt) 'i', sep,'t_i', sep, 'X(t_i)'
!            write(fileunit, fmt=csv_fmt) i, sep, i * h_op, sep, sample_path(0)
        end if
        do i=1, n_operative
            brownian_increment = brownian_path(i + 1) - brownian_path(i)
            call compute_drift(alpha, m, current_state, current_drift)
            call compute_diffusion(sigma, current_state, current_diffusion)
            call get_milstein_iteration(&
                    &h_op, &
                    &brownian_increment, &
                    &current_state ,&
                    &current_drift ,&
                    &current_diffusion, &
                    &sigma, & ! For this case diffusion is constant
                    &sigma, &
                    &x_euler, &
                    &x_next, &
                    &debug&
            &)
            euler_sample_path(i + 1) = x_euler
            sample_path(i + 1) = x_next
            if(present(debug) .AND. debug) then
                write(stdout, fmt=csv_fmt) i, sep, i * h_op, sep, x_next
            else
                if(debug) then
                    write(fileunit, fmt=csv_fmt) i, sep, i * h_op, sep, x_next
                end if
            end if
            current_state = x_next
        end do
        flush(fileunit)
    end subroutine get_sample_path
!> @brief Computes a iteration of the Milstein method
!> @details Defintion taken from the Kloedens book
!> @param [in ] delta determininstic time step-size
!> @param [in] brownian_start current observation of the Brownian process
!> the  Brownian increment.
    subroutine get_milstein_iteration(&
        &delta, &
        &brownian_increment, &
        &current_state, &
        &current_drift, &
        &current_diffusion, &
        &diffusion_derivative, &
        &sigma, &
        &x_euler, &
        &x_next, &
        &debug&
    &)
        implicit none
        real(kind=8), intent(in) :: delta, brownian_increment, &
            &current_state, current_drift, current_diffusion, &
            &diffusion_derivative, sigma
        real(kind=8), intent(out) :: x_next
        real(kind=8), optional, intent(out):: x_euler
        logical, intent(in), optional :: debug
        real(kind=8) :: euler_iteration
        euler_iteration = current_state + delta * current_drift &
            &+ sigma * current_state * brownian_increment
        if (present(x_euler)) then
            x_euler = 0.0
            x_euler = euler_iteration
        end if
        if (present(debug)) then
            if (debug) then
                print*, "X_{euler}(t_i): ", x_euler
            end if
        end if
        x_next = euler_iteration  &
            & + 0.5 * current_diffusion * diffusion_derivative * &
            & (brownian_increment * brownian_increment - delta)
    end subroutine get_milstein_iteration
end module mod_path_sampler
!ifort -c -i8 -I"${MKLROOT}/include" mod_path_sampler.f90
