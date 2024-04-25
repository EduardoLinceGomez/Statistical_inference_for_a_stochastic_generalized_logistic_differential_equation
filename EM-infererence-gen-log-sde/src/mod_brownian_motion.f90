!> @ingroup modules
!> @author E. Lince-Gomez, F. Baltazar-Larios, S. Diaz-Infante
!> @brief This module contains functions to handle various aspects of
!> one-dimensional Brownian motion.
!> @see Higham & Kloeden 2021
module mod_brownian_motion
    use mod_random_number_generator
    use iso_fortran_env, only: &
    stdin => input_unit, &
    stdout => output_unit, &
    stderr => error_unit

    implicit none
contains
!> @brief Returns a realization of a standard Brownian motion given a
!> homogenoun stencil with L operational points and a resolution of N points.
!> Such that $N = r * L$ and $h_op = r * h_res$.
!>
!> @param [in] n_resolution number of points of resoultion
!> @param [in] r_operative_factor multiplier of resoultion
!> @param [in] time_horizon final time of interval $[0, T]$
!> @param [out] brownian_path array whith the realization path using
!> the operational step-size
!> @param [optional] user_seed if given, fix the seed for the number gerneartor
!> @param [optional] debug, logical, if it is fixed as .true., then save the
!> underlying path iwth resolution step


    subroutine get_brownian_path(&
            &n_resolution, &
            &r_operative_factor, &
            &time_horizon, &
            &brownian_path, &
            &h_op, &
            &n_operative, &
            &user_seed, &
            &debug&
    &)
        implicit none
        integer, intent(in) :: n_resolution, r_operative_factor
        real(kind=8), intent(in) :: time_horizon
        real(kind=8), intent(out), allocatable :: brownian_path(:)
        real(kind=8), intent(out) :: h_op
        integer(kind=4), intent(out) :: n_operative
        integer, intent(in), optional :: user_seed
        logical, intent(in), optional :: debug
        TYPE (VSL_STREAM_STATE) :: stream
        integer(kind=4) :: i, j, idx_start, idx_end,&
                &seed_get, errcode
        integer :: brng, seed_put, method, n
        real(kind=8) h_res, xi, browinian_operative_inc, &
                browinian_resolution_inc, mean_a, std_a, w_inc, w_t_i
        real(kind=8), allocatable :: gaussian_sample(:), &
                &resolution_path(:)
        character(len=100) :: header
        character(len=100) :: filename = '../data/brownian_sample_path.csv'
        character(len=1)   :: sep = ','
        character(len=:), allocatable :: csv_fmt
        ! random number initialization
        if (present(user_seed)) then
                seed_put = user_seed
        else
                seed_put = 394569125
        end if
        brng = VSL_BRNG_MT19937
        errcode = vslNewStream(stream, brng, seed_put)
        mean_a = 0.0
        std_a = 1.0
        w_inc = 0.0
        w_t_i = 0.0
        h_res = time_horizon / n_resolution
        if (r_operative_factor > 1 ) then
            h_op = r_operative_factor * h_res
            n_operative = &
                &NINT(real(n_resolution, 4) / &
                    &real(r_operative_factor, 4))
        else if (r_operative_factor == 1) then
            h_op = h_res
            n_operative = n_resolution
        end if
        brownian_path = [(0.0, i=0, n_operative)]
        gaussian_sample = [(0.0, i=1, n_resolution)]
        resolution_path = [(0.0, i=0, n_resolution)]
!       Generating a bounch of Gaussian r.v.
        call mkl_gaussian_sampler(n_resolution, mean_a, dsqrt(h_res), &
                &seed_put, gaussian_sample)
        resolution_path(2:n_resolution + 1) = &
            &[(sum(gaussian_sample(1:i)), &
                &i = 1, &
                &size(gaussian_sample) &
            &)]
        j = 1
        header = 'i' // sep // 't_i' // sep // 'W(t_i)'

        open(unit=1, file=filename, status='replace', &
                action='write', form='formatted')
        csv_fmt = '(i4, a, f9.5, a, f9.5)'
        200 format(i4, a, f9.5, a, f9.5)
        write(1, *) header
        write(1, fmt=csv_fmt) 0, sep, 0.0, sep, 0.0
        do j=1, n_operative 
            idx_start = r_operative_factor * (j - 1) + 1
            idx_end = r_operative_factor * j
            w_inc = sum(gaussian_sample(idx_start:idx_end))
            !print *, "dWk_i: ", w_inc
            w_t_i = w_t_i + w_inc
            brownian_path(j + 1) = w_t_i
            write(1, fmt=csv_fmt) j, sep, j * h_op, sep, w_t_i
        end do
        close(unit=1)
        if (present(debug)) then
            if(debug) then
                filename = '../data/res_brownian_sample_path.csv'
                open(unit=10, file=filename, status='replace', &
                action='write', form='formatted')
                write(10, *) header
                write(10, 200) 0, sep, 0.0, sep, 0.0
                print *, "operative factor: ", r_operative_factor
                print *, "n_operative: ", n_operative
                print*, "head(B_tau_n):=", brownian_path(1:5)
                print*, '(===)random Stream state', errcode
                do i=1, n_resolution
                    write(10, 200) i, sep, &
                        &i * h_res, sep, &
                        &resolution_path(i)
                end do
                print*, "head(B_t_n):=", resolution_path(1:5)
            end if
            close(unit=10)
        end if
    end subroutine get_brownian_path
end module mod_brownian_motion
!ifort -c -i8 -I"${MKLROOT}/include" mod_brownian_motion.f90
