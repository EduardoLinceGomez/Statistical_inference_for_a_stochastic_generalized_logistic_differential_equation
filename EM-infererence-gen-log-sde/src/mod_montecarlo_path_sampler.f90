!> @ingroup modules
!> @author E. Lince-Gomez, F. Baltazar-Larios, S. Diaz-Infante
!> @brief This module enclose functions to handle sampling of realizations
!> of the process solution to the generalized logistic SDE

module mod_montecarlo_path_sampler
    use mod_brownian_motion
    use mod_path_sampler
    use ifport
    use iso_fortran_env, only: &
            stdin => input_unit, &
            stdout => output_unit, &
            stderr => error_unit
    implicit none
contains
!> @brief Returns a csv file with the sampling of
!>  n_paths trajectories. The format is optimized to read as a pandas data frame.
!> @param [in] time_horizon final time of interval $[0, T]$
!> @param [in] n_resolution number of points of resoultion
!> @param [in] r_operative_factor multiplier of resoultion
!> @param [out] n_montecarlo number of path to be sampled

!> @param [optional] user_seed if given, fix the seed for the number gerneartor
!> @param [optional] debug, logical, if it is fixed as .true., then save the
!> underlying path iwth resolution step
    
    subroutine montecarlo_path_sampler( &
        &time_horizon, &
        &n_resolution, &
        &r_operative_factor, &
        &n_montecarlo, &
        &n_operative, &
        &x_0, &
        &alpha, &
        &m, &
        &sigma, &
        &debug, &
        &file_name &
    )
    integer, intent(in) :: n_montecarlo, n_resolution, r_operative_factor
    character (len = 100), intent(in), optional :: file_name
    logical, intent(in), optional :: debug
    real(kind=8), intent(in) :: x_0, alpha, m, sigma, time_horizon
    
    real(kind=8), allocatable :: brownian_path(:), sample_path(:), mc_sample(:, :)
    real(kind=8) h_op, h_res
    integer ::  date_time(8), user_seed=394569125
    integer(kind=4) :: n_operative, i
    character*10 b(3)
    character(len=8) :: fmt
    character (len = 100) :: filename = "../data/mc_output.csv"
    
    if (present(file_name)) then
        filename = trim(file_name)
    end if
    
    call get_brownian_path(&
            &n_resolution, &
            &r_operative_factor, &
            &time_horizon, &
            &brownian_path, &
            &h_op, &
            &n_operative, &
            &user_seed, &
            &debug &
            &)
    allocate(mc_sample(n_montecarlo, n_operative + 1))
    !    mc_sample(:, :) = 0.0
    !    print *, "head(X_mc)", mc_sample(1:5, 1:10)
    call date_and_time(b(1), b(2), b(3), date_time)
    write (filename, &
            &"('../data/output_mc_', &
            &I4.4, '-', &
            &I2.2, '-', &
            &I2.2, '_', &
            &I2.2, ':', &
            &I2.2, ':'  &
            &I2.2, ':', &
            &I3.3, &
            &'.csv') &
            ") date_time(1), date_time(2), date_time(3),&
            &date_time(5), date_time(6), date_time(7), date_time(8)
    open (unit=6, carriagecontrol='fortran')
    do i=1, n_montecarlo
        user_seed = i
        call get_brownian_path(&
                &n_resolution, &
                &r_operative_factor, &
                &time_horizon, &
                &brownian_path, &
                &h_op, &
                &n_operative, &
                &user_seed, &
                &debug &
                &)
        ! To tag the current path
        !        print *, '(===)', i, ' ', trim(file_name)
        ! sample_path = [(0.0, i=0, n_operative)]
        call get_sample_path(&
                &brownian_path, &
                &n_operative, &
                &h_op, &
                &x_0, &
                &alpha, &
                &m, &
                &sigma, &
                &sample_path, &
                &debug&!, &
                ! &file_name&
                &)
        mc_sample(i, :) = sample_path
        !        print *, "(===):", i
        !        print *, "head(X_mc[i, :] )", mc_sample(i, 1:10)
        call progress(i)
    end  do
    call save_mc_csv(mc_sample, h_op, n_montecarlo, n_operative, file_name)
    ! TODO:write a function to save a csv
    ! write(1, fmt=csv_fmt) 0, sep, 0.0, sep, x_0
    end subroutine montecarlo_path_sampler
    
    subroutine progress(j)
        implicit none
        integer, intent(in) :: j
        integer::i, k
        character(len=60) :: bar
        bar="???%[....................................................]"
        i = 0
        do k=0, j
            if  (MOD(k, 20)  == 0) then
                i = i+1
                write(unit=bar(1:3),fmt="(i3)") (2 * i - 2)
                bar(5 + i: 5+i+1)= "#>"
            end if
        enddo
        ! print the progress bar.
        write(unit=6, fmt="(a1, a1, a60)") '+',char(60), bar
        return
    end subroutine progress
    
    subroutine save_mc_csv(mc_sample, h_op, n_montecarlo, n_operative, filename)
        integer(kind=4), intent(in) :: n_montecarlo, n_operative
        real(kind=8), intent(in) ::  mc_sample(n_montecarlo, n_operative), h_op
        character (len = 100), intent(in), optional :: filename
        integer(kind=4) :: i, j, fileunit
        character (len = 100) :: default_filename, header
        character(len=:), allocatable :: csv_fmt, header_fmt
        character(len=1)  :: sep = ','
        default_filename = '../data/mc_sample_.csv'
        header_fmt ='(a4, a, a4, a, a12, a, a12)'
        csv_fmt = '(i4, a, i4, a, f12.8, a, f12.8)'
        header = 'idx' // sep // 'j' // sep  // 't_j' // sep // 'X(t_{j})'
        if(present(filename)) then
            open(newunit=fileunit, file=trim(filename))
            print *, filename
        else
            open(newunit=fileunit, file=trim(default_filename))
        end if
        write(fileunit, header_fmt) header
        do i=1, n_montecarlo
            do j=1, n_operative + 1
                write(fileunit, fmt=csv_fmt) i, sep, (j - 1), sep, (j - 1) * h_op, sep, mc_sample(i, j)
            end do
        end do
    end subroutine  save_mc_csv
end module mod_montecarlo_path_sampler