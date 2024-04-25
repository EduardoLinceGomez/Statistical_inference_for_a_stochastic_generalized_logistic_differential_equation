program mod_montecarlo_path_sampler_test_dev
    use mod_brownian_motion
    use mod_path_sampler
    use mod_montecarlo_path_sampler
    use ifport
    use iso_fortran_env, only: &
            stdin => input_unit, &
            stdout => output_unit, &
            stderr => error_unit
    
    implicit none
    integer :: n_montecarlo, n_resolution, r_operative_factor, user_seed
    integer ::  date_time(8)
    character*10 b(3)
    character (len = 100) :: file_name
    character(len=8) :: fmt
    logical :: debug
    real(kind=8), allocatable :: brownian_path(:), sample_path(:), mc_sample(:, :)
    integer(kind=4) :: n_operative, i
    real(kind=8) h_op, h_res, time_horizon
    real(kind=8) :: x_0, alpha, m, sigma
    !
    time_horizon = 32.0
    n_resolution = 4096
    r_operative_factor = 64
    n_montecarlo = 1000
    n_operative = 0
   
    user_seed=123456
    debug = .false.
    x_0 = 0.1
    alpha = 0.8
    m = 1.5
    sigma = 0.1
    file_name = "../data/mod_mont_test.csv"
    
    call montecarlo_path_sampler( &
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
end program mod_montecarlo_path_sampler_test_dev
