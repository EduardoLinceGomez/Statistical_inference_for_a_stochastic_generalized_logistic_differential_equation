program mod_path_sampler_test
    use mod_brownian_motion, only: get_brownian_path
    use mod_path_sampler
    implicit none
    integer :: n_resolution, r_operative_factor, user_seed
    logical :: debug
    real(kind=8), allocatable :: brownian_path(:), sample_path(:)
    integer(kind=4) :: n_operative, i
    real(kind=8) h_op, h_res, time_horizon
    real(kind=8) :: x_0, alpha, m, sigma

    time_horizon = 32.0
    n_resolution = 4096
    r_operative_factor = 64
    n_operative = 0
    h_op = 0.0
    user_seed=123456
    debug = .false.
    x_0 = 0.1
    alpha = 0.8
    m = 1.5
    sigma = 0.1
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
    debug = .false.
    sample_path = [(0.0, i=0, n_operative)]
    call get_sample_path(&
            &brownian_path, &
            &n_operative, &
            &h_op, &
            &x_0, &
            &alpha, &
            &m, &
            &sigma, &
            &sample_path, &
            &debug&
        &)
    print*,"Hello Brownian path (n_op, h_op)= ", n_operative, h_op
    print *, "head(B_t):=", brownian_path(0:5)
    print *, "head(X_t):=", sample_path(0:5)
end program mod_path_sampler_test
