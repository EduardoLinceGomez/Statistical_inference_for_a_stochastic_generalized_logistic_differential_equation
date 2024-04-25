program mod_brownian_motion_test
    use mod_brownian_motion, only: get_brownian_path
    implicit none
    integer :: n_resolution, r_operative_factor
    real(kind=8), allocatable :: brownian_path(:)
    integer(kind=4) :: n_operative, i, seed_get, seed_put
    real(kind=8) h_op, h_res, xi, browinian_operative_inc, &
        browinian_resolution_inc, time_horizon
    time_horizon = 16.0
    n_resolution = 2048
    r_operative_factor = 64
    n_operative = 0
    h_op = 0.0
    seed_put = 123456
    call get_brownian_path(&
        &n_resolution, &
        &r_operative_factor, &
        &time_horizon, &
        &brownian_path, &
        &h_op, &
        &n_operative, &
        &user_seed=123456, &
        &debug = .true. &
    &)
    print*,"Hello Brownian path (n_op, h_op)= ", n_operative, h_op
    print *, "head(B_t):=", brownian_path(1:5)
end program mod_brownian_motion_test
