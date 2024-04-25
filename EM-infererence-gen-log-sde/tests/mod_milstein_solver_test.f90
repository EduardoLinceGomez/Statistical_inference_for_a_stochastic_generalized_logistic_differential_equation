program mod_milstein_solver_test
    use mod_random_number_generator
    use mod_stochastic_logistic_model
    use mod_milstein_solver
    implicit none
    real(kind=8) :: alpha_0, m_0, sigma_0, delta, brownian_start, brownian_end, &
        current_state, current_drift, current_diffusion, &
        diffusion_derivative, sigma, x_milstein_t
    integer seed
    parameter(seed=123)
    real(kind=8) rn_0, rn_1, rn_2
    delta = 0.001
    current_state = 0.1
    alpha_0 = 0.8
    m_0 = 1.5
    sigma_0 = 0.1
    
    brownian_start = 0.0
    !
    ! call system('clear')
    call get_brownian_increment(delta, brownian_start,  brownian_end)
    call compute_drift(alpha_0, m_0, current_state, current_drift)
    call compute_diffusion(sigma_0, current_state, current_diffusion)
    call get_milstein_iteration(delta,  &
        current_state, current_drift, current_diffusion, &
        diffusion_derivative, sigma_0, x_milstein_t)
    print*, " "
    print *,"(====) testing mod_milstein_solver"
    print *, '(====)(----) Brownian_increment: ', brownian_end
    print *, '(====)(----) : Milstein_iteration: =', x_milstein_t
    ! print *, '(====)(----) normalval: rn =', rn_2
end program

! gfortran mod_random_number_generator.f95 main_.f95 -o b.out
! gfortran mod_random_number_generator.f95 mod_stochastic_logistic_model.f95 main.f95 -o b.out
! gfortran mod_stochastic_logistic_model.f95 main.f95 -o b.out
