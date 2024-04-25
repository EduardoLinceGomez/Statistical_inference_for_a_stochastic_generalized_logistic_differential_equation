program mod_milstein_sampler_test
    use iso_fortran_env
    use mod_random_number_generator
    use mod_stochastic_logistic_model
    use mod_milstein_solver
    implicit none!
    integer(int32) grid_size, i
    real(kind=8) x_0, alpha_0, m_0, sigma_0, delta, brownian_start, t_end, &
        brownian_end, current_state, current_drift, current_diffusion, &
        diffusion_derivative, sigma, x_milstein_t
    character(len=100) :: header
    character(len=100) :: filename = './data/realization_path_gen_logistic_SDE.csv'
    character(len=1)   :: sep = ','
    parameter(grid_size=4096, x_0=0.1, t_end=16.0)
    real(kind=8) x_milstein_path(grid_size)
    !
    call random_init(repeatable=.true., image_distinct=.true.)
    delta = t_end / grid_size
    current_state = x_0
    alpha_0 = 0.8
    m_0 = 1.5
    sigma_0 = 0.1
    i = 1
    x_milstein_t = x_0
    x_milstein_path(1) = current_state
    brownian_start = 0.0
    header = 'i' // sep //'t_i' // sep //'x_milstein(t_i)'
    print*, " "
    print *,"(====) testing mod_milstein_sampler,",&
        " grid-setup and solver parameters:"
    print *, '(====)(----) [x_0, t_end, delta]:=',&
        x_milstein_path(1), t_end, delta
    ! main loop
    print 100
    100 format(10x, "i", 8x, "t_i", 7x "x_milstein(t_i)")
    200 format(8x, "---------------------------------------")
    print 200 
    print*, i, (i - 1) * delta, x_milstein_path(i)
    open(unit=1, file=filename, status='replace', &
        action='write', form='formatted')
    write(1, *) header
    300 format(i4, a, f12.4, a, f12.4)
    write(1, 300) i, sep, (i - 1) * delta, sep, x_milstein_path(i)
    do i=2, grid_size
        call get_brownian_increment(delta, brownian_start,  brownian_end)
        call compute_drift(alpha_0, m_0, current_state, current_drift)
        call compute_diffusion(sigma_0, current_state, current_diffusion)
        call get_milstein_iteration(delta,  &
            current_state, current_drift, current_diffusion, &
            diffusion_derivative, sigma_0, x_milstein_t)
        x_milstein_path(i) = x_milstein_t
        current_state = x_milstein_t
        brownian_start = brownian_end
        write(1,300) i, sep, (i - 1) * delta, sep, x_milstein_path(i)
        if(mod(i, 100) == 0) then
            print*,  i, (i-1) * delta, x_milstein_path(i)
        end if
    end do
 end program mod_milstein_sampler_test
! 
! gfortran mod_random_number_generator.f95 main_.f95 -o b.out
! gfortran mod_random_number_generator.f95 mod_stochastic_logistic_model.f95 main.f95 -o b.out
! gfortran mod_stochastic_logistic_model.f95 main.f95 -o b.out
