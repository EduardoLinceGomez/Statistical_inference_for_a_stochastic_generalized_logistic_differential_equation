program mod_random_number_generator_test
  use mod_random_number_generator
  use mod_stochastic_logistic_model
  implicit none
  real(kind=8) alpha0, m, sigma, x0, y_, y__
  real(kind=8) sigma0, m0, delta, sigmahat, gval
  integer npoints, nmc, niter, nobs, n, i, j, ntray, seed
!  parameter(sigma=0.1, alpha=0.8, m=2.0, x0=0.1)
! npoints is the number of observations
  parameter(npoints=1001, nmc=500, niter=50,&
	  seed=123, nobs=1001, ntray=1000)
  parameter(n=500, m0=1.5, x0=0.1)
  real(kind=8) rn_0, rn_1, rn_2
  real(kind=8) x(npoints), y1(npoints), y2(npoints)
  alpha0 = 0.8
  sigma0 = 0.1
  x(:) = 0.0
  call compute_drift(alpha0, m0, x0, y_)
  call compute_diffusion(sigma0, x0, y__)
  !> @todo vectorize these functions and tests
  !> @togo implement a better random number generartor  
  rn_1 = boxmuller()
  rn_2 = normalvar()
  print*, " "
  print *,"(====) testing mod_random_number_generator_test"
  print *,"(====)(----) testing, scalar f(x_0)= ", y_
  print *,"(====)(----) testing, scalar g(x_0)= ", y__
  ! print *,"(====)(++++) testing, scalar f(x[end]) ="
  ! print *,"(====)(++++) testing, scalar g(x[end]) ="
    
end program

! gfortran mod_random_number_generator.f95 main.f95 -o b.out
! gfortran mod_random_number_generator.f95 mod_stochastic_logistic_model.f95 main.f95 -o b.out
! gfortran mod_stochastic_logistic_model.f95 main.f95 -o b.out
