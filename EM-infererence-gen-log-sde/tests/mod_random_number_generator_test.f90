program mod_random_number_generator_test
  use mod_random_number_generator
!  use mkl_vsl_type
!  use mkl_vsl
!  use iso_fortran_env, only : kind=4, kind=8
! use mod_stochastic_logistic_model
  real(kind=8) :: alpha0, sigma0
  integer :: buffer_size
  integer(kind=4) :: i, j
  integer(kind=4), parameter :: seed = 123
  real(kind=8) ::  rn_0, rn_1, rn_2, mean_a
  real(kind=8), allocatable:: gaussian_sample(:)
  real(kind=8) :: alpha00

  alpha0 = 0.8
  sigma0 = 0.1
  rn_0 = 0.0
  rn_1 = 1.0
  rn_2 = 2.0
  buffer_size = 1000
  mean_a = 0.0
  alpha00 = 0.0
  gaussian_sample = [(alpha00, i=1, buffer_size)]
  call random_number(rn_0) !unif(seed)
  rn_1 = boxmuller()
  rn_2 = normalvar()
  call mkl_gaussian_sampler(buffer_size, mean_a, sigma0, 123, gaussian_sample)
  !
  call system('clear')
  print*, " "
  print *,"(====) testing mod_random_number_generator_test"
  print*, "(==) Seed:=", seed
  print *, '(====)(----) unif: rn =', rn_0
  print *, '(====)(----) boxmuller: rn =', rn_1
  print *, '(====)(----) normalval: rn =', rn_2
  print *, '(====)(----) initialized gaussian_sample =', gaussian_sample(900:buffer_size)
end program

! gfortran mod_random_number_generator.f95 main_.f95 -o b.out
! gfortran mod_random_number_generator.f95 mod_stochastic_logistic_model.f95 main.f95 -o b.out
! gfortran mod_stochastic_logistic_model.f95 main.f95 -o b.out
!fort *.o -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
