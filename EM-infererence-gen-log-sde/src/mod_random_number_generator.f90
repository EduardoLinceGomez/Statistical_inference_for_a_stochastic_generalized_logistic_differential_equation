!> @ingroup modules
!> @author E. Lince-Gomez, F. Baltazar-Larios, S. Diaz-Infante
!> @brief This module implements the Box-Muller algorithm to generate
!> random variables with standard Gaussian distribution from a
!> uniform ditributed random variable. This module enclose tree functions.
!> @see Kloeden & Platen 1992
include 'mkl_vsl.f90'
module mod_random_number_generator
    use mkl_vsl_type
    use mkl_vsl
    use iso_fortran_env, only : int32, real32
    implicit none
contains

!> @brief Returns a random variables with uniform distribution using the
!! standar gfortran random number generator, the returned value is a real 32.
!> @param [in] ix  Dummy paramer for the seed initialization for the random
!> generator
    real(real32) function unif() result(unif_)
        !real(real32), intent(out) :: u
        ! call random_seed()
        call random_number(unif_)
    end function unif
!
!> @brief Implementation of the Box-muller algorithm to genrate Gaussian random
!>  variables
!> @param [in ] seed  initial value for the uniform random generator
!> @todo: Try other random-number generators, such like mersene an others
    real(real32) function boxmuller()
        integer(int32) iset
        real fac, gset, rsq, v1, v2, seed
        save iset, gset
        data iset/0/
        1 if (iset.eq.0) then
            call random_number(seed)
            v1 = 2. * seed - 1.0
            call random_number(seed)
            v2 = 2. * seed - 1.0
            rsq = v1 ** 2 + v2 ** 2
            if (rsq.ge.1..or.rsq.eq.0)goto 1
                fac = sqrt(-2.0 * log(rsq) / rsq)
                gset = v1 * fac
                boxmuller = v2 * fac
                iset = 1
            else
                boxmuller=gset
                iset=0
            endif
          return
    end function boxmuller
!> @brief Returns a number with Gaussian distribution using the Box-Muller
!> algortihm.
!> @details This function can be ommited in the case of a step forward
!>  implelemtation
!> @param [in] seed a int32 with the initial value for the randim uniform
!> generator
!> @todo: Implement this furnction such that returns a realization path of the
!> standard Brownian motion
    real(real32) function normalvar()  result(r_x)
        r_x = boxmuller()
    end function normalvar

    subroutine mkl_gaussian_sampler(&
        &buffer_size, &
        &mean_a, &
        &std_a, &
        &user_seed, &
        &gaussian_sample, &
        &debug&
    &)
        integer, intent(in) :: buffer_size
        real(kind=8), intent(in) :: mean_a, std_a
        real(kind=8), intent(out), allocatable :: gaussian_sample(:)
        logical, intent(in), optional :: debug
        TYPE (VSL_STREAM_STATE) :: stream
        integer(kind=4) errcode
        integer(kind=4) i, j
        integer, intent(in), optional :: user_seed
        integer brng, seed, method, n
        if (present(user_seed)) then
            seed = user_seed
        else
            seed = 123
        end if
        gaussian_sample = [(0.0, i=1, buffer_size)]
    !       ***** Initializing *****
        brng = VSL_BRNG_MT19937
        method = VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2
        errcode = vslNewStream(stream, brng, seed)
        errcode = vdRngGaussian(&
            method, stream, buffer_size, gaussian_sample, mean_a, std_a&
        )
        if (present(debug)) then
            if (debug) then
              print *,"Error structure = ", stream
              print *, "Gaussian stream: ", gaussian_sample(1:5)
            end if
        end if
    !       ***** Deinitialize *****
        errcode = vslDeleteStream(stream)
    end subroutine mkl_gaussian_sampler
end module mod_random_number_generator
! ifort -c -i8 -I"${MKLROOT}/include" mod_random_number_generator.f90
