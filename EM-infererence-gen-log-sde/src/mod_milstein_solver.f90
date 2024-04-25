!> @author E. Lince-Gomez, F. Baltazar-Larios, S. Diaz-Infante
!> @brief This module compute the coeffcients of the underlyng 
!> generalized logistc SDE
!>  \f{eqnarray*}{
!>        dX(t) &=& f(X(t)) dt + g(X(t)) dB(t) 
!>        \\
!>        f(x) &:=& \alpha x \Bigg[ 1 - \Bigg(\frac{x}{K} \Bigg) ^ m \Bigg] 
!>        \\
!>        g(x) &:=& \sigma x,\qquad t>t_0 \\
!>        x_0 &=& X(0), \quad x_0 \in [0, 1] .
!> \f} 
!> 
!> @see the manuscript (...) for more details.
!> @details (details)

module mod_milstein_solver
    use iso_fortran_env, only: int32, real32
    use mod_random_number_generator
    implicit none
contains

!> @brief Simulate an increment of the related Browninar path
!> @details This can be done by simulating the whole path
!> @param [in] delta determininstic time step size
!> @param [in] brownian_start current observation of the Brownian process
!> @param [in] brownian_end observation at current time plus step-size
    subroutine get_brownian_increment(delta, brownian_start, brownian_end)
        implicit none
        real(kind=8), intent(in) :: delta, brownian_start
        real(kind=8), intent(out) :: brownian_end
        real(kind=8) var
        var = boxmuller()
        brownian_end = brownian_start + sqrt(delta) * var
    end subroutine get_brownian_increment

!> @brief Computes a iteration of the Milstein method
!> @details Defintion taken from the Kloedens book 
!> @param [in ] delta determininstic time step-size
!> @param [in] brownian_start current observation of the Brownian process
!> the  Brownian increment.
subroutine get_milstein_iteration(delta, &
    current_state, current_drift, current_diffusion, &
    diffusion_derivative, sigma, x_next)
    implicit none
    real(kind=8), intent(in) :: delta, current_state, &
        current_drift, current_diffusion, &
        diffusion_derivative, sigma 
    real(kind=8), intent(out) :: x_next
    real(kind=8) euler_iteration, xx, brownian_increment
    xx = 0.0
    call get_brownian_increment(delta, xx, brownian_increment)
    euler_iteration = current_state + delta * current_drift &
        + sigma * current_state * brownian_increment 
    
    x_next = euler_iteration  &
        + 0.5 * current_diffusion * diffusion_derivative * &
        (brownian_increment * brownian_increment - delta)
end subroutine get_milstein_iteration
end module mod_milstein_solver
