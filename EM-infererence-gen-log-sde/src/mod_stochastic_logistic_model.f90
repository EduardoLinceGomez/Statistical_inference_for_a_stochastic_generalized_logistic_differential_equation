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
!> @details 

!> @version 1.0
module mod_stochastic_logistic_model
    use iso_fortran_env
    implicit none
contains
!> @brief <BRIEF_DESCRIPTION>
!> @details <DETAILED_DESCRIPTION>
!> @param [<in or out or inout>] <PARAM1> <DESCRIPTION>
!> @param [<in or out or inout>] <PARAM2> <DESCRIPTION>
    subroutine compute_drift(alpha, m, x, y)
        implicit none
        real(kind=8), intent(in) :: alpha, m, x
        real(kind=8), intent(out) :: y
        y = alpha * x * (-1.0 * x ** m + 1.0)
        ! print*, "off-line alpha, m, X_t, f(X_t)", alpha, m, x, y
    end subroutine compute_drift
!> @brief <BRIEF_DESCRIPTION>
!> @details <DETAILED_DESCRIPTION>
!> @param [<in or out or inout>] <PARAM1> <DESCRIPTION>
!> @param [<in or out or inout>] <PARAM2> <DESCRIPTION>
    subroutine compute_diffusion(sigma, x, y)
        implicit none
        real(kind=8), intent(in) :: sigma, x
        real(kind=8), intent(out) :: y
        y = sigma * x
    end subroutine compute_diffusion
end module mod_stochastic_logistic_model
!> @brief <BRIEF_DESCRIPTION>
!> @details <DETAILED_DESCRIPTION>
!> @param [<in or out or inout>] <PARAM1> <DESCRIPTION>
!> @param [<in or out or inout>] <PARAM2> <DESCRIPTION>
