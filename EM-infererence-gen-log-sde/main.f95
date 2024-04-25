program main
  implicit none
  character(len=100) :: header
  real * 8 alpha, m, sigma, x0, alpha0, alphahat, alphahat_aux
  real * 8 sigma0, Infi, m0, delta, sigmahat, sigmahat_aux, gval
  real * 8 sum_, sum_aux, sum_sigma, sum_sigma_aux, sum_alpha, sum_alpha_aux
  real * 8 mean_m_hat, mean_m_hat_aux, b
  logical condition
  integer npoints, nmc, niter, nobs, n_sign, n, i, j, ntray, j_aux, flag
  character (len = 10) :: character_hline
  parameter(sigma=0.05, alpha=1.0, m=2, x0=0.15)
  parameter(npoints=10001, nmc=500, niter=50, nobs=2001, ntray=1000)
  parameter(n=500, m0=1.75)
  real * 8 mhat, mhat_aux, X(npoints), error, gs(51)
  real * 8 Xs(ntray, npoints), Xm(npoints), m_hat_vec(ntray), &
    m_hat_vec_aux(ntray)
  integer :: num_column
  character(len=100) :: filename = 'data/output.csv'
  character(len=1)   :: sep = ',' 
  header = 'i' // sep // 'm_hat' // sep // 'm_hat_aux' // sep // &
    'alpha_hat' // sep // 'alphahat_aux' // sep // &
    'sigma_hat' // sep // 'sigma_hat_aux'
  ! from paper
  flag = 0
  n_sign = 0
  error = 0.1
  delta = 10.0 / npoints
  !
  sum_ = 0.0
  sum_aux = 0.0
  sum_sigma = 0.0
  sum_sigma_aux = 0.0
  sum_alpha = 0.0
  sum_alpha_aux = 0.0
  
  j = 1
  j_aux = 1
  character_hline = "----------"
  call system('clear')
  print 1000, character_hline, character_hline, character_hline, &
    character_hline, character_hline, character_hline, character_hline, &
    character_hline, character_hline, character_hline, character_hline, &
    character_hline, character_hline, character_hline, character_hline, &
    character_hline, character_hline, character_hline
  1000 format(4x a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a)
  print 100 
  100 format(4x, "i", 9x,"m_hat", 10x, "m_hat_aux", &
    7x, "E(m_hat)", 8x, "E(m_hat_aux)", 3x, "sigma_hat", &
    7x, "sigma_hat_aux", 4x, "E(sigma_hat)", 4x, "E(sigma_hat_aux)", &
    7x, "alpha_hat", 7x, "alpha_hat_aux")
  call sim_gom(alpha, m, sigma, delta, x0, npoints, X)
  call Newtons2( 500, npoints, delta, m0, X, mhat)
  call Newtons2_aux( 500, npoints, delta, m0, X, mhat_aux)
  call MLE_sigma(m, X, delta, npoints, sigmahat, alphahat)
  call MLE_sigma_aux(m, X, delta, npoints, sigmahat_aux, alphahat_aux, b)
  mean_m_hat = mhat
  mean_m_hat_aux = mhat_aux
  sum_ = mhat
  sum_aux = mhat_aux
  open(unit=1, file=filename, status='replace', &
    action='write', form='formatted')
  write(1, *) header
  200 format(i3, 4X, f12.4, 4X, f12.4,&
     4X, f12.4, 4x, f12.4, 4x, f12.8, 4x, f12.8, &
      4X, f12.4, 4X, f12.4, 12x, f12.4, 4x, f12.4)
  300 format(i3, ',', f12.4, ',', f12.4, ',', &
      f12.4, ',', f12.4, ',', f12.8, ',', f12.8)
  do i=1, ntray
    call sim_gom(alpha, m, sigma, delta, x0, npoints, X)
    call Newtons2( 500, npoints, delta, m0, X, mhat)
    call Newtons2_aux( 500, npoints, delta, m0, X, mhat_aux)
    call MLE_sigma(mean_m_hat, X, delta, npoints, sigmahat, alphahat)
    call MLE_sigma_aux(mean_m_hat_aux, X, delta, npoints, sigmahat_aux, &
      alphahat_aux, b)
    m_hat_vec(i) = mhat
    if (.not.(isnan(mhat))) then
      sum_ = sum_ + mhat
      j = j + 1
      mean_m_hat = sum_ / j
    end if
    if (.not.(isnan(mhat_aux))) then
      sum_aux = sum_aux + mhat_aux
      j_aux = j_aux + 1
      mean_m_hat = sum_ / j
    end if
    condition = (sigmahat < 0.0)
    if (condition) then
      n_sign = n_sign + 1
    end if
    m_hat_vec_aux(i) = mhat_aux
    sum_sigma = sum_sigma + sigmahat
    sum_sigma_aux = sum_sigma_aux + sigmahat_aux
    sum_alpha = sum_alpha + alphahat
    sum_alpha_aux = sum_alpha_aux + alphahat_aux
    if (mod(i,100)==0) then
      print 200, i, mhat, mhat_aux,&
        sum_ / real(j), sum_aux / real(j_aux), &
        sigmahat, sigmahat_aux, &
        sum_sigma / i, sum_sigma_aux / i, sum_alpha / i, sum_alpha_aux / i 
    end if
    condition = (sigmahat_aux > 0)
    !condition = (b > 0.0)
    if(condition) then
      write(1, 300) i, mhat, mhat_aux, &    
        alphahat, alphahat_aux, &
        sigmahat, sigmahat_aux
    end if
    ! call progress(i) ! generate the progress bar.
    ! call sleep(1)
  end do
  close(unit=1)
  !
  !
  do j = 1, 51
    call g( j/10.0d0, X, delta, npoints, gval)
    gs(j) = gval
  end do
  ! initial condition for Newton method
  call Newtons2( 500, npoints, delta, m0, X, mhat)
  print*, mhat
  !
  ! data for plotting
  open(unit=10, file='datos.dat', status='replace')
  do i = 1, 51
    write(10,*) i / 10.0d0, gs(i)
  end do
  close(unit=10)
  open(unit=30, file='tray.dat', status='replace')
  do i = 1, npoints
    write(30,*) (i-1) * delta, Xm(i)
  end do
end program main
  ! Subrutines 
  ! TODO: Formated subrubtine
  subroutine MLE_B(delta,npoints,pathlam,pathlam2,sigmahat,alphahat)
    implicit none
    integer npoints,i
    real*8 delta,sigmahat,pathlam(npoints),pathlam2(npoints),alphahat,I1,I2,I3

    call ItoIntegrate(npoints,pathlam2,pathlam,I1)
    call Integrate(npoints,pathlam2,delta,I2)
    call Integrate(npoints,pathlam2(:)**2,delta,I3)

    alphahat=(sigmahat*I1+((sigmahat**2)*I2/2))/I3
    !(sigmaE*I1+((sigmaE**2)*I2/2))/I3

    return
end
! TODO: Formated subrubtine
subroutine random_stduniform(u)
   implicit none
   real*8,intent(out) :: u
   real*8 :: r
   call random_number(r)
   u = 1 - r
end subroutine random_stduniform

subroutine random_uniform(a,b,x)
   implicit none
   real*8,intent(in) :: a,b
   real*8,intent(out) :: x
   real*8 :: u
   call random_stduniform(u)
   x = (b-a)*u + a
end subroutine random_uniform

subroutine Qua_Var(npoints,path,delta,sigmahat)
  implicit none
  integer npoints,i
  real*8 sigmahat,path(npoints),delta,num,dem

    num=0.0
    dem=0.0

   do i=2,npoints
      num=num+(path(i)-path(i-1))**2
      dem=dem+path(i)**2+path(i-1)**2
  enddo

  sigmahat=SQRT(2.0*num/(dem*delta))

  return
end

subroutine LAM_GOM(npoints,sigma,m,path,pathlam,pathlam2)
   implicit none
   integer npoints,i
   real*8 sigma,path(npoints),pathlam(npoints),pathlam2(npoints),m
    pathlam(:)=log(path(:))/sigma
    do i=1,npoints
        pathlam2(i) = 1.0 - exp(m*sigma*pathlam(i))
     end do

     return
end


!> Display a porgress bar 
!! @param j: The j-iteration, j=1 ,..., N
!! @todo Handle a debuger flag to manage special cases 

subroutine progress(j) 
  implicit none
  integer(kind=4) :: j,k
  character(len=18) :: bar = "\r???% |          |"  
  ! updates the fraction of calculation done
  write(unit=bar(2:4), fmt="(i3)")  10 * j
    do k = 1, j
      bar(7 + k: 7 + k) = "*"
    enddo
  ! print the progress bar.
  write(*, '(a)', advance='no') bar 
  return
end subroutine progress

!> Find a root of the function g(m), with the Newton's method, see 
!! Burden's book to more details.
!! @param n 
!! @param npoints
!! @param delta
!! @param x0
!! @param X
!! m_root
!! @todo Handle a debuger flag to manage special cases 
subroutine Newtons2(n, npoints, delta, x0, X, xfin)
  implicit none
  integer n, i, npoints
  real * 8 puntos(n),x0,xfin,X(npoints),gval,dgval,e,delta
  real * 8 error, tol, relative_error
  logical condition1, condition2
  tol = 1e-8
  puntos(1)=x0
  do i=2,n
    call g(puntos(i-1),X,delta,npoints,gval)
    call dg(puntos(i-1),X,delta,npoints,dgval)
    puntos(i) = puntos(i-1)-gval/dgval
    error = abs(puntos(i) - puntos(i-1))
    condition1 = error < tol
    call g(puntos(i), X, delta, npoints, gval)
    condition2 = abs(gval) < tol
    if (isnan(puntos(i))) then
      print*, "Fucking NAN"
      stop
    end if
    if (condition2) then
      xfin = puntos(i)
      relative_error = (abs(2.0 - xfin) / 2.0) * 100
      !print*, "Newton converge  =)"
      !print*, "---------------------------------------------"
      !print*, "tol:", tol
      !print*, "Number of iterations, g(m_est): ", i, gval
      !print*, "m_est, relative_err(%): ", xfin, relative_error
      !print*, "============================================="
      return
    end if
    !print*,puntos(i)
  end do
  !print*, "Newton fails =("
  !print*, "---------------------------------------------"
  xfin = puntos(n)
  !relative_error = (abs(2.0 - xfin) / 2.0) * 100
  !print*, "Newton converge with tol less than", tol
  !print*, "Number of iterations: ", i, gval
  !print*, "m_estimate, abs_err: ", xfin, relative_error
  !print*, "============================================="
end

!> @brief inserts a value into an ordered array
!!
!! An array "list" consisting of n ascending ordered values. The method insert a
!! "new_entry" into the array.
!! hint: use cshift and eo-shift
!!
!! @param[in,out]   list    a real array, size: max_size
!! @param[in]       n       current values in the array
!! @param[in]       max_size    size if the array
!! @param[in]       new_entry   the value to insert
subroutine Newtons2_aux(n, npoints, delta, x0, X, xfin)
  implicit none
  integer n, i, npoints
  real * 8 puntos(n),x0,xfin,X(npoints),gval,dgval,e,delta
  real * 8 error, tol, relative_error
  logical condition1, condition2
  tol = 1e-8
  puntos(1)=x0
  do i=2,n
    call g_aux(puntos(i-1),X,delta,npoints,gval)
    call dg_aux(puntos(i-1),X,delta,npoints,dgval)
    puntos(i) = puntos(i-1)-gval/dgval
    error = abs(puntos(i) - puntos(i-1))
    condition1 = error < tol
    call g(puntos(i), X, delta, npoints, gval)
    condition2 = abs(gval) < tol
    if (isnan(puntos(i))) then
      print*, "Fucking NAN"
      stop
    end if
    if (condition2) then
      xfin = puntos(i)
      relative_error = (abs(2.0 - xfin) / 2.0) * 100
      !print*, "Newton converge  =)"
      !print*, "---------------------------------------------"
      !print*, "tol:", tol
      !print*, "Number of iterations, g(m_est): ", i, gval
      !print*, "m_est, relative_err(%): ", xfin, relative_error
      !print*, "============================================="
      return
    end if
    !print*,puntos(i)
  end do
  !print*, "Newton fails =("
  !print*, "---------------------------------------------"
  xfin = puntos(n)
  !relative_error = (abs(2.0 - xfin) / 2.0) * 100
  !print*, "Newton converge with tol less than", tol
  !print*, "Number of iterations: ", i, gval
  !print*, "m_estimate, abs_err: ", xfin, relative_error
  !print*, "============================================="
end subroutine Newtons2_aux


subroutine Newtons(lim,error,npoints,delta,X,xfin)
  implicit none
  integer npoints,flag
  real*8 puntoN,puntoA,x0,xfin,X(npoints),gval,dgval,e,delta,lim,error,dif
  10 call random_uniform(1.5d0,2.50d0,x0)
  puntoA=x0
  flag = 1
  do while(flag == 1)
    call g(puntoA,X,delta,npoints,gval)
    call dg(puntoA,X,delta,npoints,dgval)
    puntoN = puntoA-gval/dgval
    dif = abs(puntoA-puntoN)
    puntoA = puntoN
    !print*,dif,'  diferencia'
    if (dif<error) then
      call dg(puntoN,X,delta,npoints,dgval)
      !print*, dgval , " derivada"
      if (dgval< lim) then
        xfin  = puntoN
        flag = 0
      else
        goto 10
      end if
    end if
    !print*,i,gval,puntos(i-1),dgval
  end do
end

subroutine gfunc(npoints,X,m,g1,g2,g3,g4)
  implicit none
  integer npoints, i
  real*8 g1(npoints), g2(npoints), g3(npoints), g4(npoints), X(npoints), m
  do i=1, npoints
    g1(i) = (-X(i) ** m + 1.0) ** 2
    g2(i) = -(X(i) ** (m - 1.0)) * log(X(i))
    g3(i)=(-X(i)**m + 1.0 ) / X(i)
    g4(i)=(X(i) ** m)*(- X(i)**m + 1.0 ) * log(X(i))
  end do
end

subroutine g(m,X,delta,npoints,gval)
  integer npoints
  real*8 I1,I2,I3,I4,g1(npoints),g2(npoints),g3(npoints),g4(npoints),m,X(npoints),gval,delta
  call gfunc(npoints,X,m,g1,g2,g3,g4)
  call Integrate(npoints,g1,delta,I1)
  call ItoIntegrate(npoints,g2,X,I2)
  call ItoIntegrate(npoints,g3,X,I3)
  call Integrate(npoints,g4,delta,I4)
  gval = I1*I2+I3*I4
end
subroutine g_aux(m,X,delta,npoints,gval)
  integer npoints
  real*8 I1,I2,I3,I4,g1(npoints),g2(npoints),g3(npoints),g4(npoints),m,X(npoints),gval,delta
  call gfunc(npoints,X,m,g1,g2,g3,g4)
  call Integrate(npoints,g1,delta,I1)
  call ItoIntegrate_aux(npoints,g2,X,I2)
  call ItoIntegrate_aux(npoints,g3,X,I3)
  call Integrate(npoints,g4,delta,I4)
  gval = I1*I2+I3*I4
end subroutine g_aux

subroutine dgfunc(npoints, X, m, dg1, dg2, dg3, dg4, dg5, dg6, dg7, dg8)
  implicit none
  integer npoints, i
  real*8 dg1(npoints), dg2(npoints), dg3(npoints), dg4(npoints), &
         dg5(npoints), dg6(npoints), dg7(npoints), dg8(npoints), X(npoints), m
  do i=1, npoints
  !  print*,X(i)
    dg1(i) = 2.0 * ( X(i)**m) * (1.0 - X(i) ** m) * log( X(i) )
    dg2(i) = ( X(i) ** (m - 1.0)) * log(X(i))
    dg3(i) = (1.0 - X(i) ** m) ** 2.0
    dg4(i) = (X(i) ** (m - 1.0)) * (log(X(i))) ** 2
    dg5(i) = (X(i) ** (m - 1.0)) * log( X(i) )
    dg6(i) = (X(i) ** m) * (1.0 - X(i) ** m) * log(X(i))
    dg7(i) = (1.0 - X(i) ** m) / X(i)
    dg8(i) = (X(i) ** m) * (2.0 * (X(i) ** m) - 1.0) * (log(X(i))) ** 2 
  end do
end

subroutine dg(m,X,delta,npoints,dgval)
    integer npoints
    real * 8 I1, I2, I3, I4,I5, I6, I7, I8
    real * 8 dg1(npoints), dg2(npoints), dg3(npoints), dg4(npoints)
    real * 8 dg5(npoints), dg6(npoints), dg7(npoints), dg8(npoints)
    real * 8 m, X(npoints), dgval, delta
    call dgfunc(npoints,X,m,dg1,dg2,dg3,dg4,dg5,dg6,dg7,dg8)
    call Integrate(npoints,dg1,delta,I1)
    call ItoIntegrate(npoints,dg2,X,I2)
    call Integrate(npoints,dg3,delta,I3)
    call ItoIntegrate(npoints,dg4,X,I4)
    call ItoIntegrate(npoints,dg5,X,I5)
    call Integrate(npoints,dg6,delta,I6)
    call ItoIntegrate(npoints,dg7,X,I7)
    call Integrate(npoints,dg8,delta,I8)
    dgval = I1 * I2 - I3 * I4 - I5 * I6 - I7 * I8
end
subroutine dg_aux(m, X, delta, npoints, dgval)
    integer npoints
    real * 8 I1, I2, I3, I4,I5, I6, I7, I8
    real * 8 dg1(npoints), dg2(npoints), dg3(npoints), dg4(npoints)
    real * 8 dg5(npoints), dg6(npoints), dg7(npoints), dg8(npoints)
    real * 8 m, X(npoints), dgval, delta
    call dgfunc(npoints,X,m,dg1,dg2,dg3,dg4,dg5,dg6,dg7,dg8)
    call Integrate(npoints,dg1,delta,I1)
    call ItoIntegrate_aux(npoints,dg2,X,I2)
    call Integrate(npoints,dg3,delta,I3)
    call ItoIntegrate_aux(npoints,dg4,X,I4)
    call ItoIntegrate_aux(npoints,dg5,X,I5)
    call Integrate(npoints,dg6,delta,I6)
    call ItoIntegrate_aux(npoints,dg7,X,I7)
    call Integrate(npoints,dg8,delta,I8)
    dgval = I1 * I2 - I3 * I4 - I5 * I6 - I7 * I8
end

subroutine DriftParameter(alpha, m, x, y)
    implicit none
    real * 8 alpha, x, y, m
    y = alpha * x * (-1.0 * x ** m + 1)
  return
end

!  Evaluates y=difusion(x;sigma) (output) where sigma is a parameter
subroutine DiffusionParameter(sigma, x, y)
    implicit none
    real * 8 sigma, x, y
    y = sigma * x
    return
end

!  On increment in a Brownian motion.
! DIM is dimension  from
! startx is initial point
!  endx (output) is a point delta time ahead in the Brownian motion.

! TODO: Include a resolution step and used to define an operative step
subroutine BrownianStep(delta, startx, endx)
  implicit none
  real*8 delta, startx, endx, var
  call normalvar(var)
  endx = startx + sqrt(delta) * var
  return
end
! Milstein calculates for stepsize delta the next point in a diffusion.
! startx is the current location of the diffusion
! drix is the drift evaluated at startx
! drifx is the derivative of drift evaluated at startx
! sigma is the diffusion coefficient at startx.
! endx is the next point in the diffusion (output).
! W is the  Brownian increment (output).
subroutine MilsteinStep(delta,startx,drix,difx,sigma,endx)
    implicit none
    integer i
    real * 8 delta, startx, endx, drix, difx, der
    real * 8 sigma, W, xx, brow
    xx = 0.0
    call BrownianStep(delta, xx, brow)
    der = sigma
    W = brow
    endx = startx + drix * delta + difx * W & 
        + 0.5 * difx * der * (W ** 2 - delta)
    return
end
! Using the Milstein scheme we simulate a diffusion from x
! npoints steps ahead at stepsizes delta.
! brownian (output) contains the original brownian increments used to construct the diffusion.
! path (output) is the path of diffusion
subroutine sim_gom(alpha, m, sigma, delta, x, npoints, path)
  implicit none
  integer npoints, i, DIM, j
  real * 8 delta, x, path(npoints), y1 ,y2 ,z
  real * 8 alpha, m, sigma, brow, brownian(npoints)
  path(1) = x
  do i=2, npoints
    call DriftParameter(alpha, m, path(i-1), y1)
    call DiffusionParameter(sigma, path(i-1), y2)
    call MilsteinStep(delta, path(i-1), y1, y2, sigma, path(i))
  end do
  return
end
!   Generic routine to generate a uniform random number.
FUNCTION UNIF(IX)
    call random_number(x)
    UNIF = x
    RETURN
END

function boxmuller(seed)
    integer seed
    real boxmuller
    integer iset
    real fac, gset, rsq, v1, v2, unif
    save iset, gset
    data iset/0/
    1 if (iset.eq.0) then
        v1 = 2. * unif(seed) - 1.
        v2 = 2. * unif(seed) - 1.
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
end
  ! Generic routine to generate a standard normal variable.
subroutine normalvar(x)
    real * 8 x
    integer seed
    common seed
    x=boxmuller(seed)
    return
end

! TODO: Check quadratures
! Calculates the Ito integral (output) of path with respet to W from 0 to 
! npoints*delta
subroutine ItoIntegrate(npoints, path, W, integral)
    implicit none
    integer npoints, i
    real * 8 path(npoints), integral, W(npoints)
    integral = 0.0
    do i=2, npoints
     integral = integral + (W(i)-W(i-1)) * (path(i-1) + path(i)) / 2.0
   enddo
   return
end subroutine ItoIntegrate

subroutine ItoIntegrate_aux(npoints, path, W, integral)
    implicit none
    integer npoints, i
    real * 8 path(npoints), integral, W(npoints)
    integral = 0.0
    do i=2, npoints
     integral = integral + (W(i)-W(i-1)) * path(i-1) 
   enddo
   return
end subroutine ItoIntegrate_aux

! Calculates the Integral (output) of path from 0 to npoints*delta
subroutine Integrate(npoints, path, delta, integral)
    implicit none
    integer npoints,i
    real * 8 delta, path(npoints), integral
    integral = 0.0
    do i=2, npoints
     integral = integral + delta * 0.5 * (path(i-1)+path(i))
   enddo
   return
end subroutine Integrate

subroutine OU(npoints, path, pOU)
    implicit none
    integer npoints
    real*8 path(npoints), pOU(npoints)
    pOU(:)=LOG(path(:))
    return
 end

subroutine SCS(m,X,delta,npoints,CS)
  integer npoints
  real*8 CS(6),m,X(npoints),delta
  call ItoIntegrate(npoints,1.0/X,X,CS(1))
  call ItoIntegrate(npoints,(1.0-X**m)/X,X,CS(2))
  call Integrate(npoints,(1.0-X**m),delta,CS(3))
  call ItoIntegrate(npoints,LOG(X)*(X**(m-1.0)),X,CS(4))
  call Integrate(npoints,LOG(X)*(1.0-X**m)*(X**m),delta,CS(5))
  call Integrate(npoints,(1.0-X**(m))**2,delta,CS(6))
  return
end subroutine SCS

subroutine MLE_sigma(m,X,delta,npoints,sigmahat,alphahat)
  integer npoints
  real*8 CS(10),m,X(npoints),delta,sigmahat,a,b,alphahat
  call SCS(m,X,delta,npoints,CS)
  a = ( CS(2) * CS(3)) / (2.0 * CS(6)) - CS(1) / 2.0
  b = -(CS(2) * CS(4) * m) / CS(6) + ((CS(2) ** 2) * CS(5) * m) / (CS(6) ** 2) 
  ! print*, a, b
  alphahat=CS(2)/CS(6)
  sigmahat=-b/a
  return
end subroutine MLE_sigma

subroutine MLE_sigma_aux(m, X, delta, npoints, sigmahat, alphahat , b)
  integer npoints
  real*8 CS(10), m, X(npoints), delta, sigmahat, a, b, alphahat
  call SCS_aux(m,X,delta,npoints,CS)
  a = (CS(2) * CS(3)) / (2.0 * CS(6)) - CS(1) / 2.0
  b = ((CS(2) **2 ) * CS(5) * m) / (CS(6) ** 2) - (CS(2) * CS(4) * m) / CS(6)
  alphahat = CS(2) / CS(6)
  ! print*, a, b
  sigmahat=-b/a
  return
end subroutine MLE_sigma_aux

subroutine SCS_aux(m, X, delta, npoints, CS)
  integer npoints
  real * 8 CS(6), m, X(npoints), delta
  call ItoIntegrate_aux(npoints, 1.0/X, X, CS(1))
  call ItoIntegrate_aux(npoints, (1.0-X**m)/X,X,CS(2))
  call Integrate(npoints,(1.0-X**m),delta,CS(3))
  call ItoIntegrate_aux(npoints,LOG(X)*(X**(m-1.0)),X,CS(4))
  call Integrate(npoints,LOG(X)*(1.0-X**m)*(X**m),delta,CS(5))
  call Integrate(npoints,(1.0-X**(m))**2,delta,CS(6))
  return
end subroutine SCS_aux
