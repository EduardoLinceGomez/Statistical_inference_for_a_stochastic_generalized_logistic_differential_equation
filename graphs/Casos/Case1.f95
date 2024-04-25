program main
! Variables de entrada, sigma, alpha, m
! Crear una trayectoria, quitarle los datos y hacer monte carlo con los puentes
! nPuentes es la cantdad de puentes generados entre cada punto para hacer Monte Carlo
! nPuntos es la cantidad de puntos en la trayectoria original(continua)
! nObservaciones es la cantidad de puntos en la trayectoria con datos faltantes
integer nPuentes,nPuntos,nObservaciones,nObservaciones2,nPuntosPuentes,nPuntosPuentes2,nPuntosTR,nEM,i,j,k,l,n,doc1,doc2
real*8 sigma,alpha,m,delta,x0,m0,sigmaR_hat,mR_hat,alphaR_hat,tol
parameter(nPuentes=100, nPuntos=10001, nObservaciones=1001,nObservaciones2=101,nEM = 10)
parameter(nPuntosPuentes=(nPuntos-1)/(nObservaciones-1),nPuntosPuentes2=(nPuntos-1)/(nObservaciones2-1))
parameter(sigma=0.01, alpha=0.70, m=0.6 ,m0=1.5, x0=0.05,delta = 10.0 / nPuntos, tol=0.002)! sigma menor a 10^-2
real*8 T(nPuntos) ! Trayectoria continua
real*8 TR(nPuntos) ! Trayectoria reconstruida con los puentes
real*8 sigmas(nEM),ms(nEM),alphas(nEM)
real*8 TF(nObservaciones),TF2(nObservaciones2) ! Trayectorias con datos faltantes
real*8 puentesM(nPuentes,nObservaciones-1,nPuntosPuentes),puentesM2(nPuentes,nObservaciones2-1,nPuntosPuentes2)
logical flag
integer seed
character(len=100) :: header
character(len=25) :: header0
character(len=100) :: filename00 = 'EstimadoresContinuoCase1.csv'
character(len=100) :: filename11 = 'EMoutputDiezXCase1.csv'
character(len=100) :: filename12 = 'EMestimationDiezXCase1.csv'
character(len=100) :: filename21 = 'EMoutputUnoXCase1.csv'
character(len=100) :: filename22 = 'EMestimationUnoXCase1.csv'
character(len=1)   :: sep = ','
header ='sigma_hat' // sep // 'alpha_hat' // sep // &
  'm_hat' // sep // 'sigma_hatError' // sep // 'alpha_hatError' // sep // &
  'm_hatError'

header0 ='sigma_hat' // sep // 'alpha_hat' // sep //'m_hat'

seed = 2023 ! variable global

call RANDOM_SEED(seed)
open(unit=1, file=filename11, status='replace', &
  action='write', form='formatted',position="append")
close(unit=1)
open(unit=1, file=filename11, status='old', &
  action='write', form='formatted',position="append")
open(unit=2, file=filename12, status='replace', &
  action='write', form='formatted',position="append")
close(unit=2)
open(unit=2, file=filename12, status='old', &
    action='write', form='formatted',position="append")
open(unit=3, file=filename21, status='replace', &
  action='write', form='formatted',position="append")
close(unit=3)
open(unit=3, file=filename21, status='old', &
  action='write', form='formatted',position="append")
open(unit=4, file=filename22, status='replace', &
  action='write', form='formatted',position="append")
close(unit=4)
open(unit=4, file=filename22, status='old', &
  action='write', form='formatted',position="append")
open(unit=5, file=filename00, status='replace', &
  action='write', form='formatted',position="append")
close(unit=5)
open(unit=5, file=filename00, status='old', &
  action='write', form='formatted',position="append")

! write(1, *) header
write(2, *) header
! write(3, *) header
write(4, *) header
write(5, *) header0

do n=1,1000
! Creamos la trayectoria continua
12  call sim_gom(alpha, m, sigma, delta, x0, nPuntos, T)

! Creamos la trayectoria con datos faltantes
TF(1:nObservaciones)=T(1:nPuntos:nPuntosPuentes)
TF2(1:nObservaciones2)=T(1:nPuntos:nPuntosPuentes2)

! Creamos los puentes
call Qua_Var(nPuntos,T,delta,sigmaR_hat)
call Newtons2_aux( 500, nPuntos, delta, m0, T, mR_hat,flag)
if (flag .eqv. .false.) then
  go to 12
end if
call MLE_sigma_aux(mR_hat, T, delta, nPuntos, alphaR_hat)
! Checamos por valores NaN para imprimir y guardarS estimadores
if(((sigmaR_hat /= sigmaR_hat).or.(alphaR_hat /= alphaR_hat)).or.(mR_hat /= mR_hat)) then
  go to 12
else
  print*, trim("Estimaciones con trayectoria original")
  print*, trim("sigma:"),sigmaR_hat,trim("alpha:"),alphaR_hat,trim("m:"),mR_hat
  write(5, *) sigmaR_hat,sep,alphaR_hat,sep,mR_hat
end if
print*, trim("Comienzo EM"),n
write(1, *) "sigma:",sep,sigmaR_hat,sep,"alpha:",sep,alphaR_hat,sep,"m:",sep,mR_hat
write(1, *) header
write(3, *) "sigma:",sep,sigmaR_hat,sep,"alpha:",sep,alphaR_hat,sep,"m:",sep,mR_hat
write(3, *) header

print*, "10% de datos"
call EM_algorithm(1,2,flag,nEM,nPuntos,nObservaciones,nPuentes,nPuntosPuentes,puentesM,&
T,TR,TF,delta,sigmas,ms,alphas,m0)

print*, "1% de datos"
call EM_algorithm(3,4,flag,nEM,nPuntos,nObservaciones2,nPuentes,nPuntosPuentes2,puentesM2,&
T,TR,TF2,delta,sigmas,ms,alphas,m0)

end do
close(unit=1)
close(unit=2)
close(unit=3)
close(unit=4)
close(unit=5)
end program


subroutine EM_algorithm(doc1,doc2,flag,nEM,nPuntos,nObservaciones,nPuentes,nPuntosPuentes,puentesM,&
  T,TR,TF,delta,sigmas,ms,alphas,m0)
  integer doc1,doc2,nEM,i,j,k,l,nObservaciones,nPuentes,nPuntos
  real*8 delta,m0
  real*8 T(nPuntos)
  real*8 TR(nPuntos)
  real*8 puentesM(nPuentes,nObservaciones-1,nPuntosPuentes)
  real*8 TF(nObservaciones)
  real*8 sigmas(nEM),ms(nEM),alphas(nEM)
  logical flag

    do l = 1, nEM !!! 8 enero
        if (l==1) then
          ! Creamos la trayectoria con datos faltantes
          call Qua_Var(nObservaciones,TF,delta,sigmas(l))
          call Newtons2_aux( 500, nObservaciones, delta, m0, TF, ms(l),flag)
          call MLE_sigma_aux(ms(l), T, delta, nPuntos, alphas(l))
        else
          ! Generamos puentes para Monte Carlo
          11 do i =1, nPuentes
            do j =1,nObservaciones-1
              10 call puentes(TF(j),TF(j+1),sigmas(l-1),ms(l-1),alphas(l-1),delta,nPuntosPuentes,puentesM(i,j,1:nPuntosPuentes))
              do k=1,nPuntosPuentes
                if ( isnan(puentesM(i,j,k))) then  !
                  go to 10
                end if
              end do
              !print*,"puente creado en observacion",j
            end do
            !print*,"puente creado",i
          end do
          ! Reconstruimos la trayectoria haciendo Monte Carlo en los puentes
          do j=1,nObservaciones-1
            do k=1,nPuntosPuentes
              TR((j-1)*nPuntosPuentes+k)=sum(puentesM(:,j,k))/nPuentes
            end do
          end do
          TR(nPuntos)=T(nPuntos)

          call Qua_Var(nPuntos,TR,delta,sigmas(l))
          call Newtons2_aux( 500, nPuntos, delta, m0, TR, ms(l),flag)
          if (flag .eqv. .false.) then
            go to 11
          end if
          call MLE_sigma_aux(ms(l), TR, delta, nPuntos, alphas(l))
        end if
        print*,l, trim("sigma:"),sigmas(l),trim("alpha:"),alphas(l),trim("m:"),ms(l)
        write(doc1,*)sigmas(l),",",alphas(l),",",ms(l),",",sigmas(l)-sigma,",",alphas(l)-alpha,",",ms(l)-m
        !l=l+1 !!!!!!! 8 enero
        if (l==nEM) then
          write(doc2,*)sigmas(l),",",alphas(l),",",ms(l),",",sigmas(l)-sigma,",",alphas(l)-alpha,",",ms(l)-m
        end if
    end do
end subroutine

subroutine puentes(a, b, sigma, m, alpha, delta, n, puente)
  implicit none
  real(8), intent(in) :: a, b, sigma, m, alpha, delta
  integer, intent(in) :: n
  real(8), intent(out) :: puente(n)
  real(8) :: Ta(n), Tb(n)
  logical :: caso
  integer :: i

  ! Declaración de constantes
  real(8), parameter :: TOLERANCIA = 1.0E-1

  ! Llamada a sim_gom para Ta
  call sim_gom(alpha, m, sigma, delta, a, n, Ta)

  ! Bucle de cálculo
  do while (.true.)
    ! Llamada a sim_gom para Tb
    call sim_gom(alpha, m, sigma, delta, b, n, Tb)

    ! Invertir Tb
    Tb(:) = Tb(n:1:-1)

    ! Comparación de tiempos
    caso = Tb(1) > Ta(1)

    ! Construcción del puente
    do i = 1, n
      if ((caso .and. Tb(i) > Ta(i)) .or. (.not.caso .and. Tb(i) < Ta(i))) then
        puente(i) = Ta(i)
      else
        puente(i) = Tb(i)
      end if
    end do

    ! Comprobación de condición de salida
    if (abs(puente(n) - b) < TOLERANCIA) exit
  end do

end subroutine puentes


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
subroutine Newtons2_aux(n, npoints, delta, x0, X, xfin,flag)
  implicit none
  integer n, i, npoints
  real * 8 puntos(n),x0,xfin,X(npoints),gval,dgval,e,delta
  real * 8 error, tol, relative_error
  logical condition1, condition2,flag
  flag = .true.
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
      print*, "NAN value"
      flag = .false.
      return
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
subroutine MLE_sigma_aux(m, X, delta, npoints, alphahat )
  integer npoints
  real*8 CS(10), m, X(npoints), delta, a, b, alphahat
  call SCS_aux(m,X,delta,npoints,CS)
  a = (CS(2) * CS(3)) / (2.0 * CS(6)) - CS(1) / 2.0
  b = ((CS(2) **2 ) * CS(5) * m) / (CS(6) ** 2) - (CS(2) * CS(4) * m) / CS(6)
  alphahat = CS(2) / CS(6)
  ! print*, a, b
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


subroutine random_stdnormal(x)
implicit none
real*8:: x
real,parameter :: pi=3.14159265
real :: u1,u2
call random_number(u1)
call random_number(u2)
x = sqrt(-2*log(u1))*cos(2*pi*u2)
end subroutine random_stdnormal

subroutine seed(k)
implicit none
integer,DIMENSION(33):: s
integer :: k,i
do i=1,33
  s(i)=k
end do
call RANDOM_SEED(put=s)
end subroutine seed


subroutine normalvar(x,mean,s)
  implicit none
	real*8 x,s,n,mean
	call random_stdnormal(n)
	x = mean + n*sqrt(s)
	return
end subroutine normalvar

! TODO: Include a resolution step and used to define an operative step
subroutine BrownianStep(delta, startx, endx)
  implicit none
  real*8 delta, startx, endx, var
  call normalvar(var,0.0d0,1.0d0)
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
