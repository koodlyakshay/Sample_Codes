! Solve 1-D scalar advection equation
! \partial_t \phi + U \partial_x \phi - \mu \partial^2_x \phi = 0
! Use FVM integration
! \int_{\Omega} \partial_t \phi d\Omega = - (U\int_{\Omega}\partial_x \phi d\Omega - \mu\int_{\Omega}\partial^2_x \phi d\Omega)
! Choose desired numerical scheme for spatial disc term Res
! \int_{\Omega} \partial_t \phi d\Omega = Res
! Use mid point formula to approximate integral on LHS
! \partial_t \phi \Delta \Omega = Res
! Time discretization using RK4

!Ref eqn 0 <= x <= 1
! phi_t + 1.0 phi_x - 0.01 phi_xx = 0
!I.C
! phi(x,0) = exp(-((x-2)^2)/8)
!B.Cs
! phi(0,t) = 0.025/sqrt(0.000625+0.02t)exp(-(0.5-t)^2/(0.000125+0.04t))
! phi(1,t) = 0.025/sqrt(0.000625+0.02t)exp(-(1.5-t)^2/(0.000125+0.04t))
! Exact sol
! phi(x,t) = 0.025/sqrt(0.000625+0.02t)exp(-(x+0.5-t)^2/(0.000125+0.04t))
module globalvar

implicit none

integer, parameter  :: Nmax=257, scheme=4
integer             :: Nx
real                :: x(Nmax), xf(Nmax+1), phi(Nmax), phi_n(Nmax), exact(Nmax), dx(Nmax), vol(Nmax), dphi_n(Nmax)
real                :: dt, U, Lmax, Lx, Tf, t, theta, bc1, bc0, exp_ratio, init_dx 
real                :: phi_nvd, beta, gamma, eps, err, pi, mu 

! scheme - type of convective scheme
!           1 - first order upwind
!           2 - second order upwind
!           3 - QUICK
!           4 - central
!           5 - bounded central 
!           6 - FROMM
!           7 - vanAlbada
!           8 - Venkatakrishnan


end module globalvar

program onedscalar


use globalvar

implicit none

integer   :: i,j,k,stage,nstage
real      :: phi_e, phi_w, domega
real      :: resid(Nmax),rkstage(Nmax,4),rkcoeff(3)
    

Tf = 1.0      ! final time
dt = 0.0005   ! time step
beta = 0.1    ! lower limit
U = 1.0       ! advection velocity
mu = 0.01     ! diffusion coefficient




rkcoeff = (/5.d-1, 5.d-1, 1.d0/) ! RK coefficients
nstage = 4


call make_mesh
print*, 'CFL = ',U*dt/minval(dx(1:Nx-1))
domega = dx(2) ! uniform grid for now (Volume)
eps = init_dx/100.0
pi = 4.0*atan(1.0)
t = 0.0
call apply_IC


do while (t .le. Tf)
    !BC
    bc0 = 0.025/sqrt(0.000625+0.02*t)*exp(-(0.5-t)**2/(0.00125+0.04*t)) ! boundary condition
    bc1 = 0.025/sqrt(0.000625+0.02*t)*exp(-(1.5-t)**2/(0.00125+0.04*t)) ! boundary condition
    phi(1) = bc0
    phi_n(1) = bc0
    phi(Nx) = bc1
    phi_n(Nx) = bc1
    
    phi = phi_n

    !--- Time Integration (RK4) ---!
   
    !-- Stage 1 --!
    call compute_residual(resid)
    rkstage(:,1) = dt*resid(:)/(1.0*vol(:))
   
    !-- Stage 2-4 --!
    do stage = 2,nstage

      phi(:) = phi_n(:) - rkcoeff(stage-1)*rkstage(:,stage-1)

      call compute_residual(resid)

      rkstage(:,stage) = dt*resid(:)/(1.0*vol(:))
    enddo

    phi_n(:) = phi_n(:) - 1.d0/6.d0*(rkstage(:,1) + 2.0*rkstage(:,2) + 2.0*rkstage(:,3) + rkstage(:,4))
    
    call compute_gradient
    
    t = t + dt
    
enddo
print*,t,Tf/dt

call writeoutput


end program onedscalar

subroutine compute_residual(res)

! k - index of current element i.e "i" in the main loop
! res - residual of the current element

use globalvar

implicit none

integer                 :: k
real, intent (out)      :: res(Nmax)
real                    :: phi_e, phi_w, res_v(Nmax)

! Advection term
! Assuming U > 0
!1-D grid arrangement, k is given as input from the main routine, ('o' is a node, '|' is a face)
!k-1    (w)        k          (e)   k+1
! o------|---------o-----------|-----o

select case (scheme)
  case (1)
  ! First order upwind
    do k=2,Nx
      phi_e = phi_n(k)
      phi_w = phi_n(k-1)
      res(k) = phi_e - phi_w
    enddo
  case (2)
  ! Second order upwind
    do k=2,Nx
      theta = 0.0
      call general_formula(k,theta,phi_e)
      call general_formula(k-1,theta,phi_w)
      res(k) = phi_e - phi_w
    enddo
  case (3)
  ! QUICK
    do k=2,Nx
      theta = 1.0/8.0
      call general_formula(k,theta,phi_e)
      call general_formula(k-1,theta,phi_w)
      res(k) = phi_e - phi_w
    enddo
  case (4)
  ! Central
    do k=2,Nx
      theta = 1.0
      call general_formula(k,theta,phi_e)
      call general_formula(k-1,theta,phi_w)
      res(k) = phi_e - phi_w
    enddo
  case (5)
  ! Bounded central
    do k=2,Nx
      call bounded_central_scheme(k,phi_e)
      call bounded_central_scheme(k-1,phi_w)
      res(k) = phi_e - phi_w
    enddo
  case (6)
  ! FROMM
    do k=2,Nx
      theta = 0.5
      call general_formula(k,theta,phi_e)
      call general_formula(k-1,theta,phi_w)
      res(k) = phi_e - phi_w
    enddo
  case (7)
  ! Second order upwind + van albada
    do k=2,Nx
      call vanalbada(k,phi_e)
      call vanalbada(k-1,phi_w)
      res(k) = phi_e - phi_w
    enddo
  case (8)
  ! Second order upwind + venkatakrishnan
    do k=2,Nx
      call venkatakrishnan(k,phi_e)
      call venkatakrishnan(k-1,phi_w)
      res(k) = phi_e - phi_w
    enddo
end select

! Calculate diffusion contribution
call compute_diffusion_term(res_v)

! Add diffusion contribution to residual (note minus sign)
res = res - mu*res_v ! Constant mu

! Boundaries fixed
res(1) = 0.0
res(Nx) = 0.0

end subroutine compute_residual

subroutine general_formula(k,theta,phi_f)

use globalvar, only: phi_n,xf, bc0, bc1, Nx

implicit none

integer, intent (in)        :: k
integer                     :: km1,kp1
real, intent(in)            :: theta
real, intent(out)           :: phi_f
real                        :: Su,Sd,Sc,phi_c,phi_l,phi_r

km1 = k - 1
kp1 = k + 1

phi_c = phi_n(k)

if (kp1 .gt. Nx) then
    phi_r = bc1
else
    phi_r = phi_n(kp1)
endif

if (km1 .lt. 1) then
    phi_l = bc0
else
    phi_l = phi_n(km1)
endif

Su = xf(k) - xf(km1)
Sc = xf(kp1) - xf(k)
Sd = xf(kp1+1) - xf(kp1)

phi_f = theta*( (Sd/(Sc+Sd))*phi_c + (Sc/(Sc+Sd))*phi_r) + (1.0 - theta)*( ((Su+2.0*Sc)/(Su+Sc))*phi_c - (Sc/(Su+Sc))*phi_l )

end subroutine general_formula


subroutine bounded_central_scheme(k,phi_f)

use globalvar, only: phi_n, xf, x, eps, bc1, bc0, Nx, beta, gamma

implicit none

integer, intent (in)     :: k
integer                  :: km1,kp1
real, intent (out)       :: phi_f
real                     :: phi_nvd
real                     :: Su,Sd,Sc,phi_c,phi_l,phi_r

km1 = k - 1
kp1 = k + 1

phi_c = phi_n(k)

if (kp1 .gt. Nx) then
    phi_r = bc1
else
    phi_r = phi_n(kp1)
endif

if (km1 .lt. 1) then
    phi_l = bc0
else
    phi_l = phi_n(km1)
endif

Su = xf(k) - xf(km1)
Sc = xf(kp1) - xf(k)
Sd = xf(kp1+1) - xf(kp1)

phi_nvd = (phi_c - phi_l)/(phi_r - phi_l + eps)
! find blending function
gamma = phi_nvd/beta
        
if ((phi_nvd .gt. 1.0) .or. (phi_nvd .lt. 0.0)) then !not bounded, use first order upwind
    phi_f = phi_c
else if ((phi_nvd .ge. beta) .and. (phi_nvd .lt. 1.0)) then ! use central scheme
    phi_f = (Sd/(Sc+Sd))*phi_c + (Sc/(Sc+Sd))*phi_r
else ! use blended CD/UD
    phi_f = (1.0 - gamma*(1.0 - Sd/(Sc+Sd)))*phi_c + gamma*(Sc/(Sc+Sd))*phi_r
endif

end subroutine bounded_central_scheme


subroutine vanalbada(k,phi_f)

use globalvar

implicit none

integer, intent(in)  :: k
real, intent(out)    :: phi_f
real                 :: r,up1,um1,ui,psi,dm,dp

ui = phi_n(k)
um1 = phi_n(k-1)
if (k .eq. Nmax) then
  up1 = bc1
else
  up1 = phi_n(k+1)
endif

dp = up1 - ui
dm = ui-um1
r = 0.0

r = dp/(dm+1.0e-16)
psi = (r*r+r)/(r*r+1.0)

phi_f = ui + 0.5*psi*(up1-ui)

end subroutine vanalbada

subroutine venkatakrishnan(k,phi_f)

use globalvar

implicit none

integer, intent(in)   :: k
real, intent(out)     :: phi_f
real                  :: r,up1,um1,ui,psi,dm,dp
real                  :: num,den,uf,delx,epsv2,KV

delx = dx(2) ! uniform grid
ui = phi_n(k)
um1 = phi_n(k-1)
if (k .eq. Nmax) then
  up1 = bc1
else
  up1 = phi_n(k+1)
endif

uf = ui + ((up1-um1)/(2*delx))*(delx/2)
dm = uf - ui
dp = max(up1,um1,ui)-um1
KV = 0.0 !0.3
epsv2 = (K*delx)**3.0

num = (dp**2.0 + epsv2)*dm + 2*dm*dm*dp
den = dp**2.0 + 2*dm**2.0 + dp*dm + epsv2
psi = (1.0/(dm+1.0e-16))*(num/den)

phi_f = ui + ((up1-um1)/(2*delx))*(delx/2)*psi

end subroutine venkatakrishnan


subroutine make_mesh

use globalvar

implicit none

integer        :: i

open(unit = 12, file = 'diff/grid.txt',status='unknown')

Lmax = 1.0
init_dx = Lmax/(real(Nmax -1))
!init_dx = 0.02
exp_ratio = 1.0

x(1) = 0.0  ! Node
xf(1) = 0.0 ! Face
do i=2,Nmax
    dx(i-1) = init_dx*(exp_ratio)**(i-2)
    x(i) = x(i-1) + dx(i-1)
    xf(i) = 0.5*(x(i-1) + x(i))
    if (x(i) .gt. Lmax) exit
enddo
Nx = min(i,Nmax)
dx(Nx) = x(Nx) - xf(Nx)
print*,Nx
do i=1,Nx-1
 vol(i) = xf(i+1) - xf(i)
enddo
vol(Nx) = x(Nx) - xf(Nx)
do i=1,Nx
 write(12,*) i,x(i), xf(i), dx(i), vol(i)
enddo
close(12)


end subroutine make_mesh

subroutine apply_IC

use globalvar

implicit none

integer         :: i

!IC
phi = 0.0
phi_n = 0.0
open(unit = 11, file = 'diff/ic.txt',status='unknown')
do i=1,Nx
   phi(i) = exp(-((x(i)+0.5)**2)/0.00125)
   write(11,*) x(i),phi(i),(x(i)+0.5)**2/0.00125
enddo
close(11)
print*,t
bc0 = 0.025/sqrt(0.000625+0.02*t)*exp(-(0.5-t)**2/(0.00125+0.04*t)) ! boundary condition
bc1 = 0.025/sqrt(0.000625+0.02*t)*exp(-(1.5-t)**2/(0.00125+0.04*t)) ! boundary condition
phi(1) = bc0
phi_n(1) = bc0
phi(Nx) = bc1
phi_n(Nx) = bc1

phi_n = phi
call compute_gradient

end subroutine apply_IC


subroutine compute_gradient

use globalvar

implicit none

integer  :: i
! Find gradient at faces only
!1-D grid arrangement, k is given as input from the main routine, ('o' is a node, '|' is a face)
!k-1    (w)        k          (e)   k+1
! o------|---------o-----------|-----o
do i=1,Nx-1
  dphi_n(i) = (phi_n(i+1) - phi_n(i))/(x(i+1) - x(i))
enddo

dphi_n(Nx) = (phi_n(i) - phi_n(i-1))/(x(i) - x(i-1))

end subroutine compute_gradient


subroutine compute_diffusion_term(res_v)

use globalvar, only: Nmax,Nx,dphi_n

integer            :: i
real, intent (out) :: res_v(Nmax)

res_v = 0.0
do i = 2,Nx
  res_v(i) = dphi_n(i) - dphi_n(i-1)
enddo


end subroutine compute_diffusion_term

subroutine writeoutput

use globalvar

implicit none

character(len=20)    :: fname
integer              :: i

select case (scheme)
    case (1)
    ! First order upwind
    fname = 'diff/foup.txt'
    case (2)
    ! Second order upwind
    fname = 'diff/soup.txt'
    case (3)
    ! QUICK
    fname = 'diff/quik.txt'
    case (4)
    ! Central
    fname = 'diff/cent.txt'
    case (5)
    ! Bounded central
    fname = 'diff/bcen.txt'
    case (6)
    ! FROMM
    fname = 'diff/frmm.txt'
    case (7)
    ! Second order + van albada
    fname = 'diff/vasl.txt'
    case (8)
    ! Second order + venkatakrishnan
    fname = 'diff/venk.txt'
    case default
    ! Wrong optional
        print*, 'Wrong convective scheme choice, set scheme between 1 to 8'
    stop
end select

!Output files
open(unit = 10, file = fname(1:len_trim(fname)),status='unknown')
open(unit = 11, file = 'diff/exact.txt',status='unknown')
print*,'Tf',Tf
err = 0.0
do i=1,Nx
    exact(i) = 0.025/sqrt(0.000625+0.02*Tf)*exp(-(x(i)+0.5-Tf)**2/(0.00125+0.04*Tf))
    write(10,*) x(i), phi(i), dphi_n(i)
    write(11,*) x(i), exact(i)
    
    err = err + (phi(i) - exact(i))**2.0
enddo
err = sqrt(err/Nx)

print*, err
close(10)
close(11)
end subroutine writeoutput
