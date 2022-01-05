! Solve 1-D scalar advection equation
! \partial_t \phi + U \partial_x \phi = 0
! Use FVM integration
! \int_{\Omega} \partial_t \phi d\Omega = -U\int_{\Omega}\partial_x \phi d\Omega
! \int_{\Omega} \partial_t \phi d\Omega = Res
! Choose desired numerical scheme for spatial disc term Res
! \partial_t \phi \Delta \Omega = Res
! Choose desired numerical scheme for spatial disc term Res
module globalvar

implicit none

integer, parameter  :: Nmax=513, scheme=5
integer             :: Nx
real                :: x(Nmax), xf(Nmax+1), phi(Nmax), phi_n(Nmax), exact(Nmax), dx(Nmax), vol(Nmax)
real                :: dt, U, Lmax, Lx, Tf, t, theta, bc, exp_ratio, init_dx 
real                :: phi_nvd, beta, gamma, eps, err, pi

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
    

Tf = 0.4    ! final time
U = 1.0 ! advection velocity
dt = 0.00005 ! time step
beta = 0.1 ! lower limit

rkcoeff = (/5.d-1, 5.d-1, 1.d0/) ! RK coefficients
nstage = 4


call make_mesh
print*, 'CFL = ',U*dt/minval(dx(1:Nx-1))
domega = dx(2) ! uniform grid for now (Volume)
eps = init_dx/100.0
pi = 4.0*atan(1.0)

call apply_IC


do while (t .le. Tf)
    !BC
    bc = 0.0 ! boundary condition
    phi(1) = bc
    phi_n(1) = bc
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
    
    t = t + dt
    
enddo
print*,t

call writeoutput


end program onedscalar

subroutine compute_residual(res)

! k - index of current element i.e "i" in the main loop
! res - residual of the current element

use globalvar

implicit none

integer                 :: k
real, intent (out)      :: res(Nmax)
real                    :: phi_e, phi_w

! Assuming U > 0
!1-D grid arrangement, k is given as input from the main routine, ('o' is a node, '|' is a face)
!k-1    (w)        k          (e)   k+1
! o------|---------o-----------|-----o

do k=2,Nx
  
  select case (scheme)
    case (1)
    ! First order upwind
        phi_e = phi_n(k)
        phi_w = phi_n(k-1)
    case (2)
    ! Second order upwind
        theta = 0.0
        call general_formula(k,theta,phi_e)
        call general_formula(k-1,theta,phi_w)
    case (3)
    ! QUICK
        theta = 1.0/8.0
        call general_formula(k,theta,phi_e)
        call general_formula(k-1,theta,phi_w)
    case (4)
    ! Central
        theta = 1.0
        call general_formula(k,theta,phi_e)
        call general_formula(k-1,theta,phi_w)
    case (5)
    ! Bounded central
        call bounded_central_scheme(k,phi_e)
        call bounded_central_scheme(k-1,phi_w)
    case (6)
    ! FROMM
        theta = 0.5
        call general_formula(k,theta,phi_e)
        call general_formula(k-1,theta,phi_w)
    case (7)
    ! Second order upwind + van albada
        call vanalbada(k,phi_e)
        call vanalbada(k-1,phi_w)
    case (8)
    ! Second order upwind + venkatakrishnan
        call venkatakrishnan(k,phi_e)
        call venkatakrishnan(k-1,phi_w)
  end select

  res(k) = phi_e - phi_w

end do
res(1) = 0.0

end subroutine compute_residual

subroutine general_formula(k,theta,phi_f)

use globalvar, only: phi_n,xf, bc, Nx

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
    phi_r = phi_n(k)
else
    phi_r = phi_n(kp1)
endif

if (km1 .lt. 1) then
    phi_l = bc
else
    phi_l = phi_n(km1)
endif

Su = xf(k) - xf(km1)
Sc = xf(kp1) - xf(k)
Sd = xf(kp1+1) - xf(kp1)

phi_f = theta*( (Sd/(Sc+Sd))*phi_c + (Sc/(Sc+Sd))*phi_r) + (1.0 - theta)*( ((Su+2.0*Sc)/(Su+Sc))*phi_c - (Sc/(Su+Sc))*phi_l )

end subroutine general_formula


subroutine bounded_central_scheme(k,phi_f)

use globalvar, only: phi_n, xf, x, eps, bc, Nx, beta, gamma

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
    phi_r = phi_n(k)
else
    phi_r = phi_n(kp1)
endif

if (km1 .lt. 1) then
    phi_l = bc
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
  up1 = bc
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
  up1 = bc
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

open(unit = 12, file = 'grid.txt',status='unknown')

Lmax = 1.0
init_dx = Lmax/(real(Nmax -1))
!init_dx = 0.02
exp_ratio = 1.02

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

bc = 0.0 ! boundary condition

!IC
t = 0.0
phi = 0.0
phi_n = 0.0

do i=1,Nx
    if (x(i) .le. 0.25) then
        phi(i) = 1.0 !sin(4.0*pi*x(i))
    else
        phi(i) = 0.0
    endif
enddo
phi_n = phi
phi(1) = bc
phi_n(1) = bc


end subroutine apply_IC


subroutine writeoutput

use globalvar

implicit none

character(len=10)    :: fname
integer             :: i

select case (scheme)
    case (1)
    ! First order upwind
    fname = 'foup.txt'
    case (2)
    ! Second order upwind
    fname = 'soup.txt'
    case (3)
    ! QUICK
    fname = 'quik.txt'
    case (4)
    ! Central
    fname = 'cent.txt'
    case (5)
    ! Bounded central
    fname = 'bcen3.txt'
    case (6)
    ! FROMM
    fname = 'frmm.txt'
    case (7)
    ! Second order + van albada
    fname = 'vasl.txt'
    case (8)
    ! Second order + venkatakrishnan
    fname = 'venk.txt'
    case default
    ! Wrong optional
        print*, 'Wrong convective scheme choice, set scheme between 1 to 8'
    stop
end select

!Output files
open(unit = 10, file = fname(1:len_trim(fname)),status='unknown')
open(unit = 11, file = 'exact.txt',status='unknown')

err = 0.0
do i=1,Nx
    if ((x(i) .gt. U*t) .and. (x(i) .lt. U*t+0.25)) then
        exact(i) = 1.0 !sin(4.0*pi*(x(i)-U*t))
    else
        exact(i) = 0.0
    endif
    write(10,*) x(i), phi(i)
    write(11,*) x(i), exact(i)
    
    err = err + (phi(i) - exact(i))**2.0
enddo
err = sqrt(err/Nx)

print*, err
close(10)
close(11)
end subroutine writeoutput
