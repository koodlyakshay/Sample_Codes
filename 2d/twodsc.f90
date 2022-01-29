! Solve 2-D scalar advection equation
! \partial_t \phi + Ux \partial_x \phi + Uy \partial_y \phi = 0
! Use FVM integration
! \int_{\Omega} \partial_t \phi d\Omega = -\int_{\Omega}(Ux\partial_x \phi + Uy\partial_y \phi )d\Omega
! \int_{\Omega} \partial_t \phi d\Omega = Res
! Choose desired numerical scheme for spatial disc term Res
! \partial_t \phi \Delta \Omega = Res
! RK4 time discretization
!
! Exact solution phi = sin(2pi(x-Ux.t))sin(2pi(y-Uy.t))
! Periodic BC
module globalvar

implicit none

integer, parameter  :: Nmax=129, scheme=5
integer             :: Nx,Ny
real                :: x(Nmax), y(Nmax), xf(Nmax+1), yf(Nmax+1), dx(Nmax), dy(Nmax),vol(Nmax,Nmax)
real                :: phi(Nmax,Nmax), phi_n1(Nmax,Nmax), exact(Nmax,Nmax)
real                :: dphidx(Nmax,Nmax), dphidy(Nmax,Nmax)
real                :: dt, Ux, Uy, Lxmax, Lx, Lymax, Ly, exp_ratio_x, exp_ratio_y, init_dx, init_dy
real                :: phi_nvd, beta, gamma, eps, err, pi, Tf, t, theta, bc_x(Nmax), bc_y(Nmax)

! scheme - type of convective scheme
!           1 - first order upwind
!           2 - second order upwind
!           3 - QUICK
!           4 - central
!           5 - bounded central 
!           6 - FROMM
!           7 - Second order upwind + Van Albada limiter


end module globalvar

program twodscalar


use globalvar

implicit none

integer   :: i,j,k,stage,nstage
real      :: phi_e, phi_w, domega
real      :: resid(Nmax,Nmax),rkstage(Nmax,Nmax,4),rkcoeff(3)
    

Tf = 0.5     ! final time
Ux = 1.0       ! advection velocity x
Uy = 0.0       ! advection velocity y
dt = 0.00025   ! time step
beta = 0.1     ! lower limit

rkcoeff = (/5.d-1, 5.d-1, 1.d0/) ! RK coefficients
nstage = 4


call make_mesh
print*, 'CFL = ',max(Ux,Uy)*dt/minval(dx(1:Nx-1))
eps = init_dx/100.0
pi = 4.0*atan(1.0)

call apply_IC


do while (t .le. Tf)
    !BC
    call setBC(t)

    !--- Time Integration (RK4) ---!
   
    !-- Stage 1 --!
    call compute_residual(resid)
    rkstage(:,:,1) = dt*resid(:,:)/(1.0*vol(:,:))
   
    !-- Stage 2-4 --!
    do stage = 2,nstage

      phi(:,:) = phi_n1(:,:) - rkcoeff(stage-1)*rkstage(:,:,stage-1)

      call compute_residual(resid)

      rkstage(:,:,stage) = dt*resid(:,:)/(1.0*vol(:,:))
    enddo

    phi_n1(:,:) = phi_n1(:,:) - 1.d0/6.d0*(rkstage(:,:,1) + 2.0*rkstage(:,:,2) + 2.0*rkstage(:,:,3) + rkstage(:,:,4))
    
    t = t + dt
    
enddo
print*,t,Tf/dt

call writeexactsol(t)

call writeoutput

end program twodscalar

subroutine compute_residual(res)

! k - index of current element i.e "i" in the main loop
! res - residual of the current element

use globalvar

implicit none

integer                 :: k,l,dir
real, intent (out)      :: res(Nmax,Nmax)
real                    :: phi_e, phi_w, phi_nr, phi_s
real                    :: se, sw, sn, ss

! Assuming U > 0
!2-D grid arrangement, k is given as input from the main routine, ('o' is a node, '|' is a face)
!l-1,k-1(w)       l-1,k       (e)   l-1,k+1
! o------|---------o-----------|-----o
!        |                     | 
!  ------|--------(n)----------|------
!        |                     | 
!l,k-1  (w)       l,k         (e)   l,k+1
! o------|---------o-----------|-----o
!        |                     | 
!  ------|--------(s)----------|------
!        |                     | 
! o------|---------o-----------|-----o
!l+1,k-1(w)       l+1,k       (e)   l+1,k+1

select case (scheme)
  case (1)
  ! First order upwind
    do k=2,Nx
      do l=2,Ny
        phi_e = phi_n1(k,l)    ! East face      __(n)__
        phi_w = phi_n1(k-1,l)  ! West face     |       |
        phi_nr = phi_n1(k,l)   ! North face   (w)  o  (e)
        phi_s = phi_n1(k,l-1)  ! South face    |__(s)__|
        
        se = yf(l) - yf(l-1)
        sw = yf(l) - yf(l-1)
        sn = xf(k) - xf(k-1)
        ss = xf(k) - xf(k-1)
        
        res(k,l) = Ux*(phi_e*se - phi_w*sw) + Uy*(phi_nr*sn - phi_s*ss)
      enddo
    enddo
  case (2)
  ! Second order upwind
    do k=2,Nx
      do l=2,Ny
        theta = 0.0
        call general_formula(k,l,theta,phi_e,1)    ! East face
        call general_formula(k-1,l,theta,phi_w,1)  ! West face
        call general_formula(k,l,theta,phi_nr,2)   ! North face
        call general_formula(k,l-1,theta,phi_s,2)  ! South face
        
        se = yf(l) - yf(l-1)
        sw = yf(l) - yf(l-1)
        sn = xf(k) - xf(k-1)
        ss = xf(k) - xf(k-1)
        
        res(k,l) = Ux*(phi_e*se - phi_w*sw) + Uy*(phi_nr*sn - phi_s*ss)
      enddo
    enddo
  case (3)
  ! QUICK
    do k=2,Nx
      do l=2,Ny
        theta = 1.0/8.0
        call general_formula(k,l,theta,phi_e,1)
        call general_formula(k-1,l,theta,phi_w,1)
        call general_formula(k,l,theta,phi_nr,2)
        call general_formula(k,l-1,theta,phi_s,2)
        
        se = yf(l) - yf(l-1)
        sw = yf(l) - yf(l-1)
        sn = xf(k) - xf(k-1)
        ss = xf(k) - xf(k-1)
        
        res(k,l) = Ux*(phi_e*se - phi_w*sw) + Uy*(phi_nr*sn - phi_s*ss)
      enddo
    enddo
  case (4)
  ! Central
    do k=2,Nx
      do l=2,Ny
        theta = 1.0
        call general_formula(k,l,theta,phi_e,1)
        call general_formula(k-1,l,theta,phi_w,1)
        call general_formula(k,l,theta,phi_nr,2)
        call general_formula(k,l-1,theta,phi_s,2)
        
        se = yf(l) - yf(l-1)
        sw = yf(l) - yf(l-1)
        sn = xf(k) - xf(k-1)
        ss = xf(k) - xf(k-1)
        
        res(k,l) = Ux*(phi_e*se - phi_w*sw) + Uy*(phi_nr*sn - phi_s*ss)
      enddo
    enddo
  case (5)
  ! Bounded central
    do k=2,Nx
      do l=2,Ny
        call bounded_central_scheme(k,l,phi_e,1)
        call bounded_central_scheme(k-1,l,phi_w,1)
        call general_formula(k,l,theta,phi_nr,2)
        call general_formula(k,l-1,theta,phi_s,2)
        
        se = yf(l) - yf(l-1)
        sw = yf(l) - yf(l-1)
        sn = xf(k) - xf(k-1)
        ss = xf(k) - xf(k-1)
        
        res(k,l) = Ux*(phi_e*se - phi_w*sw) + Uy*(phi_nr*sn - phi_s*ss)
      enddo
    enddo
  case (6)
  ! FROMM
    do k=2,Nx
      do l=2,Ny
        theta = 0.5
        call general_formula(k,l,theta,phi_e,1)
        call general_formula(k-1,l,theta,phi_w,1)
        call general_formula(k,l,theta,phi_nr,2)
        call general_formula(k,l-1,theta,phi_s,2)
        
        se = yf(l) - yf(l-1)
        sw = yf(l) - yf(l-1)
        sn = xf(k) - xf(k-1)
        ss = xf(k) - xf(k-1)
        
        res(k,l) = Ux*(phi_e*se - phi_w*sw) + Uy*(phi_nr*sn - phi_s*ss)
      enddo
    enddo
  case (7)
  ! van albada
    do k=2,Nx
      do l=2,Ny
        call vanalbada(k,l,phi_e,1)
        call vanalbada(k-1,l,phi_w,1)
        call vanalbada(k,l,phi_nr,2)
        call vanalbada(k,l-1,phi_s,2)
         
        se = yf(l) - yf(l-1)
        sw = yf(l) - yf(l-1)
        sn = xf(k) - xf(k-1)
        ss = xf(k) - xf(k-1)
        
        res(k,l) = Ux*(phi_e*se - phi_w*sw) + Uy*(phi_nr*sn - phi_s*ss)
      enddo
    enddo
    case default
    ! Wrong option
        print*, 'Wrong convective scheme choice, set scheme between 1 to 6'
    stop
end select

res(1,:) = 0.0
res(:,1) = 0.0

end subroutine compute_residual

subroutine general_formula(k,l,theta,phi_f,dir)

use globalvar, only: phi_n1,xf,yf,Nx,Ny

implicit none

integer, intent(in)         :: k,l,dir
integer                     :: km1,kp1,lm1,lp1
real, intent(in)            :: theta
real, intent(out)           :: phi_f
real                        :: Su,Sd,Sc,phi_c,phi_l,phi_r

if (.not.((dir .eq. 1) .or. (dir .eq. 2))) then
  print*,'Wrong dir specified, should be either 1(x) or 2(y)'
  stop
endif

if (dir .eq. 1) then
  km1 = k - 1
  kp1 = k + 1
  lm1 = l
  lp1 = l
else
  km1 = k
  kp1 = k
  lm1 = l - 1
  lp1 = l + 1
endif

!Center point
phi_c = phi_n1(k,l)

!The right(x)/north(y) node, check if an index is above limits
if (kp1 .gt. Nx) kp1 = k
if (lp1 .gt. Ny) lp1 = l
phi_r = phi_n1(kp1,lp1)

!The left(x)/south(y) node, check if an index is above limits
if (km1 .lt. 1) km1 = k
if (lm1 .lt. 1) lm1 = l
phi_l = phi_n1(km1,lm1)


if (dir .eq. 1) then
  Su = xf(k) - xf(km1)
  Sc = xf(kp1) - xf(k)
  Sd = xf(kp1+1) - xf(kp1)
else 
  Su = yf(l) - yf(lm1)
  Sc = yf(lp1) - yf(l)
  Sd = yf(lp1+1) - yf(lp1)
endif
phi_f = theta*( (Sd/(Sc+Sd))*phi_c + (Sc/(Sc+Sd))*phi_r) + (1.0 - theta)*( ((Su+2.0*Sc)/(Su+Sc))*phi_c - (Sc/(Su+Sc))*phi_l )

end subroutine general_formula


subroutine bounded_central_scheme(k,l,phi_f,dir)

use globalvar, only: phi_n1, xf, x, eps, Nx, beta, gamma, Ny, yf

implicit none

integer, intent (in)     :: k,l,dir
integer                  :: km1,kp1,lm1,lp1
real, intent (out)       :: phi_f
real                     :: phi_nvd
real                     :: Su,Sd,Sc,phi_c,phi_l,phi_r

if (.not.((dir .eq. 1) .or. (dir .eq. 2))) then
  print*,'Wrong dir specified, should be either 1(x) or 2(y)'
  stop
endif

if (dir .eq. 1) then
  km1 = k - 1
  kp1 = k + 1
  lm1 = l
  lp1 = l
else
  km1 = k
  kp1 = k
  lm1 = l - 1
  lp1 = l + 1
endif

!Center point
phi_c = phi_n1(k,l)

!The right(x)/north(y) node, check if an index is above limits
if (kp1 .gt. Nx) kp1 = k
if (lp1 .gt. Ny) lp1 = l
phi_r = phi_n1(kp1,lp1)

!The left(x)/south(y) node, check if an index is above limits
if (km1 .lt. 1) km1 = k
if (lm1 .lt. 1) lm1 = l
phi_l = phi_n1(km1,lm1)

if (dir .eq. 1) then
  Su = xf(k) - xf(km1)
  Sc = xf(kp1) - xf(k)
  Sd = xf(kp1+1) - xf(kp1)
else 
  Su = yf(l) - yf(lm1)
  Sc = yf(lp1) - yf(l)
  Sd = yf(lp1+1) - yf(lp1)
endif

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

subroutine vanalbada(k,l,phi_f,dir)

use globalvar, only: phi_n1,Nx,Ny

implicit none

integer, intent(in)  :: k,l,dir
integer              :: km1,kp1,lm1,lp1
real, intent(out)    :: phi_f
real                 :: r,up1,um1,ui,psi,dm,dp

if (.not.((dir .eq. 1) .or. (dir .eq. 2))) then
  print*,'Wrong dir specified, should be either 1(x) or 2(y)'
  stop
endif

if (dir .eq. 1) then
  km1 = k - 1
  kp1 = k + 1
  lm1 = l
  lp1 = l
else
  km1 = k
  kp1 = k
  lm1 = l - 1
  lp1 = l + 1
endif

!Center point
ui = phi_n1(k,l)

!The right(x)/north(y) node, check if an index is above limits
if (kp1 .gt. Nx) kp1 = k
if (lp1 .gt. Ny) lp1 = l

!The left(x)/south(y) node, check if an index is above limits
if (km1 .lt. 1) km1 = k
if (lm1 .lt. 1) lm1 = l

if (dir .eq. 1) then
  um1 = phi_n1(km1,l)
  up1 = phi_n1(kp1,l)
else
  um1 = phi_n1(k,lm1)
  up1 = phi_n1(k,lp1)
endif

dp = up1 - ui
dm = ui-um1
r = 0.0

r = dp/(dm+1.0e-16)
psi = (r*r+r)/(r*r+1.0)

phi_f = ui + 0.5*psi*(up1-ui)

end subroutine vanalbada

!subroutine compute_gradient

!use globalvar

!implicit none

!integer  :: i,j
!! Find gradient at faces only
!
!do i=1,Nx-1
!  do j=1,Ny-1
!    dphidx(i,j) = (phi_n1(i+1,j) - phi_n1(i,j))/(x(i+1) - x(i))
!    dphidy(i,j) = (phi_n1(i,j+1) - phi_n1(i,j))/(y(j+1) - y(j))
!  enddo
!enddo

!dphidx(Nx,:) = (phi_n1(Nx,:) - phi_n1(Nx-1,:))/(x(Nx) - x(Nx-1))
!dphidy(:,Ny) = (phi_n1(:,Ny) - phi_n1(:,Ny-1))/(y(Ny) - y(Ny-1))

!end subroutine compute_gradient


!subroutine compute_diffusion_term(res_v)

!use globalvar, only: Nmax,Nx,dphi_n

!integer            :: i,j
!real, intent (out) :: res_v(Nmax,Nmax)

!res_v = 0.0
!do i=2,Nx
!  do j=2,Ny
!    res_v(i,j) = dphi_n(i) - dphi_n(i-1)
!  enddo
!enddo


!end subroutine compute_diffusion_term


subroutine make_mesh

use globalvar

implicit none

integer        :: i,j

open(unit = 12, file = 'grid.txt',status='unknown')

Lxmax = 1.0
Lymax = 1.0
init_dx = Lxmax/(real(Nmax -1))
init_dy = Lymax/(real(Nmax -1))
!init_dx = 0.02
exp_ratio_x = 1.00
exp_ratio_y = 1.00

x(1) = 0.0  ! Node
xf(1) = 0.0 ! Face
do i=2,Nmax
    dx(i-1) = init_dx*(exp_ratio_x)**(i-2)
    x(i) = x(i-1) + dx(i-1)
    xf(i) = 0.5*(x(i-1) + x(i))
    if (x(i) .gt. Lxmax) exit
enddo
Nx = min(i,Nmax)
xf(Nx+1) = x(Nx)
dx(Nx) = x(Nx) - xf(Nx)
print*,Nx

y(1) = 0.0  ! Node
yf(1) = 0.0 ! Face
do i=2,Nmax
    dy(i-1) = init_dy*(exp_ratio_y)**(i-2)
    y(i) = y(i-1) + dy(i-1)
    yf(i) = 0.5*(y(i-1) + y(i))
    if (y(i) .gt. Lymax) exit
enddo
Ny = min(i,Nmax)
yf(Ny+1) = y(Ny)
dy(Ny) = y(Ny) - yf(Ny)
print*,Ny

do i=1,Nx
  do j=1,Ny
    vol(i,j) = (xf(i+1) - xf(i))*(yf(j+1) - yf(j))
   enddo
enddo

do i=1,Nx
 write(12,*) i,x(i), xf(i), dx(i)!, vol(i)
 write(12,*) i,y(i), yf(i), dy(i)!, vol(i)
enddo
close(12)


end subroutine make_mesh

subroutine apply_IC

use globalvar, only: Nx,Ny,phi,phi_n1,x,y,pi

implicit none

real             :: t, alpha, beta, gam
integer          :: i,j

alpha =-1.0
beta = 0.5
gam = 0.5
phi = 0.0
do i=1,Nx
  do j=1,Ny
    if ((x(i) .gt. 0.1) .and. (x(i) .lt. 0.2)) phi(i,j) = 1.0
  enddo
enddo

!Apply BC
i=1 !x=0
do j=1,Ny
  phi(i,j) = 0.0 
enddo
j=1 !y=0
do i=1,Nx
  phi(i,j) = 0.0 
enddo
phi_n1 = phi


end subroutine apply_IC

subroutine setBC(time)

use globalvar, only: Nx,Ny,Nmax,phi,phi_n1,x,y

implicit none

real,intent(in)  :: time
real             :: alpha, beta, gam
integer          :: i,j

alpha =-1.0
beta = 0.5
gam = 0.5

!Apply BC 
i=1 !x=0
do j=1,Ny
  phi_n1(1,j) = phi_n1(Nx,j)
enddo
j=1 !y=0
do i=1,Nx
  phi_n1(i,1) = phi_n1(i,Ny)
enddo

phi = phi_n1

end subroutine setBC

subroutine writeexactsol(time)

use globalvar, only: Nx,Ny,Nmax,x,y,Ux,Uy,pi,phi

implicit none

real,intent(in)  :: time
real             :: exact(Nmax,Nmax),err2d(Nmax,Nmax)
real             :: norm, alpha, beta, gam
integer          :: i,j
character(len=20):: fname

alpha =-1.0
beta = 0.5
gam = 0.5
print*,Ux,Uy
print*,time,Ux*time,Ux*time+0.1
exact = 0.0
do i=1,Nx
  do j=1,Ny
    if ((x(i) .gt. 0.1+Ux*time) .and. (x(i) .lt. 0.2+Ux*time)) exact(i,j) = 1.0
  enddo
enddo

fname = 'out/exact.vtk'
call writeheader(fname)

call writevariable(fname, exact)

err2d = abs(exact - phi)
norm = 0.0
do i=1,Nx
  do j=1,Ny
    norm = norm + err2d(i,j)*err2d(i,j)
  enddo
enddo
norm = sqrt(norm)/(Nx*Ny)
print*,norm
end subroutine writeexactsol


subroutine writeoutput

use globalvar

implicit none

character(len=20)    :: fname

select case (scheme)
    case (1)
    ! First order upwind
    fname = 'out/foup.vtk'
    case (2)
    ! Second order upwind
    fname = 'out/soup.vtk'
    case (3)
    ! QUICK
    fname = 'out/quik513.vtk'
    case (4)
    ! Central
    fname = 'out/cent.vtk'
    case (5)
    ! Bounded central
    fname = 'out/bcen.vtk'
    case (6)
    ! FROMM
    fname = 'out/frmm.vtk'
    case (7)
    ! van albada
    fname = 'out/vaal.vtk'
    case default
    ! Wrong option
        print*, 'Wrong convective scheme choice, set scheme between 1 to 6'
    stop
end select

!--- Write header ---!
call writeheader(fname)

!--- Write solution data ---!
call writevariable(fname,phi)

end subroutine writeoutput

subroutine writeheader(fname)

use globalvar

implicit none

character(len=20),intent(in)    :: fname
integer                         :: i,j

!Output files
!Add routines to write output in paraview readable format here
open(unit=100,file=fname(1:len_trim(fname)),status='unknown')

!--- These three lines are compulsory. ---!
write(100,'(a)')'# vtk DataFile Version 3.0'   ! File version and identifier
write(100,'(a)')'vtk output'                   ! Header
write(100,'(a)')'ASCII'                        ! Type of data, only ascii now

!--- Write grid information. ---!
write(100,'(a)')'DATASET STRUCTURED_GRID'      ! Dataset is a keyword
                                               ! We deal with only structured grids

write(100,'(A10,A,I3,A,I3,A,I3)')'DIMENSIONS',achar(9),Nx,achar(9),Ny,achar(9),1
write(100,'(A6,A,I6,A,A6)')'POINTS',achar(9), Nx*Ny,achar(9), 'float'

!--- Write the structured grid data ---!
do i=1,Nx
  do j=1,Ny
    write(100,'(F9.6,F9.6,F9.6)') x(i), y(j), 0.0
  enddo
enddo

close(100)
end subroutine writeheader

subroutine writevariable(fname,var)

use globalvar, only: Nx,Ny,Nmax

implicit none

character(len=20),intent(in)    :: fname
real,intent(in)                 :: var(Nmax,Nmax)
integer                         :: i,j

!Output files
!Add routines to write output in paraview readable format here
open(unit=100,file=fname(1:len_trim(fname)),status='old',position='append')
write(100,'(A10,A,I6)') 'POINT_DATA',achar(9), Nx*Ny

write(100,'(a)') 'SCALARS Phi double 1'
write(100,'(a)') 'LOOKUP_TABLE default'

do i=1,Nx
  do j=1,Ny
    write(100,'(F10.6)') var(i,j)
  enddo
enddo

close(100)
end subroutine writevariable