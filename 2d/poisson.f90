!Solve 2-D Poisson equation
!\partial_x^2 \phi + \partial_y^2 \phi = Q
!Finite volume discretization
! Analytical solution, u = exp(pi*x)sin(pi*y) + (1/2)*(xy)^2
! Q = x^2 + y^2
! Uses LAPACK routine to solve resulting system of equations,
! Compile with 
module globalvar

implicit none

integer, parameter    :: Nmax=33 ! 257 gives memory issues
integer               :: Nx,Ny, Npoint, Nmatrix
real*8                :: x(Nmax), y(Nmax), xf(Nmax+1), yf(Nmax+1), dx(Nmax), dy(Nmax),vol(Nmax,Nmax)
real*8                :: phi(Nmax,Nmax), exact(Nmax,Nmax)
real*8                :: mu
real*8                :: Lxmax, Lx, Lymax, Ly, exp_ratio_x, exp_ratio_y, init_dx, init_dy
real*8                :: err, pi

end module globalvar


program poisson

use globalvar

implicit none

pi = atan(1.0)*4.0
mu = 1.0
print*,'Make grid'
call make_mesh

Npoint = Nx*Ny
print*,'Total points',Npoint
Nmatrix = ((Nx-2)*(Ny-2))
print*,'Interior points',Nmatrix

call setBC

call solve_system(Nmatrix)

call writeoutput

call writeexactsol

endprogram poisson

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

!Subroutine to set BC
! u(x,0) = 0, u(x,1) = x^2/2
! u(0,y) = sin(pi*y) = u(1,y) = exp(pi)sin(pi*y) + y^2/2
subroutine setBC

use globalvar

implicit none

integer   :: i,j
real*8    :: siny

!x = 0, i = 1
!x = 1, i = Nx
do j=1,Ny
  siny = sin(pi*y(j))
  phi(1,j) = siny
  phi(Nx,j) = exp(pi)*siny + 0.5*y(j)**2.0
enddo

!y = 0, j = 1
!y = 1, j = Ny
j=1
do i=1,Nx
  phi(i,1) = 0.0
  phi(i,Ny) = 0.5*x(i)**2.0
enddo


end subroutine setBC

! Subroutine to assemble the coefficient matrix and RHS
! Will also modify the RHS to accomodate the RHS


subroutine solve_system(Nsize)

use globalvar

implicit none

integer,intent(in)   :: Nsize
integer              :: i,j,ii,jj
integer              :: sPoint,ePoint,nrPoint,wPoint,iPoint,iMatrix
real*8               :: lCoeff(Nsize,Nsize), lRHS(Nsize)
real*8               :: ae,aw,an,as,ap
real*8               :: se,sw,sn,ss

lRHS   = 0.0
lCoeff = 0.0

! First calculate the RHS using contribution from Q = x^2 + y^2 only
! Since Dirichlet BC are applied on all boundaries, no need to solve for them
iMatrix = 1
jj = 1
do j=2,Ny-1
  ii = 1
  do i=2,Nx-1
    iPoint = ii + (jj-1)*(Ny-2)
    if (iPoint .ne. iMatrix) print*,'Wrong order'
    lRHS(iMatrix) = (x(i)**2 + y(i)**2.0)*vol(i,j)
    iMatrix = iMatrix + 1
    ii = ii + 1
  enddo
  jj = jj +1
enddo

!Check if there is an error in computing number of nodes that have to be solved
if (iMatrix-1 .ne. Nsize) then
 print*,'Nmatrix and Nsize dont match (Nsize and Nmatrix)', Nsize, iMatrix-1
 STOP
endif

print*,'Now assembling coefficient matrix'
print*,'First interior elements that dont neighbor any boundaries'
ii = 2
do i=3,Nx-2
  jj = 2
  do j=3,Ny-2
    ! Find face areas
    se = yf(j) - yf(j-1)
    sw = yf(j) - yf(j-1)
    sn = xf(i) - xf(i-1)
    ss = xf(i) - xf(i-1)    
    ! Compute coefficients
    ae = mu*se/(x(i+1) - x(i))
    aw = mu*sw/(x(i) - x(i-1))
    an = mu*sn/(y(j+1) - y(j))
    as = mu*ss/(y(j) - y(j-1))
    ap = -(ae + aw + an + as)
    ! Compute locations/indices in the matrix
    iPoint  = ii + (jj-1)*(Ny-2)
    sPoint  = ii + (jj-1-1)*(Ny-2)
    ePoint  = ii+1 + (jj-1)*(Ny-2)
    nrPoint = ii + (jj-1+1)*(Ny-2)
    wPoint  = ii-1 + (jj-1)*(Ny-2)
    ! Now assign the values
    lCoeff(iPoint,sPoint)  = as
    lCoeff(iPoint,ePoint)  = ae
    lCoeff(iPoint,nrPoint) = an
    lCoeff(iPoint,wPoint)  = aw
    lCoeff(iPoint,iPoint)  = ap
    ! No need to update RHS as all neighbors are interior
    ! Increse the index in jj
    jj = jj + 1
  enddo
  ! Increase index in ii
  ii = ii + 1
enddo

print*,'Now points that neighbor the boundaries'
print*,'South neighbor on boundaries'
j = 2
jj = 1
ii = 2
do i=3,Nx-2
! Find face areas
  se = yf(j) - yf(j-1)
  sw = yf(j) - yf(j-1)
  sn = xf(i) - xf(i-1)
  ss = xf(i) - xf(i-1)    
! Compute coefficients
  ae = mu*se/(x(i+1) - x(i))
  aw = mu*sw/(x(i) - x(i-1))
  an = mu*sn/(y(j+1) - y(j))
  as = mu*ss/(y(j) - y(j-1))
  ap = -(ae + aw + an + as)
! Compute locations/indices in the matrix
  iPoint  = ii + (jj-1)*(Ny-2)
  sPoint  = ii + (jj-1-1)*(Ny-2)
  ePoint  = ii+1 + (jj-1)*(Ny-2)
  nrPoint = ii + (jj-1+1)*(Ny-2)
  wPoint  = ii-1 + (jj-1)*(Ny-2)
! Now assign the values
  lCoeff(iPoint,ePoint)  = ae
  lCoeff(iPoint,nrPoint) = an
  lCoeff(iPoint,wPoint)  = aw
  lCoeff(iPoint,iPoint)  = ap
! Update RHS
  lRHS(iPoint) = lRHS(iPoint) - as*phi(i,j-1)
! Increse the index in ii
  ii = ii + 1
enddo

print*,'East neighbor on boundaries'
i = Nx-1
ii = Nx-2
jj = 2
do j=3,Ny-2
! Find face areas
  se = yf(j) - yf(j-1)
  sw = yf(j) - yf(j-1)
  sn = xf(i) - xf(i-1)
  ss = xf(i) - xf(i-1)    
! Compute coefficients
  ae = mu*se/(x(i+1) - x(i))
  aw = mu*sw/(x(i) - x(i-1))
  an = mu*sn/(y(j+1) - y(j))
  as = mu*ss/(y(j) - y(j-1))
  ap = -(ae + aw + an + as)
! Compute locations/indices in the matrix
  iPoint  = ii + (jj-1)*(Ny-2)
  sPoint  = ii + (jj-1-1)*(Ny-2)
  nrPoint = ii + (jj-1+1)*(Ny-2)
  wPoint  = ii-1 + (jj-1)*(Ny-2)
! Now assign the values
  lCoeff(iPoint,sPoint)  = as
  lCoeff(iPoint,nrPoint) = an
  lCoeff(iPoint,wPoint)  = aw
  lCoeff(iPoint,iPoint)  = ap
! Update RHS
  lRHS(iPoint) = lRHS(iPoint) - ae*phi(i+1,j)
! Increse the index in jj
  jj = jj + 1
enddo

print*,'North neighbor on boundaries'
j = Ny-1
jj = Ny-2
ii = 2
do i=3,Nx-2
! Find face areas
  se = yf(j) - yf(j-1)
  sw = yf(j) - yf(j-1)
  sn = xf(i) - xf(i-1)
  ss = xf(i) - xf(i-1)    
! Compute coefficients
  ae = mu*se/(x(i+1) - x(i))
  aw = mu*sw/(x(i) - x(i-1))
  an = mu*sn/(y(j+1) - y(j))
  as = mu*ss/(y(j) - y(j-1))
  ap = -(ae + aw + an + as)
! Compute locations/indices in the matrix
  iPoint  = ii + (jj-1)*(Ny-2)
  sPoint  = ii + (jj-1-1)*(Ny-2)
  ePoint  = ii+1 + (jj-1)*(Ny-2)
  nrPoint = ii + (jj-1+1)*(Ny-2)
  wPoint  = ii-1 + (jj-1)*(Ny-2)
! Now assign the values
  lCoeff(iPoint,sPoint)  = as
  lCoeff(iPoint,ePoint)  = ae
  lCoeff(iPoint,wPoint)  = aw
  lCoeff(iPoint,iPoint)  = ap
! Update RHS
  lRHS(iPoint) = lRHS(iPoint) - an*phi(i,j+1)
! Increase index in ii
  ii = ii + 1
enddo

print*,'West neighbor on boundaries'
i = 2
ii = 1
jj = 2
do j=3,Ny-2
! Find face areas
  se = yf(j) - yf(j-1)
  sw = yf(j) - yf(j-1)
  sn = xf(i) - xf(i-1)
  ss = xf(i) - xf(i-1)    
! Compute coefficients
  ae = mu*se/(x(i+1) - x(i))
  aw = mu*sw/(x(i) - x(i-1))
  an = mu*sn/(y(j+1) - y(j))
  as = mu*ss/(y(j) - y(j-1))
  ap = -(ae + aw + an + as)
! Compute locations/indices in the matrix
  iPoint  = ii + (jj-1)*(Ny-2)
  sPoint  = ii + (jj-1-1)*(Ny-2)
  ePoint  = ii+1 + (jj-1)*(Ny-2)
  nrPoint = ii + (jj-1+1)*(Ny-2)
! Now assign the values
  lCoeff(iPoint,sPoint)  = as
  lCoeff(iPoint,ePoint)  = ae
  lCoeff(iPoint,nrPoint) = an
  lCoeff(iPoint,iPoint)  = ap
! Update RHS
  lRHS(iPoint) = lRHS(iPoint) - aw*phi(i-1,j)
! Increse index in jj 
  jj = jj + 1
enddo

print*,'Corner points'
i=2
j=2
jj=1
ii=1
! Find face areas
se = yf(j) - yf(j-1)
sw = yf(j) - yf(j-1)
sn = xf(i) - xf(i-1)
ss = xf(i) - xf(i-1)    
! Compute coefficients
ae = mu*se/(x(i+1) - x(i))
aw = mu*sw/(x(i) - x(i-1))
an = mu*sn/(y(j+1) - y(j))
as = mu*ss/(y(j) - y(j-1))
ap = -(ae + aw + an + as)
! Compute locations/indices in the matrix
iPoint  = ii + (jj-1)*(Ny-2)
sPoint  = ii + (jj-1-1)*(Ny-2)
ePoint  = ii+1 + (jj-1)*(Ny-2)
nrPoint = ii + (jj-1+1)*(Ny-2)
wPoint  = ii-1 + (jj-1)*(Ny-2)
! Now assign the values
lCoeff(iPoint,ePoint)  = ae
lCoeff(iPoint,nrPoint) = an
lCoeff(iPoint,iPoint)  = ap
! Update RHS
lRHS(iPoint) = lRHS(iPoint) - as*phi(i,j-1) - aw*phi(i-1,j)

i=Nx-1
j=2
jj=1
ii=Nx-2
! Find face areas
se = yf(j) - yf(j-1)
sw = yf(j) - yf(j-1)
sn = xf(i) - xf(i-1)
ss = xf(i) - xf(i-1)    
! Compute coefficients
ae = mu*se/(x(i+1) - x(i))
aw = mu*sw/(x(i) - x(i-1))
an = mu*sn/(y(j+1) - y(j))
as = mu*ss/(y(j) - y(j-1))
ap = -(ae + aw + an + as)
! Compute locations/indices in the matrix
iPoint  = ii + (jj-1)*(Ny-2)
sPoint  = ii + (jj-1-1)*(Ny-2)
ePoint  = ii+1 + (jj-1)*(Ny-2)
nrPoint = ii + (jj-1+1)*(Ny-2)
wPoint  = ii-1 + (jj-1)*(Ny-2)
! Now assign the values
lCoeff(iPoint,wPoint)  = aw
lCoeff(iPoint,nrPoint) = an
lCoeff(iPoint,iPoint)  = ap
! Update RHS
lRHS(iPoint) = lRHS(iPoint) - as*phi(i,j-1) - ae*phi(i+1,j)

i=Nx-1
j=Ny-1
ii=Nx-2
jj=Ny-2
! Find face areas
se = yf(j) - yf(j-1)
sw = yf(j) - yf(j-1)
sn = xf(i) - xf(i-1)
ss = xf(i) - xf(i-1)    
! Compute coefficients
ae = mu*se/(x(i+1) - x(i))
aw = mu*sw/(x(i) - x(i-1))
an = mu*sn/(y(j+1) - y(j))
as = mu*ss/(y(j) - y(j-1))
ap = -(ae + aw + an + as)
! Compute locations/indices in the matrix
iPoint  = ii + (jj-1)*(Ny-2)
sPoint  = ii + (jj-1-1)*(Ny-2)
ePoint  = ii+1 + (jj-1)*(Ny-2)
nrPoint = ii + (jj-1+1)*(Ny-2)
wPoint  = ii-1 + (jj-1)*(Ny-2)
! Now assign the values
lCoeff(iPoint,sPoint)  = as
lCoeff(iPoint,wPoint)  = aw
lCoeff(iPoint,iPoint)  = ap
! Update RHS
lRHS(iPoint) = lRHS(iPoint) - an*phi(i,j+1) - ae*phi(i+1,j)

i=2
j=Ny-1
ii=1
jj=Ny-2
! Find face areas
se = yf(j) - yf(j-1)
sw = yf(j) - yf(j-1)
sn = xf(i) - xf(i-1)
ss = xf(i) - xf(i-1)    
! Compute coefficients
ae = mu*se/(x(i+1) - x(i))
aw = mu*sw/(x(i) - x(i-1))
an = mu*sn/(y(j+1) - y(j))
as = mu*ss/(y(j) - y(j-1))
ap = -(ae + aw + an + as)
! Compute locations/indices in the matrix
iPoint  = ii + (jj-1)*(Ny-2)
sPoint  = ii + (jj-1-1)*(Ny-2)
ePoint  = ii+1 + (jj-1)*(Ny-2)
nrPoint = ii + (jj-1+1)*(Ny-2)
wPoint  = ii-1 + (jj-1)*(Ny-2)
! Now assign the values
lCoeff(iPoint,sPoint)  = as
lCoeff(iPoint,ePoint)  = ae
lCoeff(iPoint,iPoint)  = ap
! Update RHS
lRHS(iPoint) = lRHS(iPoint) - an*phi(i,j+1) - aw*phi(i-1,j)


print*,'Now call the routines to solve the system'
call linsolve(Nsize,lCoeff,lRHS)

end subroutine solve_system


! Subroutine to solve the system of equations
! Will use LAPACK routines - first attempt dgesv
! SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!*  Arguments
!*  =========
!*
!*  N       (input) INTEGER
!*          The number of linear equations, i.e., the order of the
!*          matrix A.  N >= 0.
!*
!*  NRHS    (input) INTEGER
!*          The number of right hand sides, i.e., the number of columns
!*          of the matrix B.  NRHS >= 0.
!*
!*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!*          On entry, the N-by-N coefficient matrix A.
!*          On exit, the factors L and U from the factorization
!*          A = P*L*U; the unit diagonal elements of L are not stored.
!*
!*  LDA     (input) INTEGER
!*          The leading dimension of the array A.  LDA >= max(1,N).
!*
!*  IPIV    (output) INTEGER array, dimension (N)
!*          The pivot indices that define the permutation matrix P;
!*          row i of the matrix was interchanged with row IPIV(i).
!*
!*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!*          On entry, the N-by-NRHS matrix of right hand side matrix B.
!*          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
!*
!*  LDB     (input) INTEGER
!*          The leading dimension of the array B.  LDB >= max(1,N).
!*
!*  INFO    (output) INTEGER
!*          = 0:  successful exit
!*          < 0:  if INFO = -i, the i-th argument had an illegal value
!*          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
!*                has been completed, but the factor U is exactly
!*                singular, so the solution could not be computed.
subroutine linsolve(Nsize,lCoeff,lRHS)

use globalvar

implicit none

external             :: dgesv
integer,intent(in)   :: Nsize
real*8,intent(in)    :: lCoeff(Nsize,Nsize), lRHS(Nsize)
integer              :: ipiv(Nsize),info
integer              :: i,j,ii,jj, iMatrix

print*,'Call LAPACK routine to solve the matrix, solution is returned in lRHS'
call dgesv(Nsize, 1, lCoeff, Nsize, ipiv, lRHS, Nsize, info)

print*,'Assign the solution'
iMatrix=1
do j=2,Ny-1
  do i=2,Nx-1
    phi(i,j) = lRHS(iMatrix)
    iMatrix = iMatrix + 1
  enddo
enddo


end subroutine linsolve

subroutine writeexactsol

use globalvar

implicit none

integer              :: i,j
character(len=25)    :: fname

do i=1,Nx
  do j=1,Ny
    exact(i,j) = exp(pi*x(i))*sin(pi*y(j)) + 0.5*(x(i)*y(j))**2
  enddo
enddo

fname = 'out/exactpoisson.vtk'
call writeheader(fname)

call writevariable(fname,exact)

end subroutine writeexactsol


subroutine writeoutput

use globalvar

implicit none

character(len=25)    :: fname

!--- Set file name ---!
fname = 'out/poissonsol.vtk'

!--- Write header ---!
call writeheader(fname)

!--- Write solution data ---!
call writevariable(fname,phi)

end subroutine writeoutput

subroutine writeheader(fname)

use globalvar, only: x,y,Nx,Ny

implicit none

character(len=25),intent(in)    :: fname
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

character(len=25),intent(in)    :: fname
real*8,intent(in)               :: var(Nmax,Nmax)
integer                         :: i,j

!Output files
!Add routines to write output in paraview readable format here
open(unit=100,file=fname(1:len_trim(fname)),status='old',position='append')
write(100,'(A10,A,I6)') 'POINT_DATA',achar(9), Nx*Ny

write(100,'(a)') 'SCALARS Phi double 1'
write(100,'(a)') 'LOOKUP_TABLE default'

do i=1,Nx
  do j=1,Ny
    write(100,'(F15.7)') var(i,j)
  enddo
enddo

close(100)
end subroutine writevariable

